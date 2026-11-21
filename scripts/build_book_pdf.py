#!/usr/bin/env python
"""把 book/ 的 GitBook 源编译成整本 PDF（Tectonic/LaTeX 路线）。

用法（仓库根目录）：
  .venv/Scripts/python.exe scripts/build_book_pdf.py [--skip-figures]
      [--skip-compile] [--merge-cover cover.pdf] [--out 最终文件名]

流程：
  1. 调 scripts/make_figures.py 重出全部配图，写入 dist/figs/*.pdf（矢量）；
  2. 按 book/SUMMARY.md 的顺序把 Markdown 转成 LaTeX，遵守
     book/frontmatter/conventions.md 的约定：$$ 独立块配 \\tag、表格内
     竖线分隔、图号写在 alt 文本里——图号/公式号沿用书中手写编号，
     不由 LaTeX 重新编号；
  3. 调 Tectonic 编译 dist/book_body.tex；
  4. 若给定 --merge-cover，把 ReportLab 封面并成第 1 页并写元数据。

Tectonic 二进制按 PATH → 环境变量 TECTONIC → ~/.zcode-tools/tectonic/
的顺序查找。首次编译会联网下载 TeX 包（缓存于用户目录）。

额外 Python 依赖：pypdf（--merge-cover 时）；配图导出复用 make_figures
已有的 numpy/matplotlib。封面由 pdf 技能的 cover_render.py（reportlab）
单独生成，不在本脚本内。
"""

from __future__ import annotations

import argparse
import os
import pathlib
import re
import shutil
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
BOOK = ROOT / "book"
DIST = ROOT / "dist"
FIGS = DIST / "figs"
SUMMARY = BOOK / "SUMMARY.md"
FINAL_DEFAULT = "从一块电池到求解器.pdf"

# ---------------------------------------------------------------- 行内转换

CODE_SPAN = re.compile(r"`([^`]*)`")
MATH_INLINE = re.compile(r"\$([^$\n]+)\$")
BOLD = re.compile(r"\*\*([^*\n]+)\*\*")
ITALIC = re.compile(r"(?<!\*)\*([^*\n]+)\*(?!\*)")
LINK = re.compile(r"\[([^\]]+)\]\(([^)]+)\)")

_ESCAPES = {
    "\\": r"\textbackslash{}",
    "&": r"\&",
    "%": r"\%",
    "$": r"\$",
    "#": r"\#",
    "_": r"\_",
    "{": r"\{",
    "}": r"\}",
    "~": r"\textasciitilde{}",
    "^": r"\textasciicircum{}",
}


def escape_text(text: str) -> str:
    """LaTeX 特殊字符转义（仅用于普通文本段；math/code 段原样保留）。"""
    return "".join(_ESCAPES.get(ch, ch) for ch in text)


def tex_texttt(content: str) -> str:
    """行内代码：转义后在路径/标识符分隔符后给断行机会，避免顶出页边。"""
    body = escape_text(content)
    body = body.replace(r"\textbackslash{}", r"\textbackslash{}\allowbreak{}")
    body = re.sub(r"([/_:<>,()\-])", lambda m: m.group(1) + r"\allowbreak{}", body)
    return r"\texttt{" + body + "}"


class Inline:
    """行内 Markdown → LaTeX：先摘出 code 与 math，再转义并处理加粗/斜体。"""

    def __init__(self) -> None:
        self.spans: list[str] = []

    def _stash(self, tex: str) -> str:
        self.spans.append(tex)
        return f"\x00{len(self.spans) - 1}\x00"

    def convert(self, source: str) -> str:
        self.spans = []
        # Times/宋体都缺的符号换成数学命令；★☆ 由导言区指给 CJK 字体。
        source = source.replace("✓", r"\ensuremath{\checkmark}").replace(
            "✗", r"\ensuremath{\times}"
        )

        def stash_code(match: re.Match) -> str:
            return self._stash(tex_texttt(match.group(1)))

        def stash_math(match: re.Match) -> str:
            return self._stash("$" + match.group(1) + "$")

        def stash_link(match: re.Match) -> str:
            text, url = match.group(1), match.group(2)
            if url.startswith("http"):
                # \href 的 URL 参数只对 % 和 # 敏感，其余原样即可。
                safe_url = url.replace("%", r"\%").replace("#", r"\#")
                return self._stash(
                    r"\href{" + safe_url + "}{" + escape_text(text) + "}"
                )
            return text  # 仓库内部链接：只保留文字

        text = CODE_SPAN.sub(stash_code, source)
        text = LINK.sub(stash_link, text)
        text = MATH_INLINE.sub(stash_math, text)
        text = escape_text(text)
        text = BOLD.sub(lambda m: r"\textbf{" + m.group(1) + "}", text)
        text = ITALIC.sub(lambda m: r"\emph{" + m.group(1) + "}", text)
        # 放回被摘出的段（占位符里只有数字，不会被上面的步骤破坏）。
        for _ in range(3):
            replaced = re.sub(
                r"\x00(\d+)\x00",
                lambda m: self.spans[int(m.group(1))],
                text,
            )
            if replaced == text:
                break
            text = replaced
        return text


# ---------------------------------------------------------------- 块级解析

IMAGE = re.compile(r"^!\[([^\]]*)\]\(([^)]+)\)\s*$")
PART_LINE = re.compile(r"^##\s+(?!\[)(.+?)\s*$")
ENTRY_LINE = re.compile(r"^\s*\*\s+\[([^\]]+)\]\(([^)]+)\)\s*$")
UL_ITEM = re.compile(r"^(\s*)[-*]\s+(.*)$")
OL_ITEM = re.compile(r"^(\s*)(\d+)[.、]\s+(.*)$")
TABLE_SEP = re.compile(r"^\s*\|?[\s:|-]+\|?\s*$")

LST_LANGS = {
    # 只用 listings 内置语言；bat/yaml/cmake 等无定义，按纯文本排版。
    "cpp": "C++",
    "c++": "C++",
    "python": "Python",
    "bash": "bash",
    "sh": "bash",
}


def split_table_row(line: str) -> list[str]:
    """按不在 $...$ 与 `...` 内的竖线切分表格行（约定用 \\mid，无裸竖线）。"""
    stripped = line.strip()
    if stripped.startswith("|"):
        stripped = stripped[1:]
    if stripped.endswith("|"):
        stripped = stripped[:-1]
    cells: list[str] = []
    buf: list[str] = []
    in_math = in_code = False
    for ch in stripped:
        if ch == "$" and not in_code:
            in_math = not in_math
        elif ch == "`" and not in_math:
            in_code = not in_code
        if ch == "|" and not in_math and not in_code:
            cells.append("".join(buf).strip())
            buf = []
        else:
            buf.append(ch)
    cells.append("".join(buf).strip())
    return cells


def join_para(lines: list[str]) -> str:
    """合并段落内的换行：CJK 边界不加空格，拉丁边界补空格。"""
    if not lines:
        return ""
    out = lines[0]
    for line in lines[1:]:
        if out[-1].isascii() and line[0].isascii():
            out += " " + line
        else:
            out += line
    return out


def heading_level(line: str) -> int:
    match = re.match(r"^(#{1,4})\s+(.*)$", line)
    if not match:
        return 0
    return len(match.group(1))


class MarkdownToLatex:
    """单个章节文件的转换器；块级状态机见 convert()。"""

    def __init__(self) -> None:
        self.inline = Inline()
        self.figure_files: list[str] = []

    # —— 各块的输出 ——

    def _heading(self, level: int, title: str, skip_first: bool) -> list[str]:
        title_tex = self.inline.convert(title.strip())
        if level == 1:
            return [] if skip_first else ["% 未预期的 H1：" + title]
        mapping = {
            2: ("section", "section"),
            3: ("subsection", "subsection"),
            4: ("subsubsection", "subsubsection"),
        }
        if level not in mapping:
            return []
        cmd, toc = mapping[level]
        # \needspace 防孤行：标题不落在页面底部 22% 以内，避免标题与正文分家。
        return [
            "\\needspace{0.22\\textheight}",
            f"\\{cmd}*{{{title_tex}}}",
            f"\\addcontentsline{{toc}}{{{toc}}}{{{title_tex}}}",
        ]

    def _image(self, caption: str, path: str) -> list[str]:
        stem = pathlib.Path(path).stem
        self.figure_files.append(stem)
        cap = self.inline.convert(caption)
        return [
            "\\begin{figure}[H]\\centering",
            "\\includegraphics[max width=\\linewidth,max height=0.42\\textheight]{"
            + stem
            + ".pdf}",
            "\\caption*{" + cap + "}",
            "\\end{figure}",
        ]

    def _table(self, rows: list[list[str]]) -> list[str]:
        ncol = max(len(row) for row in rows)
        size = "\\small" if ncol <= 4 else "\\footnotesize"
        spec = " ".join(["X"] * ncol)
        header = rows[0]
        body = rows[1:]

        def emit_row(cells: list[str]) -> str:
            padded = cells + [""] * (ncol - len(cells))
            return " & ".join(self.inline.convert(c) for c in padded) + r" \\"

        out = [
            "{%s" % size,
            "\\begin{xltabular}{\\linewidth}{%s}" % spec,
            "\\toprule",
            emit_row(header),
            "\\midrule",
            "\\endfirsthead",
            "\\multicolumn{%d}{l}{\\footnotesize（续上表）}\\\\" % ncol,
            "\\toprule",
            emit_row(header),
            "\\midrule",
            "\\endhead",
            "\\bottomrule",
            "\\endlastfoot",
        ]
        out += [emit_row(row) for row in body]
        out += ["\\end{xltabular}}", ""]
        return out

    def _list_items(self, items: list[tuple[int, bool, str]]) -> list[str]:
        """把连续列表项分组为嵌套的 itemize/enumerate（层级按缩进//2）。"""
        out: list[str] = []
        stack: list[bool] = []  # True = enumerate

        def close_to(depth: int) -> None:
            while len(stack) > depth:
                out.append("\\end{%s}" % ("enumerate" if stack.pop() else "itemize"))

        for level, ordered, text in items:
            depth = max(1, level + 1)
            # 同深度上类型切换（少见）：先收掉旧环境再开新的。
            if len(stack) >= depth and stack[depth - 1] is not ordered:
                close_to(depth - 1)
            close_to(depth)
            while len(stack) < depth:
                out.append("\\begin{%s}" % ("enumerate" if ordered else "itemize"))
                stack.append(ordered)
            out.append("\\item " + self.inline.convert(text))
        close_to(0)
        out.append("")
        return out

    def _note(self, paras: list[str]) -> list[str]:
        out = ["\\begin{mdnote}"]
        for para in paras:
            out.append(self.inline.convert(para))
            out.append("")
        if out[-1] == "":
            out.pop()
        out.append("\\end{mdnote}")
        out.append("")
        return out

    # —— 主状态机 ——

    def convert(self, md_text: str) -> str:
        lines = md_text.splitlines()
        out: list[str] = []
        index = 0
        first_h1_seen = False

        def flush_paragraph(paras: list[str]) -> None:
            if paras:
                out.append(self.inline.convert(join_para(paras)))
                out.append("")

        while index < len(lines):
            line = lines[index]
            stripped = line.strip()

            # 代码围栏：原样进 lstlisting。
            if stripped.startswith("```"):
                lang = stripped[3:].strip().lower()
                lang_tex = LST_LANGS.get(lang, "")
                out.append(
                    "\\begin{lstlisting}"
                    + (f"[language={lang_tex}]" if lang_tex else "")
                )
                index += 1
                while index < len(lines) and not lines[index].strip().startswith("```"):
                    out.append(lines[index].rstrip())
                    index += 1
                index += 1  # 跳过收尾 ```
                out.append("\\end{lstlisting}")
                out.append("")
                continue

            # 独立公式块：$$ 独占一行，中间原样复制。
            if stripped == "$$":
                body: list[str] = []
                index += 1
                while index < len(lines) and lines[index].strip() != "$$":
                    body.append(lines[index].rstrip())
                    index += 1
                index += 1
                env = "equation" if any("\\tag{" in b for b in body) else "equation*"
                out.append(f"\\begin{{{env}}}")
                out.extend(body)
                out.append(f"\\end{{{env}}}")
                out.append("")
                continue

            # 表格：表头行 + 分隔行（含 | 且下一行是分隔线）。
            if "|" in stripped and stripped.startswith("|") and index + 1 < len(lines):
                if TABLE_SEP.match(lines[index + 1]) and "-" in lines[index + 1]:
                    rows = [split_table_row(stripped)]
                    index += 2
                    while index < len(lines):
                        candidate = lines[index].strip()
                        if not candidate or "|" not in candidate:
                            break
                        if TABLE_SEP.match(candidate):
                            index += 1
                            continue
                        rows.append(split_table_row(candidate))
                        index += 1
                    out.extend(self._table(rows))
                    continue

            # 列表项（连续收集）。
            items: list[tuple[int, bool, str]] = []
            while index < len(lines):
                candidate = lines[index]
                ul = UL_ITEM.match(candidate)
                ol = OL_ITEM.match(candidate)
                if ul:
                    items.append((len(ul.group(1)) // 2, False, ul.group(2)))
                elif ol:
                    items.append((len(ol.group(1)) // 2, True, ol.group(3)))
                elif candidate.strip() == "" and items:
                    # 列表后的空行：若下一行还是列表项则继续（宽松列表）。
                    if index + 1 < len(lines) and (
                        UL_ITEM.match(lines[index + 1]) or OL_ITEM.match(lines[index + 1])
                    ):
                        index += 1
                        continue
                    break
                else:
                    break
                index += 1
            if items:
                out.extend(self._list_items(items))
                continue

            # 引用块。
            if stripped.startswith(">"):
                quote: list[str] = []
                while index < len(lines) and lines[index].strip().startswith(">"):
                    content = lines[index].strip()[1:].strip()
                    if content:
                        quote.append(content)
                    index += 1
                # 引用内的多行合成一个段落（书里的引用块都很短）。
                out.extend(self._note([join_para(quote)] if quote else [""]))
                continue

            # 图片独占一行。
            img = IMAGE.match(stripped)
            if img:
                out.extend(self._image(img.group(1), img.group(2)))
                index += 1
                continue

            # 标题。
            level = heading_level(line)
            if level:
                if level == 1:
                    first_h1_seen = True
                    index += 1
                    continue  # 章名由 SUMMARY 提供，跳过文件内 H1
                out.append("")
                out.extend(self._heading(level, re.match(r"^#+\s+(.*)$", line).group(1), True))
                out.append("")
                index += 1
                continue

            # 水平线 / 空行。
            if stripped in ("---", "***", "___"):
                index += 1
                continue

            # 普通段落：收集到空行/特殊块为止。
            paras: list[str] = []
            while index < len(lines):
                candidate = lines[index]
                s = candidate.strip()
                if (
                    not s
                    or s.startswith("```")
                    or s == "$$"
                    or s.startswith("|")
                    or s.startswith(">")
                    or heading_level(candidate)
                    or IMAGE.match(s)
                    or UL_ITEM.match(candidate)
                    or OL_ITEM.match(candidate)
                ):
                    break
                paras.append(s)
                index += 1
            flush_paragraph(paras)
            if not paras:
                index += 1

        return "\n".join(out)


# ---------------------------------------------------------------- 文档组装

PREAMBLE = r"""\documentclass[UTF8,a4paper,oneside,openany,10pt,fontset=windows]{ctexbook}
\usepackage{amsmath,amssymb,amsbsy}
\usepackage[left=2.2cm,right=2.2cm,top=2.4cm,bottom=2.4cm]{geometry}
\usepackage{graphicx}
\usepackage[export]{adjustbox}
\usepackage{xltabular,booktabs,array}
\usepackage{float}
\usepackage{caption}
\usepackage{fancyhdr}
\usepackage{listings}
\usepackage{xcolor}
\usepackage{tcolorbox}
\usepackage{enumitem}
\usepackage{needspace}
\setmainfont{Times New Roman}
\setsansfont{Arial}
\setmonofont{Consolas}
% ★☆ 不在 Times New Roman 里，交给中文字体（宋体含这两个字形）。
\xeCJKDeclareCharClass{CJK}{"2605, "2606}

\definecolor{notebg}{HTML}{F1F5FA}
\definecolor{noteframe}{HTML}{9DB3CC}
\newtcolorbox{mdnote}{colback=notebg,colframe=noteframe,boxrule=0.5pt,arc=2pt,
  left=8pt,right=8pt,top=6pt,bottom=6pt}

\captionsetup{labelformat=empty,font=small,skip=6pt}
\lstset{
  basicstyle=\ttfamily\footnotesize,
  keywordstyle=\color{blue!55!black},
  commentstyle=\color{green!35!black}\itshape,
  stringstyle=\color{orange!75!black},
  showstringspaces=false,
  breaklines=true,
  columns=fullflexible,
  keepspaces=true,
  frame=single,
  rulecolor=\color{black!15},
  backgroundcolor=\color{black!3},
  xleftmargin=1em, xrightmargin=1em,
  aboveskip=0.8em, belowskip=0.9em,
}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small\leftmark}
\fancyhead[R]{\small\thepage}
\renewcommand{\headrulewidth}{0.3pt}
\setlength{\headheight}{14pt}

\widowpenalty=10000
\clubpenalty=10000
\displaywidowpenalty=10000
\emergencystretch=2em
\setcounter{tocdepth}{1}

\usepackage[colorlinks=true,linkcolor=blue!45!black,citecolor=blue!45!black,
  urlcolor=blue!45!black,bookmarksnumbered=true,unicode=true,
  pdftitle={从一块电池到求解器},pdfauthor={zhoukai}]{hyperref}

\graphicspath{{figs/}}

\begin{document}
"""

FRONT_HEAD = r"""\frontmatter
\tableofcontents
\clearpage
"""

MAIN_HEAD = r"""\mainmatter
"""


def parse_summary() -> list[tuple[str, list[tuple[str, str]]]]:
    """SUMMARY.md → [(部分标题, [(章标题, 文件相对路径), ...]), ...]。"""
    parts: list[tuple[str, list[tuple[str, str]]]] = []
    current = ("", [])
    for line in SUMMARY.read_text(encoding="utf-8").splitlines():
        part = PART_LINE.match(line)
        if part:
            if current[0] or current[1]:
                parts.append(current)
            current = (part.group(1).strip(), [])
            continue
        entry = ENTRY_LINE.match(line)
        if entry and not line.startswith("#"):
            title, path = entry.group(1).strip(), entry.group(2).strip()
            title = re.sub(r"（选学）$", "", title).strip() or entry.group(1).strip()
            current[1].append((title, path))
    if current[0] or current[1]:
        parts.append(current)
    return parts


def build_tex(skip_first_h1: bool = True) -> list[str]:
    converter = MarkdownToLatex()
    out = [PREAMBLE, FRONT_HEAD]
    first_part = True
    for part_title, chapters in parse_summary():
        out.append("\\part*{" + converter.inline.convert(part_title) + "}")
        out.append(
            "\\addcontentsline{toc}{part}{%s}"
            % converter.inline.convert(part_title)
        )
        if first_part and part_title.startswith("前言"):
            # 前言部分留在 frontmatter（罗马页码），正文从第 0 部分起。
            for title, path in chapters:
                out.append(_chapter(converter, title, path))
            out.append(MAIN_HEAD)
            first_part = False
            continue
        first_part = False
        for title, path in chapters:
            out.append(_chapter(converter, title, path))
    out.append("\\end{document}")
    return out


def _chapter(converter: MarkdownToLatex, title: str, path: str) -> str:
    source = (BOOK / path).read_text(encoding="utf-8")
    body = converter.convert(source)
    title_tex = converter.inline.convert(title)
    head = (
        "\\chapter*{%s}\n"
        "\\addcontentsline{toc}{chapter}{%s}\n"
        "\\markboth{%s}{}\n" % (title_tex, title_tex, title_tex)
    )
    return head + body + "\n"


# ---------------------------------------------------------------- 各阶段

def export_figures() -> None:
    """调用 make_figures 重出配图；save() 被替换为输出矢量 PDF。"""
    scripts = str(ROOT / "scripts")
    if scripts not in sys.path:
        sys.path.insert(0, scripts)
    import make_figures as mf  # noqa: PLC0415
    import matplotlib.pyplot as plt  # noqa: PLC0415

    FIGS.mkdir(parents=True, exist_ok=True)
    mf.plt.rcParams["pdf.fonttype"] = 42  # TrueType 内嵌，文本可复制

    def save_pdf(fig, name):
        fig.savefig(FIGS / (pathlib.Path(name).stem + ".pdf"))
        plt.close(fig)
        print("fig-pdf:", name)

    mf.save = save_pdf
    sys.argv = ["make_figures.py"]
    mf.main()


def find_tectonic() -> str:
    found = shutil.which("tectonic")
    if found:
        return found
    env = os.environ.get("TECTONIC_BIN")
    if env and pathlib.Path(env).exists():
        return env
    fallback = pathlib.Path.home() / ".zcode-tools" / "tectonic" / "tectonic.exe"
    if fallback.exists():
        return str(fallback)
    raise SystemExit("tectonic 未找到：请安装或设置 TECTONIC_BIN")


def compile_tex() -> None:
    subprocess.run([find_tectonic(), "book_body.tex"], cwd=DIST, check=True)


def merge_cover(cover_pdf: str, out_name: str) -> None:
    from pypdf import PdfReader, PdfWriter, Transformation  # noqa: PLC0415

    a4_w, a4_h = 595.28, 841.89

    def normalize(page):
        box = page.mediabox
        w, h = float(box.width), float(box.height)
        if abs(w - a4_w) > 2 or abs(h - a4_h) > 2:
            page.add_transformation(Transformation().scale(sx=a4_w / w, sy=a4_h / h))
            page.mediabox.lower_left = (0, 0)
            page.mediabox.upper_right = (a4_w, a4_h)
        return page

    writer = PdfWriter()
    writer.add_page(normalize(PdfReader(cover_pdf).pages[0]))
    # 整体 append（而非逐页 add_page）：连同大纲书签与内链一起导入。
    writer.append(PdfReader(DIST / "book_body.pdf"))
    writer.add_metadata(
        {
            "/Title": "从一块电池到求解器",
            "/Author": "zhoukai",
            "/Subject": "以 P2D/DFN 电池模型为例的数值计算教程与 C++ 求解器",
            "/Creator": "Z.ai",
        }
    )
    out_path = DIST / out_name
    with open(out_path, "wb") as handle:
        writer.write(handle)
    print("wrote", out_path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-figures", action="store_true",
                        help="跳过配图重出（沿用 dist/figs 已有 PDF）")
    parser.add_argument("--skip-compile", action="store_true",
                        help="只生成 .tex 不编译")
    parser.add_argument("--merge-cover", metavar="PDF",
                        help="编译后把该封面 PDF 并为第 1 页")
    parser.add_argument("--out", default=FINAL_DEFAULT,
                        help="合并封面后的最终文件名（dist/ 下）")
    args = parser.parse_args()

    DIST.mkdir(parents=True, exist_ok=True)
    if not args.skip_figures:
        export_figures()

    tex = "\n".join(build_tex())
    (DIST / "book_body.tex").write_text(tex, encoding="utf-8")
    print("wrote", DIST / "book_body.tex")

    if args.skip_compile:
        return
    compile_tex()
    if args.merge_cover:
        merge_cover(args.merge_cover, args.out)


if __name__ == "__main__":
    main()
