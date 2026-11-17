# -*- coding: utf-8 -*-
"""Build report/report.docx from report/report.md (no pandoc on this machine)."""
import os
import re
import sys

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(HERE, "report", "report.md")
OUT = os.path.join(HERE, "report", "report.docx")

BODY_FONT = "Times New Roman"
HEAD_FONT = "Calibri"
GREY = RGBColor(0x59, 0x59, 0x59)


def add_page_number(paragraph):
    run = paragraph.add_run()
    fld = OxmlElement("w:fldSimple")
    fld.set(qn("w:instr"), "PAGE")
    run._r.append(fld)


def rule(paragraph):
    pPr = paragraph._p.get_or_add_pPr()
    borders = OxmlElement("w:pBdr")
    bottom = OxmlElement("w:bottom")
    bottom.set(qn("w:val"), "single")
    bottom.set(qn("w:sz"), "6")
    bottom.set(qn("w:space"), "1")
    bottom.set(qn("w:color"), "BFBFBF")
    borders.append(bottom)
    pPr.append(borders)


def styled(par, size=11, italic=False, color=None, bold=False):
    for run in par.runs:
        run.font.name = BODY_FONT
        run.font.size = Pt(size)
        run.font.italic = italic
        run.font.bold = bold
        if color is not None:
            run.font.color.rgb = color


def add_rich(par, text, size=11, color=None):
    """Split **bold** / *italic* segments and render X_name as a real subscript."""
    for token in re.split(r"(\*\*[^*]+\*\*|\*[^*]+\*)", text):
        if not token:
            continue
        bold = token.startswith("**") and token.endswith("**")
        italic = (not bold) and token.startswith("*") and token.endswith("*") and len(token) > 2
        seg = token[2:-2] if bold else (token[1:-1] if italic else token)
        for k, piece in enumerate(re.split(r"(?<=[A-Za-z\u0391-\u03c9\u0394\u03c30-9])_([A-Za-z0-9]+)", seg)):
            if not piece:
                continue
            run = par.add_run(piece)
            run.font.name = BODY_FONT
            run.font.size = Pt(size)
            run.font.bold = bold
            run.font.italic = italic
            run.font.subscript = (k % 2 == 1)
            if color is not None:
                run.font.color.rgb = color


def count_words(text):
    body = text.split("## 1. Context", 1)[1].split("## References", 1)[0]
    body = re.sub(r"!\[\]\([^)]*\)", " ", body)
    body = re.sub(r"\*\*Figure \d+\.\*\*[^\n]*", " ", body)
    body = "\n".join(l for l in body.split("\n") if not re.match(r"^n_(band|F)\(E_F\) =", l.strip()))
    return len(re.findall(r"[A-Za-z][A-Za-z'\-]*", body))


def _add_word_count(doc, text):
    note = doc.add_paragraph()
    note.paragraph_format.space_after = Pt(6)
    run = note.add_run("Word count: %d words in the main text, excluding the title block, "
                       "displayed equations, figure captions and references." % count_words(text))
    run.font.name = BODY_FONT
    run.font.size = Pt(9)
    run.font.color.rgb = GREY


def main():
    text = open(MD, encoding="utf-8").read()
    lines = text.split("\n")

    doc = Document()
    sec = doc.sections[0]
    sec.page_height, sec.page_width = Cm(29.7), Cm(21.0)
    for attr, val in (("top_margin", 2.3), ("bottom_margin", 2.3),
                      ("left_margin", 2.4), ("right_margin", 2.4)):
        setattr(sec, attr, Cm(val))

    style = doc.styles["Normal"]
    style.font.name = BODY_FONT
    style.font.size = Pt(11)
    style.paragraph_format.space_after = Pt(6)
    style.paragraph_format.line_spacing = 1.3

    footer = sec.footer.paragraphs[0]
    footer.alignment = WD_ALIGN_PARAGRAPH.CENTER
    add_page_number(footer)
    styled(footer, size=9, color=GREY)

    pending_caption = None
    note_done = False
    i = 0
    while i < len(lines):
        line = lines[i].rstrip()
        i += 1
        if not line.strip():
            continue
        if line.startswith("# "):
            par = doc.add_paragraph()
            par.paragraph_format.space_after = Pt(4)
            run = par.add_run(line[2:])
            run.font.name = HEAD_FONT
            run.font.size = Pt(18)
            run.font.bold = True
            continue
        if line.startswith("**QXU6028") or line.startswith("Group G18"):
            par = doc.add_paragraph()
            par.paragraph_format.space_after = Pt(2)
            run = par.add_run(line.replace("**", ""))
            run.font.name = BODY_FONT
            run.font.size = Pt(10.5)
            run.font.color.rgb = GREY
            continue
        if line.strip() == "---":
            if not note_done:
                _add_word_count(doc, text)
                note_done = True
            par = doc.add_paragraph()
            par.paragraph_format.space_after = Pt(10)
            rule(par)
            continue
        if line.startswith("## "):
            par = doc.add_paragraph()
            par.paragraph_format.space_before = Pt(14)
            par.paragraph_format.space_after = Pt(5)
            par.paragraph_format.keep_with_next = True
            run = par.add_run(line[3:])
            run.font.name = HEAD_FONT
            run.font.size = Pt(13)
            run.font.bold = True
            continue
        m = re.match(r"!\[\]\((.+?)\)", line)
        if m:
            path = os.path.normpath(os.path.join(os.path.dirname(MD), m.group(1)))
            par = doc.add_paragraph()
            par.alignment = WD_ALIGN_PARAGRAPH.CENTER
            par.paragraph_format.space_before = Pt(8)
            par.paragraph_format.space_after = Pt(3)
            par.paragraph_format.keep_with_next = True
            par.add_run().add_picture(path, width=Cm(9.0))
            pending_caption = True
            continue
        if line.startswith("**Figure "):
            par = doc.add_paragraph()
            par.alignment = WD_ALIGN_PARAGRAPH.CENTER
            par.paragraph_format.space_after = Pt(10)
            add_rich(par, line, size=9, color=GREY)
            continue
        if re.match(r"^\d+\.\s", line):
            par = doc.add_paragraph()
            par.paragraph_format.left_indent = Cm(0.6)
            par.paragraph_format.first_line_indent = Cm(-0.6)
            par.paragraph_format.space_after = Pt(2)
            par.paragraph_format.line_spacing = 1.05
            add_rich(par, line, size=9.5)
            continue
        if re.match(r"^n_(band|F)\(E_F\) =", line):
            par = doc.add_paragraph()
            par.alignment = WD_ALIGN_PARAGRAPH.CENTER
            par.paragraph_format.space_before = Pt(6)
            par.paragraph_format.space_after = Pt(6)
            add_rich(par, line, size=11)
            for run in par.runs:
                run.font.italic = True
            continue
        par = doc.add_paragraph()
        par.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        add_rich(par, line)


    doc.save(OUT)
    words = 0
    body = text.split("## 1. Context", 1)[1].split("## References", 1)[0]
    body = re.sub(r"!\[\]\([^)]*\)", " ", body)
    body = re.sub(r"\*\*Figure \d+\.\*\*[^\n]*", " ", body)
    words = len(re.findall(r"[A-Za-z][A-Za-z'\-]*", body))
    print("saved", OUT)
    print("body words:", words)


if __name__ == "__main__":
    main()
