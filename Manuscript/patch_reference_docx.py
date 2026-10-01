"""Apply the house typography to the Word reference documents used by Quarto.

    python3 patch_reference_docx.py custom-reference.docx oikos-reference.docx

- Font: Palatino Linotype everywhere. The theme's major and minor Latin fonts
  are replaced, so body text, headings and captions (which refer to the theme)
  all change, and the document defaults name the font explicitly as well.
- Tables: 9 pt, set on the "Table" table style that pandoc applies to every
  table, so only text inside tables changes (tight lists use "Compact" outside
  tables and keep the body size).

Idempotent: running it again leaves an already patched file unchanged.
Line numbering and double spacing of oikos-reference.docx were added when that
file was derived from custom-reference.docx (see the commit history).
"""
import re
import shutil
import sys
import tempfile
import zipfile

FONT = "Palatino Linotype"
TABLE_HALF_POINTS = 18  # 9 pt


def patch_theme(xml):
    xml = re.sub(r'(<a:majorFont>\s*<a:latin typeface=")[^"]*(")', rf"\g<1>{FONT}\2", xml)
    return re.sub(r'(<a:minorFont>\s*<a:latin typeface=")[^"]*(")', rf"\g<1>{FONT}\2", xml)


def patch_styles(xml):
    # document defaults: explicit font next to the theme references
    def defaults(m):
        rfonts = m.group(0)
        for attr in ("w:ascii", "w:hAnsi", "w:cs"):
            if f'{attr}="' in rfonts:
                rfonts = re.sub(rf'{attr}="[^"]*"', f'{attr}="{FONT}"', rfonts)
            else:
                rfonts = rfonts.replace("<w:rFonts ", f'<w:rFonts {attr}="{FONT}" ', 1)
        return rfonts
    xml = re.sub(r"<w:rPrDefault>\s*<w:rPr>\s*<w:rFonts [^>]*/>", defaults, xml, count=1)
    # table style: 9 pt text
    m = re.search(r'<w:style w:type="table" w:default="1" w:styleId="Table">.*?</w:style>', xml, re.S)
    if m is None:
        sys.exit("no 'Table' table style found")
    style = m.group(0)
    size = f'<w:rPr><w:sz w:val="{TABLE_HALF_POINTS}" /><w:szCs w:val="{TABLE_HALF_POINTS}" /></w:rPr>'
    style = re.sub(r"\s*<w:rPr>\s*<w:sz [^>]*/>\s*<w:szCs [^>]*/>\s*</w:rPr>", "", style)
    style = style.replace("<w:tblPr>", size + "\n    <w:tblPr>", 1)
    return xml[:m.start()] + style + xml[m.end():]


for path in sys.argv[1:]:
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".docx").name
    with zipfile.ZipFile(path) as zin, zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename == "word/theme/theme1.xml":
                data = patch_theme(data.decode("utf-8")).encode("utf-8")
            elif item.filename == "word/styles.xml":
                data = patch_styles(data.decode("utf-8")).encode("utf-8")
            zout.writestr(item, data)
    shutil.move(tmp, path)
    print("patched", path)
