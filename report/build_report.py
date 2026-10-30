"""Build report.html and report.pdf from report.md.

Usage (from the project root):  .venv/bin/python report/build_report.py
Markdown -> HTML with mistune (installed with nbconvert); HTML -> PDF with headless Google Chrome.
"""
import shutil
import subprocess
import sys
from pathlib import Path

import mistune

HERE = Path(__file__).resolve().parent
CSS = """
@page { size: A4; margin: 18mm 17mm 18mm 17mm; }
body { font-family: "Charter", "Georgia", serif; font-size: 10.5pt; line-height: 1.45; color: #111; max-width: 180mm; margin: auto; }
h1 { font-size: 19pt; line-height: 1.2; margin-bottom: 4pt; }
h2 { font-size: 14pt; margin-top: 18pt; border-bottom: 1px solid #bbb; padding-bottom: 2pt; page-break-after: avoid; }
h3 { font-size: 11.5pt; margin-top: 12pt; page-break-after: avoid; }
p, li { text-align: justify; hyphens: auto; }
img { display: block; max-width: 100%; max-height: 105mm; margin: 8pt auto 2pt auto; page-break-inside: avoid; }
p:has(> img) + p > em:first-child, p > img + em { display: block; }
em { color: #222; }
table { border-collapse: collapse; margin: 6pt auto; font-size: 8.8pt; page-break-inside: avoid; }
th, td { border: 1px solid #999; padding: 2pt 5pt; text-align: left; vertical-align: top; }
th { background: #eee; }
code { font-size: 9pt; background: #f3f3f3; padding: 0 2px; }
hr { border: none; border-top: 1px solid #ccc; margin: 12pt 0; }
"""


def main():
    md = (HERE / "report.md").read_text(encoding="utf-8")
    body = mistune.create_markdown(plugins=["table", "strikethrough"])(md)
    html = f"<!doctype html><html><head><meta charset='utf-8'><title>Report</title><style>{CSS}</style></head><body>{body}</body></html>"
    out_html = HERE / "report.html"
    out_html.write_text(html, encoding="utf-8")
    chrome = shutil.which("google-chrome") or "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
    if not Path(chrome).exists():
        sys.exit(f"wrote {out_html}; Chrome not found, open the HTML and print it to PDF manually")
    subprocess.run([chrome, "--headless", "--disable-gpu", "--no-pdf-header-footer",
                    f"--print-to-pdf={HERE / 'report.pdf'}", out_html.as_uri()], check=True, capture_output=True)
    print("wrote", out_html, "and", HERE / "report.pdf")


if __name__ == "__main__":
    main()
