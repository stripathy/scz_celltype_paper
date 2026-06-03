#!/usr/bin/env python3
"""
Convert a markdown file with local image references to PDF using PyMuPDF (fitz).

Uses fitz.Story for HTML rendering with embedded images.

Usage:
    python markdown_to_pdf.py <input.md> <output.pdf>
"""

import sys
import os
import re
import base64
import markdown
import fitz  # PyMuPDF


def embed_images_in_html(html_content, base_dir):
    """Replace local image src paths with base64-embedded data URIs."""
    def replace_img(match):
        full_tag = match.group(0)
        src = match.group(1)

        # Skip if already a data URI or URL
        if src.startswith('data:') or src.startswith('http'):
            return full_tag

        # Resolve path relative to base_dir
        img_path = os.path.join(base_dir, src)
        if not os.path.exists(img_path):
            print(f"  Warning: Image not found: {img_path}")
            return full_tag

        # Read and base64 encode
        with open(img_path, 'rb') as f:
            img_data = f.read()

        ext = os.path.splitext(src)[1].lower()
        mime = {'.png': 'image/png', '.jpg': 'image/jpeg', '.jpeg': 'image/jpeg',
                '.gif': 'image/gif', '.svg': 'image/svg+xml'}.get(ext, 'image/png')

        b64 = base64.b64encode(img_data).decode('utf-8')
        data_uri = f"data:{mime};base64,{b64}"

        size_kb = len(img_data) / 1024
        print(f"  Embedded: {src} ({size_kb:.0f} KB)")

        return full_tag.replace(src, data_uri)

    return re.sub(r'<img[^>]*src="([^"]*)"[^>]*/?>',  replace_img, html_content)


def md_to_pdf(md_path, pdf_path):
    """Convert markdown file to PDF."""
    print(f"Converting: {md_path}")
    print(f"Output:     {pdf_path}")

    # Read markdown
    with open(md_path, 'r') as f:
        md_text = f.read()

    # Convert markdown to HTML
    extensions = ['tables', 'fenced_code', 'codehilite', 'toc', 'nl2br']
    html_body = markdown.markdown(md_text, extensions=extensions)

    # Embed local images
    base_dir = os.path.dirname(os.path.abspath(md_path))
    print(f"\nEmbedding images from: {base_dir}")
    html_body = embed_images_in_html(html_body, base_dir)

    # Build full HTML document with CSS styling
    css = """
    body {
        font-family: Helvetica, Arial, sans-serif;
        font-size: 10pt;
        line-height: 1.5;
        color: #1a1a1a;
        max-width: 100%;
        margin: 0;
        padding: 0;
    }
    h1 {
        font-size: 20pt;
        color: #1a365d;
        border-bottom: 2px solid #2b6cb0;
        padding-bottom: 6pt;
        margin-top: 20pt;
        margin-bottom: 10pt;
    }
    h2 {
        font-size: 16pt;
        color: #2b6cb0;
        border-bottom: 1px solid #bee3f8;
        padding-bottom: 4pt;
        margin-top: 18pt;
        margin-bottom: 8pt;
    }
    h3 {
        font-size: 13pt;
        color: #2c5282;
        margin-top: 14pt;
        margin-bottom: 6pt;
    }
    h4 {
        font-size: 11pt;
        color: #4a5568;
        margin-top: 10pt;
        margin-bottom: 4pt;
    }
    p {
        margin-top: 4pt;
        margin-bottom: 6pt;
        text-align: justify;
    }
    table {
        border-collapse: collapse;
        margin: 8pt 0;
        font-size: 9pt;
        width: 100%;
    }
    th {
        background-color: #edf2f7;
        border: 1px solid #cbd5e0;
        padding: 4pt 6pt;
        text-align: left;
        font-weight: bold;
    }
    td {
        border: 1px solid #e2e8f0;
        padding: 3pt 6pt;
    }
    tr:nth-child(even) td {
        background-color: #f7fafc;
    }
    img {
        max-width: 100%;
        height: auto;
        margin: 8pt 0;
        display: block;
    }
    code {
        background-color: #edf2f7;
        padding: 1pt 3pt;
        border-radius: 2pt;
        font-size: 9pt;
        font-family: "Courier New", monospace;
    }
    pre {
        background-color: #edf2f7;
        padding: 8pt;
        border-radius: 4pt;
        font-size: 8.5pt;
        overflow-x: auto;
    }
    strong {
        color: #1a365d;
    }
    em {
        color: #4a5568;
    }
    hr {
        border: none;
        border-top: 1px solid #e2e8f0;
        margin: 14pt 0;
    }
    ul, ol {
        margin-top: 4pt;
        margin-bottom: 6pt;
        padding-left: 20pt;
    }
    li {
        margin-bottom: 2pt;
    }
    blockquote {
        border-left: 3pt solid #bee3f8;
        margin-left: 0;
        padding-left: 10pt;
        color: #4a5568;
    }
    """

    full_html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<style>{css}</style>
</head>
<body>
{html_body}
</body>
</html>"""

    # Use fitz.Story to render HTML to PDF
    print(f"\nRendering PDF...")

    # Letter size: 8.5 x 11 inches = 612 x 792 points
    page_width = 612
    page_height = 792
    margin = 54  # 0.75 inch margins

    content_rect = fitz.Rect(margin, margin, page_width - margin, page_height - margin)

    story = fitz.Story(html=full_html)
    writer = fitz.DocumentWriter(pdf_path)

    more = True
    page_num = 0
    while more:
        page_num += 1
        dev = writer.begin_page(fitz.Rect(0, 0, page_width, page_height))
        more, _ = story.place(content_rect)
        story.draw(dev)
        writer.end_page()

    writer.close()
    print(f"\nDone! {page_num} pages written to: {pdf_path}")

    # Report file size
    size_mb = os.path.getsize(pdf_path) / (1024 * 1024)
    print(f"File size: {size_mb:.1f} MB")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python markdown_to_pdf.py <input.md> <output.pdf>")
        sys.exit(1)

    md_to_pdf(sys.argv[1], sys.argv[2])
