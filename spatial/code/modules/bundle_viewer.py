#!/usr/bin/env python3
"""
Bundle the Xenium SCZ Spatial Viewer into a single self-contained HTML file.
All JSON data is gzip-compressed, base64-encoded, and embedded as JS variables.
The browser decompresses on-the-fly using DecompressionStream API.
"""

import json
import gzip
import base64
import os
import sys

VIEWER_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'output', 'viewer')
OUTPUT_FILE = os.path.join(VIEWER_DIR, 'xenium_viewer_standalone.html')

def compress_json(filepath):
    """Read JSON file, gzip compress, base64 encode."""
    with open(filepath, 'rb') as f:
        raw = f.read()
    compressed = gzip.compress(raw, compresslevel=9)
    encoded = base64.b64encode(compressed).decode('ascii')
    return len(raw), len(compressed), encoded

def main():
    # Read the current index.html to get the full HTML/CSS/JS
    html_path = os.path.join(VIEWER_DIR, 'index.html')
    with open(html_path, 'r') as f:
        html = f.read()

    # Read and compress index.json
    index_path = os.path.join(VIEWER_DIR, 'index.json')
    with open(index_path, 'r') as f:
        index_data = json.load(f)

    raw_size, comp_size, index_b64 = compress_json(index_path)
    print(f"index.json: {raw_size/1024:.1f}KB -> {comp_size/1024:.1f}KB")

    # Compress each sample JSON
    sample_blobs = {}
    total_raw = raw_size
    total_comp = comp_size

    for sample_info in index_data['samples']:
        sid = sample_info['sample_id']
        sample_path = os.path.join(VIEWER_DIR, f'{sid}.json')
        if not os.path.exists(sample_path):
            print(f"  WARNING: {sample_path} not found, skipping")
            continue
        raw_s, comp_s, blob = compress_json(sample_path)
        sample_blobs[sid] = blob
        total_raw += raw_s
        total_comp += comp_s
        print(f"  {sid}: {raw_s/1024/1024:.1f}MB -> {comp_s/1024/1024:.1f}MB ({comp_s/raw_s*100:.0f}%)")

    print(f"\nTotal: {total_raw/1024/1024:.1f}MB raw -> {total_comp/1024/1024:.1f}MB compressed")
    # Base64 adds ~33% overhead
    b64_total = total_comp * 4 / 3
    print(f"Estimated HTML file size: ~{b64_total/1024/1024:.1f}MB")

    # Build the embedded data JS block
    data_js_lines = []
    data_js_lines.append("// ── Embedded compressed data ──")
    data_js_lines.append(f'const INDEX_B64 = "{index_b64}";')
    data_js_lines.append("const SAMPLE_BLOBS = {")
    for sid, blob in sample_blobs.items():
        data_js_lines.append(f'  "{sid}": "{blob}",')
    data_js_lines.append("};")

    # Decompression helper using DecompressionStream (modern browsers)
    # Falls back to manual inflate if not available
    decompress_js = '''
// ── Decompression utilities ──
async function decompressB64(b64str) {
  const binary = atob(b64str);
  const bytes = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);

  // Use DecompressionStream API (Chrome 80+, Firefox 113+, Safari 16.4+)
  if (typeof DecompressionStream !== 'undefined') {
    const ds = new DecompressionStream('gzip');
    const writer = ds.writable.getWriter();
    writer.write(bytes);
    writer.close();
    const reader = ds.readable.getReader();
    const chunks = [];
    while (true) {
      const { done, value } = await reader.read();
      if (done) break;
      chunks.push(value);
    }
    const totalLen = chunks.reduce((s, c) => s + c.length, 0);
    const result = new Uint8Array(totalLen);
    let offset = 0;
    for (const chunk of chunks) {
      result.set(chunk, offset);
      offset += chunk.length;
    }
    return JSON.parse(new TextDecoder().decode(result));
  } else {
    // Fallback: raw inflate (strip gzip header/trailer)
    // This is a minimal fallback - modern browsers all support DecompressionStream
    throw new Error("DecompressionStream not supported. Please use Chrome 80+, Firefox 113+, or Safari 16.4+.");
  }
}
'''

    # Now modify the HTML to use embedded data instead of fetch()
    # We need to replace:
    # 1. The fetch('index.json') call in init()
    # 2. The fetch(`${sampleId}.json`) call in loadSample()

    # Find the <script> tag and inject our data + decompression code right after it
    script_marker = '<script>'
    script_idx = html.index(script_marker) + len(script_marker)

    # Build the replacement init and loadSample functions
    # We need to replace the fetch-based versions with decompression-based ones

    # Replace the init function
    old_init = '''async function init() {
  const resp = await fetch('index.json');
  indexData = await resp.json();
  buildSampleList();
  setupEvents();
  resizeCanvas();
  // Default to Br8667 (has transcript data) or first sample
  const defaultSample = indexData.samples.find(s => s.sample_id === 'Br8667') || indexData.samples[0];
  loadSample(defaultSample.sample_id);
}'''

    new_init = '''async function init() {
  indexData = await decompressB64(INDEX_B64);
  buildSampleList();
  setupEvents();
  resizeCanvas();
  // Default to Br8667 (has transcript data) or first sample
  const defaultSample = indexData.samples.find(s => s.sample_id === 'Br8667') || indexData.samples[0];
  loadSample(defaultSample.sample_id);
}'''

    # Replace the fetch in loadSample
    old_fetch = '''  const resp = await fetch(`${sampleId}.json`);
  sampleData = await resp.json();'''

    new_fetch = '''  if (!SAMPLE_BLOBS[sampleId]) { loading.classList.remove('show'); return; }
  sampleData = await decompressB64(SAMPLE_BLOBS[sampleId]);'''

    if old_init not in html:
        print("ERROR: Could not find init() function to replace")
        sys.exit(1)
    if old_fetch not in html:
        print("ERROR: Could not find fetch() in loadSample to replace")
        sys.exit(1)

    # Make replacements
    modified_html = html.replace(old_init, new_init)
    modified_html = modified_html.replace(old_fetch, new_fetch)

    # Inject data and decompression code after <script>
    inject_block = '\n' + '\n'.join(data_js_lines) + '\n' + decompress_js + '\n'
    modified_html = modified_html[:script_idx] + inject_block + modified_html[script_idx:]

    # Write the standalone HTML
    with open(OUTPUT_FILE, 'w') as f:
        f.write(modified_html)

    file_size = os.path.getsize(OUTPUT_FILE)
    print(f"\nWrote: {OUTPUT_FILE}")
    print(f"File size: {file_size/1024/1024:.1f}MB")
    print("Done! Open this file directly in any modern browser.")

def bundle_standalone_html(viewer_dir, output_path):
    """
    Bundle viewer into a single standalone HTML file.

    Parameters
    ----------
    viewer_dir : str
        Directory containing index.html, index.json, and sample JSON files.
    output_path : str
        Path to write the standalone HTML file.
    """
    global VIEWER_DIR, OUTPUT_FILE
    VIEWER_DIR = viewer_dir
    OUTPUT_FILE = output_path
    main()


if __name__ == '__main__':
    main()
