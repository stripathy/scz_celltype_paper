#!/usr/bin/env python3
"""shared/figure_inputs.py — write and verify figure-input MANIFEST.tsv files.

Python counterpart to ``shared/figure_inputs.R``. The R module is what figure
scripts use to *refuse to draw* when a snapshot has gone stale; this module is
what the Python exporters use to *record* the provenance that makes that check
possible, and what Python-side plotting scripts use to check it themselves.

The contract is one MANIFEST.tsv beside each snapshot set, with columns:
    file, source, source_md5, snapshot_md5, source_mtime, refreshed_at, description

Large binary inputs (h5ad objects, sumstats) are fingerprinted by size+mtime
rather than md5 — hashing 1.4 GB on every run is not worth it, and any real
rewrite changes both.

Typical use in an exporter, after it has written its panel CSVs::

    from figure_inputs import write_manifest
    write_manifest(OUT, sources={
        "panel_B_genetics_vs_depletion.csv": TABLES / "gwas_vs_casecontrol_composition.csv",
        "panel_E_hcn1_vs_sag.csv":           TABLES / "sst_supertype_ephys_summary.csv",
    })
"""
from __future__ import annotations

import hashlib
import os
from datetime import datetime
from pathlib import Path

COLUMNS = ["file", "source", "source_md5", "snapshot_md5",
           "source_mtime", "refreshed_at", "description"]

# above this, fingerprint by size+mtime instead of hashing the bytes
_HASH_LIMIT = 256 * 1024 * 1024


def fingerprint(path, like: str | None = None) -> str:
    """md5 for ordinary files; a cheap size+mtime tag for very large ones.

    ``like`` pins the *kind* of fingerprint to whatever a manifest already
    recorded, so verification never depends on reader and writer agreeing about
    the size threshold — a disagreement there would make every check fail with a
    spurious "source has changed".
    """
    p = Path(os.path.expanduser(str(path)))
    if not p.exists():
        return ""
    st = p.stat()
    cheap = like.startswith("size-mtime:") if like else st.st_size > _HASH_LIMIT
    if cheap:
        return f"size-mtime:{st.st_size}:{int(st.st_mtime)}"
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(out_dir, sources: dict, descriptions: dict | None = None,
                   cheap: bool = False) -> Path:
    """Record, for each snapshot file, the source it was derived from.

    ``sources`` maps snapshot filename -> canonical source path. A snapshot may
    legitimately have several sources; pass a list and they are recorded as
    separate rows sharing the snapshot name.

    ``cheap`` forces size+mtime fingerprints regardless of file size. Use it for
    bulk inputs — the 24 per-sample Xenium h5ads are ~50 MB each, under the
    hashing threshold, but md5-ing 1.2 GB on every build is pointless when any
    real rewrite of those objects moves size and mtime anyway.
    """
    out_dir = Path(out_dir)
    descriptions = descriptions or {}
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    like = "size-mtime:" if cheap else None
    rows = []
    for name, srcs in sources.items():
        if not isinstance(srcs, (list, tuple)):
            srcs = [srcs]
        snap = out_dir / name
        snap_fp = fingerprint(snap) if snap.exists() else ""
        for src in srcs:
            sp = Path(os.path.expanduser(str(src)))
            mt = (datetime.fromtimestamp(sp.stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S")
                  if sp.exists() else "")
            rows.append([name, str(src), fingerprint(sp, like=like), snap_fp, mt, now,
                         descriptions.get(name, "")])
    path = out_dir / "MANIFEST.tsv"
    with open(path, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    return path


def read_manifest(dir_) -> list[dict]:
    """Existing manifest rows, or [] when there is no manifest yet."""
    p = Path(dir_) / "MANIFEST.tsv"
    if not p.exists():
        return []
    lines = p.read_text().splitlines()
    if not lines:
        return []
    hdr = lines[0].split("\t")
    return [dict(zip(hdr, ln.split("\t"))) for ln in lines[1:] if ln.strip()]


def merge_manifest(out_dir, sources: dict, descriptions: dict | None = None) -> Path:
    """Update the entries for ``sources`` while preserving every other row.

    Needed where two exporters write into one snapshot directory (Fig. 4's
    r_panels/): a plain write_manifest from either one would silently drop the
    other's provenance.
    """
    out_dir = Path(out_dir)
    keep = [r for r in read_manifest(out_dir) if r.get("file") not in sources]
    write_manifest(out_dir, sources, descriptions)
    fresh = read_manifest(out_dir)
    path = out_dir / "MANIFEST.tsv"
    with open(path, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in fresh + keep:
            fh.write("\t".join(str(r.get(c, "")) for c in COLUMNS) + "\n")
    return path


def problems(dir_) -> list[str]:
    """Human-readable list of stale/missing entries; empty when the set is clean."""
    d = Path(dir_)
    man = d / "MANIFEST.tsv"
    if not man.exists():
        return []
    out, seen = [], set()
    lines = man.read_text().splitlines()
    if not lines:
        return []
    hdr = lines[0].split("\t")
    for line in lines[1:]:
        f = dict(zip(hdr, line.split("\t")))
        name, src, rec = f.get("file", ""), f.get("source", ""), f.get("source_md5", "")
        if name and not (d / name).exists() and name not in seen:
            out.append(f"  {name} — snapshot file is missing")
            seen.add(name)
            continue
        sp = Path(os.path.expanduser(src))
        if not rec or not sp.exists():
            continue
        if fingerprint(sp, like=rec) != rec:
            out.append(f"  {name} — source has changed since {f.get('refreshed_at','?')}\n"
                       f"      source: {src}")
    return out


def check(dir_, strict: bool = True, refresh_cmd: str | None = None) -> bool:
    """Verify a snapshot directory. Raises when strict and something is stale."""
    probs = problems(dir_)
    if not probs:
        return True
    msg = (f"STALE FIGURE INPUTS in {dir_}\n" + "\n".join(probs) +
           "\n  The figure would render numbers that no longer match the analysis."
           f"\n  Fix: {refresh_cmd or 'run the component refresh script'}")
    if strict:
        raise SystemExit(msg)
    print(msg)
    return False
