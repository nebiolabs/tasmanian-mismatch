#!/usr/bin/env python3
"""Generate a flat coverage badge SVG from grcov markdown output.

Usage: gen-coverage-badge.py <coverage_summary.md> <output.svg>

Reads the total coverage percentage from the last '**X.XX%**' in the
grcov markdown summary (which appears in the Total row), then writes a
flat badge SVG to the output path.
"""
import re
import sys

summary_path, output_path = sys.argv[1], sys.argv[2]

text = open(summary_path).read()
# grcov markdown: "Total coverage: 96.10%"
m = re.search(r'Total coverage:\s*(\d+(?:\.\d+)?)%', text, re.IGNORECASE)
pct = float(m.group(1)) if m else 0.0
msg = f"{pct:.1f}%"

color = "green" if pct >= 80 else "goldenrod" if pct >= 60 else "firebrick"

label = "coverage"
# Approximate character widths for DejaVu Sans 11px
label_w = len(label) * 7 + 10
msg_w = len(msg) * 7 + 10
total_w = label_w + msg_w
label_cx = label_w // 2
msg_cx = label_w + msg_w // 2

svg = f"""\
<svg xmlns="http://www.w3.org/2000/svg" width="{total_w}" height="20">
  <linearGradient id="s" x2="0" y2="100%">
    <stop offset="0" stop-color="#bbb" stop-opacity=".1"/>
    <stop offset="1" stop-opacity=".1"/>
  </linearGradient>
  <rect rx="3" width="{total_w}" height="20" fill="#555"/>
  <rect rx="3" x="{label_w}" width="{msg_w}" height="20" fill="{color}"/>
  <path fill="{color}" d="M{label_w} 0h4v20h-4z"/>
  <rect rx="3" width="{total_w}" height="20" fill="url(#s)"/>
  <g fill="#fff" text-anchor="middle"
     font-family="DejaVu Sans,Verdana,Geneva,sans-serif" font-size="11">
    <text x="{label_cx}" y="15" fill="#010101" fill-opacity=".3">{label}</text>
    <text x="{label_cx}" y="14">{label}</text>
    <text x="{msg_cx}" y="15" fill="#010101" fill-opacity=".3">{msg}</text>
    <text x="{msg_cx}" y="14">{msg}</text>
  </g>
</svg>
"""

import os
os.makedirs(os.path.dirname(output_path), exist_ok=True)
open(output_path, "w").write(svg)
print(f"Coverage badge: {msg} → {output_path}")
