#!/usr/bin/env python
"""Assemble a saved Slides deck (project/deck.json + deck.css + slides/*.html) into one viewable deck.html.

Usage: python scripts/deck_to_html.py docs/artifacts/<deck-folder> [...]

The deck folders under docs/artifacts/ are copies of the decks' own content files, saved from the published
artifacts. This renders them in a plain browser: each 1920x1080 slide is scaled to the window, in deck order, with
its speaker notes (the slide's <aside>) printed underneath. Nothing is fetched except the Google Fonts the deck names.
"""
import html
import json
import re
import sys
from pathlib import Path

PAGE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title>
{fonts}
<style>
body{{margin:0;background:#2b2b2b;color:#eee;font:15px/1.5 system-ui,-apple-system,sans-serif}}
header{{padding:20px max(16px,3vw)}} header h1{{margin:0;font-size:1.3rem}} header p{{margin:4px 0 0;color:#bbb}}
.deck{{display:grid;gap:28px;padding:0 max(16px,3vw) 48px}}
.wrap{{position:relative;width:100%;aspect-ratio:16/9;overflow:hidden;background:#fff;box-shadow:0 2px 12px #0008}}
.wrap > section{{position:absolute;top:0;left:0;width:1920px;height:1080px;transform-origin:0 0;box-sizing:border-box}}
.wrap > section aside{{display:none}}
.notes{{background:#1f1f1f;border-left:3px solid #6aa;padding:10px 14px;white-space:pre-wrap;color:#ddd;font-size:14px}}
.num{{color:#9cc;font-weight:600;margin-bottom:6px}}
</style>
<style>
{css}
</style></head><body>
<header><h1>{title}</h1><p>{n} slides · speaker notes under each slide · source: <code>{src}</code></p></header>
<main class="deck">
{slides}
</main>
<script>
function fit(){{document.querySelectorAll('.wrap').forEach(w=>{{const s=w.firstElementChild;s.style.transform='scale('+(w.clientWidth/1920)+')';}});}}
addEventListener('resize',fit);fit();
</script>
</body></html>
"""


def build(folder: Path):
    proj = folder / "project"
    deck = json.loads((proj / "deck.json").read_text())
    css = (proj / "deck.css").read_text() if (proj / "deck.css").exists() else ""
    fonts = "\n".join(f'<link rel="stylesheet" href="{html.escape(f["href"])}">'
                      for f in deck.get("faces", {}).values() if f.get("href", "").startswith("https://fonts.googleapis.com/"))
    parts = []
    for i, sid in enumerate(deck["order"], 1):
        f = proj / "slides" / f"{sid}.html"
        if not f.exists():
            continue
        sec = f.read_text().strip()
        m = re.search(r"<aside>(.*?)</aside>", sec, re.S)
        notes = m.group(1).strip() if m else ""
        parts.append(f'<div><div class="num">Slide {i} · {html.escape(sid)}</div><div class="wrap">{sec}</div>'
                     + (f'<div class="notes">{notes}</div>' if notes else "") + "</div>")
    out = folder / "deck.html"
    out.write_text(PAGE.format(title=html.escape(deck.get("title", folder.name)), fonts=fonts, css=css,
                               n=len(parts), src=html.escape(str(folder)), slides="\n".join(parts)))
    print(f"{out} ({len(parts)} slides)")


if __name__ == "__main__":
    for a in sys.argv[1:]:
        build(Path(a))
