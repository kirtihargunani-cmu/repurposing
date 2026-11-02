"""Shared look for every scene: colors, fonts, timing, helpers. Change the look here, not in scenes.

Colors carry meaning and must stay consistent across the whole deck:
  UP    = regulator more active in the tumor than in healthy tissue
  DOWN  = regulator less active in the tumor
  DRUG  = anything the drug does
  GOOD  = moves toward healthy / passes a check
  BAD   = moves away from healthy / fails a check / flagged
"""
import json
from pathlib import Path

from manim import (DOWN, LEFT, RIGHT, UP, Create, FadeIn, FadeOut, Rectangle, RoundedRectangle, Text,
                   VGroup, Write, config)

# Canvas: dark, 3Blue1Brown-like.
config.background_color = "#14161B"

BG = "#14161B"
INK = "#ECEFF4"        # main text
MUTED = "#8A93A6"      # secondary text, axes
UP_C = "#E8735A"       # tumor-up (warm)
DOWN_C = "#58A6E8"     # tumor-down (cool)
DRUG_C = "#F2C14E"     # the drug (yellow, the "hero" color)
GOOD_C = "#7BC47F"
BAD_C = "#C9606B"
NEUTRAL = "#4A5163"

FONT = "DejaVu Sans"   # present on most Linux machines; swap for "Inter" if installed
TITLE_SIZE = 44
BODY_SIZE = 28
SMALL_SIZE = 20

# Timing (seconds). Slow is better than fast: one idea per beat.
T_FAST, T_MED, T_SLOW = 0.6, 1.2, 2.0

DATA = json.loads((Path(__file__).resolve().parent / "scene_data.json").read_text())


def title(text, size=TITLE_SIZE):
    return Text(text, font=FONT, font_size=size, color=INK)


def body(text, size=BODY_SIZE, color=INK):
    return Text(text, font=FONT, font_size=size, color=color)


def caption(text):
    """One-line takeaway pinned to the bottom of the frame."""
    return Text(text, font=FONT, font_size=SMALL_SIZE + 4, color=MUTED).to_edge(DOWN, buff=0.4)


def step_badge(n, label):
    """Small 'Step n · label' tag in the top-left corner, shown on every step scene."""
    t = Text(f"Step {n} · {label}", font=FONT, font_size=SMALL_SIZE, color=MUTED)
    box = RoundedRectangle(corner_radius=0.12, width=t.width + 0.4, height=t.height + 0.25,
                           stroke_color=NEUTRAL, stroke_width=1.5)
    return VGroup(box, t).to_corner(UP + LEFT, buff=0.35)


def signed_bar(value, scale=0.35, width=0.32, pos_color=UP_C, neg_color=DOWN_C):
    """A vertical bar standing on y=0: up for positive, down for negative."""
    h = max(abs(value) * scale, 0.02)
    bar = Rectangle(width=width, height=h, stroke_width=0,
                    fill_color=pos_color if value >= 0 else neg_color, fill_opacity=0.9)
    bar.shift((h / 2 if value >= 0 else -h / 2) * UP)
    return bar
