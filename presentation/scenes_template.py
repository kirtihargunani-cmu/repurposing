"""Scene template for the animated deck. Two worked scenes show the pattern; copy one per storyboard beat.

Render one scene (low quality, fast):   manim -ql presentation/scenes_template.py ReversalScore
Render final (1080p60):                 manim -qh presentation/scenes_template.py ReversalScore
As a clickable deck (manim-slides):     manim-slides render presentation/scenes_template.py Opening ReversalScore
                                        manim-slides present Opening ReversalScore
                                        manim-slides convert Opening ReversalScore deck.html   (self-contained)

Rules every scene follows:
- Text only (no MathTex / DecimalNumber) unless LaTeX is installed; formulas can be Text too.
- Numbers come from style.DATA (scene_data.json), never typed in.
- One idea per beat; end each beat with a pause (self.next_slide()) so the speaker controls pace.
- Colors keep their meaning (see style.py).
"""
from manim import (DOWN, LEFT, RIGHT, UP, Arrow, FadeIn, FadeOut, GrowFromEdge, Indicate, LaggedStart, Line,
                   ReplacementTransform, ValueTracker, VGroup, Write, always_redraw)

from style import (BAD_C, DATA, DOWN_C, DRUG_C, GOOD_C, INK, MUTED, T_FAST, T_MED, T_SLOW, UP_C, body,
                   caption, signed_bar, step_badge, title)

try:  # clickable slides if manim-slides is installed, plain video otherwise
    from manim_slides import Slide as Base
except ImportError:  # pragma: no cover
    from manim import Scene as Base

    Base.next_slide = lambda self, *a, **k: self.wait(0.8)


class Opening(Base):
    """Title card + the one-sentence idea of the whole pipeline."""

    def construct(self):
        t = title("Finding drugs that push a cancer back toward healthy")
        sub = body("41 cancers · 14 subgroups · 3,058 drugs", color=MUTED).next_to(t, DOWN, buff=0.4)
        self.play(Write(t), run_time=T_SLOW)
        self.play(FadeIn(sub, shift=UP * 0.2), run_time=T_MED)
        self.next_slide()
        self.play(FadeOut(VGroup(t, sub)), run_time=T_FAST)


class ReversalScore(Base):
    """Step 6 beat: how one score is built, regulator by regulator (LUAD + erlotinib).

    Beat 1  the tumor's regulators as bars (up = warm, down = cool)
    Beat 2  the drug's effect on the same regulators (yellow outlines)
    Beat 3  each regulator's contribution drops into a running total -> rev
    """

    N = 10  # regulators shown; the total counts all of them

    def construct(self):
        ex = DATA["worked_examples"]["TCGA-LUAD|erlotinib"]
        regs = ex["top_positive"][: self.N]
        self.add(step_badge(6, "Score each drug"))

        head = body("Lung adenocarcinoma: which regulators changed?").to_corner(UP + LEFT, buff=0.5).shift(DOWN * 0.8)
        # Chart on the left two thirds, score panel on the right third.
        axis = Line(LEFT * 6.4, RIGHT * 3.0, color=MUTED, stroke_width=2).shift(DOWN * 0.3)
        xs = [-6.0 + i * 0.95 for i in range(len(regs))]

        # Beat 1: disease bars grow from the axis.
        dis = VGroup(*[signed_bar(r["disease_nes"]).shift([x, -0.3, 0]) for r, x in zip(regs, xs)])
        names = VGroup(*[body(r["regulator"], size=14, color=MUTED).rotate(0.6).move_to([x, -2.6, 0])
                         for r, x in zip(regs, xs)])
        self.play(FadeIn(head), Write(axis), run_time=T_MED)
        self.play(LaggedStart(*[GrowFromEdge(b, DOWN if r["disease_nes"] >= 0 else UP)
                                for b, r in zip(dis, regs)], lag_ratio=0.12),
                  FadeIn(names), run_time=T_SLOW)
        cap = caption("Warm = higher in the tumor than in healthy lung · cool = lower")
        self.play(FadeIn(cap))
        self.next_slide()

        # Beat 2: the drug's effect, drawn as yellow arrows pointing the way the drug pushes.
        arrows = VGroup()
        for r, x in zip(regs, xs):
            d = 1 if r["drug_activity"] > 0 else -1
            arrows.add(Arrow([x + 0.3, -0.3, 0], [x + 0.3, -0.3 + d * min(abs(r["drug_activity"]) * 0.18, 1.8), 0],
                             color=DRUG_C, buff=0, stroke_width=5, max_tip_length_to_length_ratio=0.25))
        new_head = body("Erlotinib pushes each one the other way").move_to(head, aligned_edge=LEFT)
        self.play(ReplacementTransform(head, new_head), LaggedStart(*[GrowFromEdge(a, DOWN) for a in arrows],
                                                                     lag_ratio=0.1), run_time=T_SLOW)
        self.next_slide()

        # Beat 3: contributions fall into a running total.
        total_lbl = body("Reversal score", color=MUTED).move_to([4.8, 1.2, 0])
        # Text-based counter (DecimalNumber would need a LaTeX install).
        tracker = ValueTracker(0)
        total = always_redraw(lambda: body(f"{tracker.get_value():.1f}", size=56, color=GOOD_C)
                              .next_to(total_lbl, DOWN, buff=0.3))
        self.play(FadeIn(total_lbl), FadeIn(total))
        running = 0.0
        for b, a, r in zip(dis, arrows, regs):
            running += r["contribution"]
            self.play(Indicate(b, color=GOOD_C, scale_factor=1.15), a.animate.set_opacity(0.3),
                      tracker.animate.set_value(running), run_time=0.45)
        self.next_slide()
        all_lbl = VGroup(body(f"all {ex['n_regulators']:,} regulators", size=22),
                         body(f"{ex['n_positive']} move toward healthy", size=22, color=GOOD_C)
                         ).arrange(DOWN, buff=0.15).move_to([4.8, -1.4, 0])
        self.play(tracker.animate.set_value(ex["rev"]), FadeIn(all_lbl), run_time=T_SLOW)
        self.play(ReplacementTransform(cap, caption("No single gene carries the score: it is spread over many regulators")))
        self.next_slide()
