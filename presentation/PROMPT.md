# Prompt: build the animated pipeline deck

Paste everything below the line into Claude Code (recommended, run from `~/projects/repurposing`) or Claude
Cowork (share the `presentation/` folder with it first). Fill in the three `<...>` fields at the top.

---

## Who it's for and how it will be shown

- **Audience:** <e.g. my PI and lab meeting: biologists, comfortable with RNA-seq, not with VIPER or Celligner>
- **Length:** <e.g. 12-15 minutes, presented live, I click to advance>
- **Output:** <e.g. a clickable HTML deck plus one MP4 of the whole thing>

## What I want

An animated slide deck in the style of 3Blue1Brown. Each step of my drug-repurposing pipeline is built up visually,
one idea at a time, and each step's objects carry over into the next, so the audience watches the pipeline grow
instead of reading slides. It should explain *why* each step exists as much as *what* it does, and show real results
honestly, including where the method fails.

## Tools

- **Manim Community Edition** for the animation, **manim-slides** to turn the scenes into a deck I click through. Use
  `manim-slides convert ... deck.html` for a self-contained HTML deck. Also render one combined MP4.
- Install with conda-forge (pip often fails on cairo/pango): `conda create -p ./manim-env -c conda-forge python=3.11 manim`,
  then `pip install "manim==0.19.0" manim-slides`. **Pin Manim 0.19:** manim-slides 5.7 crashes on Manim 0.22
  (`ManimConfig has no attribute 'write_to_movie'`). Rendering with `manim-slides render` was tested with 0.19 +
  5.7; `manim-slides convert` to HTML was not checked yet. If it hangs or fails, render the MP4s and build the deck
  another way, then tell me.
- **No LaTeX is installed.** Use `Text` only. Avoid `MathTex`, `Tex`, `DecimalNumber` and `Axes` number labels: they
  need LaTeX. For changing numbers, use a `ValueTracker` + `always_redraw(Text(...))`, as in the template. If you
  want real formulas, ask me first, and I'll install a LaTeX distribution.

## Files to read first (in `presentation/`)

1. `STORYBOARD.md`: the 16 scenes, their questions, visual metaphors, beats and data keys. This is the plan. You may
   improve it, but show me the changes before animating.
2. `style.py`: colors, fonts, timing and helpers. Every scene imports it. Colors carry meaning (warm = up in tumor,
   cool = down, yellow = the drug, green = toward healthy / passes, red = flagged / fails). Never use them for anything
   else.
3. `scenes_template.py`: two finished scenes (`Opening`, `ReversalScore`) that show the pattern: data from
   `style.DATA`, one idea per beat, a `self.next_slide()` pause after each beat, a caption with the takeaway.
4. `scene_data.json`: every number to show, exported from the results by `export_scene_data.py`. **Do not type
   numbers into scenes.** If a scene needs a number that is not in the JSON, add it to `export_scene_data.py` from the
   results file and re-run the export. If you are in Cowork without the repo, tell me which number is missing and
   where you would expect it.
5. Background, if available: `DRUG_SIDE_PLAIN_LANGUAGE.md` (plain-language steps + results), `PIPELINE_DRUG_SIDE_V2.md`
   (technical), `DISEASE_SIDE_PLAIN_LANGUAGE.md` and `PIPELINE_V2.md` (disease side, Steps 1-10), all in `docs/`.

## How to work (in this order)

1. **Storyboard pass, no code.** Read the files. For each scene write the per-scene template at the bottom of
   STORYBOARD.md: opening frame, beats with one narration sentence each, closing frame, data keys, caption and honesty
   note. Decide which objects carry over from scene to scene (e.g. the regulator circles of scene 4 are the same objects
   the drug arrows hit in scene 8). Stop and show me the storyboard.
2. **One scene at a time.** After I approve, build scenes in order in `scenes.py` (keep `scenes_template.py` as is).
   For each scene:
   - render at low quality (`manim -ql`),
   - extract 3-4 frames with ffmpeg (`ffmpeg -ss <t> -i <mp4> -frames:v 1 frame.png`), one at the end of each beat,
   - **look at the frames** and fix overlaps, clipped text, crowding, unreadable labels and wrong colors before moving
     on. The template's first draft had a title colliding with the score panel and two captions on top of each other.
     Those are the kind of bugs to catch here.
3. **Self-review each scene** against this checklist before showing it:
   - Is there one idea per beat, and can a viewer say what changed after each beat?
   - Is every number on screen from `scene_data.json`?
   - Do colors keep their meaning? Is text >= 20 px at 1080p, nothing within 0.3 units of the frame edge?
   - Does the scene start from where the last one ended (continuity), not from a blank screen?
   - Is it honest? Failures and caveats from the storyboard are shown, not skipped.
4. **Final render** at `-qh` (1080p60), build the manim-slides deck and the combined MP4, then give me the commands to
   present it.

## Style (3Blue1Brown, concretely)

- Dark background, few colors, lots of empty space. One focal object at a time; everything else dims (opacity
  0.3) rather than disappearing.
- Things **transform** into each other instead of cutting: a patient dot becomes a column of gene bars, a gene bar
  becomes part of a regulator's target set. Use `ReplacementTransform`, `TransformMatchingShapes`, `LaggedStart`,
  `GrowFromEdge` and `Indicate`. Avoid `FadeOut` + `FadeIn` when a transform tells the story.
- Show the concrete case before the general rule: one gene, then all genes; one regulator, then all; one drug and one
  cancer (LUAD + erlotinib), then the whole ranking.
- Motion has meaning: bars always grow from the baseline; "toward healthy" is always motion toward the baseline.
- Text is sparse: a short headline, a caption with the takeaway, and labels only where needed. The narration carries
  the rest, so write a speaker-notes file (`NOTES.md`, one block per beat).
- Pace: about 1-2 s per motion, a pause at every beat end, about 45-90 s per scene.

## Content rules

- Running example: lung adenocarcinoma + erlotinib (works). Honest failure: anti-estrogens in luminal A breast cancer
  (the tumor-learned ESR1 map misreads them, C18 / C25). Flag example: anti-estrogens in basal breast, caught by
  `unmatched_lines` (C20).
- Lead with mechanism specificity (`spec_mech`) when showing rankings, and say why (32 drugs top >= 20 cancers by raw
  reversal; none do by spec_mech).
- Always show the evidence status: 14 of 41 cancers are matched; 27 rest on borrowed data; the stricter cell-type check
  passes 14 of 41.
- End by saying these are ranked hypotheses: the PRISM / GDSC validation has not been run on v2.
- No new claims about biology beyond what the docs say. If you want to add one, ask me.

## Deliverables

- `presentation/scenes.py` (all scenes), `presentation/NOTES.md` (speaker notes per beat), the updated
  `STORYBOARD.md`
- `presentation/deck.html` (manim-slides, self-contained) and `presentation/pipeline.mp4`
- A short list of anything you could not show honestly or that needs my decision
