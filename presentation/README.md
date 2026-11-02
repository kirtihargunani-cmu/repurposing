# presentation/: animated pipeline deck (3Blue1Brown style)

| File | What it is |
|---|---|
| `PROMPT.md` | The prompt to give Claude (Code or Cowork) to build the whole deck |
| `STORYBOARD.md` | 16 scenes: question, visual metaphor, beats, data keys, takeaway; plus a per-scene template |
| `style.py` | Colors (with fixed meanings), fonts, timing, helpers; every scene imports it |
| `scenes_template.py` | Two finished, test-rendered scenes (`Opening`, `ReversalScore`) showing the pattern |
| `export_scene_data.py` -> `scene_data.json` | Every on-screen number, exported from the results files |

`scene_data.json` is committed, so scenes render without the results; `export_scene_data.py` refreshes it from them.

Quick test (after installing Manim, see PROMPT.md):

```
cd presentation
python3 export_scene_data.py
manim -ql scenes_template.py ReversalScore      # low-quality preview in media/videos/
```

Tested 2026-11-02 without LaTeX: both template scenes render with Manim 0.19.0 + manim-slides 5.7.0
(`manim-slides render`). Manim 0.22 renders plain video but breaks manim-slides. HTML conversion (`manim-slides convert`) not yet checked.
