---
name: scientific-presentation-designer
description: Skill for designing professional, high-impact scientific and R&D presentation slide decks (.pptx) using python-pptx and clean design principles.
---

# Scientific & R&D Presentation Designer Skill

This skill provides comprehensive guidelines, reusable python-pptx templates, layout patterns, and quality rules for generating clean, editable, high-impact presentation slide decks (`.pptx`) for scientific interviews, R&D progress updates, grant proposals, and technical seminars.

---

## 🎨 Core Visual Design Philosophy

1. **Anti-AI Clutter Aesthetic**:
   - **NO** generic AI stock photos, floating icons, rounded colorful cards, rainbow gradients, or glassmorphism.
   - **DO** use clean widescreen (16:9) slides (`13.333" x 7.5"`), generous padding, high whitespace contrast, and structured container cards with subtle borders.

2. **Color Palette & Typography**:
   - **Background**: Slate Off-White (`#F8FAFC`) or Pure White (`#FFFFFF`).
   - **Text**: Dark Charcoal (`#0F172A` / `#1F2937`) for titles and body text (`#334155`).
   - **Accents**: 
     - Navy / Slate Blue (`#1E3A8A`) for primary highlights and headings.
     - Teal / Emerald (`#0D9488` / `#047857`) for secondary insights and R&D processes.
     - Amber (`#B45309`) for ongoing work & warnings.
     - Purple (`#6B21A8`) for future translational phases.
   - **Fonts**: Clean standard sans-serif (Arial, Calibri, Helvetica).

3. **Status Badges**:
   - `COMPLETED`: Soft blue background (`#EFF6FF`), navy text (`#1E40AF`).
   - `ONGOING`: Soft amber background (`#FEF3C7`), dark amber text (`#92400E`).
   - `FUTURE`: Soft purple background (`#F3E8FF`), dark purple text (`#6B21A8`).

---

## 📐 Slide Architecture & Workflow

Every scientific deck (~9–10 slides for a 20-minute talk) should follow a logical narrative:

1. **Title & Introduction**: Presenter identity, project title, central research question, project status, agenda.
2. **Problem Framing & Context**: Real-world variability vs. regulatory/public health reference points.
3. **Project Architecture**: 4-stage workflow (Process Space $\rightarrow$ Predictive Microbiology $\rightarrow$ Risk Translation $\rightarrow$ Product Translation) with explicit stakeholder mapping.
4. **Experimental Foundation**: Biological matrix, surrogate selection, dynamic heating trajectories $T(t)$, dual assay strategy (e.g., Plaque assay vs. RT-qPCR).
5. **Kinetic Modeling & Outcomes**: Candidate model comparison, selected model rationale (fit, convergence, stability, interpretability), native parameter table, major findings/matrix effects.
6. **Translational Challenges & Boundaries**: Distinguishing process artifacts (come-up time) from biological mechanisms; explicit model domain boundaries.
7. **Quantitative Risk Translation**: Multi-step risk framework (contamination $\rightarrow$ infectivity $\rightarrow$ dynamic ODE $\rightarrow$ dose $\rightarrow$ dose-response $\rightarrow$ illness probability), parameter uncertainty propagation.
8. **Safety × Quality Operating Window**: Two-objective optimization framework (Microbial safety vs. Product quality/sensory attributes).
9. **R&D Transferability & Wrap-up**: Side-by-side comparison of quantitative framework applicability to industrial R&D / fermentation (growth, acidification, post-acidification, robustness to variability).

---

## 💻 python-pptx Template Code Engine

```python
import os
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.enum.text import PP_ALIGN
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE

def create_scientific_deck(output_path):
    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)
    blank_layout = prs.slide_layouts[6]

    # Colors
    C_BG = RGBColor(248, 250, 252)
    C_CARD_BG = RGBColor(255, 255, 255)
    C_CARD_BORDER = RGBColor(226, 232, 240)
    C_TEXT_MAIN = RGBColor(15, 23, 42)
    C_TEXT_SUB = RGBColor(71, 85, 105)
    C_NAVY = RGBColor(30, 58, 138)

    # Slide Setup Helper
    def set_bg(slide):
        bg = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, 0, 0, Inches(13.333), Inches(7.5))
        bg.fill.solid(); bg.fill.fore_color.rgb = C_BG
        bg.line.fill.background()

    # Header Helper
    def add_header(slide, title_text, category_text="SCIENTIFIC PRESENTATION"):
        tb = slide.shapes.add_textbox(Inches(0.8), Inches(0.4), Inches(11.733), Inches(0.95))
        tf = tb.text_frame; tf.word_wrap = True
        tf.margin_left = tf.margin_top = tf.margin_right = tf.margin_bottom = 0
        
        p0 = tf.paragraphs[0]
        p0.text = category_text.upper()
        p0.font.name = "Arial"; p0.font.size = Pt(10); p0.font.bold = True; p0.font.color.rgb = C_NAVY
        p0.space_after = Pt(3)
        
        p1 = tf.add_paragraph()
        p1.text = title_text
        p1.font.name = "Arial"; p1.font.size = Pt(22); p1.font.bold = True; p1.font.color.rgb = C_TEXT_MAIN

    # Card Shape Helper
    def add_card(slide, left, top, width, height, bg_color=C_CARD_BG, border_color=C_CARD_BORDER):
        card = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(left), Inches(top), Inches(width), Inches(height))
        card.fill.solid(); card.fill.fore_color.rgb = bg_color
        if border_color:
            card.line.color.rgb = border_color
            card.line.width = Pt(1)
        else:
            card.line.fill.background()
        return card

    # Speaker Notes Helper
    def set_notes(slide, notes_text):
        if slide.notes_slide:
            slide.notes_slide.notes_text_frame.text = notes_text

    return prs
```

---

## 📝 Best Practices Checklist

- [x] **Native Elements**: Use native PPTX shapes, text boxes, and tables instead of pre-rendered image screenshots.
- [x] **Aspect Ratio Preservation**: When adding figures, scale width/height while preserving original aspect ratio.
- [x] **No Overlap**: Calculate container positions dynamically with safety margins (minimum `0.2"` gap).
- [x] **Speaker Notes**: Attach 150–250 word presenter scripts to every slide for seamless presentation delivery.
