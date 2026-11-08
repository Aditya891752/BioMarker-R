# Design System & UI Specification (`design.md`)

## 1. Design Philosophy: The High-Precision Laboratory Instrument
BioMarker-R is a computational genomics application, not a generic consumer SaaS product. The design philosophy is modeled after **modern clinical diagnostic workstations and peer-reviewed journals (Nature Genetics, Cell)**:
- **Maximum Contrast & Legibility**: Data points, gene identifiers, and statistical values must be instantly readable without eye strain.
- **Scientific Honesty**: Visualizations present raw statistical boundaries without distortion, glossy covers, or novelty animations.
- **Editorial Cleanliness**: Generous structural white space, strict geometric grid alignment, and typography chosen specifically for data tables and scientific figures.

---

## 2. Typography System

| Role | Font Family | Fallback | Weights | Usage |
| :--- | :--- | :--- | :--- | :--- |
| **Headings & Display** | **Plus Jakarta Sans** / **IBM Plex Sans** | `system-ui, sans-serif` | 600 (Semibold), 700 (Bold) | Section headers, metric labels, titles. |
| **Data & Scientific Codes** | **IBM Plex Mono** | `ui-monospace, monospace` | 400 (Regular), 600 (Medium) | Gene symbols, Ensembl IDs, $p$-values, $\log_2\text{FC}$, coordinate readouts. |
| **Body & Narrative** | **Source Sans 3** / **IBM Plex Sans** | `Helvetica, Arial, sans-serif` | 400 (Regular), 500 (Medium) | Biological interpretations, methodological notes. |

### Strict Typography Prohibitions
- **NO Inter as universal font**.
- **NO Space Grotesk + Instrument Serif pairing**.
- **NO faux-editorial serif italic flourishes** inside titles.
- **NO emojis in headings or navigation**.

---

## 3. Color Architecture (Scientific Chromatic Standard)

```
LIGHT THEME (Default)
├── Base Background:      #FFFFFF (Clean Pure White)
├── Surface Container:    #F8F9FA (Lab Off-White / Subdued Panel)
├── Structural Border:    #E2E8F0 (1px Solid Border - High Precision)
├── Primary Text:         #0F172A (Deep Slate 900 - High Contrast)
├── Muted / Secondary:    #475569 (Slate 600 - Minimum 4.8:1 Contrast)
│
DATA & STATISTICAL ACCENTS (Zero gradients, zero purple-blue linear fills)
├── Significant Upregulated:    #B91C1C (Crimson Red / High Expression)
├── Significant Downregulated:  #1D4ED8 (Cobalt Navy / Low Expression)
├── Non-Significant Gene:       #94A3B8 (Neutral Cool Slate)
├── Clinical Normal Cohort:     #0D9488 (Teal Diagnostic / Matched Control)
├── TNBC Tumor Cohort:          #BE185D (Deep Rose Plum / Tumor Group)
└── Key Significant Biomarker:  #047857 (Emerald Verification Marker)

HIGH-CONTRAST DARK THEME (Optional)
├── Base Background:      #0B0F17 (Deep Obsidian - Absolute contrast)
├── Surface Container:    #161F2E (Navy Charcoal Surface)
├── Structural Border:    #263346 (Crisp 1px Divider)
├── Primary Text:         #F8FAFC (Slate 50 - > 14:1 Contrast)
└── Muted Text:           #94A3B8 (Cool Gray - > 6.5:1 Contrast)
```

---

## 4. Layout Grid & Spatial Rhythm
All spacing follows a strict 8-point geometric scale:
- `4px` (micro spacing, tag padding, table cell vertical offsets)
- `8px` (compact spacing between labels and values)
- `16px` (standard component gap and input padding)
- `24px` (card internal padding, column separation)
- `32px` (section grouping breaks)
- `48px` (major page block offsets)

### Structural Layout Rules
- **No 3 Icon Boxes in a Row**: Metric statistics are laid out in a continuous diagnostic strip with direct statistical indicators (e.g. `Cohort: 154 Samples | Filtered Genes: 14,210 | Significant DEGs: 48`).
- **No Badges Above Headlines**: Titles are straightforward, authoritative, and direct (e.g., `Differential Expression: TNBC vs. Solid Tissue Normal`).
- **No Cliché Glassmorphism**: Cards use solid white backgrounds with 1px solid borders (`border border-slate-200`). No backdrop filters or blurry panels.

---

## 5. Component Specifications

### 5.1 Volcano Plot Component
- **Engine**: SVG / Canvas (Plotly.js or native D3 / SVG).
- **Axes**:
  - X-Axis: $\log_2(\text{Fold Change})$ from $-8.0$ to $+8.0$.
  - Y-Axis: $-\log_{10}(\text{Adjusted } p\text{-value})$ from $0$ to $\ge 50$.
- **Threshold Lines**: Discrete dashed lines at $\text{FDR} = 0.01$ ($-\log_{10} = 2.0$) and $|\log_2\text{FC}| = 2.0$.
- **Tooltip**: Hover card displaying:
  - Gene Symbol (e.g., `FOXA1`)
  - Ensembl ID (`ENSG00000129514`)
  - $\log_2\text{FC}$ and Adjusted $p$-value
  - Average Base Mean Expression

### 5.2 Differential Expression Data Table
- **Header**: Monospaced column names (`GENE_SYMBOL`, `ENSEMBL_ID`, `BASE_MEAN`, `LOG2_FC`, `LFC_SE`, `P_VAL`, `ADJ_P_VAL`).
- **Sorting**: Multi-column ascending/descending click sorting.
- **Numbers**: Always right-aligned with tabular numerals (`font-mono tabular-nums`).
- **Search**: Instant substring search on gene symbols.

### 5.3 Kaplan-Meier Survival Viewer
- **Curves**: Two high-contrast step functions:
  - Crimson Solid Line (`High Expression`)
  - Cobalt Solid Line (`Low Expression`)
- **Censoring Ticks**: Sharp vertical crosshairs on each step indicating patient censoring.
- **Risk Table**: Numerical table beneath the curve showing numbers at risk at 0, 12, 24, 36, 48, and 60 months.
- **Statistical Inset**: Clean box presenting Log-Rank test $p$-value and Cox Hazard Ratio ($95\%\text{ CI}$).

---

## 6. Interaction & Motion Rules
- **No Fade-in on Scroll**: Content loads immediately without layout shifts.
- **No Cursor Beams or Trail Animations**: Pointer maintains default OS behavior.
- **No Button Opacity Fades**: Buttons transition cleanly between solid resting state and 1px pressed outline or inverted fill ($< 100\text{ ms}$ crisp transition).
- **Keyboard Navigation**: Full tab navigation across all filter sliders, tables, and search inputs with a distinct 2px black focus ring (`focus:outline-none focus:ring-2 focus:ring-slate-900`).
