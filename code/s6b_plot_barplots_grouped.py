"""
Creates Fig. 2: the country decomposition of the change in AAD (left) and RP100 (right) for the Philippines, Vietnam
and Cuba (rows), panels A-F. Every panel shows the 8 decomposition terms as three bars: median (thick), p10 (thin,
above) and p90 (thin, below); bars beyond -150 % / +150 % are cut at the edge with an arrow and their value. Right of
each panel: the change in displaced people at the median (median relative change x base value).

Input:  {country}_decomposition_AAD.csv, {country}_decomposition_rp100.csv     (country decompositions)
        all_countries_base_summaries.csv                                     (base values)
Output: fig2_barplots.pdf, .svg, .png
        (folders defined in config.py)
"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory
from config import DECOMPOSITION_DIR, FIGURES_DIR, BASE_SUMMARIES_CSV, FIGURE_COUNTRIES, COUNTRY_NAMES


'''
 ----------------------------------------------------------------------------------------
     I. Setup
 ----------------------------------------------------------------------------------------
'''

# Resolution of the saved PNG
DPI = 600

COUNTRIES = FIGURE_COUNTRIES                                 # rows, top to bottom (config.py)
METRICS = {"AAD": "AAD", "rp100": "RP100"}                   # columns, left to right: metric -> label
R = 150                                                      # x axis from -R to +R (%)

OUT_DIR = FIGURES_DIR / "barplots"

BAR_TERMS = ['ALL', 'HAZ', 'EXP', 'VULN', 'HAZxEXPint', 'HAZxVULNint', 'EXPxVULNint', 'HAZxEXPxVULNint']
BAR_LABELS = ['All', 'H', 'E', 'V', 'H×E', 'H×V', 'E×V', 'H×E×V']

TERM_COLORS = {
    "ALL": "#666666",              # grey
    "HAZ": "#D9070D",              # red
    "EXP": "#3976BB",              # blue
    "VULN": "#198D56",             # green
    "HAZxEXPint": "#a8719d",       # purple
    "HAZxVULNint": "#72190e",      # brown
    "EXPxVULNint": "#298188",      # teal
    "HAZxEXPxVULNint": "k",        # black
}

# sizes in inches (final print size)
FIG_W = 7.25        # figure width: 2 columns
ROW_H = 0.24        # height of one term row (8 rows per panel)
LEFT = 0.62         # left margin: country name + term labels
PEOPLE_W = 0.45     # people column right of each panel
GAP = 0.14          # space between the two panel columns

# text sizes in points
TEXT_PT = 8         # titles of axes, rows and columns
TICK_PT = 7         # tick labels, term labels, legend
MIN_PT = 6          # people column, cut values
PANEL_PT = 9        # panel letters A-F

# bars (heights as a share of one row)
MED_H = 0.38        # median bar (thick)
Q_H = 0.15          # p10 / p90 bars (thin)
Q_GAP = 0.03        # gap between the median and the thin bars
P10_TINT, P90_TINT = 0.70, 0.35   # p10 / p90 colour: term colour mixed with this share of white


'''
 ----------------------------------------------------------------------------------------
     II. Data - Helper Functions
 ----------------------------------------------------------------------------------------
'''

def load_base_value(iso, metric):
    """Base-simulation median of the country and metric (people displaced in the base period); NaN if it is missing."""

    df = pd.read_csv(BASE_SUMMARIES_CSV)
    row = df[(df["country"] == iso) & (df["metric"] == metric) & (df["stat"] == "median")]

    if row.empty:
        return np.nan

    return float(row["value"].iloc[0])


def bar_values(iso, metric):
    """
       Extract median / p10 / p90 of decomposition terms of one country and metric, from the decomposition CSV.
       Output: {stat: array with one value per term of BAR_TERMS, in %}
    """

    table = pd.read_csv(DECOMPOSITION_DIR / iso / f"{iso}_decomposition_{metric}.csv")
    values = {}

    for stat in ("median", "p10", "p90"):
        t = table[table["stat"] == stat].set_index("decomp_term")["delta_change"]
        values[stat] = np.array([t[term] * 100.0 for term in BAR_TERMS])

    return values


def format_people(value):
    """
        Signed number of people with the most readable unit (k for thousands, M for millions);
        Values below 1 expressed with two decimals, NaN values with "n/a".
    """

    if not np.isfinite(value):
        return "n/a"
    abs_val = abs(value)
    if abs_val >= 999_950:
        return f"{value / 1_000_000:+,.2f} M"
    elif abs_val >= 1_000:
        return f"{value / 1_000:+,.1f} k"
    elif abs_val >= 1:
        return f"{value:+,.0f}"
    return f"{value:+.2f}"


def lighter(color, tint):
    """The colour mixed with the share 'tint' of white (used for the p10 / p90 bars)."""

    rgb = np.array(mcolors.to_rgb(color))
    return mcolors.to_hex((1 - tint) * rgb + tint)


'''
 ----------------------------------------------------------------------------------------
     III. Plot Configuration
 ----------------------------------------------------------------------------------------
'''

def set_style():
    """Universal sans-serif font (Arial, else Helvetica / Liberation Sans), thin axes, text kept as text in the pdf / svg."""

    installed = {f.name for f in font_manager.fontManager.ttflist}
    font = next((f for f in ["Arial", "Helvetica", "Liberation Sans"] if f in installed), "DejaVu Sans")
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": [font], "font.size": TICK_PT,
                         "axes.linewidth": 0.4, "xtick.major.width": 0.4, "xtick.major.size": 2.5,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})


def draw_bar(ax, y, v, height, color, edge_lw, where):
    """One horizontal bar. Beyond -R / +R: cut at the edge, with an arrow and its real value
       (above the bar for p10, below it for p90)."""

    if not np.isfinite(v):
        return
    # beyond the axis: bar cut just inside the edge, arrow and value
    side = 1 if v > R else (-1 if v < -R else 0)
    x = 0.975 * R * side if side else v

    ax.barh(y, x, height=height, color=color, edgecolor="#333333" if where == "center" else "#5A5A5A",
            lw=edge_lw, zorder=2)

    if side:
        halo = [pe.withStroke(linewidth=1.8, foreground="white")]
        ax.plot(x, y, marker=">" if side > 0 else "<", markersize=5, color="#222222", zorder=6, path_effects=halo)
        dy, va = {"above": (1.6, "bottom"), "below": (-1.6, "top"), "center": (0, "center")}[where]
        ax.annotate(f"{v:+,.0f}%", (x, y), xytext=(-3.0 * side, dy), textcoords="offset points",
                    ha="right" if side > 0 else "left", va=va, fontsize=MIN_PT, fontweight="bold",
                    color="#222222", zorder=7, path_effects=[pe.withStroke(linewidth=1.5, foreground="white")])


def draw_bar_panel(ax, values):
    """Draw the three bars (p10, median, p90) of every term of one country and metric. ALL on top."""

    # one row per term, ALL on top with a grey band; zero line
    y = np.arange(len(BAR_TERMS))
    ax.set_xlim(-R, R)
    ax.set_ylim(len(BAR_TERMS) - 0.5, -0.5)
    ax.axhspan(-0.5, 0.5, color="#EDEDED", lw=0, zorder=0)
    ax.axhline(0.5, color="#9A9A9A", lw=0.4, ls=(0, (1, 1.5)), zorder=1)
    ax.axvline(0, color="black", lw=0.6, zorder=3)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=3)
    ax.set_yticks(y)
    ax.set_xticks([-150, -100, -50, 0, 50, 100, 150])

    # p10 (thin, above), median (thick), p90 (thin, below)
    off = MED_H / 2 + Q_GAP + Q_H / 2
    for k, term in enumerate(BAR_TERMS):
        c = TERM_COLORS[term]
        draw_bar(ax, y[k] - off, values["p10"][k], Q_H, lighter(c, P10_TINT), 0.3, "above")
        draw_bar(ax, y[k], values["median"][k], MED_H, c, 0.4, "center")
        draw_bar(ax, y[k] + off, values["p90"][k], Q_H, lighter(c, P90_TINT), 0.3, "below")

    return y


'''
 ----------------------------------------------------------------------------------------
     IV. Plot grouped figure
 ----------------------------------------------------------------------------------------
'''

def plot_barplots_grouped():
    """Countries in rows, AAD and RP100 in columns, people column right of each panel; saved as fig2_barplots."""

    set_style()
    # layout in inches: panel size, titles above, legend below, gaps below each row
    panel_h = len(BAR_TERMS) * ROW_H
    panel_w = (FIG_W - LEFT - 2 * PEOPLE_W - GAP) / 2
    top, bottom = 0.49, 0.26
    gaps = [0.18, 0.18, 0.44]
    fig_h = top + len(COUNTRIES) * panel_h + sum(gaps) + bottom
    fig = plt.figure(figsize=(FIG_W, fig_h))

    letters = iter("ABCDEF")
    y_top = fig_h - top
    for r, iso in enumerate(COUNTRIES):
        y0 = y_top - panel_h
        for c, (metric, label) in enumerate(METRICS.items()):
            x0 = LEFT + c * (panel_w + PEOPLE_W + GAP)
            ax = fig.add_axes([x0 / FIG_W, y0 / fig_h, panel_w / FIG_W, panel_h / fig_h])
            values = bar_values(iso, metric)
            base_val = load_base_value(iso, metric)
            y = draw_bar_panel(ax, values)

            # term labels and country name only on the left
            if c == 0:
                ax.set_yticklabels(BAR_LABELS, fontsize=TICK_PT)
                ax.get_yticklabels()[0].set_fontweight("bold")
                ax.set_ylabel(COUNTRY_NAMES.get(iso, iso), fontsize=TEXT_PT, fontweight="bold", labelpad=3)
            else:
                ax.set_yticklabels([])

            # x axis only under the last row
            if r == len(COUNTRIES) - 1:
                ax.tick_params(axis="x", labelsize=TICK_PT, pad=1.5)
                fig.text((x0 + panel_w / 2) / FIG_W, (y0 - 0.20) / fig_h, f"{label} relative change (%)",
                         ha="center", va="top", fontsize=TEXT_PT)
            else:
                ax.tick_params(axis="x", labelbottom=False, length=0)

            # metric title (top row) and panel letter
            y_text = (y0 + panel_h + (0.31 if r == 0 else 0.03)) / fig_h
            if r == 0:
                fig.text((x0 + (panel_w + PEOPLE_W) / 2) / FIG_W, y_text, label, ha="center", va="bottom",
                         fontsize=TEXT_PT, fontweight="bold")
            fig.text((0.03 if c == 0 else x0 - 0.10) / FIG_W, y_text, next(letters), ha="left", va="bottom",
                     fontsize=PANEL_PT, fontweight="bold")

            # change in displaced people (median change x base value), right of the panel
            x_col = (x0 + panel_w + PEOPLE_W - 0.02) / FIG_W
            to_row = blended_transform_factory(fig.transFigure, ax.transData)
            for k, v in enumerate(values["median"]):
                people = "n/a" if base_val == 0 else format_people(v / 100.0 * base_val)
                ax.text(x_col, y[k], people, transform=to_row, ha="right", va="center", fontsize=MIN_PT,
                        fontweight="bold" if k == 0 else "normal", clip_on=False)
            if r == 0:
                ax.text(x_col, 1 + 0.03 / panel_h, "Change in\ndisplaced\npeople", ha="right", va="bottom",
                        transform=blended_transform_factory(fig.transFigure, ax.transAxes),
                        fontsize=MIN_PT, color="#444444", linespacing=1.0, clip_on=False)

        y_top = y0 - gaps[r]

    # legend below the figure
    g = "#7A7A7A"
    handles = [Patch(facecolor=lighter(g, P10_TINT), edgecolor="#5A5A5A", lw=0.3, label="p10 (thin, above)"),
               Patch(facecolor=g, edgecolor="#333333", lw=0.4, label="Median"),
               Patch(facecolor=lighter(g, P90_TINT), edgecolor="#5A5A5A", lw=0.3, label="p90 (thin, below)"),
               Line2D([], [], ls="none", marker=">", markersize=5, color="#222222", label="Beyond axis (value shown)")]
    fig.legend(handles=handles, ncol=4, frameon=False, fontsize=TICK_PT, loc="lower center",
               bbox_to_anchor=((LEFT + (FIG_W - LEFT) / 2) / FIG_W, 0.0), handlelength=1.4, columnspacing=1.4,
               borderaxespad=0.1)

    name = "fig2_barplots"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / f"{name}.pdf", facecolor="white", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(OUT_DIR / f"{name}.svg", facecolor="white", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(OUT_DIR / f"{name}.png", dpi=DPI, facecolor="white", bbox_inches="tight", pad_inches=0.02)
    print(f"    saved {OUT_DIR / name}.png")
    plt.close(fig)


def main():
    plot_barplots_grouped()
    print(f"Done -> {OUT_DIR}")

#----------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
