"""
Creates Fig. 3: the country maps of the Philippines, Vietnam and Cuba (rows). Columns: the relative change of the
exposure from 2005 to SSP3 2050, and the dominant driver of the EAD change at every exposure point at p10, median and
p90 (full colour: increase, light shade: decrease). Every grid cell (300 arcsec) is drawn as a tile of its real size.

Input:  {country}_EAD_exp_point_decomposition.csv                          (decomposition at every exposure point)
        {country}_base_2005_300as.h5, {country}_SSP3_2050_300as.h5          (base and future exposure)
Output: country_maps_grouped.png, .pdf, .svg
        (folders defined in config.py)
"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.legend as mlegend
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import SymLogNorm, BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, MultipleLocator

from climada.entity import Exposures
from climada.util.coordinates import get_country_geometries

from config import DECOMPOSITION_DIR, EXP_DIR, FIGURES_DIR, FIGURE_COUNTRIES, COUNTRY_NAMES


'''
 ----------------------------------------------------------------------------------------
     I. Setup
 ----------------------------------------------------------------------------------------
'''
 
OUT_DIR = FIGURES_DIR / "country_map"
DPI = 600                                   # resolution of the saved PNG

GRID_STEP_DEG = 300 / 3600                  # grid step of the exposures (300 arcsec): size of one tile

FIGURE_STATS = ['p10', 'median', 'p90']     # columns 2-4

# The 7 decomposition terms (without ALL) and their legend labels
SPATIAL_TERMS = ['HAZ', 'EXP', 'VULN', 'HAZxEXPint', 'HAZxVULNint', 'EXPxVULNint', 'HAZxEXPxVULNint']
SPATIAL_LABELS = ['HAZ', 'EXP', 'VULN', r'HAZxEXP$_{\mathrm{int}}$', r'HAZxVULN$_{\mathrm{int}}$', r'EXPxVULN$_{\mathrm{int}}$', r'HAZxEXPxVULN$_{\mathrm{int}}$']

# Colour of every term (full colour: increase)
TERM_COLORS = {
    "HAZ": "#D9070D",              # red
    "EXP": "#3976BB",              # blue
    "VULN": "#198D56",             # green
    "HAZxEXPint": "#a8719d",       # purple
    "HAZxVULNint": "#72190e",      # brown
    "EXPxVULNint": "#298188",      # teal
    "HAZxEXPxVULNint": "k",        # black
}
TWO_COLOUR_TINT = 0.50             # decrease: the term colour mixed 50 % with white

# Exposure maps
EXPOSURE_BASE_NAME = "base_2005_300as.h5"
EXPOSURE_FUTURE_NAME = "SSP3_2050_300as.h5"
EXPOSURE_PALETTE = "RdBu_r"                         # blue: decrease, red: increase
EXPOSURE_BINS = [-100, -50, -25, 0, 25, 50, 100, 200]   # colour classes (%); None: continuous SymLog scale

# Map style
SEA_COLOR = "#F1F5F9"
LAND_COLOR = "#FFFFFF"
NEIGHBOUR_BORDER_STYLE = dict(color="#9CA3AF", linewidth=0.35)
BORDER_PCT = 0.04                                   # margin around the country (share of its extent)

# Layout (inches) and font sizes (pt), at the final print size (2 columns)
FIG_WIDTH_IN = 7.1
LAYOUT_IN = dict(left=0.6, col_gap=0.35, shared_gap=0.17, right=0.05, top=0.17, row_gap=0.3, bottom=0.05)
LAT_TICK_COLUMNS = [0]                              # columns with latitude tick labels (0: exposure)
LON_TITLE_GAP_IN = 0.12
MAP_MAX_TICKS = 4
PANEL_LETTER_FONTSIZE = 9
TITLE_FONTSIZE = 8
AXIS_FONTSIZE = 7
TICK_FONTSIZE = 6
LEGEND_FONTSIZE = 6.5
LEGEND_NCOL = 3                                     # entries per legend line
GROUP_TITLE = "Dominant driver of the change in EAD"                                   # title over the three driver columns
GROUP_TITLE_IN = 0.2                                                                       # its height above the column titles
LON_TICK_STEP = {"VNM": 5}                                                                 # longitude tick spacing (deg) for these countries


'''
 ----------------------------------------------------------------------------------------
     II. Data
 ----------------------------------------------------------------------------------------
'''


def load_points(iso):
    """
    Extracts tables of points with a dominant driver, with the column sign (1 increase, -1 decrease), for each statistic
    """
    df = pd.read_csv(DECOMPOSITION_DIR / iso / f"{iso}_EAD_exp_point_decomposition.csv")
    df = df[df["dominant_driver"] != "NO_DATA"]
    points = {}
    for stat in FIGURE_STATS:
        p = df[df["stat"] == stat].copy()
        p["sign"] = np.where(p["dominant_value"] >= 0, 1, -1)
        points[stat] = p
    return points


def align_exposures(exp_base, exp_future):
    """
    Base and future exposure values extraction of every grid point, paired by row as in the simulations
    Relative change where the base exposure is > 0 (NaN where the base exposure is 0). Coordinates from the base exposure.
    """
    merged = pd.DataFrame({"lon": exp_base.gdf.geometry.x.values,
                           "lat": exp_base.gdf.geometry.y.values,
                           "value_base": exp_base.gdf["value"].values,
                           "value_future": exp_future.gdf["value"].values})
 
    merged["pct_diff"] = (merged["value_future"] - merged["value_base"]) / merged["value_base"] * 100.0
    merged.loc[merged["value_base"] <= 0, "pct_diff"] = np.nan
    
    return merged


def load_exposure_change(iso):
    """Uses align exposure def to calculate relative exposure change of one country, importing its base and future exposures"""
    
    exp_base = Exposures.from_hdf5(EXP_DIR / iso / f"{iso}_{EXPOSURE_BASE_NAME}")
    exp_future = Exposures.from_hdf5(EXP_DIR / iso / f"{iso}_{EXPOSURE_FUTURE_NAME}")
    
    return align_exposures(exp_base, exp_future)


def exposure_norm(exposures):
    """
    Returns a color scale for each % change magnitude, following the SymLogNorm scale, which is the same for the countries of the figure.
    Scale is computed from their exposure changes (linear up to the median of |% change|, up to the largest |% change|).
    exposures: {iso: table of load_exposure_change}
    """
    changes = []
    
    for exp in exposures.values():
        changes.append(exp["pct_diff"].dropna().values)
   
    abs_changes = np.abs(np.concatenate(changes))

    abs_changes = abs_changes[abs_changes > 0]
    linthresh = max(np.percentile(abs_changes, 50), 1.0)
    vmax = np.max(abs_changes)
    
    scale = SymLogNorm(linthresh=linthresh, vmin=-vmax, vmax=vmax, base=10)
    
    return scale


def exposure_classes():
    """
    Returns colour map and scale of the fixed exposure classes (EXPOSURE_BINS): blues for decreases, reds for increases
    (darker = larger change) and one more class, the darkest red, above the last edge.
    """
    n_neg = sum(1 for hi in EXPOSURE_BINS[1:] if hi <= 0)
    n_pos = len(EXPOSURE_BINS) - 1 - n_neg
    palette = plt.get_cmap(EXPOSURE_PALETTE)
    neg = [palette(v) for v in np.linspace(0.05, 0.40, n_neg)]
    pos = [palette(v) for v in np.linspace(0.60, 0.98, n_pos + 1)]
    cmap = ListedColormap(neg + pos[:-1])
    cmap.set_over(pos[-1])
    norm = BoundaryNorm(EXPOSURE_BINS, cmap.N)

    return cmap, norm


'''
 ----------------------------------------------------------------------------------------
     III. Drawing of one map - Helper Functions
 ----------------------------------------------------------------------------------------
'''


def set_font():
    """Universal sans-serif font (Arial, else Helvetica / Liberation Sans), embedded in the pdf, text kept as text in the svg."""
    installed = {f.name for f in font_manager.fontManager.ttflist}
    font = next((f for f in ["Arial", "Helvetica", "Liberation Sans"] if f in installed), "DejaVu Sans")
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": [font], "pdf.fonttype": 42, "svg.fonttype": "none",
                         "mathtext.fontset": "custom", "mathtext.rm": font, "mathtext.bf": f"{font}:bold"})   # also in math text


def colour(driver, sign):
    """Returns a colour in hex format for the given driver, in the shade that corresponds to its sign."""

    base = np.array(mcolors.to_rgb(TERM_COLORS[driver]))
    if sign >= 0:
        rgb = base
    else:
        rgb = (1 - TWO_COLOUR_TINT) * base + TWO_COLOUR_TINT

    return mcolors.to_hex(rgb)


def map_extent(boundary):
    """
    Area shown in the map: the rectangle around the country's outline + BORDER_PCT margin
    """

    minx, miny, maxx, maxy = boundary.total_bounds
    px = BORDER_PCT * (maxx - minx)
    py = BORDER_PCT * (maxy - miny)

    lon_range = (minx - px, maxx + px)
    lat_range = (miny - py, maxy + py)

    return lon_range, lat_range


def to_tiles(lon, lat, values):
    """
    Puts the points on the regular grid (GRID_STEP_DEG): cell centres (x, y) and a 2D array with the value of every cell.
    Column / row of a point: its distance from the westernmost / southernmost point, in grid steps. Empty cells are masked for pcolormesh (not drawn).
    """
    lon, lat, values = np.asarray(lon), np.asarray(lat), np.asarray(values, dtype=float)

    # column / row of every point on the grid
    ix = np.round((lon - lon.min()) / GRID_STEP_DEG).astype(int)
    iy = np.round((lat - lat.min()) / GRID_STEP_DEG).astype(int)

    # grid with one value per cell, empty cells NaN
    grid = np.full((iy.max() + 1, ix.max() + 1), np.nan)
    grid[iy, ix] = values

    # cell centres in degrees
    x = lon.min() + np.arange(ix.max() + 1) * GRID_STEP_DEG
    y = lat.min() + np.arange(iy.max() + 1) * GRID_STEP_DEG

    return x, y, np.ma.masked_invalid(grid)


def draw_map(ax, boundary, neighbours, xlim, ylim, lat_ticks=True, lat_title=True):
    """
    Draws elements of map without data:
    sea, land, borders, grid, lon/lat degree ticks etc.
    lat_ticks: latitude tick labels shown; lat_title: "Latitude" shown ("Longitude": once per row, plot_country_maps_grouped)
    """

    # sea and land
    ax.set_facecolor(SEA_COLOR)
    if not neighbours.empty:
        neighbours.plot(ax=ax, facecolor=LAND_COLOR, edgecolor="none", zorder=1)
        neighbours.boundary.plot(ax=ax, zorder=9, **NEIGHBOUR_BORDER_STYLE)
    boundary.plot(ax=ax, facecolor=LAND_COLOR, edgecolor="none", zorder=1)
    boundary.boundary.plot(ax=ax, color="#4B5563", linewidth=0.5, zorder=10)

    # axes, ticks and labels
    ax.set_aspect("equal")
    ax.set_xlabel("")
    ax.set_ylabel("Latitude" if lat_title else "", fontsize=AXIS_FONTSIZE, labelpad=2)

    # tick labels in degrees (E/W, N/S)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{abs(v):g}°{'E' if v >= 0 else 'W'}"))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{abs(v):g}°{'N' if v >= 0 else 'S'}"))

    # at most MAP_MAX_TICKS round tick values per axis
    ax.xaxis.set_major_locator(MaxNLocator(nbins=MAP_MAX_TICKS, steps=[1, 2, 2.5, 5, 10]))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=MAP_MAX_TICKS, steps=[1, 2, 2.5, 5, 10]))

    ax.tick_params(axis="both", labelsize=TICK_FONTSIZE, length=2, width=0.5, pad=1.5)
    ax.tick_params(axis="y", labelleft=lat_ticks)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.spines[["top", "right"]].set_visible(False)


def draw_exposure_map(ax, boundary, neighbours, exposure, cmap, norm, xlim, ylim, lat_ticks=True, lat_title=True):
    """Draws the relative exposure change (%) at every cell with base exposure, as tiles"""

    draw_map(ax, boundary, neighbours, xlim, ylim, lat_ticks, lat_title)

    changed = exposure[exposure["value_base"] > 0]
    x, y, grid = to_tiles(changed["lon"], changed["lat"], changed["pct_diff"])

    ax.pcolormesh(x, y, grid, cmap=cmap, norm=norm, alpha=0.90, shading="nearest", edgecolors="none", antialiased=False, zorder=3, rasterized=True)


def draw_driver_map(ax, boundary, neighbours, points, xlim, ylim, lat_ticks=True, lat_title=True):
    """Draws dominant driver at every cell as tiles, with the colour of the dominant driver, full for an increase, light shade for a decrease."""

    draw_map(ax, boundary, neighbours, xlim, ylim, lat_ticks, lat_title)

    # Every driver and sign is one class (0 ... 13) with its own colour: one tile layer for all cells
    classes = [(d, s) for d in SPATIAL_TERMS for s in (1, -1)]
    code = {c: i for i, c in enumerate(classes)}
    values = [code[(d, s)] for d, s in zip(points["dominant_driver"], points["sign"])]
    x, y, grid = to_tiles(points["longitude"], points["latitude"], values)

    cmap = ListedColormap([colour(d, s) for d, s in classes])
    norm = BoundaryNorm(np.arange(len(classes) + 1) - 0.5, len(classes))
    ax.pcolormesh(x, y, grid, cmap=cmap, norm=norm, shading="nearest", edgecolors="none", antialiased=False, zorder=4, rasterized=True)


def draw_exposure_colorbar(cax, norm, local_max, fontsize):
    """Creation of horizontal colour bar of the relative exposure change, from -100 % (smallest possible change) to the largest change of the figure."""
    
    cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=EXPOSURE_PALETTE), cax=cax, orientation="horizontal")
    top = min(local_max, norm.vmax) 
    cax.set_xlim(-100.0, top)
    
    # ticks
    ticks = [norm.linthresh]
    while ticks[-1] * 10 <= norm.vmax:
        ticks.append(ticks[-1] * 10)
    ticks = [-t for t in reversed(ticks)] + [0] + ticks

    start = norm(-100.0)
    span = float(norm(top) - start)

    # keep the ticks of the visible part of the bar, away from -100 %
    kept = [-100.0]
    for t in ticks:
        on_bar = -100.0 < t <= top
        safe_dist = (norm(t) - start) / span > 0.2
        if on_bar and safe_dist:
            kept.append(t)
    ticks = kept
    
    cb.set_ticks(ticks)
    cb.set_ticklabels([("0%" if t == 0 else f"{t:+.0f}%") for t in ticks])
    cb.ax.tick_params(labelsize=fontsize)
    cb.set_label("Relative exposure change (%)", fontsize=fontsize)


def draw_exposure_class_bar(cax, cmap, norm, fontsize):
    """Creation of horizontal colour bar of the exposure classes: one block per class, the arrow for the class above the last edge."""

    cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal", extend="max", spacing="uniform")
    cb.set_ticks(EXPOSURE_BINS)
    cb.set_ticklabels([f"{t:.0f}" for t in EXPOSURE_BINS])
    cb.ax.tick_params(labelsize=TICK_FONTSIZE, length=2, width=0.5, pad=1.5)
    cb.outline.set_linewidth(0.5)
    cb.set_label("Relative exposure change (%)", fontsize=fontsize)


def draw_legend(ax, all_points, fontsize):
    """Designs the driver's Legend with colour shades for each driver and sign"""

    # frame without ticks
    ax.set_xticks([])
    ax.set_yticks([])
 
    for spine in ax.spines.values():                     
        spine.set_edgecolor("#9CA3AF")
        spine.set_linewidth(0.5)
 
    present = set(zip(all_points["dominant_driver"], all_points["sign"]))

    # driver / sign combinations shown in the maps
    increase = []
    decrease = []
    for d in SPATIAL_TERMS:
        if (d, 1) in present:
            increase.append(d)
        if (d, -1) in present:
            decrease.append(d)

    rows = []
 
    if increase:
        rows.append((1, "$\\bf{Increase}$", increase))
    if decrease:
        rows.append((-1, "$\\bf{Decrease}$", decrease))

    # separator between Increase and Decrease
    if len(rows) == 2:                                                  
        ax.axvline(0.5, ymin=0.12, ymax=0.88, color="#9CA3AF", linewidth=0.5)

    for r, (sign, header, drivers) in enumerate(rows):

        handles = []
        labels = []
        for d in drivers:
            dot = Line2D([0], [0], marker="s", linestyle="none", markersize=5, markeredgewidth=0, color=colour(d, sign))
            handles.append(dot)
            labels.append(SPATIAL_LABELS[SPATIAL_TERMS.index(d)])

        # one legend group per sign, side by side
        leg = mlegend.Legend(ax, handles, labels, title=header, handlelength=1.0, ncol=min(len(drivers), LEGEND_NCOL), loc="center",
                             bbox_to_anchor=((r + 0.5) / len(rows), 0.5), fontsize=fontsize,
                             title_fontsize=fontsize, frameon=False, columnspacing=0.8, handletextpad=0.3, borderpad=0.2, labelspacing=0.3)

        ax.add_artist(leg)


'''
 ----------------------------------------------------------------------------------------
     IV. Figure with all countries
 ----------------------------------------------------------------------------------------
'''


def shared_extents(boundaries):
    """
    Extracts map area (lon_range, lat_range) of every country of the figure, to maintain ratios and panel size among countries.
    All panels have the same width, so each row's height follows from its map's shape. 
    Countries taller than wide get the same shape, that of the widest of them, so their rows have the same height.
    Their map area is widened in terms of degrees to the sides around its centre. 
    """

    extents = {}
    ratio = {}
    for iso, boundary in boundaries.items(): 
        lon_range, lat_range = map_extent(boundary)
        extents[iso] = (lon_range, lat_range)         
        ratio[iso] = (lon_range[1] - lon_range[0]) / (lat_range[1] - lat_range[0])
    
    upright = []
    for iso in boundaries:
        if ratio[iso] <= 1:
            upright.append(iso)    
    
    target = 0
    for iso in upright:
        if ratio[iso] > target:
            target = ratio[iso]

    # Widens the narrower countries to the target shape, around their centre (the height stays)
    for iso in upright:
        
        (x0, x1), (y0, y1) = extents[iso]

        cx = (x0 + x1) / 2
        cy = (y0 + y1) / 2
        
        w = x1 - x0
        h = y1 - y0                                                   
        
        if w / h < target:
            w = h * target
    
        extents[iso] = ((cx - w / 2, cx + w / 2), (cy - h / 2, cy + h / 2))
    
    return extents


def plot_country_maps_grouped():
    """
    Builds final plot, with one row per country of FIGURE_COUNTRIES: exposure | p10 | median | p90
    with the exposure colour bar and the driver legend below"""

    set_font()

    # ----Data----
    
    boundaries, neighbours, points, exposures = {}, {}, {}, {}

    for iso in FIGURE_COUNTRIES:
        boundaries[iso] = get_country_geometries([iso])
        points[iso] = load_points(iso)
        exposures[iso] = load_exposure_change(iso)

    # map area of every country
    extents = shared_extents(boundaries) 

    # exposure colour scale: fixed classes, or SymLog scale from the countries of the figure
    if EXPOSURE_BINS:
        exp_cmap, norm = exposure_classes()
    else:
        exp_cmap, norm = EXPOSURE_PALETTE, exposure_norm(exposures)

    # neighbouring countries in every map area
    for iso in FIGURE_COUNTRIES: 
        (x0, x1), (y0, y1) = extents[iso]
        gdf = get_country_geometries(extent=(x0, x1, y0, y1))
        neighbours[iso] = gdf[gdf["ISO_A3"] != iso]

    # ----Layout----

    L = LAYOUT_IN
    FW = FIG_WIDTH_IN

    n_panels = 1 + len(FIGURE_STATS)
    gaps = [L["col_gap"] if k in LAT_TICK_COLUMNS else L["shared_gap"] for k in range(1, n_panels)]
    pw = (FW - L["left"] - L["right"] - sum(gaps)) / n_panels

    # row height of every country, from its map shape
    row_h = {}
    for iso in FIGURE_COUNTRIES:
        (x0, x1), (y0, y1) = extents[iso]
        row_h[iso] = pw * (y1 - y0) / (x1 - x0)

    # all driver points (for the legend)
    tables = []
    
    for iso in FIGURE_COUNTRIES:
        for stat in FIGURE_STATS:
            tables.append(points[iso][stat])
        
    all_points = pd.concat(tables)
    
    n_increase = all_points.loc[all_points["sign"] == 1, "dominant_driver"].nunique()
    n_decrease = all_points.loc[all_points["sign"] == -1, "dominant_driver"].nunique()

    # legend lines needed for the drivers shown
    legend_lines = max(int(np.ceil(n_increase / LEGEND_NCOL)), int(np.ceil(n_decrease / LEGEND_NCOL)), 1)
    
    # height of the legend box and colour bar
    bottom_h = max(0.2 + 0.13 * legend_lines, 0.45)

    FH = L["top"] + GROUP_TITLE_IN + sum(row_h.values()) + len(FIGURE_COUNTRIES) * L["row_gap"] + bottom_h + L["bottom"]
    fig = plt.figure(figsize=(FW, FH))

    def box(x_in, y_top_in, w_in, h_in):
        
        """Rectangle given in inches from the top-left corner -> [left, bottom, width, height] as figure fractions."""
        
        left = x_in / FW
        bottom = 1 - (y_top_in + h_in) / FH
        width = w_in / FW
        height = h_in / FH
        
        return [left, bottom, width, height]
        
    # ----Maps Creation----

    # x position of every column: exposure, p10, median, p90
    col_x = [L["left"]]
    for gap in gaps:
        col_x.append(col_x[-1] + pw + gap)

    y = L["top"] + GROUP_TITLE_IN
    letters = ["A", "B", "C"]
    
    for r, iso in enumerate(FIGURE_COUNTRIES):
        
        xlim, ylim = extents[iso]
        row_height = row_h[iso]

        # row label: panel letter and country name
        fig.text(0.07 / FW, 1 - (y + row_height / 2) / FH, f"({letters[r]}) {COUNTRY_NAMES[iso]}",
                 rotation=90, ha="center", va="center", fontsize=PANEL_LETTER_FONTSIZE, fontweight="bold")

        # column 1: exposure change
        ax = fig.add_axes(box(col_x[0], y, pw, row_height))
        draw_exposure_map(ax, boundaries[iso], neighbours[iso], exposures[iso], exp_cmap, norm, xlim, ylim,
                          lat_ticks=0 in LAT_TICK_COLUMNS, lat_title=True)
        if iso in LON_TICK_STEP:
            ax.xaxis.set_major_locator(MultipleLocator(LON_TICK_STEP[iso]))
        
        if r == 0:
            ax.set_title("Exposure rel. change", fontsize=TITLE_FONTSIZE, fontweight="bold", pad=3)

        # columns 2-4: dominant driver at p10, median, p90
        for k, stat in enumerate(FIGURE_STATS, start=1):
            
            ax = fig.add_axes(box(col_x[k], y, pw, row_height))
            draw_driver_map(ax, boundaries[iso], neighbours[iso], points[iso][stat], xlim, ylim,
                            lat_ticks=k in LAT_TICK_COLUMNS, lat_title=False)
            if iso in LON_TICK_STEP:
                ax.xaxis.set_major_locator(MultipleLocator(LON_TICK_STEP[iso]))
            if r == 0:
                ax.set_title(stat, fontsize=TITLE_FONTSIZE, fontweight="bold", pad=3)

        # One "Longitude" title per row, centred below the four maps (below their tick labels)
        row_centre = (col_x[0] + col_x[-1] + pw) / 2
        fig.text(row_centre / FW, 1 - (y + row_height + LON_TITLE_GAP_IN) / FH, "Longitude", ha="center", va="top", fontsize=AXIS_FONTSIZE)

        y += row_height + L["row_gap"]

    # Group title over the three driver columns, with a line under it spanning these columns
    x_left, x_right = col_x[1], col_x[-1] + pw
    fig.text((x_left + x_right) / 2 / FW, 1 - 0.07 / FH, GROUP_TITLE, ha="center", va="center",
             fontsize=TITLE_FONTSIZE, fontweight="bold")
    fig.add_artist(Line2D([x_left / FW, x_right / FW], [1 - 0.17 / FH, 1 - 0.17 / FH],
                          transform=fig.transFigure, color="#4B5563", linewidth=0.6))

    #----Colour Bar and Legend creation----

    # largest exposure increase (end of the colour bar)
    largest_increase = 0
    for iso in FIGURE_COUNTRIES:
        country_max = exposures[iso]["pct_diff"].max()
        if country_max > largest_increase:
            largest_increase = country_max

    bar_top = y + (bottom_h - 0.3) / 2
    bar_ax = fig.add_axes(box(col_x[0], bar_top, pw, 0.08))
    if EXPOSURE_BINS:
        draw_exposure_class_bar(bar_ax, exp_cmap, norm, LEGEND_FONTSIZE)
    else:
        draw_exposure_colorbar(bar_ax, norm, largest_increase, LEGEND_FONTSIZE)

    # driver legend
    legend_left = col_x[1]
    legend_right = col_x[3] + pw
    legend_ax = fig.add_axes(box(legend_left, y, legend_right - legend_left, bottom_h))
    draw_legend(legend_ax, all_points, LEGEND_FONTSIZE)
    
    #----Save----
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / "country_maps_grouped.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(path.with_suffix(".pdf"), dpi=DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(path.with_suffix(".svg"), dpi=DPI, bbox_inches="tight", facecolor="white")
    print(f"Done -> {path}")
    plt.close(fig)

#----------------------------------------------------------------------------------

if __name__ == "__main__":
    plot_country_maps_grouped()