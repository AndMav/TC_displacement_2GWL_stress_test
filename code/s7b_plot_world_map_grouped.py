"""
Creates Fig. 4: the dominant driver of the EAD change at every exposure point of all countries, at (A) p10,
(B) median and (C) p90 (full colour: increase, light shade: decrease), on Eckert IV world maps.
Every grid cell (300 arcsec) is drawn as a tile of its real size.

Input:  world_points.csv                                  (dominant driver at every exposure point)
Output: world_maps_p10_median_p90.png, .pdf, .svg
        (folders defined in config.py)
"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.legend as mlegend
import matplotlib.path as mpath
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from config import WORLD_DATA_DIR, FIGURES_DIR

'''
 ----------------------------------------------------------------------------------------
     I. Setup
 ----------------------------------------------------------------------------------------
'''

OUT_DIR = FIGURES_DIR / "world_map"
DPI = 600                                                      # resolution of the saved PNG

FIGURE_STATS = ['p10', 'median', 'p90']                        # maps from top to bottom

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

# Eckert IV world maps at the final print size: 4.6 in maps -> figure ~5.0 in wide, ~7.8 in deep
MAP_WIDTH_IN = 4.6
LAT_RANGE = None                                               # whole world, rounded outline; (south, north): crop

GRID_STEP_DEG = 300 / 3600                                     # grid step of the exposures (300 arcsec): size of one tile
BORDER_LINE_WIDTH_PT = 0.3                                     # country borders and coastlines

# Layout (inches) and font sizes (pt)
LAYOUT_IN = dict(left=0.2, right=0.02, top=0.02, map_gap=0.06, legend_gap=0.08, bottom=0.02)
LEGEND_HEIGHT_IN = 0.5
STAT_FONTSIZE = 9                                              # panel letter and statistic
LEGEND_FONTSIZE = 6.5
LEGEND_MARKER_SIZE_PT = 5

'''
 ----------------------------------------------------------------------------------------
     II. Drawing of one map and of the legend
 ----------------------------------------------------------------------------------------
'''

def set_font():
    """Universal sans-serif font (Arial, else Helvetica / Liberation Sans), embedded in the pdf, text kept as text in the svg."""
    installed = {f.name for f in font_manager.fontManager.ttflist}
    font = next((f for f in ["Arial", "Helvetica", "Liberation Sans"] if f in installed), "DejaVu Sans")
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": [font], "pdf.fonttype": 42, "svg.fonttype": "none",
                         "mathtext.fontset": "custom", "mathtext.rm": font, "mathtext.bf": f"{font}:bold"})   # also in math text


def colour(driver, sign):
    """ Provides color shade (decrease, increase) for the given driver"""
    
    base = np.array(mcolors.to_rgb(TERM_COLORS[driver]))
    if sign >= 0:
        return mcolors.to_hex(base)
    shaded = mcolors.to_hex((1 - TWO_COLOUR_TINT) * base + TWO_COLOUR_TINT)
    
    return shaded


def map_height(width_in):
    """Calculation of figure height in inches for a given map width"""
    
    proj = ccrs.EckertIV()                                                 
    x_east = proj.transform_point(180, 0, ccrs.PlateCarree())[0]
    x_west = proj.transform_point(-180, 0, ccrs.PlateCarree())[0]
    
    south, north = LAT_RANGE if LAT_RANGE else (-90, 90)
    
    y_north = proj.transform_point(0, north, ccrs.PlateCarree())[1]
    y_south = proj.transform_point(0, south, ccrs.PlateCarree())[1]

    height = width_in * (y_north - y_south) / (x_east - x_west)
    
    return height


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


def draw_world_map(ax, points):
    """World Map creation. Can be cropped to LAT_RANGE with a rectangular frame if set"""

    # land and ocean
    ax.add_feature(cfeature.LAND, facecolor="#F8F9FA", edgecolor="none", zorder=0)
    ax.add_feature(cfeature.OCEAN, facecolor="#EAF3FB", edgecolor="none", zorder=0)

    # country borders and coastlines, above the tiles
    ax.add_feature(cfeature.BORDERS, linewidth=BORDER_LINE_WIDTH_PT, edgecolor="#374151", zorder=20)
    ax.add_feature(cfeature.COASTLINE, linewidth=BORDER_LINE_WIDTH_PT, edgecolor="#374151", zorder=20)

    ax.set_global()
    ax.spines["geo"].set_linewidth(0.5)                       # thin grey outline of the map
    ax.spines["geo"].set_edgecolor("#374151")

    # crop to LAT_RANGE, with a rectangular frame
    if LAT_RANGE:                                              
        south, north = LAT_RANGE
        ax.set_ylim(ax.projection.transform_point(0, south, ccrs.PlateCarree())[1],
                    ax.projection.transform_point(0, north, ccrs.PlateCarree())[1])
        ax.set_boundary(mpath.Path([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]]), transform=ax.transAxes)
        ax.set_facecolor("#EAF3FB")  

    # Tiles: every driver and sign is one class (0 ... 13) with its own colour, one tile layer for all cells

    classes = []
    for driver in SPATIAL_TERMS:                     
        for sign in (1, -1):
            classes.append((driver, sign))

    code = {}
    for number in range(len(classes)):                 
        code[classes[number]] = number

    # class of every point
    drivers = points["dominant_driver"].values
    signs = points["sign"].values
    values = []
    
    for i in range(len(drivers)):                     
        values.append(code[(drivers[i], signs[i])])   

    x, y, grid = to_tiles(points["longitude"], points["latitude"], values)              

    # colour of every class
    class_colours = []
    for driver, sign in classes:
        class_colours.append(colour(driver, sign))
    
    cmap = mcolors.ListedColormap(class_colours)

    boundaries = []
    for number in range(len(classes) + 1):
        boundaries.append(number)
    
    norm = mcolors.BoundaryNorm(boundaries, len(classes))
    
    ax.pcolormesh(x, y, grid, cmap=cmap, norm=norm, alpha=0.90, shading="nearest", edgecolors="none", antialiased=False,
                  zorder=4, transform=ccrs.PlateCarree(), rasterized=True)


def draw_legend(ax, all_points, fontsize):
    """Framed legend filling the axes: one row per sign that occurs (Increase above, Decrease below), separated by a
    line. Every driver has its own column, the same in both rows, spread over the whole width; a driver that does
    not occur with a sign leaves its place empty."""

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
 
    # separator between Increase (above) and Decrease (below)
    if len(rows) == 2:
        ax.axhline(0.5, xmin=0.02, xmax=0.98, color="#9CA3AF", linewidth=0.5)
 
    for r, (sign, header, drivers) in enumerate(rows):
 
        # one legend row per sign, the same 7 columns in both rows
        handles = []
        labels = []
        for d in SPATIAL_TERMS:
            dot = Line2D([0], [0], marker="s", linestyle="none", markersize=LEGEND_MARKER_SIZE_PT, markeredgewidth=0, color=colour(d, sign))
            handles.append(dot)
            labels.append(SPATIAL_LABELS[SPATIAL_TERMS.index(d)])
 
        leg = mlegend.Legend(ax, handles, labels, title=header, handlelength=1.0, ncol=len(SPATIAL_TERMS), loc="center",
                             mode="expand", bbox_to_anchor=(0.01, 1 - (r + 1) / len(rows), 0.98, 1 / len(rows)),
                             fontsize=fontsize, title_fontsize=fontsize, frameon=False, handletextpad=0.3, borderpad=0.2, columnspacing=0.5)
 
        # drivers that do not occur with this sign: hidden, their place stays empty
        for d, handle, text in zip(SPATIAL_TERMS, leg.legend_handles, leg.get_texts()):
            if d not in drivers:
                handle.set_visible(False)
                text.set_visible(False)
 
        ax.add_artist(leg)


'''
 ----------------------------------------------------------------------------------------
     III. Plot creation
 ----------------------------------------------------------------------------------------
'''
 
def plot_world_maps_grouped():

    set_font()
    world = pd.read_csv(WORLD_DATA_DIR / "world_points.csv")
    world["sign"] = np.where(world["dominant_value"] >= 0, 1, -1)
 
    #----Layout----
    L = LAYOUT_IN
    map_h = map_height(MAP_WIDTH_IN)
    n = len(FIGURE_STATS)
    FW = L["left"] + MAP_WIDTH_IN + L["right"]
    FH = (L["top"] + n * map_h + (n - 1) * L["map_gap"]+ L["legend_gap"] + LEGEND_HEIGHT_IN + L["bottom"])
 
    def box(left, top, width, height):

        """Rectangle given in inches from the top-left corner -> [left, bottom, width, height] as figure fractions."""

        left_frac = left / FW
        bottom_frac = 1 - (top + height) / FH
        width_frac = width / FW
        height_frac = height / FH

        return [left_frac, bottom_frac, width_frac, height_frac]
 
    fig = plt.figure(figsize=(FW, FH))
 
    # ----Maps Creation----
    y = L["top"]
    
    for i, stat in enumerate(FIGURE_STATS):

        ax = fig.add_axes(box(L["left"], y, MAP_WIDTH_IN, map_h), projection=ccrs.EckertIV())
        draw_world_map(ax, world[world["stat"] == stat])

        # panel letter and statistic, left of the map
        ax.text(-0.01, 0.5, f"({'ABC'[i]}) {stat}", transform=ax.transAxes, rotation=90, ha="right", va="center",
                fontsize=STAT_FONTSIZE, fontweight="bold")

        y += map_h + L["map_gap"]


    # ----Legend Creation---- below the maps, same left edge and width
    legend_top = FH - L["bottom"] - LEGEND_HEIGHT_IN
    legend_ax = fig.add_axes(box(L["left"], legend_top, MAP_WIDTH_IN, LEGEND_HEIGHT_IN))
    draw_legend(legend_ax, world, LEGEND_FONTSIZE)
 
    stem = "world_maps_" + "_".join(FIGURE_STATS)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_DIR / f"{stem}.pdf", dpi=DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_DIR / f"{stem}.svg", dpi=DPI, bbox_inches="tight", facecolor="white")
    print(f"Done: {OUT_DIR / stem}.png")
    plt.close(fig)

#----------------------------------------------------------------------------------

if __name__ == "__main__":
    plot_world_maps_grouped()