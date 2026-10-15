"""
Publication-quality grayscale maps of the 1949 San Juan and 1971 Teneguía
lava flows on La Palma, Canary Islands.

Generates THREE maps:
  1. 1949 San Juan close-up with hillshade DEM + numbered SRA stations
  2. 1971 Teneguía close-up with hillshade DEM + numbered SRA stations
  3. La Palma island overview showing both flow locations

Each close-up includes:
  - DEM hillshade background (grayscale)
  - Numbered SRA field sample stations
  - Flow-path line connecting stations in order
  - Vent location marker
  - UTM coordinate grid with tick labels
  - Metric scale bar
  - North arrow
  - CRS annotation

The overview map includes:
  - Both DEM hillshades mosaicked
  - Flow paths highlighted in red
  - Vent locations marked with red stars
  - Scale bar, north arrow, coordinate grid

Outputs 300 dpi PNG + PDF.

Author: Auto-generated for Joshua's La Palma methodology
"""

import os
from pathlib import Path
import csv
import numpy as np
import rasterio
from rasterio.windows import from_bounds
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle
from matplotlib.lines import Line2D
from matplotlib import ticker
from pyproj import Transformer

# ── Paths ────────────────────────────────────────────────────────────────────
BASE      = Path(__file__).resolve().parent
SRA_DIR   = os.path.join(BASE, "sra")
RES_DIR   = os.path.join(BASE, "Results Figures")
CSV_PATH  = os.path.join(SRA_DIR, "sheet2_out.csv")

DEM_1085  = os.path.join(BASE, "MDT02-REGCAN95-HU28-1085-3-4-COB2.tif")  # northern (1949)
DEM_1087  = os.path.join(BASE, "MDT02-REGCAN95-HU28-1087-1-2-COB2.tif")  # southern (1971)

# ── Known vent locations (approximate, UTM REGCAN95 Zone 28N) ────────────────
VENT_1949 = (219650, 3167450)
VENT_1971 = (220900, 3153250)

# ── UTM to WGS84 transformer ────────────────────────────────────────────────
transformer = Transformer.from_crs("EPSG:4083", "EPSG:4326", always_xy=True)


# ══════════════════════════════════════════════════════════════════════════════
#  HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def compute_hillshade(dem, azimuth=315, altitude=45, z_factor=1.0):
    """Compute a hillshade raster from a DEM array."""
    az_rad  = np.radians(azimuth)
    alt_rad = np.radians(altitude)
    dy, dx = np.gradient(dem * z_factor)
    slope  = np.arctan(np.sqrt(dx**2 + dy**2))
    aspect = np.arctan2(-dy, dx)
    hs = (np.sin(alt_rad) * np.cos(slope) +
          np.cos(alt_rad) * np.sin(slope) * np.cos(az_rad - aspect))
    return np.clip(hs, 0, 1)


def read_dem_window(tif_path, x_min, x_max, y_min, y_max, pad=200):
    """Read a windowed subset of a GeoTIFF and return (array, extent)."""
    with rasterio.open(tif_path) as ds:
        window = from_bounds(
            x_min - pad, y_min - pad, x_max + pad, y_max + pad,
            ds.transform
        )
        data = ds.read(1, window=window).astype(np.float32)
        win_transform = ds.window_transform(window)
        nrows, ncols = data.shape
        left   = win_transform.c
        top    = win_transform.f
        right  = left + ncols * win_transform.a
        bottom = top  + nrows * win_transform.e
        extent = (left, right, bottom, top)
        if ds.nodata is not None:
            data[data == ds.nodata] = np.nan
        data[data < -1000] = np.nan
    return data, extent


def read_dem_full(tif_path):
    """Read the full GeoTIFF and return (array, extent)."""
    with rasterio.open(tif_path) as ds:
        data = ds.read(1).astype(np.float32)
        b = ds.bounds
        extent = (b.left, b.right, b.bottom, b.top)
        if ds.nodata is not None:
            data[data == ds.nodata] = np.nan
        data[data < -1000] = np.nan
    return data, extent


def add_scale_bar(ax, extent, bar_length_m, y_frac=0.06, x_frac=0.05):
    """Add a neat alternating black/white scale bar."""
    x_min, x_max, y_min, y_max = extent
    width  = x_max - x_min
    height = y_max - y_min
    bar_x = x_min + width * x_frac
    bar_y = y_min + height * y_frac
    bar_h = height * 0.008

    n_segments = 4
    seg_len = bar_length_m / n_segments
    for i in range(n_segments):
        color = "k" if i % 2 == 0 else "w"
        rect = plt.Rectangle(
            (bar_x + i * seg_len, bar_y), seg_len, bar_h,
            facecolor=color, edgecolor="k", linewidth=0.6, zorder=20,
            clip_on=False
        )
        ax.add_patch(rect)

    # Tick labels
    if bar_length_m >= 1000:
        labels = [(0, bar_x),
                  (bar_length_m / 2000, bar_x + bar_length_m / 2),
                  (bar_length_m / 1000, bar_x + bar_length_m)]
        fmt = lambda v: f"{v:.0f} km" if v == bar_length_m / 1000 else f"{v:.0f}"
    else:
        labels = [(0, bar_x),
                  (bar_length_m / 2, bar_x + bar_length_m / 2),
                  (bar_length_m, bar_x + bar_length_m)]
        fmt = lambda v: f"{int(v)} m" if v == bar_length_m else f"{int(v)}"

    for val, x_pos in labels:
        ax.text(x_pos, bar_y - height * 0.012,
                fmt(val),
                ha="center", va="top", fontsize=7, fontweight="bold",
                color="k", zorder=20,
                path_effects=[pe.withStroke(linewidth=2, foreground="white")])


def add_north_arrow(ax, extent, x_frac=0.93, y_frac=0.92, size=0.06):
    """Add a north arrow."""
    x_min, x_max, y_min, y_max = extent
    width  = x_max - x_min
    height = y_max - y_min
    cx = x_min + width * x_frac
    cy = y_min + height * y_frac
    arrow_len = height * size

    ax.annotate("",
                xy=(cx, cy + arrow_len * 0.5),
                xytext=(cx, cy - arrow_len * 0.5),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=1.8,
                                mutation_scale=14),
                zorder=25)
    ax.text(cx, cy + arrow_len * 0.58, "N",
            ha="center", va="bottom", fontsize=10, fontweight="bold",
            color="k", zorder=25,
            path_effects=[pe.withStroke(linewidth=2.5, foreground="white")])


def add_utm_grid(ax, extent, max_ticks=8):
    """Add UTM coordinate ticks and grid to axes."""
    x_min, x_max, y_min, y_max = extent
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))

    x_range = x_max - x_min
    y_range = y_max - y_min

    for step in [100, 200, 250, 500, 1000, 2000, 5000]:
        if x_range / step <= max_ticks:
            x_step = step
            break
    else:
        x_step = 1000

    for step in [100, 200, 250, 500, 1000, 2000, 5000]:
        if y_range / step <= max_ticks:
            y_step = step
            break
    else:
        y_step = 1000

    x_start = int(np.ceil(x_min / x_step) * x_step)
    y_start = int(np.ceil(y_min / y_step) * y_step)
    ax.set_xticks(np.arange(x_start, x_max, x_step))
    ax.set_yticks(np.arange(y_start, y_max, y_step))

    ax.tick_params(axis="both", labelsize=7, direction="in", length=4,
                   width=0.6, color="k", pad=3)
    ax.set_xlabel("Easting (m)", fontsize=8, labelpad=6)
    ax.set_ylabel("Northing (m)", fontsize=8, labelpad=6)
    ax.grid(True, linewidth=0.3, alpha=0.4, color="k", linestyle="--")


# ══════════════════════════════════════════════════════════════════════════════
#  CLOSE-UP MAP BUILDER (for each individual flow)
# ══════════════════════════════════════════════════════════════════════════════

def make_flow_map(
    flow_name, short_name, dem_path, sites, vent_utm,
    scale_bar_m, pad_m=300,
):
    """Generate a single publication-quality close-up map figure."""
    print(f"\n{'='*60}")
    print(f"  Generating: {flow_name}")
    print(f"{'='*60}")

    eastings  = [s["easting"] for s in sites] + [vent_utm[0]]
    northings = [s["northing"] for s in sites] + [vent_utm[1]]
    x_min, x_max = min(eastings) - pad_m, max(eastings) + pad_m
    y_min, y_max = min(northings) - pad_m, max(northings) + pad_m

    x_range = x_max - x_min
    y_range = y_max - y_min
    if x_range > y_range * 1.6:
        extra = (x_range / 1.3 - y_range) / 2
        y_min -= extra; y_max += extra
    elif y_range > x_range * 1.6:
        extra = (y_range / 1.3 - x_range) / 2
        x_min -= extra; x_max += extra

    y_range = y_max - y_min
    x_range = x_max - x_min

    print(f"  Extent: E [{x_min:.0f}, {x_max:.0f}], N [{y_min:.0f}, {y_max:.0f}]")

    dem, extent = read_dem_window(dem_path, x_min, x_max, y_min, y_max, pad=100)
    hs = compute_hillshade(dem, azimuth=315, altitude=45)

    fig_w = 7.5
    fig_h = fig_w * (y_range / x_range)
    fig_h = max(min(fig_h, 11.0), 5.5)

    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h), dpi=300)
    fig.patch.set_facecolor("white")

    # Hillshade
    ax.imshow(hs, cmap="gray", vmin=0, vmax=1,
              extent=extent, origin="upper", aspect="equal", zorder=1)
    ax.imshow(dem, cmap="gray", alpha=0.15,
              extent=extent, origin="upper", aspect="equal", zorder=2)

    # Flow path (Bright yellow with solid black outline for maximum visibility)
    sorted_sites = sorted(sites, key=lambda s: s["site"])
    path_x = [s["easting"] for s in sorted_sites]
    path_y = [s["northing"] for s in sorted_sites]
    
    # Outer black stroke
    ax.plot(path_x, path_y, color="k", linewidth=4.5, zorder=7, solid_capstyle="round")
    # Inner bright yellow path
    ax.plot(path_x, path_y, color="#FFD700", linewidth=2.5, zorder=8, solid_capstyle="round")

    # Vent
    # Removed per user request

    # Numbered sample sites
    for s in sorted_sites:
        ex, ny = s["easting"], s["northing"]
        ax.plot(ex, ny, "o", markersize=8, markerfacecolor="white",
                markeredgecolor="k", markeredgewidth=1.2, zorder=12)
        ax.text(ex, ny, str(s["site"]),
                ha="center", va="center", fontsize=5.5, fontweight="bold",
                color="k", zorder=13)

        # Elevation label — alternate left/right to reduce overlaps
        if s["site"] % 2 == 0:
            offset_x, ha = x_range * 0.025, "left"
        else:
            offset_x, ha = -x_range * 0.025, "right"
        ax.text(ex + offset_x, ny + y_range * 0.01,
                f'{s["elevation"]:.0f} m',
                ha=ha, va="center", fontsize=5.5, color="k", zorder=14,
                fontstyle="italic",
                path_effects=[pe.withStroke(linewidth=2, foreground="white")])

    final_extent = (x_min, x_max, y_min, y_max)
    add_utm_grid(ax, final_extent)
    add_scale_bar(ax, final_extent, bar_length_m=scale_bar_m)
    add_north_arrow(ax, final_extent)

    ax.set_title(flow_name, fontsize=12, fontweight="bold", pad=12)
    ax.text(0.99, 0.01, "CRS: REGCAN95 / UTM zone 28N",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=5.5, color="k", alpha=0.6, zorder=20,
            path_effects=[pe.withStroke(linewidth=1.5, foreground="white")])

    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    ax.set_aspect("equal")
    plt.tight_layout(pad=1.5)

    # Save
    for d in [SRA_DIR, RES_DIR]:
        for ext in ["png", "pdf"]:
            path = os.path.join(d, f"map_{short_name}_pub.{ext}")
            fig.savefig(path, dpi=300, bbox_inches="tight",
                        facecolor="white", edgecolor="none")
            print(f"  ✓ {path}")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
#  ISLAND OVERVIEW MAP (Map 3)
# ══════════════════════════════════════════════════════════════════════════════

def make_overview_map(sites_1949, sites_1971):
    """
    Generate a La Palma island overview map showing both DEM hillshades
    and both flow paths marked in red.
    """
    print(f"\n{'='*60}")
    print(f"  Generating: La Palma Island Overview")
    print(f"{'='*60}")

    # Read full extents of both DEMs
    print("  Reading DEM 1085 (northern)...")
    dem_n, ext_n = read_dem_full(DEM_1085)
    print(f"    Shape: {dem_n.shape}, extent: E [{ext_n[0]:.0f},{ext_n[1]:.0f}], N [{ext_n[2]:.0f},{ext_n[3]:.0f}]")

    print("  Reading DEM 1087 (southern)...")
    dem_s, ext_s = read_dem_full(DEM_1087)
    print(f"    Shape: {dem_s.shape}, extent: E [{ext_s[0]:.0f},{ext_s[1]:.0f}], N [{ext_s[2]:.0f},{ext_s[3]:.0f}]")

    # Compute hillshades
    print("  Computing hillshades...")
    hs_n = compute_hillshade(dem_n, azimuth=315, altitude=45)
    hs_s = compute_hillshade(dem_s, azimuth=315, altitude=45)

    # Combined extent
    x_min = min(ext_n[0], ext_s[0])
    x_max = max(ext_n[1], ext_s[1])
    y_min = min(ext_n[2], ext_s[2])
    y_max = max(ext_n[3], ext_s[3])

    # Add some padding
    pad = 500
    x_min -= pad; x_max += pad
    y_min -= pad; y_max += pad

    x_range = x_max - x_min
    y_range = y_max - y_min

    print(f"  Combined extent: E [{x_min:.0f},{x_max:.0f}], N [{y_min:.0f},{y_max:.0f}]")
    print(f"  Range: {x_range:.0f} × {y_range:.0f} m")

    fig_w = 7.5
    fig_h = fig_w * (y_range / x_range)
    fig_h = max(min(fig_h, 14.0), 6.0)

    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h), dpi=300)
    fig.patch.set_facecolor("white")

    # Set background to light gray for ocean
    ax.set_facecolor("#e8e8e8")

    # Plot both hillshades
    ax.imshow(hs_n, cmap="gray", vmin=0, vmax=1,
              extent=ext_n, origin="upper", aspect="equal", zorder=1)
    ax.imshow(hs_s, cmap="gray", vmin=0, vmax=1,
              extent=ext_s, origin="upper", aspect="equal", zorder=1)

    # Elevation tinting
    ax.imshow(dem_n, cmap="gray", alpha=0.15,
              extent=ext_n, origin="upper", aspect="equal", zorder=2)
    ax.imshow(dem_s, cmap="gray", alpha=0.15,
              extent=ext_s, origin="upper", aspect="equal", zorder=2)

    # ── 1949 San Juan flow path — BRIGHT RED ────────────────────────────────
    sorted_1949 = sorted(sites_1949, key=lambda s: s["site"])
    path_x_49 = [s["easting"] for s in sorted_1949]
    path_y_49 = [s["northing"] for s in sorted_1949]

    # Thick white outline for contrast
    ax.plot(path_x_49, path_y_49, color="white", linewidth=4.5, zorder=7, solid_capstyle="round")
    # Bright red inner line
    ax.plot(path_x_49, path_y_49, color="#FF1493", linewidth=2.5, zorder=8, solid_capstyle="round")  # DeepPink/Neon Red

    # Vent star
    # Removed per user request

    # Label
    # Removed per user request

    # ── 1971 Teneguía flow path — BRIGHT RED ────────────────────────────────
    sorted_1971 = sorted(sites_1971, key=lambda s: s["site"])
    path_x_71 = [s["easting"] for s in sorted_1971]
    path_y_71 = [s["northing"] for s in sorted_1971]

    # Thick white outline for contrast
    ax.plot(path_x_71, path_y_71, color="white", linewidth=4.5, zorder=7, solid_capstyle="round")
    # Bright red inner line
    ax.plot(path_x_71, path_y_71, color="#FF1493", linewidth=2.5, zorder=8, solid_capstyle="round")  # DeepPink/Neon Red

    # Vent star
    # Removed per user request

    # Label
    # Removed per user request

    # ── DEM tile boundary outlines (dashed black) ────────────────────
    # Removed per user request

    # ── Legend ───────────────────────────────────────────────────────
    legend_elements = [
        Line2D([0], [0], color="#FF1493", linewidth=2.5, 
               path_effects=[pe.withStroke(linewidth=4.5, foreground="white")],
               label="Lava flow path"),
    ]
    leg = ax.legend(handles=legend_elements, loc="lower left",
                    fontsize=7, framealpha=0.9, edgecolor="k",
                    frameon=True, fancybox=False)
    leg.get_frame().set_linewidth(0.6)
    leg.set_zorder(30)

    # ── Cartographic elements ────────────────────────────────────────
    final_extent = (x_min, x_max, y_min, y_max)
    add_utm_grid(ax, final_extent, max_ticks=6)
    add_scale_bar(ax, final_extent, bar_length_m=2000, y_frac=0.04, x_frac=0.65)
    add_north_arrow(ax, final_extent, x_frac=0.93, y_frac=0.93)

    ax.set_title("La Palma — Study Area Overview", fontsize=13,
                 fontweight="bold", pad=12)
    ax.text(0.99, 0.01, "CRS: REGCAN95 / UTM zone 28N",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=5.5, color="k", alpha=0.6, zorder=20,
            path_effects=[pe.withStroke(linewidth=1.5, foreground="white")])

    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    ax.set_aspect("equal")
    plt.tight_layout(pad=1.5)

    # Save
    for d in [SRA_DIR, RES_DIR]:
        for ext in ["png", "pdf"]:
            path = os.path.join(d, f"map_lapalma_overview_pub.{ext}")
            fig.savefig(path, dpi=300, bbox_inches="tight",
                        facecolor="white", edgecolor="none")
            print(f"  ✓ {path}")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
#  LOAD DATA & RUN
# ══════════════════════════════════════════════════════════════════════════════

print("Loading field sample data...")
sites_1949 = []
sites_1971 = []

with open(CSV_PATH, "r", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    for row in reader:
        eruption = row["Eruption"].strip()
        entry = {
            "site":     int(row["Site"]),
            "easting":  float(row["Easting (X)"]),
            "northing": float(row["Northing (Y)"]),
            "elevation": float(row["Elevation (m)"]),
        }
        if "1949" in eruption:
            sites_1949.append(entry)
        elif "1971" in eruption:
            sites_1971.append(entry)

print(f"  1949 San Juan:  {len(sites_1949)} sites")
print(f"  1971 Teneguía: {len(sites_1971)} sites")

# ── Map 1: 1949 San Juan close-up ────────────────────────────────────────────
make_flow_map(
    flow_name="1949 San Juan (Llano del Banco) Lava Flow — SRA Stations",
    short_name="1949_sanjuan",
    dem_path=DEM_1085,
    sites=sites_1949,
    vent_utm=VENT_1949,
    scale_bar_m=200,
    pad_m=200,
)

# ── Map 2: 1971 Teneguía close-up ────────────────────────────────────────────
make_flow_map(
    flow_name="1971 Teneguía Lava Flow — SRA Stations",
    short_name="1971_teneguia",
    dem_path=DEM_1087,
    sites=sites_1971,
    vent_utm=VENT_1971,
    scale_bar_m=500,
    pad_m=350,
)

# ── Map 3: Island overview ──────────────────────────────────────────────────
make_overview_map(sites_1949, sites_1971)

print(f"\n{'='*60}")
print("  All three maps generated successfully!")
print(f"{'='*60}")
