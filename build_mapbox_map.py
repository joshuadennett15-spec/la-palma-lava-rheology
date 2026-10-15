"""
Build Mapbox GL JS interactive map for La Palma methodology figure.

Displays:
  - Two REGCAN95 DEM tile coverage polygons (HU28-1087-1-2, HU28-1085-3-4)
  - 1949 San Juan lava flow field sample sites (10 stations)
  - 1971 Teneguía lava flow field sample sites (20 stations)

Reads GeoTIFF bounds via rasterio, converts all UTM coordinates
(REGCAN95 / UTM Zone 28N) to WGS84 lat/lon via pyproj, and outputs
a self-contained HTML file with Mapbox GL JS.

Author: Auto-generated for Joshua's La Palma methodology
"""

import os
from pathlib import Path
import json
import csv
import rasterio
from pyproj import Transformer

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent
DEM_1087 = os.path.join(BASE, "MDT02-REGCAN95-HU28-1087-1-2-COB2.tif")
DEM_1085 = os.path.join(BASE, "MDT02-REGCAN95-HU28-1085-3-4-COB2.tif")
CSV_PATH = os.path.join(BASE, "sra", "sheet2_out.csv")
OUT_HTML = os.path.join(BASE, "la_palma_methodology_map.html")

MAPBOX_TOKEN = os.environ.get("MAPBOX_TOKEN")
if not MAPBOX_TOKEN:
    raise SystemExit(
        "MAPBOX_TOKEN is not set. Create a public token at https://account.mapbox.com/ "
        "and export it (export MAPBOX_TOKEN=pk...) before running this script."
    )


# ── Coordinate transformer (REGCAN95 / UTM 28N  →  WGS 84) ──────────────────
# EPSG:4083 = REGCAN95 / UTM zone 28N
# EPSG:4326 = WGS 84  (lon, lat)
transformer = Transformer.from_crs("EPSG:4083", "EPSG:4326", always_xy=True)


def utm_to_lonlat(easting, northing):
    """Convert REGCAN95 UTM 28N coords to WGS84 [lon, lat]."""
    lon, lat = transformer.transform(easting, northing)
    return [round(lon, 6), round(lat, 6)]


# ── Read DEM bounds ──────────────────────────────────────────────────────────
def get_dem_bounds_lonlat(tif_path):
    """Return [[lon,lat], ...] polygon corners for a GeoTIFF in WGS84."""
    with rasterio.open(tif_path) as ds:
        b = ds.bounds  # left, bottom, right, top in native CRS
        crs = ds.crs
        print(f"  {os.path.basename(tif_path)}: CRS={crs}, bounds={b}")
        
        # If it's already in a projected CRS close to REGCAN95 UTM28N,
        # use the transformer. Otherwise try auto-detecting.
        t = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
        
        # Corners: BL, BR, TR, TL, BL (closed ring)
        corners_xy = [
            (b.left,  b.bottom),
            (b.right, b.bottom),
            (b.right, b.top),
            (b.left,  b.top),
            (b.left,  b.bottom),  # close ring
        ]
        corners_lonlat = []
        for x, y in corners_xy:
            lon, lat = t.transform(x, y)
            corners_lonlat.append([round(lon, 6), round(lat, 6)])
        return corners_lonlat


print("Reading DEM bounds...")
dem_1087_poly = get_dem_bounds_lonlat(DEM_1087)
dem_1085_poly = get_dem_bounds_lonlat(DEM_1085)
print()

# ── Read field sample coordinates ────────────────────────────────────────────
print("Reading field sample data...")
flow_1949 = []
flow_1971 = []

with open(CSV_PATH, "r", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    for row in reader:
        eruption = row["Eruption"].strip()
        site = row["Site"].strip()
        easting = float(row["Easting (X)"])
        northing = float(row["Northing (Y)"])
        elev = float(row["Elevation (m)"])
        lonlat = utm_to_lonlat(easting, northing)

        entry = {
            "site": int(site),
            "lonlat": lonlat,
            "elevation": elev,
            "easting": easting,
            "northing": northing,
        }
        if "1949" in eruption:
            flow_1949.append(entry)
        elif "1971" in eruption:
            flow_1971.append(entry)

print(f"  1949 San Juan: {len(flow_1949)} sites")
print(f"  1971 Teneguía: {len(flow_1971)} sites")

# ── Build GeoJSON features ──────────────────────────────────────────────────
def make_point_features(sites, eruption_label, color):
    features = []
    for s in sites:
        features.append({
            "type": "Feature",
            "properties": {
                "eruption": eruption_label,
                "site": s["site"],
                "elevation": s["elevation"],
                "color": color,
                "label": f"Site {s['site']} ({s['elevation']:.0f} m)",
            },
            "geometry": {
                "type": "Point",
                "coordinates": s["lonlat"],
            }
        })
    return features


def make_line_feature(sites, eruption_label, color):
    """Create a line connecting all sites in order."""
    coords = [s["lonlat"] for s in sorted(sites, key=lambda x: x["site"])]
    return {
        "type": "Feature",
        "properties": {"eruption": eruption_label, "color": color},
        "geometry": {"type": "LineString", "coordinates": coords}
    }


def make_polygon_feature(corners, label, color):
    return {
        "type": "Feature",
        "properties": {"label": label, "color": color},
        "geometry": {"type": "Polygon", "coordinates": [corners]}
    }


# Points
features_1949_pts = make_point_features(flow_1949, "1949 San Juan", "#e67e22")
features_1971_pts = make_point_features(flow_1971, "1971 Teneguía", "#3498db")

# Lines connecting sites
line_1949 = make_line_feature(flow_1949, "1949 San Juan", "#e67e22")
line_1971 = make_line_feature(flow_1971, "1971 Teneguía", "#3498db")

# DEM polygons
dem_poly_1087 = make_polygon_feature(dem_1087_poly, "DEM HU28-1087-1-2 (2 m)", "#2ecc71")
dem_poly_1085 = make_polygon_feature(dem_1085_poly, "DEM HU28-1085-3-4 (2 m)", "#9b59b6")

# Combine into GeoJSON collections
geojson_points = {
    "type": "FeatureCollection",
    "features": features_1949_pts + features_1971_pts
}

geojson_lines = {
    "type": "FeatureCollection",
    "features": [line_1949, line_1971]
}

geojson_dems = {
    "type": "FeatureCollection",
    "features": [dem_poly_1087, dem_poly_1085]
}

# ── Compute map center ───────────────────────────────────────────────────────
all_lons = [s["lonlat"][0] for s in flow_1949 + flow_1971]
all_lats = [s["lonlat"][1] for s in flow_1949 + flow_1971]
center_lon = (min(all_lons) + max(all_lons)) / 2
center_lat = (min(all_lats) + max(all_lats)) / 2

print(f"\nMap center: [{center_lon:.6f}, {center_lat:.6f}]")
print(f"1949 flow lon range: [{min(l[0] for l in [s['lonlat'] for s in flow_1949]):.6f}, {max(l[0] for l in [s['lonlat'] for s in flow_1949]):.6f}]")
print(f"1971 flow lon range: [{min(l[0] for l in [s['lonlat'] for s in flow_1971]):.6f}, {max(l[0] for l in [s['lonlat'] for s in flow_1971]):.6f}]")

# ── Generate HTML ────────────────────────────────────────────────────────────
html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>La Palma — DEM Coverage &amp; Lava Flow Sample Sites</title>
    <script src="https://api.mapbox.com/mapbox-gl-js/v3.3.0/mapbox-gl.js"></script>
    <link href="https://api.mapbox.com/mapbox-gl-js/v3.3.0/mapbox-gl.css" rel="stylesheet" />
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{ font-family: 'Inter', 'Segoe UI', system-ui, sans-serif; }}
        #map {{ position: absolute; top: 0; bottom: 0; width: 100%; }}

        /* ── Floating legend panel ───────────────────────────────── */
        .legend {{
            position: absolute;
            bottom: 40px;
            left: 16px;
            background: rgba(15, 15, 25, 0.92);
            backdrop-filter: blur(12px);
            -webkit-backdrop-filter: blur(12px);
            border: 1px solid rgba(255,255,255,0.08);
            border-radius: 14px;
            padding: 18px 22px;
            color: #e8e8ec;
            font-size: 13px;
            line-height: 1.6;
            z-index: 10;
            min-width: 260px;
            box-shadow: 0 8px 32px rgba(0,0,0,0.45);
        }}
        .legend h3 {{
            margin: 0 0 12px;
            font-size: 15px;
            font-weight: 700;
            letter-spacing: 0.3px;
            color: #fff;
            border-bottom: 1px solid rgba(255,255,255,0.1);
            padding-bottom: 8px;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 10px;
            margin-bottom: 6px;
            cursor: pointer;
            padding: 3px 4px;
            border-radius: 6px;
            transition: background 0.2s;
        }}
        .legend-item:hover {{
            background: rgba(255,255,255,0.06);
        }}
        .legend-swatch {{
            width: 16px;
            height: 16px;
            border-radius: 50%;
            flex-shrink: 0;
            border: 2px solid rgba(255,255,255,0.15);
        }}
        .legend-swatch.rect {{
            border-radius: 3px;
            width: 18px;
            height: 12px;
            opacity: 0.65;
        }}
        .legend-label {{
            flex: 1;
        }}
        .legend-section {{
            margin-top: 10px;
            font-size: 11px;
            text-transform: uppercase;
            letter-spacing: 1px;
            color: rgba(255,255,255,0.4);
            margin-bottom: 4px;
        }}

        /* ── Title bar ───────────────────────────────────────────── */
        .title-bar {{
            position: absolute;
            top: 16px;
            left: 50%;
            transform: translateX(-50%);
            background: rgba(15, 15, 25, 0.92);
            backdrop-filter: blur(12px);
            -webkit-backdrop-filter: blur(12px);
            border: 1px solid rgba(255,255,255,0.08);
            border-radius: 12px;
            padding: 12px 28px;
            color: #fff;
            font-size: 16px;
            font-weight: 600;
            z-index: 10;
            text-align: center;
            box-shadow: 0 4px 24px rgba(0,0,0,0.35);
            letter-spacing: 0.3px;
        }}
        .title-bar small {{
            display: block;
            font-size: 11px;
            font-weight: 400;
            color: rgba(255,255,255,0.5);
            margin-top: 3px;
        }}

        /* ── Popup styling ───────────────────────────────────────── */
        .mapboxgl-popup-content {{
            background: rgba(20, 20, 35, 0.95) !important;
            color: #e8e8ec !important;
            border-radius: 10px !important;
            padding: 14px 18px !important;
            font-family: 'Inter', 'Segoe UI', system-ui, sans-serif !important;
            font-size: 13px !important;
            box-shadow: 0 6px 24px rgba(0,0,0,0.5) !important;
            border: 1px solid rgba(255,255,255,0.08) !important;
        }}
        .mapboxgl-popup-anchor-bottom .mapboxgl-popup-tip {{
            border-top-color: rgba(20, 20, 35, 0.95) !important;
        }}
        .mapboxgl-popup-anchor-top .mapboxgl-popup-tip {{
            border-bottom-color: rgba(20, 20, 35, 0.95) !important;
        }}
        .mapboxgl-popup-close-button {{
            color: rgba(255,255,255,0.5) !important;
            font-size: 18px !important;
        }}
        .popup-title {{
            font-weight: 700;
            font-size: 14px;
            margin-bottom: 6px;
        }}
        .popup-row {{
            display: flex;
            justify-content: space-between;
            gap: 16px;
            padding: 2px 0;
        }}
        .popup-key {{
            color: rgba(255,255,255,0.5);
        }}
    </style>
</head>
<body>

<div id="map"></div>

<div class="title-bar">
    La Palma — DEM Coverage &amp; Lava Flow Sample Sites
    <small>REGCAN95 2 m DSM tiles &bull; SRA field stations</small>
</div>

<div class="legend">
    <h3>Map Legend</h3>

    <div class="legend-section">Lava Flows</div>
    <div class="legend-item" onclick="toggleLayer('flow-points-1949','flow-line-1949')">
        <span class="legend-swatch" style="background:#e67e22;"></span>
        <span class="legend-label">1949 San Juan (10 sites)</span>
    </div>
    <div class="legend-item" onclick="toggleLayer('flow-points-1971','flow-line-1971')">
        <span class="legend-swatch" style="background:#3498db;"></span>
        <span class="legend-label">1971 Teneguía (20 sites)</span>
    </div>

    <div class="legend-section">DEM Tiles (2 m)</div>
    <div class="legend-item" onclick="toggleLayer('dem-fill','dem-outline')">
        <span class="legend-swatch rect" style="background:#2ecc71;"></span>
        <span class="legend-label">HU28-1087-1-2</span>
    </div>
    <div class="legend-item" onclick="toggleLayer('dem-fill','dem-outline')">
        <span class="legend-swatch rect" style="background:#9b59b6;"></span>
        <span class="legend-label">HU28-1085-3-4</span>
    </div>
</div>

<script>
mapboxgl.accessToken = '{MAPBOX_TOKEN}';

const map = new mapboxgl.Map({{
    container: 'map',
    style: 'mapbox://styles/mapbox/satellite-streets-v12',
    center: [{center_lon}, {center_lat}],
    zoom: 11.5,
    pitch: 0,
    bearing: 0,
}});

map.addControl(new mapboxgl.NavigationControl(), 'top-right');
map.addControl(new mapboxgl.ScaleControl({{ unit: 'metric' }}), 'bottom-right');

// ── GeoJSON data (embedded) ─────────────────────────────────────────────
const pointsData = {json.dumps(geojson_points)};
const linesData  = {json.dumps(geojson_lines)};
const demsData   = {json.dumps(geojson_dems)};

map.on('load', () => {{

    // ── Terrain (3D hillshade) ──────────────────────────────────────
    map.addSource('mapbox-dem', {{
        type: 'raster-dem',
        url: 'mapbox://mapbox.mapbox-terrain-dem-v1',
        tileSize: 512,
        maxzoom: 14
    }});
    map.setTerrain({{ source: 'mapbox-dem', exaggeration: 1.3 }});

    // ── DEM tile polygons ───────────────────────────────────────────
    map.addSource('dems', {{ type: 'geojson', data: demsData }});
    map.addLayer({{
        id: 'dem-fill',
        type: 'fill',
        source: 'dems',
        paint: {{
            'fill-color': ['get', 'color'],
            'fill-opacity': 0.12,
        }}
    }});
    map.addLayer({{
        id: 'dem-outline',
        type: 'line',
        source: 'dems',
        paint: {{
            'line-color': ['get', 'color'],
            'line-width': 2.5,
            'line-dasharray': [4, 3],
            'line-opacity': 0.75,
        }}
    }});

    // ── Flow path lines ─────────────────────────────────────────────
    map.addSource('flow-lines', {{ type: 'geojson', data: linesData }});
    map.addLayer({{
        id: 'flow-line-1949',
        type: 'line',
        source: 'flow-lines',
        filter: ['==', ['get', 'eruption'], '1949 San Juan'],
        paint: {{
            'line-color': '#e67e22',
            'line-width': 3,
            'line-opacity': 0.8,
        }}
    }});
    map.addLayer({{
        id: 'flow-line-1971',
        type: 'line',
        source: 'flow-lines',
        filter: ['==', ['get', 'eruption'], '1971 Teneguía'],
        paint: {{
            'line-color': '#3498db',
            'line-width': 3,
            'line-opacity': 0.8,
        }}
    }});

    // ── Sample site points ──────────────────────────────────────────
    map.addSource('flow-points', {{ type: 'geojson', data: pointsData }});

    // Outer glow
    map.addLayer({{
        id: 'flow-points-glow',
        type: 'circle',
        source: 'flow-points',
        paint: {{
            'circle-radius': 10,
            'circle-color': ['get', 'color'],
            'circle-opacity': 0.2,
            'circle-blur': 0.6,
        }}
    }});

    // Main circle
    map.addLayer({{
        id: 'flow-points-1949',
        type: 'circle',
        source: 'flow-points',
        filter: ['==', ['get', 'eruption'], '1949 San Juan'],
        paint: {{
            'circle-radius': 6,
            'circle-color': '#e67e22',
            'circle-stroke-color': '#fff',
            'circle-stroke-width': 1.5,
            'circle-opacity': 0.95,
        }}
    }});
    map.addLayer({{
        id: 'flow-points-1971',
        type: 'circle',
        source: 'flow-points',
        filter: ['==', ['get', 'eruption'], '1971 Teneguía'],
        paint: {{
            'circle-radius': 6,
            'circle-color': '#3498db',
            'circle-stroke-color': '#fff',
            'circle-stroke-width': 1.5,
            'circle-opacity': 0.95,
        }}
    }});

    // Site number labels
    map.addLayer({{
        id: 'flow-labels',
        type: 'symbol',
        source: 'flow-points',
        layout: {{
            'text-field': ['to-string', ['get', 'site']],
            'text-size': 11,
            'text-offset': [0, -1.4],
            'text-font': ['DIN Pro Medium', 'Arial Unicode MS Bold'],
            'text-allow-overlap': false,
        }},
        paint: {{
            'text-color': '#fff',
            'text-halo-color': 'rgba(0,0,0,0.7)',
            'text-halo-width': 1.5,
        }}
    }});

    // DEM label annotations
    map.addLayer({{
        id: 'dem-labels',
        type: 'symbol',
        source: 'dems',
        layout: {{
            'text-field': ['get', 'label'],
            'text-size': 12,
            'text-font': ['DIN Pro Medium', 'Arial Unicode MS Bold'],
        }},
        paint: {{
            'text-color': ['get', 'color'],
            'text-halo-color': 'rgba(0,0,0,0.8)',
            'text-halo-width': 1.5,
        }}
    }});

    // ── Popups ──────────────────────────────────────────────────────
    const popup = new mapboxgl.Popup({{ closeButton: true, closeOnClick: true, maxWidth: '280px' }});

    ['flow-points-1949', 'flow-points-1971'].forEach(layerId => {{
        map.on('click', layerId, (e) => {{
            const p = e.features[0].properties;
            popup.setLngLat(e.lngLat)
                .setHTML(`
                    <div class="popup-title" style="color:${{p.color}}">${{p.eruption}} — Site ${{p.site}}</div>
                    <div class="popup-row"><span class="popup-key">Elevation</span><span>${{p.elevation}} m</span></div>
                    <div class="popup-row"><span class="popup-key">Coordinates</span><span>${{e.lngLat.lat.toFixed(5)}}°N, ${{Math.abs(e.lngLat.lng).toFixed(5)}}°W</span></div>
                `)
                .addTo(map);
        }});
        map.on('mouseenter', layerId, () => {{ map.getCanvas().style.cursor = 'pointer'; }});
        map.on('mouseleave', layerId, () => {{ map.getCanvas().style.cursor = ''; }});
    }});
}});

// ── Layer toggle helper ─────────────────────────────────────────────────
function toggleLayer(...ids) {{
    ids.forEach(id => {{
        const v = map.getLayoutProperty(id, 'visibility');
        map.setLayoutProperty(id, 'visibility', v === 'none' ? 'visible' : 'none');
    }});
}}
</script>
</body>
</html>"""

# ── Write HTML ───────────────────────────────────────────────────────────────
with open(OUT_HTML, "w", encoding="utf-8") as f:
    f.write(html)

print(f"\n✅ Map written to: {OUT_HTML}")
print("Open this file in a browser to view the interactive map.")
