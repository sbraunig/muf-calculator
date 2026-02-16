#!/usr/bin/env python3
"""
MUF Calculator Dashboard — Stage 2 Streamlit App
Run with:  streamlit run muf_dashboard.py

Requires muf_calculator_v2.py in the same directory.
"""

import streamlit as st
import math
import sys
from datetime import datetime, timezone

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

# ---------------------------------------------------------------------------
# Import the Stage 1 calculation engine
# ---------------------------------------------------------------------------
sys.path.insert(0, ".")
from muf_calculator_v2 import (
    grid_to_latlon,
    latlon_to_grid,
    haversine_distance,
    calculate_bearing,
    grid_distance,
    calculate_muf,
    estimate_foF2_from_sfi,
    get_space_weather_summary,
    HF_BANDS,
    DEFAULT_F2_HEIGHT_KM,
)

# ---------------------------------------------------------------------------
# Page configuration — must be the FIRST Streamlit command
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="MUF Calculator",
    page_icon="📡",
    layout="wide",
)

# ---------------------------------------------------------------------------
# Cached data fetching (prevents API spam on every widget interaction)
# ---------------------------------------------------------------------------
@st.cache_data(ttl=600)  # Re-fetch from NOAA every 10 minutes
def fetch_space_weather():
    """Fetch current space weather with caching."""
    try:
        return get_space_weather_summary()
    except Exception as e:
        return {
            "solar_flux": None,
            "kp_index": None,
            "overall": f"Could not fetch data: {e}",
            "errors": [str(e)],
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }


# ---------------------------------------------------------------------------
# Helper: 24-hour foF2 simulation
# ---------------------------------------------------------------------------
def simulated_foF2_hourly(hour_utc, sfi_val):
    """Simple sinusoidal model of foF2 diurnal variation."""
    base = 0.034 * sfi_val + 1.0
    phase = (hour_utc - 14) * (2 * math.pi / 24)
    variation = 0.5 * (1 + math.cos(phase))
    factor = 0.4 + 0.6 * variation
    return round(base * factor, 2)


# ===========================================================================
# SIDEBAR — Space Weather & Settings
# ===========================================================================
with st.sidebar:
    st.header("🌤️ Space Weather")
    st.caption(f"UTC: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M')}")

    weather = fetch_space_weather()

    # --- Solar Flux ---
    if weather.get("solar_flux"):
        sfi = weather["solar_flux"]["value"]
        st.metric(
            label="Solar Flux Index (SFI)",
            value=f"{sfi:.0f}",
            help="10.7 cm radio flux. Higher = more ionization = higher MUFs.",
        )
        st.caption(weather["solar_flux"]["description"])
    else:
        sfi = 120.0
        st.warning("SFI unavailable — using default 120")

    # --- Kp Index ---
    if weather.get("kp_index"):
        kp = weather["kp_index"]["value"]
        if kp < 4:
            st.metric("Kp Index", f"{kp:.1f} ✅")
        elif kp < 5:
            st.metric("Kp Index", f"{kp:.1f} ⚠️")
        else:
            st.metric("Kp Index", f"{kp:.1f} 🔴")
        st.caption(weather["kp_index"]["description"])
    else:
        kp = 2.0
        st.warning("Kp unavailable — using default 2.0")

    st.divider()
    st.caption(f"**Overall:** {weather.get('overall', 'N/A')}")

    if weather.get("errors"):
        for err in weather["errors"]:
            st.caption(f"⚠️ {err}")

    st.divider()
    st.header("⚙️ Settings")
    layer_height = st.slider(
        "F2 Layer Height (km)", 200, 400, 300, step=10,
        help="Typical F2 layer height: 250-350 km. Higher during solar max.",
    )

    st.divider()
    st.caption("Data: [NOAA SWPC](https://www.swpc.noaa.gov/)")
    st.caption("foF2 maps: [prop.kc2g.com](https://prop.kc2g.com/)")
    if st.button("🔄 Refresh NOAA Data"):
        st.cache_data.clear()
        st.rerun()


# ===========================================================================
# MAIN PANEL — Title & Inputs
# ===========================================================================
st.title("📡 MUF Calculator Dashboard")
st.markdown("*Real-time HF propagation analysis for amateur radio operators*")

# --- Grid Locator Inputs ---
st.header("Path Configuration")

col_g1, col_g2 = st.columns(2)
with col_g1:
    my_grid = st.text_input(
        "🏠 Your Grid Locator",
        value="CN88",
        max_chars=8,
        help="4 or 6 character Maidenhead grid (e.g., CN88 or CN88mm)",
    )
with col_g2:
    dx_grid = st.text_input(
        "🎯 Target Grid Locator",
        value="IO91wm",
        max_chars=8,
        help="Grid locator of the station you want to work",
    )

# --- Preset Quick Paths ---
presets = {
    "Custom (use inputs above)": None,
    "Canada BC → UK London": ("CN88", "IO91wm"),
    "Canada BC → Japan Tokyo": ("CN88", "PM95"),
    "Canada BC → Australia Sydney": ("CN88", "QF56"),
    "Canada BC → Germany Munich": ("CN88", "JN58"),
    "US East Coast → UK": ("FN31", "IO91wm"),
    "US East Coast → Brazil": ("FN31", "GG87"),
}
preset_choice = st.selectbox("Quick Presets", list(presets.keys()))
if presets[preset_choice] is not None:
    my_grid, dx_grid = presets[preset_choice]

# --- foF2 Source ---
foF2_mode = st.radio(
    "foF2 Source",
    options=["Estimate from SFI", "Manual Entry"],
    horizontal=True,
    help="'Estimate' uses current solar flux; 'Manual' lets you enter ionosonde data from prop.kc2g.com",
)

if foF2_mode == "Manual Entry":
    foF2_value = st.slider("foF2 (MHz)", 1.0, 15.0, 7.0, step=0.1)
    foF2_source_label = "manual"
else:
    estimate = estimate_foF2_from_sfi(sfi)
    foF2_value = estimate["estimated_foF2_mhz"]
    st.info(
        f"Estimated foF2: **{foF2_value} MHz** from SFI={sfi:.0f} — "
        f"_{estimate['notes']}_"
    )
    st.caption(
        "⚠️ Rough estimate. For accurate foF2, check "
        "[prop.kc2g.com](https://prop.kc2g.com/) and use Manual Entry."
    )
    foF2_source_label = "estimated from SFI"


# ===========================================================================
# CALCULATIONS
# ===========================================================================
try:
    # Validate grids and compute path
    path = grid_distance(my_grid, dx_grid)
    muf_result = calculate_muf(foF2_value, path["distance_km"], layer_height)

    if muf_result.get("error"):
        st.error(f"Propagation error: {muf_result['error']}")
        st.stop()

    # =======================================================================
    # RESULTS — Metrics Strip
    # =======================================================================
    st.header("📊 Results")

    col_m1, col_m2, col_m3, col_m4, col_m5 = st.columns(5)
    col_m1.metric("Distance", f"{path['distance_km']:,.0f} km")
    col_m2.metric("Bearing", f"{path['bearing']:.1f}°")
    col_m3.metric("Hops", f"{muf_result['num_hops']}")
    col_m4.metric("MUF", f"{muf_result['muf_mhz']:.1f} MHz")
    col_m5.metric("FOT (85%)", f"{muf_result['fot_mhz']:.1f} MHz")

    st.caption(
        f"Path: {path['grid1']} ({path['lat1']:.2f}°N, {path['lon1']:.2f}°E) → "
        f"{path['grid2']} ({path['lat2']:.2f}°N, {path['lon2']:.2f}°E) | "
        f"Return bearing: {path['reverse_bearing']:.1f}° | "
        f"foF2: {foF2_value:.1f} MHz ({foF2_source_label}) | "
        f"Inc. angle: {muf_result['incidence_angle_deg']:.1f}°"
    )

    # =======================================================================
    # BAND RECOMMENDATION TABLE
    # =======================================================================
    st.subheader("Band Availability")

    band_data = []
    for freq, band_name in HF_BANDS:
        if freq <= muf_result["fot_mhz"]:
            status = "✅ Excellent"
            color = "green"
        elif freq <= muf_result["muf_mhz"]:
            status = "⚠️ Marginal"
            color = "orange"
        else:
            status = "❌ Closed"
            color = "red"
        band_data.append({
            "Band": band_name,
            "Freq (MHz)": freq,
            "Status": status,
        })

    st.dataframe(
        pd.DataFrame(band_data),
        use_container_width=True,
        hide_index=True,
    )

    # =======================================================================
    # CHART 1: MUF vs Distance
    # =======================================================================
    st.subheader("MUF vs. Distance")

    fig1, ax1 = plt.subplots(figsize=(10, 5))
    distances = list(range(500, 10001, 100))

    for foF2 in [4.0, 6.0, 8.0, 10.0, 12.0]:
        mufs = []
        for d in distances:
            r = calculate_muf(foF2, d, layer_height)
            mufs.append(r["muf_mhz"] if r and "muf_mhz" in r else 0)
        is_current = abs(foF2 - foF2_value) < 1.0
        lw = 3 if is_current else 1.5
        alpha = 1.0 if is_current else 0.5
        label = f"foF2={foF2}" + (" ◀ current" if is_current else "")
        ax1.plot(distances, mufs, linewidth=lw, alpha=alpha, label=label)

    ax1.axvline(x=path["distance_km"], color="red", linestyle="--", alpha=0.7,
                label=f"Path: {path['distance_km']:.0f} km")

    band_colors = ['#ff9999', '#ffcc99', '#ffff99', '#ccff99', '#99ff99',
                   '#99ffcc', '#99ccff', '#cc99ff', '#ff99cc']
    for i, (freq, band) in enumerate(HF_BANDS):
        ax1.axhline(y=freq, color=band_colors[i], alpha=0.3, linestyle="--")
        ax1.text(10100, freq + 0.3, band, fontsize=7, color="gray")

    ax1.set_xlabel("Distance (km)", fontsize=11)
    ax1.set_ylabel("MUF (MHz)", fontsize=11)
    ax1.set_title("Maximum Usable Frequency vs. Distance", fontsize=13)
    ax1.legend(fontsize=8, loc="upper left")
    ax1.grid(True, alpha=0.2)
    ax1.set_xlim(500, 10000)
    ax1.set_ylim(0, 50)
    plt.tight_layout()
    st.pyplot(fig1)
    plt.close(fig1)

    # =======================================================================
    # CHART 2: 24-Hour Propagation Forecast
    # =======================================================================
    st.subheader("24-Hour Propagation Forecast")

    hours = list(range(24))
    foF2_by_hour = [simulated_foF2_hourly(h, sfi) for h in hours]
    muf_by_hour, fot_by_hour = [], []

    for fh in foF2_by_hour:
        r = calculate_muf(fh, path["distance_km"], layer_height)
        if r and "muf_mhz" in r:
            muf_by_hour.append(r["muf_mhz"])
            fot_by_hour.append(r["fot_mhz"])
        else:
            muf_by_hour.append(0)
            fot_by_hour.append(0)

    fig2, ax2 = plt.subplots(figsize=(10, 5))
    ax2.plot(hours, muf_by_hour, "r-o", linewidth=2, markersize=4, label="MUF")
    ax2.plot(hours, fot_by_hour, "g--", linewidth=2, label="FOT (85%)")
    ax2.fill_between(hours, fot_by_hour, muf_by_hour, alpha=0.15, color="orange",
                     label="Marginal zone")
    ax2.fill_between(hours, 0, fot_by_hour, alpha=0.1, color="green",
                     label="Reliable zone")

    for freq, band in HF_BANDS:
        if freq < max(muf_by_hour) * 1.1:
            ax2.axhline(y=freq, color="gray", alpha=0.3, linestyle=":")
            ax2.text(23.5, freq, f" {band}", fontsize=7, va="center", color="gray")

    current_hour = datetime.now(timezone.utc).hour
    ax2.axvline(x=current_hour, color="blue", linewidth=2, alpha=0.6,
                label=f"Now ({current_hour:02d}Z)")

    ax2.set_xlabel("Hour (UTC)", fontsize=11)
    ax2.set_ylabel("Frequency (MHz)", fontsize=11)
    ax2.set_title(f"24-Hour Forecast: {path['grid1']} → {path['grid2']}  "
                  f"({path['distance_km']:,.0f} km)", fontsize=13)
    ax2.set_xlim(0, 23)
    ax2.xaxis.set_major_locator(ticker.MultipleLocator(2))
    ax2.grid(True, alpha=0.2)
    ax2.legend(fontsize=8, loc="upper right")
    plt.tight_layout()
    st.pyplot(fig2)
    plt.close(fig2)

    # =======================================================================
    # CHART 3: Band Availability Heatmap
    # =======================================================================
    st.subheader("Band Availability Heatmap")

    fig3, ax3 = plt.subplots(figsize=(12, 4.5))

    availability = []
    for freq, band in HF_BANDS:
        row = []
        for muf_h, fot_h in zip(muf_by_hour, fot_by_hour):
            if freq <= fot_h:
                row.append(2)
            elif freq <= muf_h:
                row.append(1)
            else:
                row.append(0)
        availability.append(row)

    cmap = ListedColormap(["#ff6b6b", "#ffd93d", "#6bcf7f"])
    ax3.imshow(availability, aspect="auto", cmap=cmap, vmin=0, vmax=2,
               interpolation="nearest")

    ax3.set_xticks(range(24))
    ax3.set_xticklabels([f"{h:02d}" for h in range(24)], fontsize=8)
    ax3.set_yticks(range(len(HF_BANDS)))
    ax3.set_yticklabels([f"{b} ({f})" for f, b in HF_BANDS], fontsize=9)
    ax3.set_xlabel("Hour (UTC)", fontsize=11)
    ax3.set_title(f"Band Availability: {path['grid1']} → {path['grid2']}  |  "
                  f"SFI: {sfi:.0f}  |  Kp: {kp:.1f}", fontsize=13)

    ax3.axvline(x=current_hour, color="blue", linewidth=2, alpha=0.7)

    legend_elements = [
        Patch(facecolor="#6bcf7f", label="Excellent"),
        Patch(facecolor="#ffd93d", label="Marginal"),
        Patch(facecolor="#ff6b6b", label="Closed"),
    ]
    ax3.legend(handles=legend_elements, loc="upper right", fontsize=9)
    plt.tight_layout()
    st.pyplot(fig3)
    plt.close(fig3)

    # =======================================================================
    # MAP
    # =======================================================================
    st.subheader("🗺️ Signal Path")

    map_data = pd.DataFrame({
        "lat": [path["lat1"], path["lat2"]],
        "lon": [path["lon1"], path["lon2"]],
    })
    st.map(map_data, zoom=1, use_container_width=True)

    st.caption(
        f"📍 **{path['grid1']}** ({path['lat1']:.2f}°N, {path['lon1']:.2f}°E) → "
        f"📍 **{path['grid2']}** ({path['lat2']:.2f}°N, {path['lon2']:.2f}°E) | "
        f"Bearing: {path['bearing']:.1f}° | Return: {path['reverse_bearing']:.1f}°"
    )

except ValueError as e:
    st.error(f"❌ Invalid input: {e}")
    st.info("Please enter valid 4 or 6 character Maidenhead grid locators (e.g., CN88 or FN31pr).")
except Exception as e:
    st.error(f"❌ Unexpected error: {e}")
    st.info("Check that muf_calculator_v2.py is in the same directory as this file.")

# =========================================================================
# FOOTER
# =========================================================================
st.divider()
st.caption(
    "MUF Calculator v2 — Stage 2 Streamlit Dashboard | "
    "Space weather data from [NOAA SWPC](https://www.swpc.noaa.gov/) | "
    "foF2 maps: [prop.kc2g.com](https://prop.kc2g.com/) | "
    "73 & Good DX!"
)
