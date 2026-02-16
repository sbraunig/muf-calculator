#!/usr/bin/env python3
"""
Enhanced MUF (Maximum Usable Frequency) Calculator for Ham Radio - v2
Now with:
  - Maidenhead grid locator support
  - Great-circle distance & bearing calculations
  - Live NOAA Space Weather data (solar flux, Kp index)
  - Improved propagation estimates

Usage:
    Import as module:  from muf_calculator_v2 import *
    Interactive:       python muf_calculator_v2.py
    Command line:      python muf_calculator_v2.py <grid1> <grid2>

Examples:
    python muf_calculator_v2.py FN31pr IO91wm
    python muf_calculator_v2.py CN88 FN20
"""

import math
import sys
import json
from datetime import datetime, timezone

# We use try/except so the module loads even without 'requests' installed
try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

# =============================================================================
# CONSTANTS
# =============================================================================

EARTH_RADIUS_KM = 6371.0       # Earth's mean radius in km
DEFAULT_F2_HEIGHT_KM = 300.0    # Typical F2 layer height in km

# Amateur HF bands: (lower edge MHz, band name)
HF_BANDS = [
    (1.800,  "160m"),
    (3.500,  "80m"),
    (7.000,  "40m"),
    (10.100, "30m"),
    (14.000, "20m"),
    (18.068, "17m"),
    (21.000, "15m"),
    (24.890, "12m"),
    (28.000, "10m"),
]

# NOAA Space Weather API endpoints (free, no key required)
NOAA_URLS = {
    "solar_wind_mag":   "https://services.swpc.noaa.gov/products/solar-wind/mag-1-day.json",
    "kp_index":         "https://services.swpc.noaa.gov/products/noaa-planetary-k-index.json",
    "solar_flux":       "https://services.swpc.noaa.gov/json/f107_cm_flux.json",
    "geomag_forecast":  "https://services.swpc.noaa.gov/products/noaa-planetary-k-index-forecast.json",
    "sunspot_number":   "https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json",
}


# =============================================================================
# SECTION 1: MAIDENHEAD GRID LOCATOR
# =============================================================================

def grid_to_latlon(grid_square):
    """
    Convert a Maidenhead grid locator to latitude/longitude.

    Maidenhead locators use a hierarchical system:
      - 2 letters (field):     18 x 18 grid, each 20° lon x 10° lat
      - 2 digits (square):     10 x 10 sub-grid, each 2° lon x 1° lat
      - 2 letters (subsquare): 24 x 24 sub-grid, each 5' lon x 2.5' lat

    Examples:
        'FN31'    -> 4-character (center of square)
        'FN31pr'  -> 6-character (center of subsquare)

    Returns:
        tuple: (latitude, longitude) in decimal degrees
               Returns the CENTER of the specified grid area.

    Raises:
        ValueError: If the grid locator format is invalid.
    """
    grid = grid_square.strip()

    if len(grid) not in (4, 6, 8):
        raise ValueError(
            f"Grid locator must be 4, 6, or 8 characters, got {len(grid)}: '{grid}'"
        )

    # --- Field (first 2 letters: A-R) ---
    field_lon_letter = grid[0].upper()
    field_lat_letter = grid[1].upper()

    if not ('A' <= field_lon_letter <= 'R') or not ('A' <= field_lat_letter <= 'R'):
        raise ValueError(f"Field letters must be A-R, got '{grid[0:2]}'")

    lon = (ord(field_lon_letter) - ord('A')) * 20 - 180
    lat = (ord(field_lat_letter) - ord('A')) * 10 - 90

    # --- Square (2 digits: 0-9) ---
    if not grid[2].isdigit() or not grid[3].isdigit():
        raise ValueError(f"Square characters must be digits 0-9, got '{grid[2:4]}'")

    lon += int(grid[2]) * 2
    lat += int(grid[3]) * 1

    if len(grid) >= 6:
        # --- Subsquare (2 letters: a-x) ---
        sub_lon = grid[4].lower()
        sub_lat = grid[5].lower()

        if not ('a' <= sub_lon <= 'x') or not ('a' <= sub_lat <= 'x'):
            raise ValueError(f"Subsquare letters must be a-x, got '{grid[4:6]}'")

        lon += (ord(sub_lon) - ord('a')) * (2 / 24)
        lat += (ord(sub_lat) - ord('a')) * (1 / 24)

        if len(grid) == 8:
            # Extended square (2 more digits: 0-9)
            lon += int(grid[6]) * (2 / 240)
            lat += int(grid[7]) * (1 / 240)
            # Center of extended square
            lon += (2 / 240) / 2
            lat += (1 / 240) / 2
        else:
            # Center of subsquare
            lon += (2 / 24) / 2
            lat += (1 / 24) / 2
    else:
        # Center of square (4-char locator)
        lon += 1.0    # half of 2°
        lat += 0.5    # half of 1°

    return (round(lat, 6), round(lon, 6))


def latlon_to_grid(lat, lon, precision=3):
    """
    Convert latitude/longitude to a Maidenhead grid locator.

    Args:
        lat: Latitude in decimal degrees (-90 to 90)
        lon: Longitude in decimal degrees (-180 to 180)
        precision: 1=field(2 char), 2=square(4 char), 3=subsquare(6 char)

    Returns:
        str: Maidenhead grid locator string

    Raises:
        ValueError: If coordinates are out of range.
    """
    if not (-90 <= lat <= 90):
        raise ValueError(f"Latitude must be -90 to 90, got {lat}")
    if not (-180 <= lon <= 180):
        raise ValueError(f"Longitude must be -180 to 180, got {lon}")

    # Shift to positive range
    lon += 180
    lat += 90

    grid = ""

    # Field
    grid += chr(int(lon / 20) + ord('A'))
    grid += chr(int(lat / 10) + ord('A'))

    if precision >= 2:
        # Square
        lon_remainder = lon % 20
        lat_remainder = lat % 10
        grid += str(int(lon_remainder / 2))
        grid += str(int(lat_remainder / 1))

    if precision >= 3:
        # Subsquare
        lon_sub = (lon % 2) * 12
        lat_sub = (lat % 1) * 24
        grid += chr(int(lon_sub) + ord('a'))
        grid += chr(int(lat_sub) + ord('a'))

    return grid


# =============================================================================
# SECTION 2: GREAT-CIRCLE DISTANCE AND BEARING
# =============================================================================

def haversine_distance(lat1, lon1, lat2, lon2):
    """
    Calculate the great-circle distance between two points using the
    Haversine formula.

    This is the standard method for computing shortest-path distance on a
    sphere — important for radio propagation since signals follow great-circle
    paths.

    Args:
        lat1, lon1: First point (decimal degrees)
        lat2, lon2: Second point (decimal degrees)

    Returns:
        float: Distance in kilometers
    """
    # Convert to radians
    lat1_r = math.radians(lat1)
    lat2_r = math.radians(lat2)
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)

    # Haversine formula
    a = (math.sin(dlat / 2) ** 2 +
         math.cos(lat1_r) * math.cos(lat2_r) * math.sin(dlon / 2) ** 2)

    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    return EARTH_RADIUS_KM * c


def calculate_bearing(lat1, lon1, lat2, lon2):
    """
    Calculate the initial bearing (azimuth) from point 1 to point 2.

    This gives the compass direction you'd point your beam antenna for the
    short-path to the target station.

    Args:
        lat1, lon1: Origin point (decimal degrees)
        lat2, lon2: Destination point (decimal degrees)

    Returns:
        float: Bearing in degrees (0-360, where 0=North)
    """
    lat1_r = math.radians(lat1)
    lat2_r = math.radians(lat2)
    dlon_r = math.radians(lon2 - lon1)

    x = math.sin(dlon_r) * math.cos(lat2_r)
    y = (math.cos(lat1_r) * math.sin(lat2_r) -
         math.sin(lat1_r) * math.cos(lat2_r) * math.cos(dlon_r))

    bearing = math.degrees(math.atan2(x, y))
    return bearing % 360  # Normalize to 0-360


def grid_distance(grid1, grid2):
    """
    Calculate distance and bearing between two Maidenhead grid locators.

    Returns:
        dict with keys:
            distance_km, distance_mi, bearing, reverse_bearing,
            lat1, lon1, lat2, lon2
    """
    lat1, lon1 = grid_to_latlon(grid1)
    lat2, lon2 = grid_to_latlon(grid2)

    dist_km = haversine_distance(lat1, lon1, lat2, lon2)
    bearing = calculate_bearing(lat1, lon1, lat2, lon2)
    reverse_bearing = calculate_bearing(lat2, lon2, lat1, lon1)

    return {
        "grid1": grid1.upper()[:4] + grid1[4:].lower() if len(grid1) > 4 else grid1.upper(),
        "grid2": grid2.upper()[:4] + grid2[4:].lower() if len(grid2) > 4 else grid2.upper(),
        "lat1": lat1,
        "lon1": lon1,
        "lat2": lat2,
        "lon2": lon2,
        "distance_km": round(dist_km, 1),
        "distance_mi": round(dist_km * 0.621371, 1),
        "bearing": round(bearing, 1),
        "reverse_bearing": round(reverse_bearing, 1),
    }


# =============================================================================
# SECTION 3: NOAA SPACE WEATHER API
# =============================================================================

def fetch_noaa_data(data_type, timeout=10):
    """
    Fetch live space weather data from NOAA SWPC.

    These are free, public JSON endpoints — no API key required.
    NOAA updates them every 1–15 minutes depending on the data type.

    Args:
        data_type: One of 'kp_index', 'solar_flux', 'sunspot_number',
                   'geomag_forecast', 'solar_wind_mag'
        timeout: Request timeout in seconds

    Returns:
        dict or list: Parsed JSON response from NOAA

    Raises:
        ConnectionError: If the request fails
        ImportError: If 'requests' library is not installed
    """
    if not HAS_REQUESTS:
        raise ImportError(
            "The 'requests' library is required for live data.\n"
            "Install it with: pip install requests"
        )

    if data_type not in NOAA_URLS:
        raise ValueError(
            f"Unknown data type '{data_type}'. "
            f"Available: {list(NOAA_URLS.keys())}"
        )

    url = NOAA_URLS[data_type]

    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        raise ConnectionError(f"Failed to fetch {data_type} from NOAA: {e}")


def get_current_kp():
    """
    Get the most recent Kp index value.

    The Kp index (0-9) measures geomagnetic disturbance:
      0-1: Quiet        — Excellent HF propagation
      2-3: Unsettled    — Good propagation
      4:   Active       — Some degradation, especially at high latitudes
      5+:  Storm        — Significant HF disruption

    Returns:
        dict: {value, timestamp, description}
    """
    data = fetch_noaa_data("kp_index")

    # NOAA returns a list of lists: [time_tag, Kp, Kp_fraction, ...]
    # Skip the header row
    if len(data) < 2:
        raise ValueError("No Kp data available from NOAA")

    latest = data[-1]  # Most recent entry
    kp_value = float(latest[1])

    # Classify the Kp value
    if kp_value < 2:
        desc = "Quiet — Excellent HF conditions"
    elif kp_value < 4:
        desc = "Unsettled — Good HF conditions"
    elif kp_value < 5:
        desc = "Active — Some HF degradation"
    elif kp_value < 7:
        desc = "Minor/Moderate Storm — HF disrupted"
    else:
        desc = "Major Storm — Severe HF blackout likely"

    return {
        "value": kp_value,
        "timestamp": latest[0],
        "description": desc,
    }


def get_current_solar_flux():
    """
    Get the most recent 10.7 cm solar flux index (SFI).

    SFI correlates strongly with ionospheric ionization and therefore MUF:
      < 70:   Solar minimum — Low MUFs, poor DX on higher bands
      70-100: Low activity  — 20m usually open, 15m sometimes
      100-150: Moderate     — 15m regularly open, 10m sometimes
      150-200: High         — All bands frequently open
      > 200:  Very high     — Excellent DX on all bands

    Returns:
        dict: {value, timestamp, description}
    """
    data = fetch_noaa_data("solar_flux")

    if not data:
        raise ValueError("No solar flux data available from NOAA")

    latest = data[-1]  # Most recent entry
    flux_value = float(latest.get("flux", latest.get("f107", 0)))
    time_tag = latest.get("time_tag", "unknown")

    if flux_value < 70:
        desc = "Solar minimum — Low MUFs expected"
    elif flux_value < 100:
        desc = "Low activity — 20m likely open"
    elif flux_value < 150:
        desc = "Moderate — 15m/20m should be open"
    elif flux_value < 200:
        desc = "High activity — All HF bands active"
    else:
        desc = "Very high — Excellent DX conditions"

    return {
        "value": flux_value,
        "timestamp": time_tag,
        "description": desc,
    }


def get_space_weather_summary():
    """
    Get a combined summary of current space weather conditions.

    Returns:
        dict with solar_flux, kp_index, and overall assessment
    """
    summary = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "solar_flux": None,
        "kp_index": None,
        "overall": "Unknown — Could not retrieve data",
        "errors": [],
    }

    # Fetch solar flux
    try:
        summary["solar_flux"] = get_current_solar_flux()
    except Exception as e:
        summary["errors"].append(f"Solar flux: {e}")

    # Fetch Kp index
    try:
        summary["kp_index"] = get_current_kp()
    except Exception as e:
        summary["errors"].append(f"Kp index: {e}")

    # Generate overall assessment
    if summary["solar_flux"] and summary["kp_index"]:
        sfi = summary["solar_flux"]["value"]
        kp = summary["kp_index"]["value"]

        if sfi >= 100 and kp < 4:
            summary["overall"] = "Good — Favorable HF propagation expected"
        elif sfi >= 100 and kp >= 4:
            summary["overall"] = "Mixed — High solar flux but disturbed geomagnetic field"
        elif sfi < 100 and kp < 3:
            summary["overall"] = "Fair — Lower bands (40m/80m) should work, higher bands limited"
        else:
            summary["overall"] = "Poor — Low ionization and/or disturbed conditions"

    return summary


# =============================================================================
# SECTION 4: MUF CALCULATION (improved from v1)
# =============================================================================

def calculate_hops(distance_km, layer_height_km=DEFAULT_F2_HEIGHT_KM):
    """
    Estimate the number of ionospheric hops needed for a given distance.
    Single hop max is typically ~3500 km for F2 layer.
    """
    max_single_hop = 2 * math.sqrt(
        (EARTH_RADIUS_KM + layer_height_km) ** 2 - EARTH_RADIUS_KM ** 2
    )
    return max(1, math.ceil(distance_km / max_single_hop))


def calculate_incidence_angle(distance_km, layer_height_km=DEFAULT_F2_HEIGHT_KM):
    """
    Calculate the angle of incidence at the ionosphere using spherical Earth
    geometry.

    Returns:
        tuple: (angle_radians, num_hops)
    """
    num_hops = calculate_hops(distance_km, layer_height_km)
    hop_distance = distance_km / num_hops
    half_distance = hop_distance / 2

    # Central angle subtended at Earth's center
    central_angle = half_distance / EARTH_RADIUS_KM
    reflection_radius = EARTH_RADIUS_KM + layer_height_km

    # Law of sines: sin(incidence) / R_earth = sin(central_angle) / R_reflection
    sin_incidence = (EARTH_RADIUS_KM * math.sin(central_angle)) / reflection_radius
    sin_incidence = max(-1.0, min(1.0, sin_incidence))

    return math.asin(sin_incidence), num_hops


def calculate_muf(critical_freq_mhz, distance_km, layer_height_km=DEFAULT_F2_HEIGHT_KM):
    """
    Calculate the Maximum Usable Frequency (MUF).

    MUF = foF2 × sec(θ)

    Where foF2 is the critical frequency and θ is the incidence angle.

    Args:
        critical_freq_mhz: F2 layer critical frequency (foF2) in MHz
        distance_km: Great-circle distance in km
        layer_height_km: Height of F2 layer (default 300 km)

    Returns:
        dict with muf, fot, num_hops, incidence_angle_deg, and band info
    """
    incidence_angle, num_hops = calculate_incidence_angle(distance_km, layer_height_km)
    angle_deg = math.degrees(incidence_angle)

    cos_angle = math.cos(incidence_angle)
    if cos_angle <= 0:
        return {
            "muf_mhz": None,
            "fot_mhz": None,
            "num_hops": num_hops,
            "incidence_angle_deg": angle_deg,
            "error": "Path geometry does not support F2 propagation",
        }

    muf = critical_freq_mhz / cos_angle
    fot = muf * 0.85  # Frequency of Optimum Traffic

    # Determine usable bands
    usable_bands = []
    for freq, band_name in HF_BANDS:
        if freq <= muf:
            if freq <= fot:
                status = "Excellent"
            else:
                status = "Marginal"
            usable_bands.append({
                "band": band_name,
                "freq_mhz": freq,
                "status": status,
            })

    return {
        "muf_mhz": round(muf, 2),
        "fot_mhz": round(fot, 2),
        "num_hops": num_hops,
        "incidence_angle_deg": round(angle_deg, 1),
        "critical_freq_mhz": critical_freq_mhz,
        "distance_km": distance_km,
        "layer_height_km": layer_height_km,
        "usable_bands": usable_bands,
    }


def estimate_foF2_from_sfi(sfi, hour_utc=None, lat=None):
    """
    Rough estimate of foF2 based on solar flux index.

    THIS IS A VERY SIMPLIFIED MODEL — real foF2 depends on location,
    season, time of day, and many other factors. For accurate values,
    use ionosonde data (GIRO) or the IRI model.

    This is included as a learning tool to show the relationship between
    solar activity and ionospheric critical frequency.

    Args:
        sfi: Solar flux index (10.7 cm flux)
        hour_utc: Hour of day in UTC (0-23). If None, uses current time.
        lat: Latitude in degrees (optional, for day/night estimate)

    Returns:
        dict: {estimated_foF2, confidence, notes}
    """
    if hour_utc is None:
        hour_utc = datetime.now(timezone.utc).hour

    # Base foF2 estimate from empirical SFI relationship
    # foF2 roughly follows: foF2 ≈ 0.03 * SFI + some_offset
    # This is a gross simplification!
    base_foF2 = 0.034 * sfi + 1.0

    # Very rough day/night factor
    # Real ionosphere has complex diurnal variation
    if 6 <= hour_utc <= 18:
        time_factor = 1.0  # Daytime
        time_note = "Daytime estimate (higher ionization)"
    elif hour_utc <= 4 or hour_utc >= 22:
        time_factor = 0.5  # Nighttime — F2 weakens significantly
        time_note = "Nighttime estimate (reduced ionization)"
    else:
        time_factor = 0.75  # Transition
        time_note = "Twilight estimate (transitional)"

    estimated = round(base_foF2 * time_factor, 1)

    return {
        "estimated_foF2_mhz": estimated,
        "sfi_used": sfi,
        "hour_utc": hour_utc,
        "confidence": "Low — Use real ionosonde data when available",
        "notes": time_note,
        "suggestion": "For accurate foF2: check https://prop.kc2g.com/ or GIRO DIDBase",
    }


# =============================================================================
# SECTION 5: CONVENIENCE / COMBINED FUNCTIONS
# =============================================================================

def full_path_analysis(grid1, grid2, critical_freq_mhz=None,
                       layer_height_km=DEFAULT_F2_HEIGHT_KM):
    """
    Complete path analysis between two grid locators.

    If critical_freq_mhz is not provided, attempts to estimate it
    from live NOAA solar flux data.

    Args:
        grid1: Origin Maidenhead grid locator
        grid2: Destination Maidenhead grid locator
        critical_freq_mhz: Manual foF2 value, or None to estimate
        layer_height_km: F2 layer height (default 300 km)

    Returns:
        dict: Complete analysis results
    """
    # Step 1: Path geometry
    path = grid_distance(grid1, grid2)

    # Step 2: Space weather (if available)
    weather = None
    try:
        weather = get_space_weather_summary()
    except Exception:
        pass

    # Step 3: Determine foF2
    foF2_source = "manual"
    if critical_freq_mhz is None:
        if weather and weather["solar_flux"]:
            sfi = weather["solar_flux"]["value"]
            estimate = estimate_foF2_from_sfi(sfi)
            critical_freq_mhz = estimate["estimated_foF2_mhz"]
            foF2_source = "estimated from SFI"
        else:
            critical_freq_mhz = 7.0  # Fallback default
            foF2_source = "default (no live data available)"

    # Step 4: Calculate MUF
    muf_result = calculate_muf(critical_freq_mhz, path["distance_km"], layer_height_km)

    return {
        "path": path,
        "space_weather": weather,
        "foF2_source": foF2_source,
        "muf": muf_result,
    }


# =============================================================================
# SECTION 6: DISPLAY / CLI
# =============================================================================

def print_path_info(path):
    """Pretty-print path information."""
    print(f"\n{'='*55}")
    print(f"  PATH: {path['grid1']}  →  {path['grid2']}")
    print(f"{'='*55}")
    print(f"  Distance:     {path['distance_km']:>8.1f} km  ({path['distance_mi']:.1f} mi)")
    print(f"  Bearing:      {path['bearing']:>8.1f}°")
    print(f"  Return:       {path['reverse_bearing']:>8.1f}°")
    print(f"  Origin:       {path['lat1']:.4f}°N, {path['lon1']:.4f}°E")
    print(f"  Destination:  {path['lat2']:.4f}°N, {path['lon2']:.4f}°E")


def print_weather_info(weather):
    """Pretty-print space weather data."""
    print(f"\n{'-'*55}")
    print(f"  SPACE WEATHER (live from NOAA)")
    print(f"{'-'*55}")

    if weather.get("solar_flux"):
        sf = weather["solar_flux"]
        print(f"  Solar Flux (SFI):  {sf['value']:.0f}")
        print(f"    → {sf['description']}")

    if weather.get("kp_index"):
        kp = weather["kp_index"]
        print(f"  Kp Index:          {kp['value']:.1f}")
        print(f"    → {kp['description']}")

    print(f"  Overall:           {weather['overall']}")

    if weather.get("errors"):
        for err in weather["errors"]:
            print(f"  ⚠ {err}")


def print_muf_results(muf_result, foF2_source=""):
    """Pretty-print MUF calculation results."""
    print(f"\n{'-'*55}")
    print(f"  MUF CALCULATION")
    print(f"{'-'*55}")

    if muf_result.get("error"):
        print(f"  Error: {muf_result['error']}")
        return

    print(f"  foF2:             {muf_result['critical_freq_mhz']:.1f} MHz ({foF2_source})")
    print(f"  Layer height:     {muf_result['layer_height_km']:.0f} km")
    print(f"  Incidence angle:  {muf_result['incidence_angle_deg']:.1f}°")
    print(f"  Number of hops:   {muf_result['num_hops']}")
    print()
    print(f"  >>> MUF:          {muf_result['muf_mhz']:.2f} MHz")
    print(f"  >>> FOT (85%):    {muf_result['fot_mhz']:.2f} MHz")

    bands = muf_result.get("usable_bands", [])
    if bands:
        print(f"\n  {'Band':<8} {'Freq':<10} {'Status'}")
        print(f"  {'─'*8} {'─'*10} {'─'*12}")
        for b in bands:
            print(f"  {b['band']:<8} {b['freq_mhz']:<10.3f} {b['status']}")
    else:
        print("\n  No amateur HF bands usable at this MUF.")

    print(f"{'='*55}")


def main():
    """Command-line interface."""
    print("=" * 55)
    print("   MUF Calculator v2 — Enhanced for Ham Radio")
    print("=" * 55)

    if len(sys.argv) >= 3:
        # Command-line mode
        grid1 = sys.argv[1]
        grid2 = sys.argv[2]
        foF2 = float(sys.argv[3]) if len(sys.argv) > 3 else None

        try:
            result = full_path_analysis(grid1, grid2, foF2)
            print_path_info(result["path"])
            if result["space_weather"]:
                print_weather_info(result["space_weather"])
            print_muf_results(result["muf"], result["foF2_source"])
        except Exception as e:
            print(f"Error: {e}")
            sys.exit(1)
    else:
        # Interactive mode
        print("\nEnter your grid locator (e.g., FN31pr):")
        grid1 = input("  Your grid:   ").strip()
        grid2 = input("  Target grid: ").strip()

        foF2_input = input("  foF2 MHz (Enter for auto-estimate): ").strip()
        foF2 = float(foF2_input) if foF2_input else None

        try:
            result = full_path_analysis(grid1, grid2, foF2)
            print_path_info(result["path"])
            if result["space_weather"]:
                print_weather_info(result["space_weather"])
            print_muf_results(result["muf"], result["foF2_source"])
        except ValueError as e:
            print(f"Input error: {e}")
        except Exception as e:
            print(f"Error: {e}")

    print("\n73! Good DX!")


if __name__ == "__main__":
    main()
