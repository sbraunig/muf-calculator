# 📡 MUF Calculator Dashboard

Real-time HF propagation dashboard for amateur radio operators. Calculates Maximum Usable Frequency (MUF) between any two Maidenhead grid locators using live NOAA space weather data. Features band availability heatmaps, 24-hour forecasts, and interactive signal path maps. Built with Python and Streamlit.

## Features

- **Live Space Weather** — Pulls real-time Solar Flux Index (SFI) and Kp geomagnetic index from NOAA SWPC
- **Maidenhead Grid Locator Support** — Enter 4 or 6 character grid squares (e.g., CN88, FN31pr)
- **Great-Circle Path Analysis** — Distance, bearing, and return bearing using the Haversine formula
- **MUF & FOT Calculation** — Maximum Usable Frequency and Frequency of Optimum Traffic with multi-hop support
- **Band Recommendations** — Shows which amateur HF bands (160m–10m) are usable for your path
- **24-Hour Propagation Forecast** — Simulated MUF variation over a full day with current-time marker
- **Band Availability Heatmap** — Visual grid showing which bands are open at which hours
- **Interactive Signal Path Map** — Map displaying both station locations
- **Preset Quick Paths** — Common DX paths pre-configured for quick selection

## Screenshot

After launching, the dashboard displays a sidebar with live space weather conditions and a main panel with path configuration, results, charts, and map.

## Installation

### Requirements
- Python 3.10 or higher
- pip

### Dependencies
- `streamlit~=1.30.0` - Web framework
- `pandas~=1.5.0` - Data processing
- `numpy~=1.24.0` - Scientific calculations
- `matplotlib~=3.6.0` - Visualization
- `requests~=2.28.0` - HTTP requests (NOAA API)

All dependencies are listed in `requirements.txt`.

### Setup

```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/muf-calculator.git
cd muf-calculator

# Create a virtual environment (recommended)
python -m venv venv

# Activate the virtual environment
# On Linux/Mac:
source venv/bin/activate
# On Windows:
# venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Run

```bash
streamlit run muf_dashboard.py
```

The app will open automatically in your browser at `http://localhost:8501`.

## Usage

1. Enter your Maidenhead grid locator in the **Your Grid Locator** field (e.g., `CN88`)
2. Enter the target station's grid in the **Target Grid Locator** field (e.g., `IO91wm`)
3. Choose an **foF2 source**: estimate from live SFI data or enter a manual value from [prop.kc2g.com](https://prop.kc2g.com/)
4. View results: MUF, FOT, band availability, charts, and signal path map
5. Adjust the **F2 Layer Height** slider in the sidebar for different ionospheric conditions

## Project Structure

```
muf-calculator/
├── muf_dashboard.py        # Streamlit web dashboard
├── muf_calculator_v2.py    # Calculation engine (grid conversion, distance,
│                           #   NOAA API, MUF math)
├── requirements.txt        # Python dependencies
└── README.md
```

The project separates the calculation engine (`muf_calculator_v2.py`) from the user interface (`muf_dashboard.py`). The engine can also be used standalone from the command line:

```bash
python muf_calculator_v2.py CN88 IO91wm
```

## Data Sources

- **Solar Flux Index & Kp Index** — [NOAA Space Weather Prediction Center](https://www.swpc.noaa.gov/) (free, no API key required)
- **foF2 Reference Maps** — [prop.kc2g.com](https://prop.kc2g.com/)

## How It Works

The MUF is calculated using the secant law:

**MUF = foF2 × sec(θ)**

Where **foF2** is the critical frequency of the ionospheric F2 layer and **θ** is the angle of incidence, derived from the great-circle distance between stations and the F2 layer height. For paths exceeding single-hop range (~3500 km), the calculator automatically determines the number of hops required.

The **Frequency of Optimum Traffic (FOT)** is 85% of MUF — frequencies below FOT are the most reliable for communication.

## Limitations

- The foF2 estimate from SFI is a rough approximation. For accurate values, use real ionosonde data from [prop.kc2g.com](https://prop.kc2g.com/) or GIRO DIDBase.
- The 24-hour forecast uses a simplified sinusoidal model of diurnal ionospheric variation. Real propagation depends on location, season, and solar cycle position.
- The calculator models F2 layer propagation only. Sporadic-E and other propagation modes are not included.

## Contributing

Contributions are welcome. Some ideas for future enhancements:

- Real ionosonde data integration from GIRO DIDBase
- Sporadic-E layer support
- Day/night terminator awareness
- VOACAP integration for professional-grade predictions
- Great-circle arc line on the map using pydeck

## License

This project is open source. See [LICENSE](LICENSE) for details.

## Acknowledgments

- [NOAA Space Weather Prediction Center](https://www.swpc.noaa.gov/) for free space weather data APIs
- [prop.kc2g.com](https://prop.kc2g.com/) for real-time foF2 maps

---

*73 & Good DX!*
