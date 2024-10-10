# OMI SO2 Analysis

This project processes and analyzes sulfur dioxide (SO₂) data from OMI (Ozone Monitoring Instrument) Aura satellite's Level 2 swath files (`.nc4` format). The project specifically aims to locate and visualize sulfur column amounts over a specific geographic region using latitude and longitude data. The code employs Haversine distance calculations to find the closest matching lat/lon pairs and visualizes the differences in a heatmap.

## Dependencies

The project uses the following Python libraries:
- `netCDF4`
- `numpy`
- `scipy`
- `xarray`
- `tqdm`
- `matplotlib`

Install dependencies via:
```bash
pip install -r requirements.txt
```

## Currently being worked on:
- Optimize speed of stepping algorithm
- Fix issue with the distance threshold being useless
- Generalize the code to also work for other agents other than Sulphur Column Amount
- Generalize the code to be able to run through multiple nc4 files generating a 3d time plot showing the changes in pollutant concentration over time

