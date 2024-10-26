# OMI SO2 Analysis

This project processes and analyzes sulfur dioxide (SO₂) data from OMI (Ozone Monitoring Instrument) Aura satellite's Level 2 swath files (`.nc4` format). The project specifically aims to process low level nxtrack and ntime data to into readable latitude and longitude data by employing Haversine distance calculations to find the closest matching lat/lon pairs. It locates sulfur column amounts over a specific geographic region using latitude and longitude data downloaded from the [NASA OMI dataset](https://disc.gsfc.nasa.gov/datasets/OMSO2_003/summary?keywords=aura%20sulfur) and visualizes the data and differences in heatmaps.

## Dependencies

The project uses the following Python libraries:
- `netCDF4`
- `numpy`
- `xarray`
- `matplotlib`
- `haversine`
- `imageio`
- `h5py`
- `os`
- `tqdm` *(optional, for progress bars)*
- `haversine` *(optional, for faster haversine function)*



Install dependencies via:
```bash
pip install -r requirements.txt
```

## Currently being worked on:
- Fix issue with the distance threshold being useless
- Generalize the code to also work for other agents other than Sulphur Column Amount

