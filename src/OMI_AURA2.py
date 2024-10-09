from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import griddata
import xarray as xr
from tqdm import tqdm
from math import radians, cos, sin, asin, sqrt

def find_closest_index(value, array):
    """
    Find the index of the value in the array that is closest to the input value.
    """
    array = np.array(array)
    return np.abs(array - value).argmin()

def get_sulfur_column_amount(ntime_input, nxtrack_input, data, ntime_range=(594.0, 649.0), nxtrack_range=(1, 10)):
    """
    Retrieve the sulfur column amount for the nearest `ntime` and `nxtrack` in the nc4 file.
    
    Parameters:
    - ntime_input: Input `ntime` float value.
    - nxtrack_input: Input `nxtrack` float value.
    - data: 2D array from nc4 file.
    - ntime_range: Tuple representing the start and end of the `ntime` range.
    - nxtrack_range: Tuple representing the start and end of the `nxtrack` range.
    
    Returns:
    - The sulfur column amount at the nearest (ntime, nxtrack) position.
    """

    ntime_indices = np.linspace(ntime_range[0], ntime_range[1], data.shape[0])
    nxtrack_indices = np.linspace(nxtrack_range[0], nxtrack_range[1], data.shape[1])

    # Find the closest index for ntime and nxtrack
    closest_ntime_idx = find_closest_index(ntime_input, ntime_indices)
    closest_nxtrack_idx = find_closest_index(nxtrack_input, nxtrack_indices)

    # Retrieve the sulfur column amount from the 2D data array
    sulfur_column_amount = data[closest_ntime_idx, closest_nxtrack_idx]

    return sulfur_column_amount
# Example usage:
# Assuming `sulfur_data` is the 2D array read from your nc4 file
# sulfur_amount = get_sulfur_column_amount(612.5, 1.5, sulfur_data)
#from: https://stackoverflow.com/questions/29545704/fast-haversine-approximation-python-pandas
def haversine_np(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    
    All args must be of equal length.    
    
    """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    
    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2
    
    c = 2 * np.arcsin(np.sqrt(a))
    km = 6378.137 * c
    return km
def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance in kilometers between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles. Determines return value units.
    return c * r

def get_lat_lon(Latitude, Longitude, latlon, ntimes_nxtrack, threshold):

    ntimes, nxtrack = ntimes_nxtrack
    current_ntime = ntimes // 2
    current_nxtrack = nxtrack // 2
    target_lat, target_lon = latlon
    tries = 0
    while True:
        # Get the current lat/lon at (current_ntime, current_nxtrack)
        current_lat = Latitude[current_ntime][current_nxtrack]
        current_lon = Longitude[current_ntime][current_nxtrack]
        
        # Calculate the current haversine difference
        distance = haversine_np(current_lon, current_lat, target_lon, target_lat)
        
        # If the current difference meets the threshold, break the loop
        if distance <= threshold:
            # print("broke")
            break
        
        # Initialize variables to track the best move
        min_diff = distance
        best_move = None
        
        # Iterate over the 8 possible directions
        for dntime, dnxtrack in directions:
            new_ntime = current_ntime + dntime
            new_nxtrack = current_nxtrack + dnxtrack
            
            # Ensure the new indices are within bounds
            if 0 <= new_ntime < ntimes and 0 <= new_nxtrack < nxtrack:
                new_lat = Latitude[new_ntime][new_nxtrack]
                new_lon = Longitude[new_ntime][new_nxtrack]
                
                # Skip NaN values
                if np.isnan(new_lat) or np.isnan(new_lon):
                    continue
                
                # Calculate the difference for the new point
                new_diff = haversine_np(new_lon,new_lat, target_lon, target_lat)
                
                
                # If this move results in a smaller difference, update the best move
                if new_diff < min_diff:
                    min_diff = new_diff
                    # print(min_diff)
                    best_move = (new_ntime, new_nxtrack)
        
        # If no better move is found, return the current lat/lon and exit
        if best_move is None:
            # print(min_diff, end = " ")
            # print(f"No better move found. Closest lat/lon pair: {current_lat}, {current_lon} at ({current_ntime}, {current_nxtrack})")
            #return current_lat, current_lon, best_move[0], best_move[1], min_diff

            return current_lat,current_lon,current_ntime,current_nxtrack, min_diff
        
        # Otherwise, step to the best move
        current_ntime, current_nxtrack = best_move
        tries +=1
        # print(f"Stepped to ({current_ntime}, {current_nxtrack}) with lat/lon: {Latitude[current_ntime][current_nxtrack]}, {Longitude[current_ntime][current_nxtrack]}")
    # return current_lat,current_lon,current_ntime,current_nxtrack, new_diff
    return current_lat, current_lon, best_move[0], best_move[1], min_diff


# Open the .nc4 file
file_path = r"C:\Users\erich\Downloads\OMI-Aura_L2-OMSO2_2024m0910t1735-o107221_v003-2024m0911t122019.SUB.nc4"

nc_file = Dataset(file_path, mode='r')
#setup data
DataFields = nc_file["HDFEOS"]["SWATHS"]["OMI Total Column Amount SO2"]["Data Fields"]
GeoLocationFields = nc_file["HDFEOS"]["SWATHS"]["OMI Total Column Amount SO2"]["Geolocation Fields"]

Sulphur = DataFields.variables['ColumnAmountSO2'][:]
Latitude, Longitude = GeoLocationFields["Latitude"][:], GeoLocationFields["Longitude"][:]

min_lat = np.nanmin(Latitude)
max_lat = np.nanmax(Latitude)
min_lon = np.nanmin(Longitude)
max_lon = np.nanmax(Longitude)

#PARAMETERS
THRESHOLD = 0.1 
directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
ntimes = len(Latitude)
nxtrack = len(Latitude[0])

disections = 100
array_600x600 = np.full((disections, disections), np.nan)  
lon_steps = np.linspace(min_lon, max_lon, disections)
lat_steps = np.linspace(min_lat, max_lat, disections)

#min/max of lat/lon of the downloaded omi aura data
min_lat, min_lon = -18.788, -74.855
max_lat, max_lon = -12.788, -68.855

#min/max of lat/lon of the geo2d lat/lon file
data_min_lat = np.nanmin(Latitude)
data_max_lat = np.nanmax(Latitude)
data_min_lon = np.nanmin(Longitude)
data_max_lon = np.nanmax(Longitude)

for i in tqdm(range(disections)):
    
    for j in tqdm(range(disections)):
        current_lon = lon_steps[j]
        current_lat = lat_steps[i]


        lat, lon, time, track, diff = get_lat_lon(Latitude = Latitude,
                                            Longitude= Longitude,
                                            latlon = (current_lat, current_lon),
                                            ntimes_nxtrack= (ntimes, nxtrack),
                                            threshold = THRESHOLD)
        # print(lat, lan, end=" ")
        # if j == 8:
        #     print(lat,lon,time,track)
        #     print(Sulphur[time][track])
        array_600x600[i][j] = diff
    
lat, lon, time, track, diff = get_lat_lon(Latitude=Latitude,
                                          Longitude=Longitude,
                                        #   latlon=(-17.8713188,-70.1685),
                                          latlon = (min_lat,min_lon),
                                          ntimes_nxtrack=(ntimes,nxtrack),
                                          threshold=THRESHOLD)
print(lat,lon)

print(min_lat,min_lon)

print(diff)
# print(haversine_np(lon,lat,-70.1685,-17.8713188))
print(haversine_np(lat,lon,min_lat,min_lon))
print(haversine(lat,lon,min_lat,min_lon))






# print(Sulphur[time][track])
# print(len(array_600x600[0]))

# print(f" min : {np.min(array_600x600)} \n max : {np.max(array_600x600)}")

import matplotlib.pyplot as plt
plt.figure(figsize=(8, 8))  # Adjust the size for better clarity
plt.imshow(array_600x600, cmap='viridis', interpolation='nearest', aspect='auto')
plt.colorbar(label='Value')

# Show the heatmap
plt.title('600x600 Array Heatmap')
plt.show()