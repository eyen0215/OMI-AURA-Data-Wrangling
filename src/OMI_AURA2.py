from netCDF4 import Dataset
import numpy as np
import xarray as xr
from tqdm import tqdm
from haversine import haversine

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

def plot_array(arrays, labels = "Value", cmap = 'viridis'):
    """
    Plot an array heatmap on the given axis

    Parameters:
    - arrays: Mandatory, list  of 2D numpy array
    - labels: Optional, list of  strings for heatmap label
    - cmap: Optional, input plt color map (default 'viridis')

    Returns:
    - None
    """
    import matplotlib.pyplot as plt
    num_arrays = len(arrays)

    #if labels not provided, generate default labels
    if labels is None:
        labels = [f"Array {i+1}" for i in range(num_arrays)]

    _, axes = plt.subplots(1, num_arrays, figsize=(8 * num_arrays, 8))

    #if only one array is provided, convert axes to list
    if num_arrays ==1:
        axes = [axes] 
    
    for i, (ax, array) in enumerate(zip(axes, arrays)):
        cax = ax.imshow(array, cmap=cmap, interpolation='nearest', aspect='auto')
        ax.set_title(f"Array {i+1} Heatmap")
        plt.colorbar(cax, ax=ax, label=labels[i])
    
    plt.tight_layout()
    plt.show()



def get_lat_lon(Latitude, Longitude, latlon, ntimes_nxtrack, threshold, directions):
    """
    Finds the closest latitude and longitude pair to a target latitude and longitude by iteratively searching
    through a Longitude and Latitude grid. The search begins at the midpoint of the grid and "steps" to neighboring points
    in one of eight directions (up, down, left, right, and diagonals) until either a point within a specified
    Haversine distance threshold is found or no closer point exists.

    Parameters:
    Latitude : 2D array-like
        A 2D array of latitude values corresponding to the geographic grid.
    
    Longitude : 2D array-like
        A 2D array of longitude values corresponding to the geographic grid.
    
    latlon : tuple
        A tuple (target_lat, target_lon) containing the target latitude and longitude to search for.
        
    ntimes_nxtrack : tuple
        A tuple containing the dimensions (ntimes, nxtrack) of the geographic grid.
    
    threshold : float
        The maximum allowable Haversine distance (in kilometers) between the target and a point for the search
        to terminate.
    
    directions : list of tuples
        A list of tuples representing the possible directions to "step" to neighboring points in the grid.
        Each tuple is of the form (dntime, dnxtrack), where `dntime` and `dnxtrack` are the relative changes
        in the ntimes and nxtrack dimensions.

    Returns:
    --------
    current_lat : float
        The latitude of the closest point found during the search.
    
    current_lon : float
        The longitude of the closest point found during the search.
    
    current_ntime : int
        The `ntime` index of the closest point found during the search.
    
    current_nxtrack : int
        The `nxtrack` index of the closest point found during the search.
    
    min_diff : float
        The Haversine distance (in kilometers) between the closest point found and the target latitude/longitude.
    """
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
        distance = haversine((current_lat, current_lon), (target_lat, target_lon))
        
        # If the current difference meets the threshold, break the loop
        if distance <= threshold:
            break
        
        # Initialize variables to track the best move
        min_diff = distance
        best_move = None
        
        # Iterate over the possible directions
        for dntime, dnxtrack in directions:
            new_ntime = current_ntime + dntime
            new_nxtrack = current_nxtrack + dnxtrack
            
            # Ensure the new indices are within bounds
            if 0 <= new_ntime < ntimes and 0 <= new_nxtrack < nxtrack:
                new_lat = Latitude[new_ntime][new_nxtrack]
                new_lon = Longitude[new_ntime][new_nxtrack]
                
                # Skip NaN values
                if str(new_lat) == "--" or str(new_lon) == "--":
                    continue
                
                # Calculate the difference for the new point
                
                new_diff = haversine((new_lat, new_lon), (target_lat, target_lon))

                # new_diff = haversine((new_lon,new_lat), (target_lon, target_lat))
                
                # If this move results in a smaller difference, update the best move
                if new_diff < min_diff:
                    min_diff = new_diff
                    best_move = (new_ntime, new_nxtrack)
        
        # If no better move is found, return the current lat/lon and exit
        if best_move is None:
            
            return current_lat,current_lon,current_ntime,current_nxtrack, min_diff
        
        # Otherwise, step to the best move
        current_ntime, current_nxtrack = best_move
        tries +=1
    return current_lat, current_lon, best_move[0], best_move[1], min_diff


# Open the .nc4 file
file_path = "src\OMI-Aura_L2-OMSO2_2024m0910t1735-o107221_v003-2024m0911t122019.SUB.nc4"
nc_file = Dataset(file_path, mode='r')

#setup data
DataFields = nc_file["HDFEOS"]["SWATHS"]["OMI Total Column Amount SO2"]["Data Fields"]
GeoLocationFields = nc_file["HDFEOS"]["SWATHS"]["OMI Total Column Amount SO2"]["Geolocation Fields"]

Sulfur = DataFields.variables['ColumnAmountSO2'][:]
Latitude, Longitude = GeoLocationFields["Latitude"][:], GeoLocationFields["Longitude"][:]

#PARAMETERS

THRESHOLD = 0.1 #Minimum loss when finding closest (ntime,nxtrack). Can be adjusted for speed
DIRECTIONS = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)] #Directions that can be stepped. Can be adjusted for interesting effects
DISECTIONS = 100 #dimensions of the array. ex: DISECTIONS = 100 generates 100x100 resultant arrays. Can be adjusted for speed

ntimes = len(Latitude)
nxtrack = len(Latitude[0])


diff_array = np.full((DISECTIONS, DISECTIONS), np.nan)  
sulf_array = np.full((DISECTIONS, DISECTIONS), np.nan)

#setup lat/lon steps through area
lon_steps = np.linspace(np.nanmin(Longitude), np.nanmax(Longitude), DISECTIONS)
lat_steps = np.linspace(np.nanmin(Latitude), np.nanmax(Latitude), DISECTIONS)

#min/max of lat/lon of the specified omi aura data
#used for testing the results of corners and sides of data range
min_lat, min_lon = -18.788, -74.855
max_lat, max_lon = -12.788, -68.855

#min/max of lat/lon of the geo2d lat/lon file
#used to find difference between the dimensions of the actual data downloaded and the specified data downloaded
'''
this is relevant because when downloading the satellite from the NASA OMI AURA website, we specified a certain
data range (seen above in min_lat, min_lon, max_lat, max_lon) but most of the data actually downloaded only 
corresponds to a subsection of this area. Thus there are lat/lon points that don't have any corresponding data
and we should decide whether we should keep them as NaN or smooth them out
'''
data_min_lat = np.nanmin(Latitude)
data_max_lat = np.nanmax(Latitude)
data_min_lon = np.nanmin(Longitude)
data_max_lon = np.nanmax(Longitude)

for i in tqdm(range(DISECTIONS)):
    
    for j in (range(DISECTIONS)):
        current_lon = lon_steps[j]
        current_lat = lat_steps[i]


        lat, lon, time, track, diff = get_lat_lon(Latitude = Latitude,
                                            Longitude= Longitude,
                                            latlon = (current_lat, current_lon),
                                            ntimes_nxtrack= (ntimes, nxtrack),
                                            threshold = THRESHOLD,
                                            directions = DIRECTIONS)
        
        diff_array[i][j] = diff
        sulf_array[i][j] = Sulfur[time][track]


plot_array([diff_array, sulf_array],
           ["Loss", "Vertical Sulfur Column Amount"])
