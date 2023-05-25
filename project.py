import pandas as pd
import glob
import os
import folium
import contextily as ctx
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

capacitydf=pd.read_excel('C:\\Users\\illes\\Documents\\repcsiadat_beviteé.xlsx')
capacitydf=capacitydf[['ICAO Code','Capacity merged']]
capacitydf=capacitydf.dropna()
capacitydf=capacitydf.sort_values(by=['Capacity merged'], ascending=False)
capacitydf=capacitydf.drop_duplicates(subset=['ICAO Code'], keep="first")

#%%for ciklusos betöltés
#create list of .csv.gz files in the folder
all_files = glob.glob("C:/Users/illes/Documents/*.csv.gz")

df_list = []

#read all the files
for filename in all_files:
    df = pd.read_csv(filename, compression='gzip')
    df=df[df['on_ground']==False]
    df=df.round({'latitude': 2, 'longitude': 2})
    df=df.drop_duplicates(subset=['latitude', 'longitude'], keep="last")
    df_list.append(df)
    
#concat all dataframes into one from the list
df = pd.concat(df_list, ignore_index=True)
df = df.drop_duplicates(subset=['latitude', 'longitude'], keep="last")
df = pd.merge(capacitydf, df, left_on='ICAO Code',right_on='icao_actype',how='right')
df = df[df['Capacity merged']>0]
#%%filterezés
df['timestamp'] = pd.to_datetime(df['timestamp'])
df.set_index("timestamp",inplace=True)
df_sampled=df.groupby('icao_address').resample('10min').mean()
df_sampled=df_sampled.dropna()
#%%most popular plane types
# Group the data by tail number and aircraft type, and count the occurrences
popular_aircraft = df.groupby(['icao_address', 'icao_actype']).size().reset_index()

# Group the data by aircraft type and sum the counts to get the popularity
popular_aircraft = popular_aircraft['icao_actype'].value_counts()

# Sort the data by popularity in descending order
popular_aircraft = popular_aircraft.sort_values(ascending=False)

# Optionally, limit the number of aircraft types to plot (e.g., top 10)
popular_aircraft = popular_aircraft.head(10)

# Plot the data using a bar plot
plt.figure(figsize=(10, 6))
popular_aircraft.plot(kind='bar')
plt.xlabel('Aircraft Type')
plt.ylabel('Count')
plt.title('Popular Aircraft Types')
plt.show()

#%%ábrázolás időről

import matplotlib.pyplot as plt
df_sampled=df_sampled.reset_index()
# convert the timestamp column to a datetime object
df_sampled['timestamp'] = pd.to_datetime(df_sampled['timestamp'])

# create a histogram of the timestamp column
df_sampled['timestamp'].hist(bins='auto',figsize=(10, 5))

# set the x-axis label
plt.xlabel('Timestamp')

# set the y-axis label
plt.ylabel('Frequency')

# show the plot
plt.show()
#%%ábrázolás a magasság eloszlásáról
plt.hist(df_sampled['altitude_baro'], bins=50, orientation='horizontal')

# Set the labels for the x-axis and y-axis
plt.ylabel('Barometric altitude')
plt.xlabel('Frequency')

# Set the y-axis limits
plt.ylim(0, 45000)

# Show the plot
plt.show()

#%%kaliningrad
from mpl_toolkits.basemap import Basemap
kaliningrad=df[(df['origin_airport_icao']=='UMKK')|(df['destination_airport_icao']=='UMKK')]
# Set up the map
# Create a scatter plot of the points in the dataframe
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_aspect('equal')
ax.scatter(kaliningrad['longitude'], kaliningrad['latitude'], color='blue',s=2)
# Add map background
ctx.add_basemap(ax, crs='EPSG:4326', source=ctx.providers.Stamen.TerrainBackground)

plt.show()
#%%skandinav orszagok

# Define the boundaries of Scandinavian countries
scandinavia_boundaries = {
    'norway': {'lon_min': 4.0, 'lon_max': 31.5, 'lat_min': 57.5, 'lat_max': 71.5},
    'sweden': {'lon_min': 11.0, 'lon_max': 24.5, 'lat_min': 55.0, 'lat_max': 69.0},
    'finland': {'lon_min': 19.0, 'lon_max': 31.5, 'lat_min': 59.5, 'lat_max': 70.0},
    'denmark': {'lon_min': 7.5, 'lon_max': 15.0, 'lat_min': 54.5, 'lat_max': 57.5},
}

# Filter the dataframe based on the coordinates
scandinavian_countries = df[
    (df['longitude'] >= scandinavia_boundaries['norway']['lon_min']) &
    (df['longitude'] <= scandinavia_boundaries['norway']['lon_max']) &
    (df['latitude'] >= scandinavia_boundaries['norway']['lat_min']) &
    (df['latitude'] <= scandinavia_boundaries['norway']['lat_max'])
    |
    (df['longitude'] >= scandinavia_boundaries['sweden']['lon_min']) &
    (df['longitude'] <= scandinavia_boundaries['sweden']['lon_max']) &
    (df['latitude'] >= scandinavia_boundaries['sweden']['lat_min']) &
    (df['latitude'] <= scandinavia_boundaries['sweden']['lat_max'])
    |
    (df['longitude'] >= scandinavia_boundaries['finland']['lon_min']) &
    (df['longitude'] <= scandinavia_boundaries['finland']['lon_max']) &
    (df['latitude'] >= scandinavia_boundaries['finland']['lat_min']) &
    (df['latitude'] <= scandinavia_boundaries['finland']['lat_max'])
    |
    (df['longitude'] >= scandinavia_boundaries['denmark']['lon_min']) &
    (df['longitude'] <= scandinavia_boundaries['denmark']['lon_max']) &
    (df['latitude'] >= scandinavia_boundaries['denmark']['lat_min']) &
    (df['latitude'] <= scandinavia_boundaries['denmark']['lat_max'])
]


fig, ax = plt.subplots(figsize=(10, 10))
ax.set_aspect('equal')
ax.scatter(scandinavian_countries['longitude'], scandinavian_countries['latitude'], color='blue',s=2,alpha=0.01)
# Add map background
ctx.add_basemap(ax, crs='EPSG:4326', source=ctx.providers.Stamen.TerrainBackground)

plt.show()
# The resulting dataframe will contain the data for Scandinavian countries
#%%
atlantic=df[(df['longitude']<-9)]
# Extract the latitude and longitude data
latitudes = atlantic['latitude']
longitudes = atlantic['longitude']

# Set up the map
map = Basemap(projection='merc', lat_0=0, lon_0=0,
              resolution='h', area_thresh=0.1,
              llcrnrlon=min(longitudes)-10, llcrnrlat=min(latitudes)-10,
              urcrnrlon=max(longitudes)+10, urcrnrlat=max(latitudes)+10)

# Draw coastlines and countries
map.drawcoastlines()
map.drawcountries()

# Convert latitude and longitude to map coordinates
x, y = map(longitudes, latitudes)

# Plot the data points
map.plot(x, y, 'bo', markersize=0.1)

# Show the map
plt.show()
#%%plotting samples
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_aspect('equal')
ax.scatter(df_sampled.iloc[::10]['longitude'], df_sampled.iloc[::10]['latitude'], color='blue',s=15,alpha=0.1)
# Add map background
ctx.add_basemap(ax, crs='EPSG:4326', source=ctx.providers.Stamen.TerrainBackground)

plt.show()
#%%algoritmus
from shapely.geometry import Point
from shapely.ops import cascaded_union
import numpy as np

def find_optimal_circles(df, circles_dict):
    # Create shapely points from x and y columns in dataframe
    points = [Point(x, y) for x, y in zip(df['latitude'], df['longitude'])]
    
    # Create a list to hold the circle objects
    circles = []
    covered=[]
    radius_list=[]
    
    # Iterate over the radii and numbers of circles in the dictionary
    for radius, num_circles in circles_dict.items():
        
        # Iterate through the specified number of circles
        for i in range(num_circles):
            
            # Initialize the maximum coverage and optimal circle for this iteration
            max_coverage = 0
            optimal_circle = None
            
            # Iterate through all points to find the circle with the most coverage
            for j, point in enumerate(points):
                
                # Create a circle around the current point
                circle = point.buffer(radius)
                
                # Calculate the coverage of the circle
                coverage = sum([circle.contains(p) for p in points])
                
                # If the coverage is greater than the current maximum, update the maximum
                if coverage > max_coverage:
                    max_coverage = coverage
                    optimal_circle = circle
            
            # Remove the points covered by the optimal circle from the list of points
            points = [p for p in points if not optimal_circle.contains(p)]
            
            # Add the optimal circle to the list of circles
            circles.append(optimal_circle)
            covered.append(max_coverage)
            radius_list.append(radius)
            
            #if number of covered datapoints is less than 5, break out of loop
            #if max_coverage<5:
            #    break
            
            # If there are no more points left, break out of the loop
            if not points:
                break
    
    return circles, covered, radius_list

#%%futtatás
#middles,covered_points=find_optimal_circles(df.head(100),1,5)
#%%földdel a háttérben
def plot_optimal_circles(df, circles_dict):
    # Find the optimal circles
    circles, covered, radius = find_optimal_circles(df, circles_dict)
    
    # Create a dictionary from the two lists
    data = {'circles': circles, 'covered': covered,'radius':radius}
    
    # Create a dataframe from the dictionary
    circlesdf = pd.DataFrame(data)
    
    # Create a scatter plot of the points in the dataframe
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect('equal')
    ax.scatter(df['longitude'], df['latitude'], color='blue',s=2,alpha=0.1)
    
    # Plot the optimal circles
    for index, row in circlesdf.iterrows():
        circle = row['circles']
        radius = row['radius']
        ax.add_patch(Circle((circle.centroid.y, circle.centroid.x), radius, color='red', fill=False))
    
    # Add map background
    ctx.add_basemap(ax, crs='EPSG:4326', source=ctx.providers.Stamen.TerrainBackground)
    
    plt.show()
    return circlesdf

#dictionary with radius and number of circles
circles_dict = {3: 50, 2: 40}
coveragedf=plot_optimal_circles(df_sampled.iloc[::10],circles_dict)
coveragedf.to_excel('50db3radius_40db2radius.xlsx')
#%%folium map print

fig, ax = plt.subplots(figsize=(10, 10))
ax.set_aspect('equal')
ax.scatter(df_sampled.iloc[::10]['longitude'], df_sampled.iloc[::10]['latitude'], color='blue',s=2,alpha=0.)

# Plot the optimal circles
for index, row in coveragedf.iterrows():
    circle = row['circles']
    radius = row['radius']
    ax.add_patch(Circle((circle.centroid.y, circle.centroid.x), radius, color='red', fill=False))

# Add map background
ctx.add_basemap(ax, crs='EPSG:4326', source=ctx.providers.Stamen.TerrainBackground)
plt.show()
#%%checkolni hogy pont benne van e
#circles az elso oszlop a df-ben
def is_coordinate_covered(coordinate, circles):
    point = Point(coordinate[0], coordinate[1])
    
    for circle in circles:
        if circle.contains(point):
            return True
    
    return False
#%%covered area
from shapely.ops import cascaded_union

def calculate_covered_area(circles):
    # Perform union operation to combine all circles into a single polygon
    coverage = cascaded_union(circles)

    # Calculate the total covered area
    total_area = coverage.area

    # Calculate the overlap area
    overlap_area = sum(circle.intersection(coverage).area for circle in circles)

    # Subtract the overlap area from the total covered area
    covered_area = total_area - overlap_area

    return covered_area
asd=calculate_covered_area(coveragedf['circles'])

print(sum(coveragedf['covered'])/len(df_sampled.iloc[::10])* 100)
