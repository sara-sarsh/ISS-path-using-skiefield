from datetime import datetime
import pytz
# from skyfield.plotting import plot_sky
import matplotlib.pyplot as plt
# from geopy import Nominatim
# from tzwhere import tzwhere
# from pytz import timezone, utc
# matplotlib to help display our star map
# skyfield for star data 
from skyfield.api import Star, Topos,load, wgs84
from skyfield.data import hipparcos,stellarium
from matplotlib.collections import LineCollection
from skyfield.projections import build_stereographic_projection
import numpy as np
# Load data for satellites
stations_url = 'http://celestrak.com/NORAD/elements/stations.txt'
satellites = load.tle_file(stations_url)
by_name = {sat.name: sat for sat in satellites}
satellite = by_name['ISS (ZARYA)']


planets = load('de421.bsp')
earth = planets['earth']

# Set up location for Isfahan
isfahan = Topos('32.6539 N', '51.6660 E')

# Load timescale and get a range of dates
ts = load.timescale()
t0 = ts.utc(2024, 5, 21)
t1 = ts.utc(2024, 6, 20)

# Use the satellite's `find_events` method
t, events = satellite.find_events(isfahan, t0, t1, altitude_degrees=30.0)

# Calculate the altitude of the ISS at each event
ephemeris = load('de421.bsp')
highest_altitude = 0
highest_event_time = None

# Define event names for the indices returned by find_events
event_names = ['rise', 'culminate', 'set']

# Initialize variables to store the rise time and the highest altitude event time
highest_altitude = 0
highest_event_time = None
rise_time = None
set_time=None
# Convert UTC to Isfahan's local time
isfahan_timezone = pytz.timezone('Asia/Tehran')

for ti, event in zip(t, events):
    # Convert to local time
    local_time = ti.astimezone(isfahan_timezone)
    # Save the rise time if the event is a rise or set
    if event_names[event] == 'set':
        set_time = local_time
        
    if event_names[event] == 'rise':
        rise_time = local_time
       
    # Check if it's night time
    if local_time.hour >= 18 or local_time.hour <= 6:
        difference = satellite - isfahan
        topocentric = difference.at(ti)
        alt, az, distance = topocentric.altaz()
            
        # Check if this event has the highest altitude so far
        if alt.degrees > highest_altitude:
            highest_altitude = alt.degrees
            highest_event_time = local_time

# Print the result
if highest_event_time:
    print(f"The ISS will transit over Isfahan at its highest altitude of {highest_altitude:.2f}° on {highest_event_time.strftime('%Y-%m-%d %H:%M:%S')} local time.")
else:
    print("No high-altitude night-time transits were found for the given date range.")




# de421 shows position of earth and sun in space
eph = load('de421.bsp')
# hipparcos dataset contains star location data
with load.open(hipparcos.URL) as f:
    stars = hipparcos.load_dataframe(f)
# constellation dataset
url = ('https://raw.githubusercontent.com/Stellarium/stellarium/master'
       '/skycultures/modern_st/constellationship.fab')

with load.open(url) as f:
    constellations = stellarium.parse_constellations(f)
# Set up location for Isfahan
isfahan = Topos('32.6539 N', '51.6660 E')
lat=32.6539
long=51.6660
# Get the current time in UTC
ts = load.timescale()
# utc_now = datetime.utcnow().replace(tzinfo=pytz.utc)
t = ts.from_datetime(rise_time)
# # define datetime and convert to utc based on our timezone
# timezone_str = tzwhere.tzwhere().tzNameAt(lat, long)
# local = timezone(timezone_str)
# # get UTC from local timezone and datetime
# local_dt = local.localize(dt, is_dst=None)
# utc_dt = local_dt.astimezone(utc)
# find location of earth and sun and set the observer position
sun = eph['sun']
earth = eph['earth']
# And the constellation outlines list.
edges = [edge for name, edges in constellations for edge in edges]
edges_star1 = [star1 for star1, star2 in edges]
edges_star2 = [star2 for star1, star2 in edges]
# # define observation time from our UTC datetime
# ts = load.timescale()
# t = ts.from_datetime(utc_dt)
# define an observer using the world geodetic system data
observer = wgs84.latlon(latitude_degrees=lat, longitude_degrees=long).at(t)
position = observer.from_altaz(alt_degrees=90, az_degrees=0)
ra, dec, distance = observer.radec()
center_object = Star(ra=ra, dec=dec)
center = earth.at(t).observe(center_object)
projection = build_stereographic_projection(center)
star_positions = earth.at(t).observe(Star.from_dataframe(stars))
stars['x'], stars['y'] = projection(star_positions)
chart_size = 10
max_star_size = 100
limiting_magnitude = 10
bright_stars = (stars.magnitude <= limiting_magnitude)
magnitude = stars['magnitude'][bright_stars]

# Define a function to convert satellite altaz to x, y for plotting
# Create a fake star at the given alt, az
# earth = wgs84.latlon(latitude_degrees=alt, longitude_degrees=long)
t = ts.from_datetime(rise_time)

# Observe the fake star from the observer's position
ISS_positions  = earth.at(t).observe(Star.from_dataframe(stars))
# Initialize a list to store the ISS path coordinates
iss_path={'x': [], 'y': []}
# Calculate the path of the ISS on the night of the highest event
if highest_event_time and rise_time and set_time:
    # Define the time range for the path: from rise to set
    t_path = ts.utc(rise_time.year, rise_time.month, rise_time.day, range(24))  # Generate times for the entire day
    for ti in t_path:
        # Calculate the topocentric position of the ISS
        difference = satellite - isfahan
        topocentric = difference.at(ti)
        alt, az, distance = topocentric.altaz()
        if alt.degrees > 0:  # Only plot when the ISS is above the horizon
            # Convert altaz to projection coordinates
            x, y = projection(topocentric)
            iss_path['x'].append(x)
            iss_path['y'].append(y)

# iss_path_x = []
# iss_path_y = []
# t_path = ts.utc(rise_time.year, rise_time.month, rise_time.day, range(24))  # Generate times for the entire day
# for ti in t_path:
#     # Project the apparent position onto the stereographic projection
#     projection = build_stereographic_projection(star_positions)
#     iss_path_x['x'], iss_path_y['y'] = projection(star_positions)
    
# Project the apparent position onto the stereographic projection
projection = build_stereographic_projection(star_positions)

# iss_path_x['x'], iss_path_y['y'] = projection(star_positions)
 

# Calculate the path of the ISS on the night of the highest event
# if highest_event_time and rise_time and set_time:
    # Define the time range for the path: from rise to set
    # t_path = ts.utc(rise_time.year, rise_time.month, rise_time.day, range(24))  # Generate times for the entire day
    # for ti in t_path:
    #     # if ti >= ts.from_datetime(rise_time) and ti <= ts.from_datetime(set_time):
          
        # # Calculate the topocentric position of the ISS
        # difference = satellite - isfahan
        # topocentric = difference.at(ti)
        # alt, az, distance = topocentric.altaz()
        # # if alt.degrees > 0:  # Only plot when the ISS is above the horizon
        # # Convert altaz to projection coordinates
        # x, y = altaz_to_projection(alt, az, projection)
        # iss_path_x.append(x.radians)
        # iss_path_y.append(y.radians)


fig, ax = plt.subplots(figsize=(chart_size, chart_size))
# Draw the constellation lines.
   
border = plt.Circle((0, 0), 1, color='navy', fill=True)
ax.add_patch(border)
marker_size = max_star_size * 10 ** (magnitude / -2.5)
# Calculate the constellation lines
xy1 = stars[['x', 'y']].loc[edges_star1].values
xy2 = stars[['x', 'y']].loc[edges_star2].values
lines_xy = np.rollaxis(np.array([xy1, xy2]), 1)
ax.add_collection(LineCollection(lines_xy, colors='#ffff', linewidths=0.15))

ax.scatter(stars['x'][bright_stars], stars['y'][bright_stars],
 s=marker_size, color='white', marker='.', linewidths=0, 
 zorder=2)
ax.add_collection(LineCollection(lines_xy, colors='#ffff', linewidths=0.15))
# Plot the ISS path
ax.plot(iss_path['x'], iss_path['y'], color='yellow', linewidth=2, linestyle='-', label='ISS Path')

ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
plt.axis('off')
plt.show()