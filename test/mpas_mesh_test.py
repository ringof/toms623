import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np
from toms623 import trintrp

mesh_filename = 'x1.2562.grid.nc' # 480 km mesh
mesh_nc = Dataset(mesh_filename)
lats = mesh_nc.variables['latCell'][:]
print lats.min(), lats.max()
lons = mesh_nc.variables['lonCell'][:]
print lons.min(), lons.max()
# fake test data.
nexp = 8
icos_data = np.cos(nexp*lons)*np.sin(0.5*lons)**nexp*np.cos(lats)**nexp +\
np.sin(lats)**nexp

import time
t1 = time.clock()
print 'triangulation', lons.shape, lats.shape
tri = trintrp(lons, lats)
print 'done triangulation',time.clock()-t1
nlons = 3072; nlats = nlons/2 # 1x1 mesh
delta = 360./nlons
olons = delta*np.arange(nlons)
olats = -90.0 + 0.5*delta + delta*np.arange(nlats)
olons = np.radians(olons);  olats = np.radians(olats)
print olons.shape, olats.shape
olons, olats = np.meshgrid(olons, olats)
print olons.shape, olats.shape

print 'linear interp'
t1 = time.clock()
latlon_data = tri.interp_linear(olons,olats,icos_data)
print 'done linear interp',time.clock()-t1

latlon_datax = np.cos(nexp*olons)*np.sin(0.5*olons)**nexp*np.cos(olats)**nexp +\
np.sin(olats)**nexp
print (np.abs(latlon_datax-latlon_data)).max()
print latlon_data.min(), latlon_data.max()
m = Basemap(projection='ortho',lon_0=180,lat_0=40)
x,y = m(np.degrees(olons), np.degrees(olats))
m.drawcoastlines()
m.drawmapboundary()
m.contourf(x,y,latlon_data,15)
plt.show()
