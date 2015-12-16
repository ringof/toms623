import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np
from toms623 import trintrp
import time

# MPAS model mesh files from
# https://mpas-dev.github.io/atmosphere/atmosphere_meshes.html
mesh_filename = 'x1.2562.grid.nc' # 480 km mesh
#mesh_filename = 'x1.40962.grid.nc' # 120 km mesh
#mesh_filename = 'x1.2621442.grid.nc' # 15 km mesh
mesh_nc = Dataset(mesh_filename)
lats = mesh_nc.variables['latCell'][:]
print 'min/max lats:',lats.min(), lats.max()
lons = mesh_nc.variables['lonCell'][:]
print 'min/max lons:',lons.min(), lons.max()

# fake test data.
nexp = 8
icos_data = np.cos(nexp*lons)*np.sin(0.5*lons)**nexp*np.cos(lats)**nexp +\
np.sin(lats)**nexp

t1 = time.clock()
print 'triangulation of', len(lons),' points'
tri = trintrp(lons, lats, reorder='lat')
print 'triangulation took',time.clock()-t1,' secs'

nlons = 360; nlats = nlons/2 # 1 degree output mesh
delta = 360./nlons
olons = delta*np.arange(nlons)
olats = -90.0 + 0.5*delta + delta*np.arange(nlats)
olons = np.radians(olons);  olats = np.radians(olats)
olons, olats = np.meshgrid(olons, olats)

t1 = time.clock()
latlon_data = tri.interp_linear(olons,olats,icos_data)
print 'linear interp took',time.clock()-t1,' secs'

latlon_datax = np.cos(nexp*olons)*np.sin(0.5*olons)**nexp*np.cos(olats)**nexp +\
np.sin(olats)**nexp
print 'max abs error:',(np.abs(latlon_datax-latlon_data)).max()
print 'min/max field:',latlon_data.min(), latlon_datax.max()

# make plot on output mesh
m = Basemap(projection='ortho',lon_0=180,lat_0=20)
x,y = m(np.degrees(olons), np.degrees(olats))
m.drawcoastlines()
m.drawmapboundary()
m.contourf(x,y,latlon_data,15)
plt.show()
