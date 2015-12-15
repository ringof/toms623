import _toms623
import numpy as np

# Performs a triangulation, then interpolates within
# each triangle.
# Algorithm:
# R. J. Renka, "ALGORITHM 623:  Interpolation on the Surface of a
# Sphere", ACM Trans. Math. Software, Vol. 10, No. 4, December 1984,
# pp. 437-439.

class trintrp(object):
    def __init__(self, lons, lats,reorder = None):
        # given mesh points (lons,lats in radians)
        # define triangulation.  Must be called before interpolation.
        # n is size of input mesh (length of 1-d arrays lons and lats)
        if len(lons.shape) != 1 or len(lats.shape) != 1:
            raise ValueError('lons and lats must be 1d')
        npts = len(lons)
        if len(lats) != npts:
            raise ValueError('lons and lats must have same length')
        # compute cartesian coords on unit sphere.
        x,y,z = _toms623.trans(lats.astype(np.float64),\
                lons.astype(np.float64),npts)
        # reorder the coordinates by latitude.
        if reorder == 'lat':
            dummy = lats.astype(np.float64)
        elif reorder == 'lon':
            dummy = lons.astype(np.float64)
        elif reorder is None:
            dummy = None
        else:
            raise ValueError('illegal value for reorder kwarg')
        if dummy is not None:
            self.ind = _toms623.reordr(4,dummy,x,y,z,npts)
        else:
            self.ind = np.arange(1,npts+1,1)
        self.ind = self.ind - 1 # change to zero based indexing
        # create the triangulation.
        iadj,iend,ierr = _toms623.trmesh(x,y,z,npts)
        if ierr != 0:
            raise ValueError('warning: ierr = %s in trmesh' % ierr)
        self.lons = lons; self.lats = lats; self.npts = npts
        self.x = x; self.y = y; self.z = z
        self.iadj = iadj; self.iend = iend
    def interp_linear(self,olons,olats,data):
        # given a triangulation, perform interpolation on
        # points olons,olats (in radians), return result in odata.
        # nptso is number of points to interpolate to
        # (size of 1-d arrays olons,olats,odata).
        # 1-d input array data (length npts) contains values on
        # triangulated mesh to interpolate.
        # This version uses linear interpolation within each triangle.
        shapeout = olons.shape
        if len(shapeout) not in [1,2]:
            raise ValueError('olons,olats must be 1d or 2d')
        olons1 = olons.ravel(); olats1 = olats.ravel()
        nptso = len(olons1)
        if len(olats1) != nptso:
            raise ValueError('lons and lats must have same length')
        if len(data) != self.npts:
            raise ValueError('input data wrong size')
        # reorder input data based on sorting of nodes.
        data_reordered = data[self.ind].astype(np.float64)
        odata,ierr = \
        _toms623.intrpc0_n(olats1.astype(np.float64),olons1.astype(np.float64),\
                     self.x, self.y, self.z, data_reordered,\
                     self.iadj,self.iend,self.npts,nptso)
        if ierr != 0:
            raise ValueError('ierr = %s in trmesh' % ierr)
        return odata.reshape(shapeout)
    def interp_nn(self,olons,olats,data):
        # given a triangulation, perform interpolation on
        # points olons,olats (in radians), return result in odata.
        # nptso is number of points to interpolate to
        # (size of 1-d arrays olons,olats,odata).
        # 1-d input array data (length npts) contains values on
        # triangulated mesh to interpolate.
        # This version uses linear interpolation within each triangle.
        shapeout = olons.shape
        if len(shapeout) not in [1,2]:
            raise ValueError('olons,olats must be 1d or 2d')
        olons1 = olons.flatten(); olats1 = olats.flatten()
        nptso = len(olons1)
        if len(olats1) != nptso:
            raise ValueError('lons and lats must have same length')
        if len(data) != self.npts:
            raise ValueError('input data wrong size')
        data_reordered = data[self.ind].astype(np.float64)
        odata,ierr = \
        _toms623.intrpnn_n(olats1.astype(np.float64),olons1.astype(np.float64),\
                     self.x, self.y, self.z, data_reordered,\
                     self.iadj,self.iend,self.npts,nptso)
        if ierr != 0:
            raise ValueError('ierr = %s in trmesh' % ierr)
        return odata.reshape(shapeout)
