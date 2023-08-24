# Copyright (C) 2013-2016 Martin Vejmelka, UC Denver
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR
# A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import sys
import logging
import os.path as osp
from utils import inq, ensure_dir
from ingest.rtma_source import RTMA
from wrf.wps_format import WPSFormat
from geo.write_geogrid import write_geogrid_var
from geo.var_wisdom import get_wisdom
from six.moves import range

class FuelMoistureModel:
    """
    A pure-python implementation of independent moisture models running on a grid.
    """

    def __init__(self, m0, Tk=None, P0=None):
        """
        Initialize the model with given position and moisture levels.

        :param m0: the initial condition for the entire grid (shape: grid0 x grid1 x num_fuels)
        :param Tk: drying/wetting time constants of simulated fuels (one per fuel), default [1 10 100 1000]
        :param P0: initial state error covariance
        """
        logging.info('FuelMoistureModel.__init__ m0 %s, Tk %s, P0 %s' % (inq(m0),inq(Tk), inq(P0)))
        self.Tk = np.array([1.0, 10.0, 100.0]) * 3600    # nominal fuel delays
        self.r0 = 0.05                                   # threshold rainfall [mm/h]
        self.rk = 8.0                                    # saturation rain intensity [mm/h]
        self.Trk = 14.0 * 3600                           # time constant for wetting model [s]
        self.S = 2.5                                     # saturation intensity [dimensionless]

        s0, s1, k = m0.shape
        if Tk is not None:
            self.Tk = Tk
        dim = k + 2
        assert k == len(self.Tk)
        self.m_ext = np.zeros((s0, s1, dim))
        self.m_ext[:, :, :k] = m0

        # note: the moisture advance will proceed by fuel moisture types
        # thus we only need space for one class at a time
        self.m_i = np.zeros((s0, s1))
        #self.mn_i = np.zeros((s0,s1))
        self.rlag = np.zeros((s0, s1))
        self.equi = np.zeros((s0, s1))
        self.model_ids = np.zeros((s0, s1))

        self.EdA = np.zeros((s0, s1))
        self.EwA = np.zeros((s0, s1))

        # state covariance current and forecasted
        self.P = np.zeros((s0, s1, dim, dim))
        self.P2 = np.zeros((dim, dim))
        for s in np.ndindex((s0, s1)):
            self.P[s[0], s[1], :, :] = P0

        # fill out the fixed parts of the jacobian
        self.J = np.zeros((s0, s1, dim, dim))
        self.Jii = np.zeros((s0, s1))

        # note: the observation operator H is common for all dims
        self.H = np.zeros((k, dim))


    def advance_model(self, Ed, Ew, r, dt, mQ=None):
        """
        Advance the model by one time step of length dt

        :param Ed: drying equilibrium [grid shape]
        :param Ew: wetting equilibrium [grid shape]
        :param r: rain intensity during the dt time unit [grid shape] [mm/h]
        :param dt: integration step in seconds
        :param mQ: the proess noise matrix [shape dim x dim, common across grid points]
        """

        logging.info('fuel_moisture_model: advance_model')
        Tk = self.Tk
        k = Tk.shape[0]                 # number of fuel classes
        m_ext = self.m_ext
        equi = self.equi
        rlag = self.rlag
        model_ids = self.model_ids
        m_i = self.m_i
        #mn_i = self.mn_i
        J = self.J
        Jii = self.Jii
        P2 = self.P2
        P = self.P
        EdA, EwA = self.EdA, self.EwA

        # add assimilated equilibrium difference, which is shared across spatial locations
        EdA[:, :] = Ed
        EdA += m_ext[:, :, k]
        EwA[:, :] = Ew
        EwA += m_ext[:, :, k]

        dS = m_ext[:, :, k + 1]

        logging.info('GMM: (min, mean, max) EdA (%g,%g,%g)  EwA (%g,%g,%g)' % (np.amin(EdA),np.mean(EdA),np.amax(EdA), np.amin(EwA),np.mean(EdA),np.amax(EdA)))

        assert np.all(EdA >= EwA)

        # re-initialize the Jacobian (in case we must recompute it)
        J[:] = 0.0
        J[:, :, k, k] = 1.0
        J[:, :, k + 1, k + 1] = 1.0

        # where rainfall is above threshold (spatially different), apply
        # saturation model, equi and rlag are specific to fuel type and
        # location
        for i in range(k):

            # initialize equilibrium with current moisture
            m_i[:, :] = m_ext[:, :, i]
            equi[:, :] = m_i
            rlag[:, :] = 0.0
            model_ids[:] = 4

            # on grid locations where there is rain, modify equilibrium
            has_rain = r > self.r0
            no_rain = np.logical_not(has_rain)

            equi[has_rain] = self.S + dS[has_rain]
            rlag[has_rain] = 1.0 / self.Trk * (1.0 - np.exp(- (r[has_rain] - self.r0) / self.rk))
            model_ids[has_rain] = 3

            # equilibrium is selected according to current moisture level
            is_drying = np.logical_and(no_rain, m_i > EdA)
            equi[is_drying] = EdA[is_drying]
            model_ids[is_drying] = 1

            is_wetting = np.logical_and(no_rain, m_i < EwA)
            equi[is_wetting] = EwA[is_wetting]
            model_ids[is_wetting] = 2

            rlag[no_rain] = 1.0 / Tk[i]

            dead_zone = (model_ids == 4)

            # optional inspection of changes at given position
            #dg_pos = (86,205)
            #mi_old = m_i[dg_pos]

            # select appropriate integration method according to change for each fuel
            # and location
            rlag *= dt
            change = rlag
            big_change = change > 0.01
            m_i[big_change] += (equi[big_change] - m_i[big_change]) * (1.0 - np.exp(-change[big_change]))
            small_change = np.logical_not(big_change)
            m_i[small_change] += (equi[small_change] - m_i[small_change]) * change[small_change] * (
                1.0 - 0.5 * change[small_change])

            #print('Diag at 86,205 fuel %d: model_id %d equi %g rlag %g change %g value %g -> %g' % 
            #        (i, model_ids[dg_pos], equi[dg_pos], rlag[dg_pos], change[dg_pos], mi_old, m_i[dg_pos]))

            # store in state matrix
            m_ext[:, :, i] = m_i

            # update model state covariance if requested using the old state (the jacobian must be computed as well)
            if mQ is not None:

                # partial m_i/partial m_i
                Jii[big_change] = np.exp(-change[big_change])
                Jii[small_change] = 1.0 - change[small_change] * (1.0 - 0.5 * change[small_change])
                Jii[dead_zone] = 1.0
                J[:, :, i, i] = Jii

                # partial m_i/partial deltaE
                Jii[:] = 0.0
                norain_big = np.logical_and(no_rain, big_change)
                norain_small = np.logical_and(no_rain, small_change)
                Jii[norain_big] = 1.0 - np.exp(-change[norain_big])
                Jii[norain_small] = change[norain_small] * (1.0 - 0.5 * change[norain_small])
                Jii[dead_zone] = 0.0
                J[:, :, i, k] = Jii

                # partial m_i / partial deltaS
                Jii[:] = 0.0
                rain_big = np.logical_and(has_rain, big_change)
                rain_small = np.logical_and(has_rain, small_change)
                Jii[rain_big] = 1.0 - np.exp(-change[rain_big])
                Jii[rain_small] = change[rain_small] * (1.0 - 0.5 * change[rain_small])
                Jii[dead_zone] = 0.0
                J[:, :, i, k+1] = Jii


        # transformed to run in-place with one pre-allocated temporary
        dom_shape = P[:, :, 0, 0].shape
        for i in range(dom_shape[0]):
            for j in range(dom_shape[1]):
                # Ps = Js P[s,:,:] Js^T
                Ps = P[i, j, :, :]
                Js = J[i, j, :, :]
                np.dot(Js, Ps, P2)
                np.dot(P2, Js.T, Ps)

                # P[s,:,:] = Js P[s,:,:] Js^T + mQ
                Ps += mQ
                P[i, j, :, :] = Ps


    def get_state(self):
        """
        Return the current state for all grid points. READ-ONLY under normal circumstances.
        
        :return: the extended state with dimensions dom[0] x dom[1] x (k+2)
        """
        return self.m_ext


    def get_state_covar(self):
        """
        Return the state covariance for all grid points. READ-ONLY under normal circumstances.
        
        :return: the state covariance with dimensions dom[0] x dom[1] x (k+2) x (k+2)
        """
        return self.P


    def kalman_update(self, O, V, fuel_types, Kg=None):
        """
        Updates the state using the observation at the grid point.

        :param O: the observations [shape grid_size x len(fuel_types)]
        :param V: the (diagonal) variance matrix of the measurements [shape grid_size x len(fuel_types) x len(fuel_types)]
        :param fuel_types: the fuel types for which the observations exist (used to construct observation vector)
        """
        m_ext, P, H, P2 = self.m_ext, self.P, self.H, self.P2
        H[:] = 0.0

        # construct an observation operator tailored to observed fuel types
        for i in fuel_types:
            H[i, i] = 1.0
        Ho = H[fuel_types, :]

        for i, j in np.ndindex(P[:, :, 0, 0].shape):
            # re-use P2 to store forecast covariance for position s
            P2[:, :] = P[i, j, :, :]

            if np.any(np.diag(P2) < 0.0):
                print(("ERROR: negative diagonal %s" % (str(np.diag(P2)))))

            # compute Kalman gain
            I = np.dot(np.dot(Ho, P2), Ho.T) + V[i, j, :, :]
            K = np.dot(np.dot(P2, Ho.T), np.linalg.inv(I))

            if K[1, 0] > 1.0 or K[1, 0] < 0.0:
                print(("ERROR at %s, Kalman gain %g out of bounds! I=%g V=%g" % (
                    str((i, j)), K[fuel_type, 0], I, V[s[0], s[1], 0, 0])))

            # update state and state covariance
            m_ext[i, j, :] += np.dot(K, O[i, j, :] - m_ext[i, j, fuel_types])
            P[i, j, :, :] -= np.dot(np.dot(K, Ho), P2)

            if Kg is not None:
                Kg[i, j, :] = K[:, 0]


    def kalman_update_single(self, O, V, fuel_type, Kg=None):
        """
        Perform an optimized Kalman update assuming there is only one observed fuel type.
        This function is more efficient.

        :param O: the observations [shape grid_size x len(fuel_types)]
        :param V: the (diagonal) variance matrix of the measurements [shape grid_size x len(fuel_types) x len(fuel_types)]
        :param fuel_type: the fuel type for which the observations exist (used to construct observation operator)
        :param Kg: the output Kalman gain
        """
        m_ext, P, P2 = self.m_ext, self.P, self.P2

        #print("******* GRID MODEL ********")
        #print("grid: O = %.18e  V = %.18e   fuel_type = %d" % (O[0,0], V[0,0,0,0], fuel_type))

        dom_shape = P[:, :, 0, 0].shape
        for i in range(dom_shape[0]):
            for j in range(dom_shape[1]):
                P2[:, :] = P[i, j, :, :]
                #print(P2)

                if np.any(np.diag(P2) < 0.0):
                    print(("ERROR: negative diagonal %s" % (str(np.diag(P2)))))

                I = P2[fuel_type, fuel_type] + V[i, j, 0, 0]
                K = P2[:, fuel_type] / I

                if K[fuel_type] > 1.0 or K[fuel_type] < 0.0:
                    print(("WARN: at %s, Kalman gain %g out of bounds! I=%g V=%g" % (
                        str((i, j)), K[fuel_type], I, V[i, j, 0, 0])))

                m_ext[i, j, :] += K * (O[i, j, 0] - m_ext[i, j, fuel_type])

                # NOTE: this update may produce matrices that are not symmetric positive definite
                P2 -= np.dot(K[:, np.newaxis], P2[fuel_type:fuel_type + 1, :])
                P[i, j, :, :] = P2

                if Kg is not None:
                    Kg[i, j, :] = K


    def kalman_update_single2(self, O, V, fuel_type, Kg):
        """
        Perform an optimized Kalman update assuming there is only one observed fuel type.
        This function is more efficient again than kalman_update_single, since it's vectorized.

        :param O: the observations [shape grid_size x len(fuel_types)]
        :param V: the (diagonal) variance matrix of the measurements [shape grid_size x len(fuel_types) x len(fuel_types)]
        :param fuel_type: the fuel type for which the observations exist (used to construct observation operator)
        :param Kg: the output Kalman gain
        """

        m_ext, P = self.m_ext, self.P

        if np.any(O < 0.0):
            print("ERROR: found negative observations!")
            sys.exit(5)

        # I is dom_shape
        I = P[:, :, fuel_type, fuel_type] + V[:, :, 0, 0]
        # K is dom_shape x dim
        Kg[:, :, :] = P[:, :, :, fuel_type] / I[:, :, np.newaxis]

        if np.any(np.logical_or(Kg[:, :, fuel_type] > 1.0, Kg[:, :, fuel_type] < 0.0)):
            print("ERROR: some Kalman gains for 10-hr fuel are out of bounds.")

        # m_ext is dom_shape x dim, O[:,:,0] and m_ext[:,:,fuel_type] are dom_shape
        m_ext += Kg * (O - m_ext[:, :, fuel_type:fuel_type + 1])

        # try to do a tensor product here
        # Kg[:,:,:,newaxis] is shape dom_shape x dim x 1
        # P[:,:,fuel_type:fuel_type+1,:] is dom_shape x 1 x dim
        dom_shape = P[:, :, 0, 0].shape
        for i in range(dom_shape[0]):
            for j in range(dom_shape[1]):
                P[i, j, :, :] -= np.dot(Kg[i, j, :, np.newaxis], P[i, j, fuel_type:fuel_type + 1, :])


    def to_pickle(self, path):
        """
        Store the model in a cPickle file.
        
        :param path: the path where to store the model
        """
        cPickle.dump(self, path)
                

    def to_netcdf(self, path, data_vars):
        """
        Store the model in a netCDF file that attempts to be displayable
        using standard tools and loosely follows the WRF 'standard'.
        
        :param path: the path where to store the model
        :param data_vars: dictionary of additional variables to store for visualization
        """
        import netCDF4
        

        d = netCDF4.Dataset(path, 'w', format='NETCDF4')
        
        d0, d1, k = self.m_ext.shape

        d.createDimension('fuel_moisture_classes_stag', k)
        d.createDimension('south_north', d0)
        d.createDimension('west_east', d1)
        ncfmc = d.createVariable('FMC_GC', 'f4', ('south_north', 'west_east','fuel_moisture_classes_stag'))
        ncfmc[:,:,:] = self.m_ext
        logging.info('fuel_moisture_model.to_netcdf: writing extended state as FMC_GC %s covariance as FMC_COV %s' % (inq(self.m_ext),inq(self.P)))
        ncfmc_cov = d.createVariable('FMC_COV', 'f4', ('south_north', 'west_east','fuel_moisture_classes_stag', 'fuel_moisture_classes_stag'))
        ncfmc_cov[:,:,:,:] = self.P
        for v in data_vars:
            d.createVariable(v, 'f4', ('south_north', 'west_east'))[:,:]=data_vars[v]
       
        
        d.close()

    def to_geogrid(self, path, index, lats=[], lons=[]):
        """
        Store model to geogrid files
        """
        not_coord = len(lats) == 0 or len(lons) == 0
        test_latslons=True
        if not_coord:
            test_latslons=False

        logging.info("fmda.fuel_moisture_model.to_geogrid path=%s lats %s lons %s" % (path, inq(lats), inq(lons)))
        logging.info("fmda.fuel_moisture_model.to_geogrid: geogrid_index="+str(index))
        ensure_dir(path)
        xsize, ysize, n = self.m_ext.shape
        if n != 5:
            logging.error('wrong number of extended state fields, expecting 5')

        if not not_coord:
            x=int(xsize*0.5)
            y=int(ysize*0.5)
            index.update({'known_x':float(y),'known_y':float(x),'known_lat':float(lats[x-1,y-1]),'known_lon':float(lons[x-1,y-1])})
            logging.info("fmda.fuel_moisture_model.to_geogrid: geogrid updated="+str(index))

        FMC_GC = np.zeros((xsize, ysize, 5))
        FMC_GC[:,:,:3] = self.m_ext[:,:,:3]
        if test_latslons:
            logging.info("fmda.fuel_moisture_model.to_geogrid: storing lons lats to FMC_GC(:,:,4:5) to test in WRF against XLONG and XLAT")
            FMC_GC[:,:,3] = lons
            FMC_GC[:,:,4] = lats
        FMEP = self.m_ext[:,:,3:]
        
        write_geogrid_var(path,'FMC_GC',FMC_GC,index,bits=32)
        write_geogrid_var(path,'FMEP',FMEP,index,bits=32)
        
    def to_wps_format(self, path, index, lats, lons, time_tag):
        """
        Store model to wps format files
        """
        test_latslons=True

        logging.info("fmda.fuel_moisture_model.to_wps_format path=%s lats %s lons %s" % (path, inq(lats), inq(lons)))
        logging.info("fmda.fuel_moisture_model.to_wps_format: index="+str(index))
        ny, nx, n = self.m_ext.shape
        if n != 5:
            logging.error('wrong number of extended state fields, expecting 5')

        m = 7
        var = np.zeros((ny, nx, m))
        var[:,:,:3] = self.m_ext[:,:,:3]
        var[:,:,3] = lons
        var[:,:,4] = lats
        var[:,:,5:] = self.m_ext[:,:,3:]
        arrs = [var[:,:,k].swapaxes(0,1) for k in range(m)]
        startloc = "SWCORNER"
        startlat = lats[0,0]
        startlon = lons[0,0]
        params = {
            'ifv': 5, 'hdate': "{}:00:00".format(time_tag), 'xfcst': 0.,
            'map_source': "WRF-SFIRE Wildland Fire Information and Forecasting System",
            'field': ["FM1H", "FM10H", "FM100H", "FMXLON", "FMXLAT", "FMEP0", "FMEP1"], 
            'units': ["1", "1", "1", "degrees", "degrees", "1", "1"],
            'desc': ["1h Fuel Moisture Content","10h Fuel Moisture Content","100h Fuel Moisture Content",
                     "Longitude for testing", "Latitude for testing", "Drying/Wetting Equilibrium Adjustment",
                     "Rain Equilibrium Adjustment"],
            'xlvl': 200100., 'nx': nx, 'ny': ny, 'iproj': 3, 'startloc': startloc, 
            'startlat': startlat, 'startlon': startlon, 'dx': index['dx'], 'dy': index['dy'], 
            'xlonc': index['stdlon'], 'truelat1': index['truelat1'], 'truelat2': index['truelat2'], 
            'earth_radius': index['radius'], 'is_wind_earth_rel': 0, 'slab': arrs
        }
        WPSFormat.from_params(**params).to_file(osp.join(path,'FMDA:{}'.format(time_tag)))

    @classmethod
    def from_netcdf(cls, path):
        """
        Construct a fuel moisture model from data stored in a netCDF file.
        
        :param path: the path to the netCDF4 file
        """
        import netCDF4

        logging.info("reading from netCDF file " + path)
        d = netCDF4.Dataset(path)
        ncfmc = d.variables['FMC_GC'][:,:,:]

        d0, d1, k = ncfmc.shape
        P = d.variables['FMC_COV'][:,:,:,:]

        logging.info('fuel_moisture_model.from_netcdf: reading FMC_GC %s FMC_COV %s' % (inq(ncfmc),inq(P)))
        
        Tk = np.array([1.0, 10.0, 100.0]) * 3600
        
        fm = cls(ncfmc[:,:,:k-2], Tk)
        
        fm.m_ext[:,:,k-2:] = ncfmc[:,:,k-2:]
        fm.P[:,:,:,:] = P

        logging.info('fuel_moisture_model.from_netcdf: err %s' % np.max(np.abs(ncfmc[:,:,:k-2]-fm.m_ext[:,:,:k-2])))
        logging.info('fuel_moisture_model.from_netcdf: extended state fmc %d fields + 2 parameters' % (k - 2))

        return fm
        
