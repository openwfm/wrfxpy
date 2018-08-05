import numpy as np
import logging

def interpolate2height_old(var,height,level):
      """
      Interpolate 3d variable to a given height
      :param var: the 3d array to be interpolated, 1st axis is vertical
      :param height: the 3d array of heights of the nodes on which v is defined
      :param level: the target height to interpolate to
      :return: interpolated value, or
      """
      maxlayer = var.shape[0]-1
      z = np.zeros([var.shape[1],var.shape[2]])
      r = z
      for i in range(0, var.shape[1]):
          for j in range(0, var.shape[2]):
             k = np.searchsorted(height[:,i,j],level)
             if k==0 or k>maxlayer:
                 r[i,j] = np.nan
                 #logging.error("Need height[0,%s,%s]=%s < level=%s <= height[%s,%s,%s]=%s" \
                 #   % (i,j,height[0,i,j],level,maxlayer,i,j,height[maxlayer,i,j]))
                 #return z
             else:
                 # interpolate in the interval height[k-1,i,j] to height[k,i,j]
                 r[i,j]=var[k-1,i,j]+(var[k,i,j]-var[k-1,i,j]) \
                     * (level - height[k-1,i,j])/(height[k,i,j] - height[k-1,i,j])
      return r

def interpolate2height(var,height,level):
      """
      Interpolate 3d variable to a given height
      :param var: the 3d array to be interpolated, 1st axis is vertical
      :param height: the 3d array of heights of the nodes on which v is defined
      :param level: the target height to interpolate to
      :return: interpolated value, or
      """
      ix, tx = index8height(height,level)
      r = np.zeros([var.shape[1],var.shape[2]])
      maxlayer = var.shape[0]-1
      for i in range(0, var.shape[1]):
          for j in range(0, var.shape[2]):
             k = ix[i,j]
             t = tx[i,j]
             if k==0 or k>maxlayer:
                 #r[i,j] = np.nan
                 r[i,j] = 0
             else:
                 #r[i,j]=var[k,i,j]+(var[k+1,i,j]-var[k,i,j])*tx[i,j] 
                 r[i,j] = var[k,i,j]*(1.0-t) + var[k+1,i,j]*t 
      return r

def index8height(height,level):
      """
      Find index and fraction at given height
      :param height: the 3d array of heights of the nodes on which v is defined
      :param level: the target height to interpolate to
      :return: index of level as integer part and fractional part
      If level is <= smallest height or > largest g=height, return all zeros
      """
      maxlayer = height.shape[0]-1
      iz = np.zeros([height.shape[1],height.shape[2]],dtype=np.int_)
      tz = np.zeros([height.shape[1],height.shape[2]])
      ix = np.copy(iz)
      tx = np.copy(tz)
      for i in range(0, height.shape[1]):
          for j in range(0, height.shape[2]):
             k = np.searchsorted(height[:,i,j],level)
             if k==0 or k>maxlayer:
                 logging.error("Need height[0,%s,%s]=%s < level=%s <= height[%s,%s,%s]=%s" \
                    % (i,j,height[0,i,j],level,maxlayer,i,j,height[maxlayer,i,j]))
                 return iz, tz
             ix[i,j]=k-1    # integer part
             # interpolation in the interval height[k-1,i,j] to height[k,i,j]
             tx[i,j]= (level - height[k-1,i,j])/(height[k,i,j] - height[k-1,i,j])
      return ix,tx 

def pressure(d,t):
      """
      Compute pressure height of mesh centers
      :param d: open NetCDF4 dataset
      :param t: number of timestep
      """
      return d.variables['P'][t,:,:,:] + d.variables['PB'][t,:,:,:]

def pressure8w(d,t):
      """
      Compute pressure height at mesh cell bottoms a.k.a. w-points 
      :param d: open NetCDF4 dataset
      :param t: number of timestep
      """
      ph = pressure(d,t)
      ph8w = ph[0:ph.shape[0]-1,:,:]  
      # average from 2nd layer up 
      ph8w[1:,:,:] = 0.5*(ph[1:,:,:] + ph[0:ph.shape[0]-1,:,:])
      return ph8w 

def u8p(d,t):
      """
      Compute horizontal wind u at mesh cell centers a.k.a. p-points
      :param d: open NetCDF4 dataset
      :param t: number of timestep
      """
      u = d.variables['U'][t,:,:,:]
      return 0.5*(u[:,:,0:u.shape[2]-1]+u[:,:,1:])

def v8p(d,t):
      """
      Compute horizontal wind v at mesh cell centers a.k.a. p-points
      :param d: open NetCDF4 dataset
      :param t: number of timestep
      """
      v = d.variables['V'][t,:,:,:]
      return 0.5*(v[:,0:v.shape[1]-1,:]+v[:,1:,:])

def w8p(d,t):
      """
      Compute vertical wind w at mesh cell centers a.k.a. p-points
      :param d: open NetCDF4 dataset
      :param t: number of timestep
      """
      w = d.variables['W'][t,:,:,:]
      return 0.5*(w[0:w.shape[0]-1,:,:]+w[1:,:,:])

def height8w(d,t):
      """
      Compute height at mesh bottom a.k.a. w-points 
      :param d: open NetCDF4 dataset
      :param t: number of timestep
      """
      ph = d.variables['PH'][t,:,:,:]  
      phb = d.variables['PHB'][t,:,:,:]
      return (phb + ph)/9.81 # geopotential height at W points

def height8p(d,t):
      """
      Compute height of mesh centers (p-points)
      :param d: open NetCDF4 dataset
      :param t: number of timestep
      """
      z8w = height8w(d,t)
      return 0.5*(z8w[0:z8w.shape[0]-1,:,:]+z8w[1:,:,:])

def height8p_terrain(d,t):
      """
      Compute height of mesh centers (p-points) above terrain
      :param d: open NetCDF4 dataset
      :param t: number of timestep
      """
      z8w = height8w(d,t)
      h =  0.5*(z8w[0:z8w.shape[0]-1,:,:]+z8w[1:,:,:])
      for i in range(0, h.shape[2]):
          for j in range(0, h.shape[1]):
              h[:,i,j] -= z8w[0,i,j]
      return h

def wcloud(d,t):
      """
      Compute the cloud water mass intensity in mesh cells (kg/m^2)
      :param d: open NetCDF4 dataset
      :param t: number of timestep
      :return: cloud water mass in mesh cells per surface area (kg/m^2)
      """
      P = d.variables['P'][t,:,:,:]   # dry air pressure (Pa )
      T = d.variables['T'][t,:,:,:] + d.variables['T00'][t]  # temperature (K)
      r_d = 287                       # specific gas constant (J/kg/K)
      rho = P/(r_d * T)               # dry air density  (kg/m^3)
      z8w = height8w(d,t)             # height of mesh bottom (m)
      dz = z8w[1:,:,:]-z8w[0:z8w.shape[0]-1,:,:] # mesh cell heights (m)
      qcloud = d.variables['QCLOUD'][t,:,:,:] # cloud water mixing ratio (kg water/kg dry air)
      return rho * qcloud * dz        # cloud water mass intensity kg/m^2
      
 

