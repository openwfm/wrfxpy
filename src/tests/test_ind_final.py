import netCDF4 as nc4
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import numpy as np
from clamp2mesh import interpolate_coords
import sys

plot = True

nc = nc4.Dataset('../WRF4/WPS-subgrid_shift/geo_em_ind/geo_em.d01_final.nc')
srx = nc.getncattr('sr_x')
sry = nc.getncattr('sr_y')
X = np.array(nc['XI'][0])
Y = np.array(nc['YI'][0])
n,m = X.shape
X[X > m+1] = np.nan
Y[Y > n+1] = np.nan
Xs = np.array(nc['XIS'][0])
Ys = np.array(nc['YIS'][0])
Xs[Xs > m+1] = np.nan
Ys[Ys > n+1] = np.nan

ll = 10
label_ticks = np.arange(-1,ll,2)
fig, axs = plt.subplots(1,2,constrained_layout=True,figsize=(15,15))
im = axs[0].imshow(X[:ll,:ll],origin='lower',interpolation='None')
for i in range(ll):
     for j in range(ll):
              text = axs[0].text(i, j, round(X[j,i],4),
                	ha="center", va="center", color="r", fontsize=8,fontweight='bold')
axs[0].set_title('Indexes X')
axs[0].set_xlabel("i")
axs[0].set_ylabel("j")
axs[0].set_xticklabels(label_ticks)
axs[0].set_yticklabels(label_ticks)
im = axs[1].imshow(Y[:ll,:ll],origin='lower',interpolation='None')
for i in range(ll):
     for j in range(ll):
              text = axs[1].text(i, j, round(Y[j,i],4),
                        ha="center", va="center", color="r", fontsize=8,fontweight='bold')
axs[1].set_title('Indexes Y')
axs[1].set_xlabel("i")
axs[1].set_ylabel("j")
axs[1].set_xticklabels(label_ticks)
axs[1].set_yticklabels(label_ticks)
fig.suptitle('WPS first %dx%d array' % (ll,ll))
fig, axs = plt.subplots(1,2,constrained_layout=True,figsize=(15,15))
im = axs[0].imshow(Xs[:ll,:ll],origin='lower',interpolation='None')
for i in range(ll):
     for j in range(ll):
              text = axs[0].text(i, j, round(Xs[j,i],4),
                	ha="center", va="center", color="r", fontsize=8,fontweight='bold')
axs[0].set_title('Indexes X')
axs[0].set_xlabel("i")
axs[0].set_ylabel("j")
axs[0].set_xticklabels(label_ticks)
axs[0].set_yticklabels(label_ticks)
im = axs[1].imshow(Ys[:ll,:ll],origin='lower',interpolation='None')
for i in range(ll):
     for j in range(ll):
              text = axs[1].text(i, j, round(Ys[j, i],4),
                        ha="center", va="center", color="r", fontsize=8,fontweight='bold')
axs[1].set_title('Indexes Y')
axs[1].set_xlabel("i")
axs[1].set_ylabel("j")
axs[1].set_xticklabels(label_ticks)
axs[1].set_yticklabels(label_ticks)
fig.suptitle('WPS sugbrid first %dx%d array' % (ll,ll))
plt.show()

Xo,Yo = np.meshgrid(range(1,n+1),range(1,m+1))
Xd = X-Xo
Yd = Y-Yo
print('Differences x-direction min=%f max=%f abs_max=%f' % (np.nanmin(Xd),np.nanmax(Xd),np.nanmax(abs(Xd))))
print('Differences y-direction min=%f max=%f abs_max=%f' % (np.nanmin(Yd),np.nanmax(Yd),np.nanmax(abs(Yd))))

Xsi,Ysi = interpolate_coords(X,Y,srx,sry)

kx = int(np.ceil(srx/2.))
ky = int(np.ceil(sry/2.))
Xds = np.ma.masked_invalid(Xs[ky:-4*ky,kx:-4*kx]-Xsi[ky:-4*ky,kx:-4*kx])
Yds = np.ma.masked_invalid(Ys[ky:-4*ky,kx:-4*kx]-Ysi[ky:-4*ky,kx:-4*kx])
print('Differences subgrid x-direction min=%f max=%f abs_max=%f' % (np.nanmin(Xds),np.nanmax(Xds),np.nanmax(abs(Xds))))
print('Differences subgrid y-direction min=%f max=%f abs_max=%f' % (np.nanmin(Yds),np.nanmax(Yds),np.nanmax(abs(Yds))))

if plot:
    cmap = cm.coolwarm
    cmap.set_bad('white',1.)
    n,m = Xd.shape
    x,y = np.meshgrid(range(1,n+1),range(1,m+1))
    fig = plt.figure(figsize=(15,15))
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax.plot_surface(x,y,Xd,cmap=cmap)
    ax.set_xlabel("i")
    ax.set_ylabel("j")
    ax.set_zlabel("Differences")
    ax.set_title("X indexes")
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    ax.plot_surface(x,y,Yd,cmap=cmap)
    ax.set_xlabel("i")
    ax.set_ylabel("j")
    ax.set_zlabel("Differences")
    ax.set_title("Y indexes")
    fig.suptitle("Difference in the atmospheric grid indexes")

    ns,ms = Xds.shape
    xs,ys = np.meshgrid(range(1,ns+1),range(1,ms+1))
    fig = plt.figure(figsize=(15,15))
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax.plot_surface(xs,ys,Xds,cmap=cmap)
    ax.set_xlabel("i")
    ax.set_ylabel("j")
    ax.set_zlabel("Differences")
    ax.set_title("X indexes")
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    ax.plot_surface(xs,ys,Yds,cmap=cmap)
    ax.set_xlabel("i")
    ax.set_ylabel("j")
    ax.set_zlabel("Differences")
    ax.set_title("Y indexes")
    fig.suptitle("Difference in the subgrid indexes")
    
    plt.show()
