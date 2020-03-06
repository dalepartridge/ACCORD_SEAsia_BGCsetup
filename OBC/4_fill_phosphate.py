import netCDF4
import glob
from seapy.progressbar import progress
import numpy as np

files = glob.glob('hadgem-data/accord_bdytrc*.nc')

for f in progress(files):
    nc = netCDF4.Dataset(f,'a')
    v = nc.variables['nitrate'][:]
    for i in range(75):
        x = np.percentile(v[0,i,0,:],5)
        v[0,i,0,:][v[0,i,0,:]<x] = x
    nc.variables['nitrate'][:] = v
    nc.variables['phosphate'][:] = v/16.0
    v = nc.variables['silicate'][:]
    for i in range(75):
        x = np.percentile(v[0,i,0,:],5)
        v[0,i,0,:][v[0,i,0,:]<x] = x
    nc.variables['silicate'][:] = v
    v = nc.variables['DIC'][:]
    for i in range(75):
        x = np.percentile(v[0,i,0,:],5)
        v[0,i,0,:][v[0,i,0,:]<x] = x
    nc.variables['DIC'][:] = v
    v = nc.variables['oxygen'][:]
    for i in range(75):
        x = np.percentile(v[0,i,0,:],5)
        v[0,i,0,:][v[0,i,0,:]<x] = x
    nc.variables['oxygen'][:] = v
    v = nc.variables['TA'][:]
    for i in range(75):
        x = np.percentile(v[0,i,0,:],5)
        v[0,i,0,:][v[0,i,0,:]<x] = x
    nc.variables['TA'][:] = v
    nc.close()

