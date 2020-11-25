import sys
import netCDF4
import numpy as np

sfile = sys.argv[1]
v = sys.argv[2]
maskfile = sys.argv[3]
maskvar='mask'

ncm = netCDF4.Dataset(maskfile)
mask = np.abs(np.squeeze(ncm.variables[maskvar][:])-1)
ncm.close()

nc = netCDF4.Dataset(sfile,'a')
for i in range(nc.variables[v].shape[0]):
    print('Filling mask for record {} of {}'.format(i+1,nc.variables[v].shape[0]))
    d = np.ma.masked_array(np.squeeze(nc.variables[v][i,:]),mask=mask)
    d[d.mask] = 0
    nc.variables[v][i,:] = d[np.newaxis,:]
nc.close()
