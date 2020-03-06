import xarray as xr
import pandas as pd
import datetime as dt
import numpy as np
import netCDF4
import os

dat = pd.read_csv('oxidized_reduced_Ndeposition.csv')
lon = np.unique(dat['longitude'][:])
lat = np.unique(dat['latitude'][:])

#Reshape columns into 2D fields
v = {}
for i in dat:
    v[i] = (('lat','lon'),dat[i].values.reshape(len(lon),len(lat)).T)
Ndep = xr.Dataset(data_vars=v, coords={'lon':lon,'lat':lat})
Ndep = Ndep.set_coords(['latitude','longitude','pixel_area_km2'])

#Find (unique) years there is data
t = []
for i in Ndep.keys():
    if int(i[-4:]) not in t:
        t.extend([int(i[-4:])])
time_units = 'days since 1984-01-01'
time = netCDF4.date2num([dt.datetime(x,1,1) for x in t],time_units)

#Combine years for each variable
v = {}
for s in set([x[:-5] for x in Ndep.keys()]):
    var = np.zeros((len(time),len(lat),len(lon)))
    for i,k in enumerate(t):
        var[i,:,:] = Ndep[s+'_'+str(k)].values
    v[s] = (('time','lat','lon'),var,{'units':'kg / km^2 / y^1'})
v['area'] = (('lat','lon'),Ndep['pixel_area_km2'].values,{'units': 'km^2'})
Ndep = xr.Dataset(data_vars=v, coords={'time':time,'lon':lon,'lat':lat})

#Mean in time
Ndep = Ndep.mean('time')

#Interpolate onto grid
nc = netCDF4.Dataset('domain_cfg.nc')
Ndep = Ndep.interp(lat=np.mean(nc.variables['nav_lat'][:],axis=1),
                   lon=np.mean(nc.variables['nav_lon'][:],axis=0))
nc.close()

# Create NETCDF file
dims = {'y': 554, 'x': 684, 'time_counter': None}
vars = [{'name':'N3_flux','type':'double',
           'dims': ('time_counter','y','x'),
           'units':"mmol*m-2*s-1"},
        {'name':'N4_flux','type':'double',
           'dims': ('time_counter','y','x'),
           'units':"mmol*m-2*s-1"},
        {'name':'time_counter','type':'double',
           'dims': ('time_counter'),
           'units':"unit"}]
nco = netCDF4.Dataset('Ndep.nc','w')
for d in dims:
    nco.createDimension(d,dims[d])
for v in vars:
    dat=nco.createVariable(v['name'],v['type'],v['dims'])
    dat.units = v['units']

nco.author = os.getenv('USER')
nco.history = dt.datetime.now().strftime("Created on %a, %B %d, %Y at %H:%M")
nco.title = 'Nitrogen Deposition from University of Minnesota'

# Set variables
nco.variables['time_counter'][:] = 1
nco.variables['N3_flux'][:] = Ndep['total_oxidized'].values[np.newaxis,:,:]/(14*365*86400)
nco.variables['N4_flux'][:] = Ndep['total_reduced'].values[np.newaxis,:,:]/(14*365*86400)
nco.close()

