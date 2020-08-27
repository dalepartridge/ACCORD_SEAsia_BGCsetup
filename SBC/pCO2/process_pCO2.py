'''
Script to convert the global RCP8.5 CO2 combined hindcast/forecast file into
a NEMO surface forcing file

Inputs
------
nx/ny - integer, domain size
ys/e - integer, desired start and end year. Note the script adds one
                year to the start and end to ensure full coverage
domnam - str, short tag name for the domain
f_in - str, file path to the RCP8.5 dat file
f_out - str, output filename

'''

import pandas as pd
import datetime as dt
import numpy as np
import netCDF4
import os

#######################################
# USER INPUTS
ny=324
nx=454
ys=1990
ye=2100
domnam='ACCORD'
f_in='RCP85_MIDYR_CONC.DAT'
f_out='pCO2a.nc'
#######################################

# Read Data
dat = pd.read_csv(f_in, sep='\s+', skiprows=38, usecols=[0,3])
sl = np.s_[ys-1-1765:ye+2-1765] #Find year indexes
time = np.array([dt.datetime(y,7,1) for y in dat['YEARS'][sl]])
co2 = np.ones((ye-ys+3,ny,nx))
for i,k in enumerate(np.array(dat['CO2'][sl])):
    co2[i,:] *= k

# Create NETCDF file
dims = {'y': ny, 'x': nx, 't': None}
vars = [{'name':'pCO2a','type':'double',
           'dims': ('t','y','x'),
           'units':"ppm"},
        {'name':'t','type':'double',
           'dims': ('t'),
           'units':"days since 1970-01-01 00:00:00"}]
nco = netCDF4.Dataset(f_out,'w')
for d in dims:
    nco.createDimension(d,dims[d])
for v in vars:
    dat=nco.createVariable(v['name'],v['type'],v['dims'])
    dat.units = v['units']

nco.author = os.getenv('USER')
nco.history = dt.datetime.now().strftime("Created on %a, %B %d, %Y at %H:%M")
nco.title = 'Atmospheric partial pressure of CO2 for '+domnam+' domain'
nco.comment = 'Global annual values taken from climate scenaris RCP8.5'

# Set variables
nco.variables['t'][:] = netCDF4.date2num(time,nco.variables['t'].units)
nco.variables['pCO2a'][:] = co2
nco.close()

