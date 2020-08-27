'''
Script to calculate and fill river files with oxygen, DIC and total alkalinity

Oxygen - Values calculated at saturation as function of temperature based on 
         formula from https://www.waterontheweb.org/under/waterquality/oxygen.html
DIC - Values calculated using temperature, zero salinity, alkalinity and pCO2 
      based upon calculation in ERSEM carbonate module (ICALC=4)
TA - Set to a constant mmol / m3 value of 1500, based on GEMS/GLORI project
'''

import netCDF4
import numpy as np
import datetime as dt

ncd = netCDF4.Dataset('domain_cfg.nc')
SA = ncd.variables['e1t'][:]*ncd.variables['e2t'][:] 
lon = ncd.variables['nav_lon'][:]
lat = ncd.variables['nav_lat'][:]
ncd.close()

# Define year to calculate for / begin loop over years
year=2018 
#for year in np.arange(syear,eyear+1): #Loop until end of script

#################### Load air temperature ##############################
nct = netCDF4.Dataset('/work/n01/n01/annkat/SEAsia_HadGEM_R12/SURFACE_FORCING/HAD_INT_6hr_T150_y{}.nc'.format(year))
#Ttime = netCDF4.num2date(nct.variables['time'][:],nct.variables['time'].units) # For some reason this doesnt work
Ttime = np.array([dt.datetime(year,1,1) + dt.timedelta(days=i) for i in nct.variables['time'][:]])
Tfull = nct.variables['T150'][:]
Tlon,Tlat = np.meshgrid(nct.variables['lon'][:],nct.variables['lat'][:])
nct.close()

# Convert Temperatures to monthly averages
l = np.array([i.month for i in Ttime])
T = np.zeros((12,Tfull.shape[1],Tfull.shape[2]))
for k in range(12):
    T[k,:,:] = np.mean(np.squeeze(Tfull[np.where(l==k+1),:,:]),axis=0)

#################### Load river runoff #############################
ncr = netCDF4.Dataset('rivers_y{}.nc'.format(year),'a')
r = ncr.variables['rorunoff'][:]

# Create new variables
rdim = ncr.variables['rorunoff'].dimensions
ncr.createVariable('DIOrunoff','double',rdim)
ncr.createVariable('DICrunoff','double',rdim)
ncr.createVariable('TArunoff','double',rdim)

#################### Dissolved Oxygen ##############################
dat = np.zeros(r.shape)
for j in range(r.shape[1]):
    for i in range(r.shape[2]):
        if r[0,j,i] > 0:
            c = np.maximum(np.abs(Tlon-lon[j,i]),np.abs(Tlat-lat[j,i]))
            ll_idx = np.where(c == np.min(c))
            Tr = T[:,ll_idx[0][0],ll_idx[1][0]].ravel()
            dat[:,j,i] = ((np.exp(7.7117-1.31403*np.log(Tr+45.93))) * (1-np.exp(11.8571-(3840.7/(Tr+273.15)) - \
                          (216961/((Tr+273.15)**2)))) * (1-(0.000975-(0.00001426*Tr)+0.00000006436*(Tr**2)))) / \
                          (1-np.exp(11.8571-(3840.7/(Tr+273.15))-(216961/((Tr+273.15)**2)))) / \
                          (1-(0.000975-(0.00001426*Tr)+(0.00000006436*(Tr**2)))) # Dissolved oxygen at saturation (mg/L)

dat *= (1000/32) # Convert to mmol/m3
dat *= (r/1000) * SA # Convert to mmol/s

ncr.variables['DIOrunoff'][:] = dat


###################### DIC ########################################
def K (T):
    '''
    Function to calculate carbonate equilibirum constants at surface with S=0
    '''
    K0 = np.exp(93.4517*(100/T) - 60.2409 + 23.3585*np.log(T/100))
    K1 = 10**(126.34048 - 6320.813*(1/T) - 19.568224*np.log(T))  
    K2 = 10**(90.18333 - 5143.692*(1/T) - 14.613358*np.log(T))
    return K0, K1, K2

#Load pCO2
ncp = netCDF4.Dataset('/work/n01/n01/dapa/ACCORD/BGC-INPUTS/SBC/pCO2_y2018.nc')
pco2 = ncp.variables['pCO2a'][year-1979,:,:]
ncp.close()

TA = 1500

dat = np.zeros(r.shape)
for j in range(r.shape[1]):
    for i in range(r.shape[2]):
        if r[0,j,i] > 0:
            c = np.maximum(np.abs(Tlon-lon[j,i]),np.abs(Tlat-lat[j,i]))
            ll_idx = np.where(c == np.min(c))
            Tr = T[:,ll_idx[0][0],ll_idx[1][0]].ravel()
            k0,k1,k2 = K(Tr)
            Pr = (k1/k2)*k0*pco2[j,i]
            x = np.sqrt((8*TA+Pr)*Pr)
            dat[:,j,i] = TA/2 - Pr/8 + x/8 + pco2[j,i]*k0 

dat *= (r/1000) * SA # Convert to mmol/s
ncr.variables['DICrunoff'][:] = dat

################ Total Alkalinity ##################################
# Alkalinity set to constant value

dat = TA # mmol/m3
dat *= (r/1000) * SA # Convert to mmol/s
ncr.variables['TArunoff'][:] = dat

ncr.close()

