'''
Copyright 2018 Yuri Artioli, Plymouth Marine Laboratory 

Permission is hereby granted, free of gridge, to any person obtaining a copy of 
this software and associated documentation files (the "Software"), to deal in 
the Software without restriction, including without limitation the rights to 
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
of the Software, and to permit persons to whom the Software is furnished to do 
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
SOFTWARE.
'''


import matplotlib
matplotlib.use('agg')

from netCDF4 import Dataset as DS
from numpy import zeros, arange, array,  where, abs, mod, minimum, linspace, dot, maximum, minimum,savetxt,load,save,float32
from scipy.interpolate import interp2d,interp1d,griddata
import seawater as SEA
from numpy.linalg import inv as INVMAT
from scipy.sparse import diags as DIAG
from matplotlib.pyplot import figure, plot, close,title,suptitle

try:
    from joblib import delayed,Parallel
    force_serial=False
except:
    print('joblib module not present, turning off parallel parts of the code')
    force_serial=True
import yaml as Y
import sys

def define_bdy(coordsfile,meshfile,first_bdylayer_only,grid):
    # this function initialises the BDY and store all important info in a dictionnary
    BDY={}
    
    # gather the initial cell centre depth of the domain grid points
    fmesh=DS(meshfile)
    depth_0=fmesh.variables['gdept'][:].squeeze()
    depth_levels=depth_0.shape[0]
    cell_thickness=fmesh.variables['e3'+grid.lower()][:].squeeze()
    fmesh.close()
        
    # extract the lat, lon coordinates of the boundary from the bdy_coordinate file
    # bdy_i, bdy_j are the grid indices of the boundary points
    # bdy_lat, bdy_lon are the coorindates of the boundary points
    # bdy_r is the "r" index in the bdy, i.e. the ordinal number of the bdy layer (generally the baoundary is imposed on a band along the bpoundary that is R cells thick
    fcoords=DS(coordsfile)
    BDY['nbj']=fcoords.variables['nbj'+grid.lower()][:].squeeze()-1
    BDY['nbi']=fcoords.variables['nbi'+grid.lower()][:].squeeze()-1
    BDY['lat']=fcoords.variables['gphi'+grid.lower()][:].squeeze()
    BDY['lon']=fcoords.variables['glam'+grid.lower()][:].squeeze()
    BDY['nbr']=fcoords.variables['nbr'+grid.lower()][:].squeeze()

    if first_bdylayer_only:
        end_of_first_bdylayer=(BDY['nbr']==1).sum()
        BDY['nbj']=BDY['nbj'][:end_of_first_bdylayer]
        BDY['nbi']=BDY['nbi'][:end_of_first_bdylayer]
        BDY['lat']=BDY['lat'][:end_of_first_bdylayer]
        BDY['lon']=BDY['lon'][:end_of_first_bdylayer]
    fcoords.close()
    
    #extract depth for all bdy points:
    #bdy_depth=zeros((depth_levels,BDY['nbj'].size))
    #bdy_cell_thickness=zeros((depth_levels,BDY['nbj'].size))
    #for position in range(BDY['nbj'].size):
    #    bdy_depth[:,position]=depth_0[:,BDY['nbj'][position],BDY['nbi'][position]]
    #    bdy_cell_thickness[:,position]=cell_thickness[:,BDY['nbj'][position],BDY['nbi'][position]]
    BDY['depth']=depth_0.data#bdy_depth
    BDY['e3']=cell_thickness.data#bdy_cell_thickness
    return BDY

def parallel_extraction(ncfilename,varname,BDY,scaling,bias,surface,ncfilename_surface,time):
    ncfile=DS(ncfilename)
    tmp_var=ncfile.variables[varname][time]
    var=tmp_var[:,BDY['nbj'],BDY['nbi']]*scaling-bias
    ncfile.close()
    print (time)
    if surface:
        ncfile_surface=DS(ncfilename_surface)
        tmp_svar=ncfile_surface.variables[varname][time]
        surface_var=tmp_svar[BDY['nbj'],BDY['nbi']]*scaling-bias
        ncfile_surface.close()
        return time,var,surface_var
    else:
        return time,var

def extract_values(BDY,varfile,varname,surface=False,surface_varfile='',scaling=1.,bias=0.,fill_mask=False):
    # this function extracte the variable "varname" from the "varfile" netcdf file uding the BDY metadata in the BDY dictionnary
    # if the surface monthly values need to be extracted then the logical switch needs to be set to true and the name of the netcdf file with the surface monthly value need to be read.
    # it also extract the original depth information
    
    # open the files and check what is the name of the depth variable
    fvar=DS(varfile)
    if 'lev' in fvar.variables.keys():
        zvar='lev'
    elif 'depth' in fvar.variables.keys():
        zvar='depth'
    else: 
        print ('No recognised depth array in the file')
        zvar=''
    if surface:
        surface_fvar=DS(surface_varfile)
    else:
        #this is needed because in case of parallel extraction surface_fvar is called, but not used
        surface_fvar=None
    
    #extract variable
    # if the length of the time dimension is big (i.e. monthly), then extract data one timstep at time to save memory
    print ('extracting values for '+varname)
    if len(fvar.dimensions['time'])<300:
        tmp_var=fvar.variables[varname][:].squeeze()
        if len(tmp_var.shape)==3: tmp_var=array([tmp_var]) #this is to add a fictional temporal dimension in case of the variable has only one time record (long term climatological mean)
        if zvar!='': var=tmp_var[:,:,BDY['nbj'],BDY['nbi']]*scaling-bias  #this is to extract main variables that are not depth resolved, like ssh
        else: var=tmp_var[:,BDY['nbj'],BDY['nbi']]*scaling-bias
            
        if surface:
            tmp_svar=surface_fvar.variables[varname][:].squeeze()
            surface_var=tmp_svar[:,BDY['nbj'],BDY['nbi']]*scaling-bias
            surface_fvar.close()
    else:
        timelen=len(fvar.dimensions['time'])
        if zvar!='':
            nlev=len(fvar.dimensions[zvar])
            var=zeros((timelen,nlev,BDY['nbj'].size))
        else:
            var=zeros((timelen,BDY['nbj'].size))
        if surface: surface_var=zeros((timelen,BDY['nbj'].size))
        print ('given the big size the extraction is one timestep at time')
        print ('this might take few minutes')
        for time in range(timelen):
            tmp_var=fvar.variables[varname][time]
            if zvar!='':
                var[time]=tmp_var[:,BDY['nbj'],BDY['nbi']]*scaling-bias
            else:
                var[time]=tmp_var[BDY['nbj'],BDY['nbi']]*scaling-bias
            if surface:
               tmp_svar=surface_varfile.variables[varname][time]
               surface_var[time]=tmp_svar[BDY['nbj'],BDY['nbi']]*scaling-bias
            if mod(time,50)==0:
                print ('%4i out of %4i'%(time,timelen))
        if surface: surface_fvar.close()
        # TODO this need to be coded to speed up
        #else:
        #    print ('given the big size the extraction is one timestep at time in parallel')
        #    messy_out=Parallel(n_jobs=njobs,verbose=1)#(delayed(parallel_extraction)(varfile,varname,BDY,scaling,bias,surface,surface_varfile,time) for time in range(timelen))
        #    for mess in messy_out:
        #        var[mess[0]]=mess[1]
        #    if surface:
        #        surface_var[mess[0]]=mess[2]          
            
    # extract depth information
    if zvar!='':
        original_depth=fvar.variables[zvar][:]
    else:
        original_depth=0.
    fvar.close()
    #try:
        #bdy_is_masked=var.mask.any()
    #except:
        #bdy_is_masked=False
    #if bdy_is_masked:
        #print ('filling regional seapoint in the ESM land, CAREFULE does not work from nbr>1')
        #for time in range(var.shape[0]):
            #if zvar!='':
                #for z in range(var.shape[1]):
                    #if ~(var[time,z].mask.all()):
                             #xin=BDY['nbi'][~var[time,z].mask.ravel()]
                             #yin=BDY['nbj'][~var[time,z].mask.ravel()]
                             #vin=var[time,z][~var[time,z].mask.ravel()]
                             ##print z ,xin.shape,yin.shape,vin.shape
                             #var[time,z]=griddata((xin,yin),vin,(BDY['nbi'],BDY['nbj']),method='nearest')
                #if surface:
                    #if ~(surface_var[time].mask.all()):
                             #xin=BDY['nbi'][~surface_var[time].mask.ravel()]
                             #yin=BDY['nbj'][~surface_var[time].mask.ravel()]
                             #vin=surface_var[time][~surface_var[time].mask.ravel()]
                             #surface_var[time]=griddata((xin,yin),vin,(BDY['nbi'],BDY['nbj']),method='nearest')               
            #else:
                #if ~(var[time].mask.all()):
                         #xin=BDY['nbi'][~var[time].mask.ravel()]
                         #yin=BDY['nbj'][~var[time].mask.ravel()]
                         #vin=var[time][~var[time].mask.ravel()]
                         #var[time]=griddata((xin,yin),vin,(BDY['nbi'],BDY['nbj']),method='nearest')               
        #var.mask=False
        #if surface: surface_var.mask=False

    if surface:
        return var,original_depth,surface_var
    else:
        return var,original_depth
    
def depth_interpolation(variable,original_depth,new_depth,remove_seasonal_cycle_at_depth=False):
    # this function linearly interpolates all layer in the depth point of the AMM boundary
    # variable is the original 3D array (time,depth, position)
    # original_depth is the depth in the original variable
    # new depths are taken from the class attribute  bdy_depth
    var_interpolated=zeros((variable.shape[0],new_depth.shape[0],variable.shape[2]))
    #this put a flag where the script extrapolate at depth 
    if remove_seasonal_cycle_at_depth:
      for position in range(variable.shape[2]):
        interpolator=interp1d(original_depth,variable[:,:,position],axis=1,bounds_error=False,fill_value=(variable[:,0,position],variable[:,-1,position].mean(0)))
        var_interpolated[:,:,position]=interpolator(new_depth[:,position])
    else:
      for position in range(variable.shape[2]):
        interpolator=interp1d(original_depth,variable[:,:,position],axis=1,bounds_error=False,fill_value=(variable[:,0,position],variable[:,-1,position]))
        var_interpolated[:,:,position]=interpolator(new_depth[:,position])
    return var_interpolated

def calculate_MLD(S,T,z,DT=0.2,ref=10.,njobs=0):
    # function to identify the index of the cell along the depth axis that is closer to the MLD calculated following Kara definition
    
    # read the time and depth dimension and create array of outputs
    dim0=T.shape[0]
    dim2=T.shape[2]
    mld_index=zeros((dim0,dim2),dtype='int64')
    
    if njobs<=0:
        for position in range(dim2):
            # calculate_single_MLD returns two value, here we only need the first because it's being called serially
            mld_index[:,position]=calculate_single_MLD(position,S[:,:,position],T[:,:,position],z,DT=DT,ref=ref)[0]
    else:
        messy_out=Parallel(n_jobs=njobs,verbose=1)(delayed(calculate_single_MLD)(position,S[:,:,position],T[:,:,position],z,DT=DT,ref=ref) for position in range(dim2))
        for mess in messy_out:
            mld_index[:,mess[1]]=mess[0]

    return  mld_index

def calculate_single_MLD(position,S,T,z,DT=0.2,ref=10.):
    # index of the reference depth
    iref=abs(z[:,position]-ref).argmin(0) 
    # this calculate the MLD in a specific BDY position during the whole period
    # calculate density
    dens=SEA.dens(S.ravel(),T.ravel(),zeros(T.size)).reshape(S.shape)
    
    #calculate the two possible threshold density by increasing (plus) or decreasing (minus) 
    dens_minus=SEA.dens(S[:,iref].ravel(),(T[:,iref]-DT).ravel(),zeros(T[:,iref].size))
    dens_plus=SEA.dens(S[:,iref].ravel(),(T[:,iref]+ DT).ravel(),zeros(T[:,iref].size))
    
    # pick the density treshold that is higher of the density at reference depth, because at low temperature the sign flips
    dens_delta=where(dens[:,iref]>dens_minus,dens_plus,dens_minus)
    mld_index=abs(dens.T-dens_delta[:]).argmin(0)
    
    # set the lower limit to the reference depth
    mld_index=maximum(mld_index,iref)
    
    if mod(position,100)==0:
        print ('%4i out of %4i'%(position,z.shape[1]))
        
    return mld_index,position


def calculate_surface_anomalies(surface_monthly_var,method='stretch'):
    nyears=int(surface_monthly_var.shape[0]/12)
    npositions=surface_monthly_var.shape[1]
    annual_surface_mean=surface_monthly_var.reshape(nyears,12,npositions).mean(1)
    years=arange(0.5,nyears)
    months=arange(0,nyears,1/12.)
    if method=='interp':
        time_interp=interp1d(years,annual_surface_mean,axis=0,bounds_error=False,fill_value=(annual_surface_mean[0],annual_surface_mean[-1]))
        annual_mean_at_monthly_resolution=time_interp(months)
    elif method=='stretch':
        annual_mean_at_monthly_resolution=annual_surface_mean.repeat(12,axis=0)
    elif method=='killworth':
        print ('method to be re-tested')
        #annual_mean_at_monthly_resolution=killworth_interpolation(arange(0.5,nyears),annual_mean,arange(0,nyears,1/12.),z=arange(nhor))
    else:
        print ('Method %s is not recognised'%method)
        annual_mean_at_monthly_resolution=zeros(months.size)
    surface_anom=(surface_monthly_var-annual_mean_at_monthly_resolution)/annual_mean_at_monthly_resolution
    # enforce limits for stability
    if (surface_anom.min()<-.99):
        print ('clipping %i value below -.99 out of %i'%((surface_anom<-.99).sum(),surface_anom.size))
        print ('minimum anomaly: %f'%surface_anom.min())
        surface_anom=maximum(-0.99,surface_anom)
    if (surface_anom.max()>4.):
        print ('clipping %i values above 4 out of %i'%((surface_anom>4).sum(),surface_anom.size))
        print ('maximum anomaly: %f'%surface_anom.max())
    surface_anom=minimum(4.,surface_anom)
    return surface_anom

def seasonal_amplitude_with_depth(var_in):
    # this method calculates how the amplitude of the seasonal cycle of the variable scale with depth in the reference dataset (WOA)
    # this is done over a climatological year (usually not neeed as WOA is already a climatology of 12 months
    
    if mod(var_in.shape[0],12)==0:
        bgc_shape=var_in.shape
        tmp=var_in.reshape(bgc_shape[0]/12,12,bgc_shape[1],bgc_shape[2]).mean(0)   # this calculates climatological year
        amplitude=(tmp.max(0)-tmp.min(0))/tmp.mean(0)
        rel_amplitude=amplitude/amplitude[0]
    else:
        print ('time dimension is not multiple of 12: are you sure it contains seasonal cycle?')
    return rel_amplitude

def stretch_profile(var_to_stretch,MLD_in,z_in,MLD_out,z_out,id=-1):
    # this method stretches the vertical profile of var_stretch such as the associated MLD_in moves to coincide with MLD_out
    # var_stretch is a 2D (depth,position) array
    # MLD is a 1D array (position)

    
    # initialise stretched variable
    dims=var_to_stretch.shape
    stretched_var=zeros(var_to_stretch.shape)
    
    for position in range(dims[1]):
        z_top_in = z_in[:MLD_in[position]+1,position]/float(z_in[MLD_in[position],position])
        z_top_out= z_out[:MLD_out[position]+1,position]/float(z_out[MLD_out[position],position])
        interp_top=interp1d(z_top_in,var_to_stretch[:MLD_in[position]+1,position],bounds_error=False,fill_value=(var_to_stretch[0,position],var_to_stretch[MLD_in[position],position]))
        stretched_var[:MLD_out[position]+1,position]=interp_top(z_top_out)
        if (MLD_out[position]==var_to_stretch.shape[0]-1):
           #if the MLD in the data is at the seafloor, no need to interpolate below the MLD
           continue
        elif (MLD_in[position]==var_to_stretch.shape[0]-1):
            # if the reference MLD reaches the seafloor and the model doesn't, extend the last value to the seafloor
            stretched_var[MLD_out[position]+1:,position]=stretched_var[MLD_out[position],position]
        else:
            # if neither MLD reaches the seafloor interpolates the reference profile below the MLD
            z_bottom_in = (z_in[MLD_in[position]:,position]-z_in[MLD_in[position],position])/float(z_in[-1,position]-z_in[MLD_in[position],position])
            #print mon,pos,MLDi_ref[m,pos],z_ref.size
            z_bottom_out = (z_out[MLD_out[position]:,position]-z_out[MLD_out[position],position])/float(z_out[-1,position]-z_out[MLD_out[position],position])
            int_bottom=interp1d(z_bottom_in,var_to_stretch[MLD_in[position]:,position],bounds_error=False,fill_value=(var_to_stretch[MLD_in[position]+1,position],var_to_stretch[-1,position]))
            stretched_var[MLD_out[position]+1:,position]=int_bottom(z_bottom_out)[1:]

    # if the function is run by Parallel, need o return the timeID
    if id==-1:
        return stretched_var
    else:
        return stretched_var,id

def apply_monthly_cycle(annual_mean,surface_anomalies,rel_amplitude_ref_profile,
                        MLD_in,z_in,MLD_out,z_out,method='stretch',njobs=0,time_shift=0):
    # this method extends the surface seasonal cycle at depth
    # to do this it takes the profile of the rel_amplitude from the input (WOA) dataset and stretches it to adjust to the model MLD (out)
    # MLD_in and z_in refer to the input/reference dataset (e.g. WOA)
    # MLD_out, z_out are the ESM model values
    # the default method to apply the seasonal cycle to the annual mean is stretching, if a more precise conservation of the annual mean is required than the Killworth method is suggested, however this introduce spurious seasonal cycle
    # njobs is the number of processors in case you want to run the function in parallel with Parallel. 0 is to run serial
    # time_shift is the number of month that need to be shifted in case the ESM annual mean is not calculated from Jan to Dec.
    nyears=annual_mean.shape[0]
    nz=annual_mean.shape[1]
    npositions=annual_mean.shape[2]
    
    # here the vertical profiles are stretched repeated for the whole timeseries
    amplitude_profile=rel_amplitude_ref_profile.reshape(1,nz,npositions).repeat(surface_anomalies.shape[0],axis=0)
    
    # check if the "WOA" MLD has the same time extension as the data one or if it is just a climatology
    if MLD_in.shape[0]==MLD_out.shape[0]:
        MLDin_is_climatology=False
    elif MLD_in.shape[0]==12:
        MLDin_is_climatology=True
    else:
        print ('the MLD of the reference profiles has not a correct time dimension')
        print(' this must be either identical to the data MLD, or a climatological one (i.e. 12)')
        return

    # the 3D anomalies are then corrected so the relative vertical profile matches the MLD of the input data
    # if time_shift!=0 then the amplitude_profile need to be shifted (because the first month in the ESM anomalies won't be Jan)
    # when amplitude_profiles or MLD are climatology it can be wrapped around
    # while when full time series are provided better to limit to 0,int(min(max(0,time+time_shift),MLD_out.shape[0])) in order to avoid mixing future and past 
    if njobs>0:
        # in case the stretching is run in parallel the first dump from parallel is stored in "messy_out" and then ordered
        if MLDin_is_climatology:
            messy_out=Parallel(n_jobs=njobs,verbose=1)(delayed(stretch_profile)(amplitude_profile[time+time_shift],MLD_in[mod(time,12)+time_shift],z_in,MLD_out[int(min(max(0,time+time_shift),MLD_out.shape[0]))],z_out,time) for time in range(amplitude_profile.shape[0]))
            for mess in messy_out:
                amplitude_profile[mess[1]]=mess[0]
        else:
            if time_shift!=0:
                print ('WARNING: MLD_ref and MLD_ESM have the same number of time points, therefore they should be aligned')
                print ('a time_shift of %i month(s) has been specified in the yaml file, therefore the MLD_ref will be wrapped around')
            messy_out=Parallel(n_jobs=njobs,verbose=1)(delayed(stretch_profile)(amplitude_profile[int(min(max(0,time+time_shift),MLD_in.shape[0]))],MLD_in[int(min(max(0,time+time_shift),MLD_in.shape[0]))],z_in,MLD_out[int(min(max(0,time+time_shift),MLD_out.shape[0]))],z_out,time) for time in range(amplitude_profile.shape[0]))
            for mess in messy_out:
                amplitude_profile[mess[1]]=mess[0]
    else:
        for time in range(amplitude_profile.shape[0]):
            out_time=int(min(max(0,time+time_shift),MLD_out.shape[0]))
            if mod(time,25)==0: print ('%4i out of %4i'%(time,amplitude_profile.shape[0]))
            if MLDin_is_climatology:
                time_in=mod(time,12)+time_shift
            else:
                if time_shift!=0:
                    print ('WARNING: MLD_ref and MLD_ESM have the same number of time points, therefore they should be aligned')
                    print ('a time_shift of %i month(s) has been specified in the yaml file, therefore the MLD_ref will be wrapped around')
                time_in=int(min(max(0,time+time_shift),MLD_in.shape[0]))
            time_out=int(min(max(0,time+time_shift),MLD_out.shape[0]))
            amplitude_profile[time]=stretch_profile(amplitude_profile[time_in],MLD_in[time_in],z_in,MLD_out[time_out],z_out)
        
    # here the surface anomalies are brought at depth
    # the profile need to be transposed (depth,time,posistion) to multiply
    # then the array is transposed back
    anomalies_3D=surface_anomalies*amplitude_profile.transpose(1,0,2)
    anomalies_3D=anomalies_3D.transpose(1,0,2)

    #here the 3D anomalies are applied to the 3D annual mean
    years=arange(0.5,nyears)
    months=arange(0,nyears,1/12.)
    if method=='stretch':
        repeated_annual_mean=annual_mean.repeat(12,axis=0)
        monthly_3D=repeated_annual_mean*(1+anomalies_3D)
    elif method=='interp':
        print ('this method has to be recoded')
        #time_interp=interp1d(years,annual_mean,axis=0,bounds_error=Fal#se,fill_value=(annual_mean[0],annual_mean[-1]))
        #interpolated_annual_mean=time_interp(months)
        #time_interp=interp1d(years,annual_factor_3D_bgc,axis=0,bounds_error=False,fill_value=(annual_factor_3D_bgc[0],annual_factor_3D_bgc[-1]))
        #tmp=time_interp(months)
    elif method=='killworth':
        print ('this method has to be recoded')
        #interpolated_annual_mean=zeros((months.size,nz,npositions))
        #tmp=zeros((months.size,nz,npositions))
        #for pos in range(npositions):
        #    inteprolated_annual_mean[:,:,pos]= killworth_interpolation(years,annual_mean[:,:,pos],months,z=z_interp)
        #    tmp[:,:,pos]= killworth_interpolation(years,annual_factor_3D_bgc[:,:,pos],months,z=.z_interp)
        #monthly_3D_anomalies=(monthly_factor_3D_bgc-tmp)/tmp            
    return monthly_3D

def calculate_scale_factor(var_input,var_ref,y0=20,y1=51,
                           monthly=False,zlim=0,z_values=0):
        
    # this calculate the 2D scaling factor need to bring the input variable to the same average  value of the reference.
    
    # var_input is the variable to bescaled (usually the ESM)
    # var_ref is the reference variable (usually WOA)
    # y0, y1 are the first and last year that the reference period is representative of (expressed as indices from the beginning of the time series)
    #monthly is a logical switch that indicates if var_input contains monthly values
    # zlim is the depth of deepest data in the reference dataset
    # z_values is the depth 2D (depth, position) array of the BDY points
    
    if monthly:
       data=var_input[y0*12:y1*12].mean(0)
    else:
       data=var_input[y0:y1].mean(0)
    ref=var_ref.mean(0)
    scaling_factor=ref/data

    if zlim>0:
        scaling_factor=where(z_values>zlim,1,scaling_factor)
       #for position in range(scaling_factor.shape[-1]):
       #    above=(z_values[:,position]<=zlim).sum()
       #    scaling_factor[above:,position]=scaling_factor[above-1,position]
    
    return scaling_factor

def killworth_interpolation(x,y,xnew,ext=3,z=0):
    # this function interpolate monthly anual means into monthly time series by "conserving" the annual mean
    #however it introduces spurious seasonal cycle
    x_ext=zeros(x.size+ext*2)
    x_ext[ext:-ext]=x
    x_ext[:ext]=x_ext[ext]+arange(-ext,0)
    x_ext[-ext:]=x_ext[-ext-1]+arange(1,ext+1)
    if (len(y.shape)==1):
        y_ext=zeros(y.size+ext*2)
    elif (len(y.shape)==2):
        y_ext=zeros((y.shape[0]+ext*2,y.shape[1]))
    y_ext[ext:-ext]=y
    y_ext[:ext]=y_ext[ext]
    y_ext[-ext:]=y_ext[-ext-1]
    MAT=DIAG([array([.75,]).repeat(x_ext.size),array([.125,]).repeat(x_ext.size-1),array([.125,]).repeat(x_ext.size-1)],[0,1,-1]).toarray()
    MAT1=INVMAT(MAT)
    y_ext_corrected=dot(MAT1,y_ext)
    if (len(y_ext_corrected.shape)==1):
        int_killworth=interp1d(x_ext,y_ext_corrected)
        ynew=int_killworth(xnew)
    elif (len(y_ext_corrected.shape)==2):
        int_killworth=interp2d(x_ext,z,y_ext_corrected.T)
        ynew=int_killworth(xnew,z).T
    return ynew

def write_bdy(var_to_write,ncname,varname,y0=1960,y1=2099):
    # this write the BDY in a NEMO netcdf BDY file
    count=0
    if var_to_write.shape[0]!=(y1-y0+1)*12:
        print ('time size of the array does not correspondes with the dates provided')
        print ('size: ',var_to_write.shape[0],'dates: %4i-%4i'%(y0,y1))
        return
    if len(var_to_write.shape)==3:
        #save 3D variable
        for y in range(y0,y1+1):
            if mod(y,10)==0: print('saving year %3i'%y)
            for m in range(1,13):
                ncfile=DS(ncname+'_y%4im%02i.nc'%(y,m),'a')
                ncfile.variables[varname][0,:,0,:]=var_to_write[count]
                count+=1
                ncfile.close()
    elif len(var_to_write.shape)==2:
        # save 2D variables like barotropic velocities and ssh
        for y in range(y0,y1+1):
            if mod(y,10)==0: print('saving year %3i'%y)
            for m in range(1,13):
                ncfile=DS(ncname+'_y%4im%02i.nc'%(y,m),'a')
                ncfile.variables[varname][0,0,:]=var_to_write[count]
                count+=1
                ncfile.close()
    return

def read_TS_from_file(BDY,bdyfile,varname,ystart,yend):
    timelen=(yend-ystart+1)*12 # this assumes monthly T,S
    positions=BDY['depth'].shape[1]
    depth_levels=BDY['depth'].shape[0]
    TSbdy=zeros((timelen,depth_levels,positions))
    for y in range(ystart,yend+1):
        if mod(y,10)==0: print ('extracting year: ',y)
        for m in range(1,13):
            time=(y-ystart)*12+m-1
            ncbdyfile=DS(bdyfile+'_y%4im%02i.nc'%(y,m))
            TSbdy[time]=ncbdyfile.variables[varname][:].squeeze()
            ncbdyfile.close()
    return TSbdy

def plot_check(bdyvar,refvar,depbdy,REFy0,REFy1,varname):
    positions=(100,200,331)
    dims=bdyvar.shape
    totyears=int(dims[0]/12.)
    bdyvar=bdyvar.reshape(totyears,12,dims[1],dims[2])
    months=('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
    for pos in positions:
        fig=figure(figsize=(9,12))
        for m in range(12):
            fig.add_subplot(4,3,m+1)
            plot(bdyvar[:REFy0,m,:,pos].T,-depbdy[:,pos],c='C0')
            plot(bdyvar[REFy1:,m,:,pos].T,-depbdy[:,pos],c='C2')
            plot(bdyvar[REFy0:REFy1,m,:,pos].T,-depbdy[:,pos],c='C1')
            plot(refvar[m,:,pos],-depbdy[:,pos],'--',lw=2,c='r')
            title(months[m])
        suptitle('position %s\nblue="past", green="future", orange="present", red=REFERENCE'%pos)
        fig.savefig(varname+'_BDY_%i.png'%pos)
        close(fig)
    return

if __name__=='__main__':
    BGC_WOA_vars=['NO3','PO4','Si','O2','NH4'] #TODO NH4 not actually a WOA variable, but for now need to be here
    BGC_GLODAP_vars=['DIC','bioalk','TA']
    physics_vars=['uo','vo','ssh'] #T and S are not liste here because they are special cases
    if len(sys.argv)>1:
        conf_filename=sys.argv[1]
    else:
        conf_filename='3_BDY_extract_config.yaml'
    conf_file=open(conf_filename)
    print ('reading configuration from '+conf_filename)
    Yconfiguration=Y.load(conf_file)
    ystart=Yconfiguration['y0']
    yend=Yconfiguration['y1']
    
    meshfile=Yconfiguration['meshfile']
    coordsfile=Yconfiguration['coordsfile']
    grid=Yconfiguration['grid']
    
    #check if only the first layer is required (this is needed in case we're extracting depth integrated barotropic velocities or ssh)
    try:
        first_bdylayer_only=Yconfiguration['first_bdylayer_only']
    except:
        first_bdylayer_only=False
    try:
        njobs=Yconfiguration['njobs']
    except:
        njobs=0
    if force_serial: njobs=0
    try:
        ESM_dict=Yconfiguration['ESM_dict']
    except:
        ESM_dict={'NO3':'DIN','Si':'SIL','O2':'OXY','DIC':'DIC','TA':'ALK'}
    try:
        out_dict=Yconfiguration['out_dict']
    except:
        out_dict={'TA': 'TA', 'NH4':'ammonium','NO3':'nitrate','PO4':'phosphate','Si':'silicate','temp':'votemper','sal':'vosaline','O2':'oxygen','DIC':'DIC','bioalk':'bioalk','uo':'vozocrtx','vo':'vomecrty','ssh':'sossheig','small_pon':'small_pon','large_poc':'large_poc','calcite_c':'calcite_c','R3_c':'R3_c','nanophytoplankton_Chl':'nanophytoplankton_Chl','small_pop':'small_pop','R2_c':'R2_c','picophytoplankton_n':'picophytoplankton_n','microphytoplankton_p':'microphytoplankton_p','picophytoplankton_c':'picophytoplankton_c','microphytoplankton_c':'microphytoplankton_c','medium_pos': 'medium_pos','microphytoplankton_n':'microphytoplankton_n','picophytoplankton_p':'picophytoplankton_p','nanophytoplankton_p':'nanophytoplankton_p','large_pop':'large_pop','large_pos':'large_pos','diatoms_p':'diatoms_p','diatoms_s':'diatoms_s','diatoms_n':'diatoms_n','microphytoplankton_Chl':'microphytoplankton_Chl','large_pon':'large_pon','small_poc':'small_poc','nanophytoplankton_n':'nanophytoplankton_n','medium_poc':'medium_poc','medium_pon':'medium_pon','nanophytoplankton_c':'nanophytoplankton_c','picophytoplankton_Chl':'picophytoplankton_Chl','diatoms_c':'diatoms_c','diatoms_Chl':'diatoms_Chl','medium_pop': 'medium_pop'}
    try:
        data_dict=Yconfiguration['data_dict']
    except:
        data_dict={'NO3':'n_an','PO4':'p_an','Si':'i_an','temp':'t_an','sal':'s_an','O2':'o_an','DIC':'TCO2','bioalk':'TAlk','uo':'uo','vo':'vo'}
    variables=Yconfiguration['variables']
    
    #check if at least one variable is a BGC variable, so the script knows it needs to calculate the MLD
    BGC=True
    for varname in variables.keys():
        if (varname in BGC_WOA_vars) or (varname in BGC_GLODAP_vars): BGC=True

    #here the information about themetadata are read and sotred tin the BDY dictionnary
    BDY=define_bdy(coordsfile,meshfile,first_bdylayer_only,grid)

    # check that temperature and salinity are in the variable list
    # if not then make sure that the Tbdy and Sbdy field are defined so the MLD can be calculated
    
    #### TODO this can go after new method being implemented
    #if not (('temp' in variables.keys()) and ('sal' in variables.keys())):
        #try:
            #print ('case not coded yet')
            ##need to think how to merge all files
            #Tesm_file=Yconfiguration['Tbdy']
            #Sesm_file=Yconfiguration['Sbdy']
            #Tref_file=Yconfiguration['Tref']
            #Sref_file=Yconfiguration['Sref']
        #except:
            #print (' the configuration Yaml file is missing important information')
            #print (' information on T and S need to be passed either as variable to extract (add to the "variables" list),')
            #print (' either as bdy file already extracted (add specific Tbdy and Sbdy keys')
    
    #here the temperature and salinity field are either extracted from the ESM file and interpolated or loaded directly from the appropriate bdyfile
    # T & S are converted to float32 to guarantee that MLD calculations are identical regardless if T&S are extracted directly from ESM files or read from bdy files (where they are generally stored as flaot32)
    if 'temp' in variables.keys():
        #extract the temperatureand salinity field from the ESM file and interpolate along the vertical on the BDY grid
        print ('extracting temperature')
        var_keys=Yconfiguration['variables']['temp']

        if 'read_from' in var_keys.keys():
            print ('reading temperature from bdy files '+var_keys['read_from'])
            Tesm_zinterp=float32(read_TS_from_file(BDY,var_keys['read_from'],out_dict['temp'],ystart,yend))
        elif 'ESM_file' in var_keys.keys():
            print ('extracting ESM temperature from  '+var_keys['ESM_file'])
            try:
                ESMbias=var_keys['ESM_bias']
            except:
                ESMbias=0.
            try:
                ESMscale=var_keys['ESM_scale']
            except:
                ESMscale=1.
            Tesm,Tesm_original_depth=extract_values(BDY,var_keys['ESM_file'],'thetao',bias=ESMbias,scaling=ESMscale)
            Tesm_zinterp=float32(depth_interpolation(Tesm,Tesm_original_depth,BDY['depth']))
        else:
            print ('Temperature needs to be extracted by ESM specifying the "ESM_file" key, OR')
            print ('by passing direct link to alreaqdy extracted bdy file through the "read_from" key')
            print ('neither of the two keys is present')
            Tesm_zinterp=None #this would force a crash later

    if 'sal' in variables.keys():
        print ('extracting salinity')
        var_keys=Yconfiguration['variables']['sal']
        if 'read_from' in var_keys.keys():
            print ('reading salinity from bdy files '+var_keys['read_from'])
            Sesm_zinterp=float32(read_TS_from_file(BDY,var_keys['read_from'],out_dict['sal'],ystart,yend))
        elif 'ESM_file' in var_keys.keys():
            print ('extracting ESM salinity from  '+var_keys['ESM_file'])
            try:
                ESMbias=var_keys['ESM_bias']
            except:
                ESMbias=0.
            try:
                ESMscale=var_keys['ESM_scale']
            except:
                ESMscale=1.
            Sesm,Sesm_original_depth=extract_values(BDY,var_keys['ESM_file'],'so',bias=ESMbias,scaling=ESMscale)
            Sesm_zinterp=float32(depth_interpolation(Sesm,Sesm_original_depth,BDY['depth']))
        else:
            print ('Salinity needs to be extracted by ESM specifying the "ESM_file" key, OR')
            print ('by passing direct link to alreaqdy extracted bdy file through the "read_from" key')
            print ('neither of the two keys is present')
            Sesm_zinterp=None #this would force a crash later

    #if there is BGC variable to extract and manipulate, then read the reference T and S data as well and calculate MLD
    if BGC:
        '''
        print ('extracting temperature from reference file')
        var_keys=Yconfiguration['variables']['temp']
        try:
            REFbias=var_keys['REF_bias']
        except:
            REFbias=0.
        try:
            REFscale=var_keys['REF_scale']
        except:
            REFscale=1.
        Tref,Tref_original_depth=extract_values(BDY,var_keys['REF_file'],data_dict['temp'],bias=REFbias,scaling=REFscale)
        Tref_zinterp=depth_interpolation(Tref,Tref_original_depth,BDY['depth'])
        
        print ('extracting salinity from reference file')
        var_keys=Yconfiguration['variables']['sal']
        try:
            REFbias=var_keys['REF_bias']
        except:
            REFbias=0.
        try:
            REFscale=var_keys['REF_scale']
        except:
            REFscale=1.
        Sref,Sref_original_depth=extract_values(BDY,var_keys['REF_file'],'s_an',bias=REFbias,scaling=REFscale)
        Sref_zinterp=depth_interpolation(Sref,Sref_original_depth,BDY['depth'])
        print ('calculating MLD in ESM')
        MLDesm=calculate_MLD(Sesm_zinterp,Tesm_zinterp,BDY['depth'],njobs=njobs)
        print ('calculating MLD in reference dataset')
        MLDref=calculate_MLD(Sref_zinterp,Tref_zinterp,BDY['depth'],njobs=njobs)
        savetxt('MLD_ESM.txt',MLDesm)
        savetxt('MLD_ref.txt',MLDref)
    '''    
    for varname in variables.keys():
        var_keys=Yconfiguration['variables'][varname]
        
        if varname=='temp':
            # temperature values have been already extracted, just need to save them if needed
            if 'output_file' in var_keys.keys():
                print('saving temperature bdy')
                outfile=var_keys['output_file']
                outvar=out_dict[varname]            
                write_bdy(Tesm_zinterp,outfile,outvar,y0=ystart,y1=yend)
            continue
        if varname=='sal':
            # salinity values have been already extracted, just need to save them if needed
            if 'output_file' in var_keys.keys():
                print('saving salinity bdy')
                outfile=var_keys['output_file']
                outvar=out_dict[varname]            
                write_bdy(Sesm_zinterp,outfile,outvar,y0=ystart,y1=yend)
            continue
        
        if 'fixed_value' in var_keys.keys():
            outfile=var_keys['output_file']
            outvar=out_dict[varname]
            value=var_keys['fixed_value']
            tmp_file=DS(outfile+'_y%4im%02i.nc'%(ystart,1))
            tmp_dims=tmp_file.variables[outvar][:].shape
            #monthly_3D=tmp_file.variables[outvar][:]
            tmp_file.close()
            monthly_3D=zeros(((yend-ystart+1)*12,tmp_dims[1],tmp_dims[3]))
            monthly_3D[:]=value
            write_bdy(monthly_3D,outfile,outvar,y0=ystart,y1=yend)
            continue

        # set the bias and scaling factor before extraction
        try:
            ESMbias=var_keys['ESM_bias']
        except:
            ESMbias=0.
        try:
            ESMscale=var_keys['ESM_scale']
        except:
            ESMscale=1.
        try:
            REFbias=var_keys['REF_bias']
        except:
            REFbias=0.
        try:
            REFscale=var_keys['REF_scale']
        except:
            REFscale=1.

        # read the reference period for the scaling factor
        try:
            REFy0=var_keys['REF_y0']
        except:
            REFy0=1980
        try:
            REFy1=var_keys['REF_y1']
        except:
            REFy1=2012
            
        #read the switch that impose to strecth annual means into annualy constant monthly means
        try:
            annual_to_monthly=var_keys['annual_to_monthly']
        except:
            annual_to_monthly=False
        
        # check if the present day mean needs to be scaled to REFERENCE values - default = False
        try:
            scale_to_REF=var_keys['scale_to_REF']
        except:
            scale_to_REF=False
        
        # check if tehre is the need to acocunt for temporal shift in comparing MLD_ESM and MLD_ref
        # this might be needed in model like HadGEM2-ES where annual means are calculated from December of the year before to November
        # this might be supersseded if the values of the time variable are read from the files instead of just assuming they are in the right order
        try:
            time_shift=var_keys['time_shift']
        except:
            time_shift=0
        
        if (varname in ('uo','vo')):
            
            # f the variable is a physical one (othern than T and S) then extract and save straight away
            PHYSesm,PHYSesm_original_depth=extract_values(BDY,var_keys['ESM_file'],ESM_dict[varname],bias=ESMbias,scaling=ESMscale)
            PHYSesm_zinterp=depth_interpolation(PHYSesm,PHYSesm_original_depth,BDY['depth'])
            outfile=var_keys['output_file']
            outvar=out_dict[varname]
            if var_keys['save_3D']:
                write_bdy(PHYSesm_zinterp,outfile,outvar,y0=ystart,y1=yend)
            if var_keys['save_averaged']:
                print('calculating depth averaged integrated barotropic velocities')
                PHYSesm_depth_average=(PHYSesm_zinterp*BDY['e3']).sum(1)/(BDY['e3'].sum(0))
                if 'zo' in outvar:
                    outvar=outvar.replace('zo','bt')
                elif 'me' in outvar:
                    outvar=outvar.replace('me','bt')
                write_bdy(PHYSesm_depth_average,outfile,outvar,y0=ystart,y1=yend)
        
        if varname=='ssh':
            ssh_esm,dummy=extract_values(BDY,var_keys['ESM_file'],ESM_dict[varname],bias=ESMbias,scaling=ESMscale)
            outfile=var_keys['output_file']
            outvar=out_dict[varname]
            write_bdy(ssh_esm,outfile,outvar,y0=ystart,y1=yend)
            
        
        if (varname in BGC_WOA_vars) or (varname in BGC_GLODAP_vars):
            # extract the BGC variable from the reference dataset
            print ('Extract %s from the reference dataset'%varname)
            if 'REF_file' in var_keys.keys():
                BGCref,BGCref_original_depth=extract_values(BDY,var_keys['REF_file'],data_dict[varname],bias=REFbias,scaling=REFscale)
                BGCref_zinterp=depth_interpolation(BGCref,BGCref_original_depth,BDY['depth'])

            #read from the configuration file if the reconstruction of 3D monthly field is needed
            try:
                reconstruct3D=var_keys['reconstruct3D']
                ESM_surface_file=var_keys['ESM_surface_file']
            except:
                reconstruct3D=False

            if reconstruct3D:
                # if the reconstrution of the 3D monthly field is needed than check if this has to be done mimicking the depth profile of the seasonality of another variable
                try:
                    reconstruct_from_variable=var_keys['reconstruct_from_variable']
                except:
                    reconstruct_from_variable=''
                if reconstruct_from_variable!='':
                    # if the profile need to be mimicked, than read further metadata on the profile to be mimicked
                    reconstruct_from_variable_file=var_keys['reconstruct_from_variable_file']
                    try:
                        reconstruct_scale=var_keys['reconstruct_scale']
                    except:
                        reconstruct_scale=1.
                    try:
                        reconstruct_bias=var_keys['reconstruct_bias']
                    except:
                        reconstruct_bias=0.
                        

            print ('Extract %s from the ESM '%varname)
            if reconstruct3D:
                print ('reconstruct')
                BGCesm,BGCesm_original_depth,BGCesm_surface=extract_values(BDY,var_keys['ESM_file'],ESM_dict[varname],bias=ESMbias,scaling=ESMscale,surface=True,surface_varfile=ESM_surface_file)
                BGCesm_zinterp=depth_interpolation(BGCesm,BGCesm_original_depth,BDY['depth'])
                
                # extract the vertical profile of the amplitude of the                seasonal cycle with depth in the reference dataset
                if reconstruct_from_variable=='':
                   rel_amplitude_ref_profile=seasonal_amplitude_with_depth(BGCref_zinterp)
                else:
                    # if the profile of the amplitude need to be mimicked, then extract the variable to be mimicked and calculate the profile for that variable
                    BGCrecon,BGCrecon_original_depth=extract_values(BDY,reconstruct_from_variable_file,data_dict[reconstruct_from_variable],bias=reconstruct_bias,scaling=reconstruct_scale)
                    BGCrecon_zinterp=depth_interpolation(BGCrecon,BGCrecon_original_depth,BDY['depth'])
                    rel_amplitude_ref_profile=seasonal_amplitude_with_depth(BGCrecon_zinterp)
                ESM_surface_anomalies=calculate_surface_anomalies(BGCesm_surface)
                print ('Reconstructing seasonal cycle at depth')
                monthly_3D=apply_monthly_cycle(BGCesm_zinterp,ESM_surface_anomalies,rel_amplitude_ref_profile,MLDref,BDY['depth'],MLDesm,BDY['depth'],njobs=njobs,time_shift=time_shift)
            else:
                print ('NO reconstruct')
                BGCesm,BGCesm_original_depth=extract_values(BDY,var_keys['ESM_file'],ESM_dict[varname],bias=ESMbias,scaling=ESMscale)
                BGCesm_zinterp=depth_interpolation(BGCesm,BGCesm_original_depth,BDY['depth'])
                if annual_to_monthly:
                    monthly_3D=BGCesm_zinterp.repeat(12,axis=0)
                else:
                    monthly_3D=BGCesm_zinterp
                    
            # calculate the scale factor need to bring the long term present day mean of the ESM model to the reference value
            # here the long term mean is calculated between 1980 and 2012
            # 2012 is the last year in WOA v2
            # 1980 is a semi-arbitrary year to account for period with more data
            # to delete - scale_factor=calculate_scale_factor(monthly_3D,BGCref_zinterp,y0=20,y1=51,zlim=500,z_values=BDY['depth'])
            if scale_to_REF:
                REFy0=REFy0-ystart
                REFy1=REFy1-ystart+1
                scale_factor=calculate_scale_factor(monthly_3D,BGCref_zinterp,y0=REFy0,y1=REFy1,monthly=True)
            elif 'scale_as' in var_keys.keys():
                scale_as_varname=var_keys['scale_as']
                try:
                    #try read the scaling factor from file
                    scale_factor=loadtxt(scale_as_varname+'_scaling_factor.dat')
                    print ('scale factor loaded from %s_scaling_factor.dat'%scale_as_varname)
                except:
                    #otherwise calculate from variable
                    print ('calculating scaling factor for %s from ESM and REF values of %s'%(varname,scale_as_varname))
                    try:
                        scale_as_ESMbias=var_keys['ESM_bias']
                    except:
                        scale_as_ESMbias=0.
                    try:
                        scale_as_ESMscale=var_keys['ESM_scale']
                    except:
                        scale_as_ESMscale=1.
                    try:
                        scale_as_REFbias=var_keys['REF_bias']
                    except:
                        scale_as_REFbias=0.
                    try:
                        scale_as_REFscale=var_keys['REF_scale']
                    except:
                        scale_as_REFscale=1.                    
                    scale_as_BGCesm,scale_as_BGCesm_original_depth=extract_values(BDY,var_keys['scale_as_ESM_file'],ESM_dict[scale_as_varname],bias=scale_as_ESMbias,scaling=scale_as_ESMscale,surface=False)
                    scale_as_BGCesm_zinterp=depth_interpolation(scale_as_BGCesm,scale_as_BGCesm_original_depth,BDY['depth'])
                    scale_as_BGCref,scale_as_BGCref_original_depth=extract_values(BDY,var_keys['scale_as_REF_file'],data_dict[scale_as_varname],bias=scale_as_REFbias,scaling=scale_as_REFscale)
                    scale_as_BGCref_zinterp=depth_interpolation(scale_as_BGCref,scale_as_BGCref_original_depth,BDY['depth'])
                    REFy0=REFy0-ystart
                    REFy1=REFy1-ystart+1
                    scale_factor=calculate_scale_factor(monthly_3D,scale_as_BGCref_zinterp,y0=REFy0,y1=REFy1,monthly=True)
                    
            else:
                scale_factor=1.
                
            
            try:
                if var_keys['save_scale']: savetxt(varname+'_scaling_factor.dat',scale_factor)
            except:
                pass                    

            # apply the scale factor
            monthly_3D=monthly_3D*scale_factor

            
            #save the boundary in the appropriate files
            outfile=var_keys['output_file']
            outvar=out_dict[varname]
            write_bdy(monthly_3D,outfile,outvar,y0=ystart,y1=yend)
            try:
                if var_keys['plotcheck']: 
                    plot_check(monthly_3D,BGCref_zinterp,BDY['depth'],REFy0,REFy1,varname)
            except:
                print ('no plotting for',varname)
                pass                    
        
