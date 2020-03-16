from netCDF4 import Dataset as DS
from netCDF4 import default_fillvals as FV
from scipy.interpolate import interp2d,griddata
from numpy import meshgrid,where,zeros,isnan,mod,arange
import sys
from numpy.ma import masked_where as MW
from matplotlib.mlab import find
import yaml as Y
from matplotlib.pyplot import figure,pcolormesh,colorbar,clim,title,close,suptitle

try:
    from joblib import Parallel, delayed
    FORCE_serial=False
except:
    print ('joblib module not present, forcing serial run')
    print ('the execution might take a while')
    FORCE_serial=True


def create_file(input_filename,output_filename,varname,reg_lat,reg_lon):
    # this function create the output netcdf files
    # it requires the file names and the size of the dimensions of the regional model
    # the input filename is used to create the netCDF structure
    # for now it work only on regularly gridded models
    
    if len(reg_lat.shape)==2:
        Ylen=reg_lat.shape[0]
        lat_values=reg_lat[:,0]
    else:
        Ylen=reg_lat.size
        lat_values=reg_lat
    if len(reg_lon.shape)==2:
        Xlen=reg_lon.shape[1]
        lon_values=reg_lon[0]
    else:
        Xlen=reg_lon.size
        lon_values=reg_lon
        
    ncin=DS(input_filename)
    ncout=DS(output_filename,'w')
    #create dimensions and variables
    ncout.createDimension('time',None)
    ncout.createDimension('lat',Ylen)
    ncout.createDimension('lon',Xlen)
    
    ncout.comment='The data contained in this file are an interpolation of ESM output\nEach Z level has been interpolated on the regional model grid using a cubic spline'
    ncout.institution='Plymouth Marine Laboratory'
    ncout.contact='Yuri Artioli - yuti@pml.ac.uk'
    ncout.license='These files are generally for internal use, please write to the "contact" before using the file'
    
    if 'lev' in ncin.dimensions.keys():
        # depth dimensions is not present in surface only outputs
        ncout.createDimension('lev',len(ncin.dimensions['lev']))
        ncout.createVariable('lev','f8',('lev',))
        ncout.variables['lev'][:]=ncin.variables['lev'][:]
        for at in ncin.variables['lev'].ncattrs():
           ncout.variables['lev'].setncattr(at,ncin.variables['lev'].getncattr(at))
    elif 'deptht' in ncin.dimensions.keys():
        # this is for GLODAP format
        ncout.createDimension('lev',len(ncin.dimensions['deptht']))
        ncout.createVariable('lev','f8',('lev',))
        ncout.variables['lev'][:]=ncin.variables['deptht'][:]
        for at in ncin.variables['deptht'].ncattrs():
           ncout.variables['lev'].setncattr(at,ncin.variables['deptht'].getncattr(at))
    ncout.createVariable('time','f8',('time',))
    ncout.variables['time'][:]=ncin.variables['time_average_1mo'][:]
    for at in ncin.variables['time_average_1mo'].ncattrs():
        ncout.variables['time'].setncattr(at,ncin.variables['time_average_1mo'].getncattr(at))
    ncout.createVariable('lat','f8',('lat',))
    ncout.variables['lat'][:]=lat_values
    ncout.variables['lat'].units='degrees_north'
    ncout.variables['lat'].long_name = 'latitude'
    ncout.variables['lat'].standard_name = 'latitude'
    ncout.createVariable('lon','f8',('lon',))
    ncout.variables['lon'][:]=lon_values
    ncout.variables['lon'].units='degrees_east'
    ncout.variables['lon'].long_name = 'longitude'
    ncout.variables['lon'].standard_name = 'longitude'
    if 'lev' in ncout.dimensions.keys():
        ncout.createVariable(varname,'f4',('time','lev','lat','lon'),fill_value=FV['f4'])
    else:
        ncout.createVariable(varname,'f4',('time','lat','lon'),fill_value=FV['f8'])
    ncout.variables[varname].units=ncin.variables[varname].getncattr('units')
    ncout.variables[varname].long_name=ncin.variables[varname].getncattr('long_name')
    try:
        ncout.variables[varname].standard_name=ncin.variables[varname].getncattr('standard_name')
    except:
        pass
    ncout.close()
    ncin.close()
    return

def interpolate_layer(x_in,y_in,v_in,x_out,y_out,Ylen,Xlen,z,method):
    #this function interpolate a 2D layer (v_in) on the regional grid x_out,y_out
    # x_in, y_in and x_out,y_out are 2D arrays with coordinates of all grid points
    # Y,X are the size of the regional domain dimensions

    xin=x_in.ravel()[~v_in.mask.ravel()]
    yin=y_in.ravel()[~v_in.mask.ravel()]
    vin=v_in.ravel()[~v_in.mask.ravel()]
    try:
        var_out=griddata((xin,yin),vin,(x_out.ravel(),y_out.ravel()),method=method).reshape(Ylen,Xlen)
    except:
        #if an exception is raised (usually because v_in is fully masked) the output layer is set to the bottom left corner value of v_in
        var_out=(zeros(x_out.shape)+v_in.data.ravel()[0]).reshape(Ylen,Xlen)
    return z,var_out

def interpolate_on_regional_grid(nfile,outname,varname,reg_lat,reg_lon,reg_mask,njobs=-1,method='cubic',verbose=0):
    ncin=DS(nfile)
    ncout=DS(outname,'a')
    v_in=ncin.variables[varname][:]
    if len(v_in.shape)==3:
        depth=False
    if len(v_in.shape)==4:
        depth=True
    lat_in=ncin.variables['nav_lat'][:]
    lon_in=ncin.variables['nav_lon'][:]
    ncin.close()
    timelen=v_in.shape[0]
    # if longitude is not centered on 0, then correct longitude
    if (lon_in>180).any(): lon_in=where(lon_in>180,lon_in-360,lon_in)
    
    # if coordinates are single coordinates on regular grid, then creates 2D arrys of coordinates for all gridpoints
    if len(lat_in.shape)==1: 
        x_in,y_in=meshgrid(lon_in,lat_in)
    else: 
        x_in=lon_in.copy()
        y_in=lat_in.copy()
    if len(reg_lat.shape)==1: 
        x_out,y_out=meshgrid(reg_lon,reg_lat)
    else: 
        x_out=reg_lon.copy()
        y_out=reg_lat.copy()
    Ylen=y_out.shape[0]
    Xlen=x_out.shape[1]
    if method=='NN':
       #proj=Basemap(llcrnrlon=reg_lon.min(),llcrnrlat=reg_lat.min(),urcrnrlon=reg_lon.max(),urcrnrlat=reg_lat.max(),resolution='l',projection='lcc',lat_1=reg_lat.min(),lat_2=reg_lat.max(),lat_0=(reg_lat.min()+reg_lat.max())*.5,lon_0=0.5*(reg_lon.min()+reg_lon.max()))
       x_in,y_in=proj(x_in,y_in)
       x_out,y_out=proj(x_out,y_out)
    
    if depth:  # the input varaible has depth dimension
        depthlen=v_in.shape[1]
        for time in range(timelen):
            if mod(time,10)==0: 
                print ('interpolating time point # %4i out of %4i'%(time,timelen))
            var_out=zeros((depthlen,Ylen,Xlen))
            if njobs>0:
                messy_out=Parallel(n_jobs=njobs,verbose=verbose)(delayed(interpolate_layer)(x_in,y_in,v_in[time,z],reg_lon,reg_lat,Ylen,Xlen,z,method) for z in range(depthlen))
                for mess in messy_out:
                    var_out[mess[0]]=mess[1]
            else:
                for z in range(depthlen):
                    xin=x_in.ravel()[~v_in[time,z].mask.ravel()]
                    yin=y_in.ravel()[~v_in[time,z].mask.ravel()]
                    vin=v_in[time,z].ravel()[~v_in[time,z].mask.ravel()]
                    var_out[z]=griddata((xin,yin),vin,(reg_lon.ravel(),reg_lat.ravel()),method=method).reshape(Ylen,Xlen)
            var_out=MW(reg_mask,var_out)
            
            #the variable is saved at each time point to overcome potential memory constraint
            ncout.variables[varname][time]=var_out
    else:
        var_out=zeros((v_in.shape[0],Ylen,Xlen))
        if njobs>0:
            messy_out=Parallel(n_jobs=njobs,verbose=verbose)(delayed(interpolate_layer)(x_in,y_in,v_in[time],reg_lon,reg_lat,Ylen,Xlen,time,method) for time in range(timelen))
            for mess in messy_out:
                var_out[mess[0]]=mess[1]
            #to mask the full array need to extend the surface regional mask in time
            mask_to_use=reg_mask.reshape(1,reg_mask.shape[0],reg_mask.shape[1]).repeat(timelen,axis=0)
        var_out=MW(mask_to_use,var_out)
        ncout.variables[varname][:]=var_out
    ncout.close()
    return

def horizontal_NN_interpolation(lon,lat,v_in,z):
    # this function fills the land points with the nearest neighbour in the same Z depth level
    #if mod(z,10)==0: print ('%4i layers done'%z)
    Ylen=lat.shape[0]
    Xlen=lat.shape[1]
    xin=lon.ravel()[~v_in.mask.ravel()]
    yin=lat.ravel()[~v_in.mask.ravel()]
    vin=v_in.ravel()[~v_in.mask.ravel()]
    vout=griddata((xin,yin),vin,(lon.ravel(),lat.ravel()),method='nearest').reshape(Ylen,Xlen)
    return [z,vout]

def fill_vertical(v_in):
    # this function fills the deep masked region of sea areas with the last valid point on the vertical
    #v_in is a 3D array (depth,lat,lon)
    #if mod(x,10000)==0: print ('%6i grid points done'%x)
    v_in=MW(isnan(v_in),v_in)
    seafloor_index=(~v_in[:].mask).sum(0)-1
    dims=v_in.shape
    # these ii and jj indices are needed to select the values of v_in at seafloor using fancy indexing in a fast way
    ii=arange(dims[2]).reshape(1,dims[2]).repeat(dims[1],axis=0)
    jj=arange(dims[1]).reshape(dims[1],1).repeat(dims[2],axis=1)
    seafloor=v_in[seafloor_index,jj,ii]
    vout=where(v_in.mask,seafloor,v_in)
    return vout

def fill_masked_values(outname,varname,tmask,verbose=0):
    ncout=DS(outname,'a')
    lat=ncout.variables['lat'][:]
    lon=ncout.variables['lon'][:]
    # if longitude is not centered on 0, then correct longitude
    if (lon>180).any(): lon=where(lon>180,lon-360,lon)
    # if coordinates are single coordinates on regular grid, then creates 2D arrys of coordinates for all gridpoints
    if len(lat.shape)==1: lon,lat=meshgrid(lon,lat)
    Ylen=lat.shape[0]
    Xlen=lon.shape[1]
    ndims=len(ncout.variables[varname].dimensions)
    timelen=len(ncout.dimensions['time'])
    if ndims==4: 
        #3D variables are first fill vertically with last admissable value
        # then inteprolated with NN on Z levels
        # the I/O is done inside the loop because for big monthly array the memory requirement can be too big
        if njobs>0:
            for time in range(timelen): 
                v_in=ncout.variables[varname][time]
                # mask potential nans generated by th previous interpolation if gridpointsds were outside of the convex hull of the input data
                v_in=MW(isnan(v_in),v_in)
                if mod(time,10)==0: print('filling timestep #%3i out of %3i'%(time,timelen))
                var_out=fill_vertical(v_in)  # this could speed up by parallel looping over the flattened horizontal dimension
                var_out=MW(var_out>1e10,var_out)
                messy_out=Parallel(n_jobs=njobs,verbose=verbose)(delayed(horizontal_NN_interpolation)(lon,lat,var_out[z],z) for z in range(v_in.shape[0]))
                for mess in messy_out:
                    var_out[mess[0]]=mess[1]
                #mask on the regional model Tmask
                var_out=MW(tmask.reshape(1,Ylen,Xlen).repeat(var_out.shape[0],axis=0)==0,var_out)
                ncout.variables[varname][time,:,:,:]=var_out

            #TODO this part of the code can be recovedred to speed up
            #flatten the array to optimize parallelisation for vertical gap filling
            #print ('filling deep Z layers: %i to go'%(dims[2]*dims[3]))
            #tmp_v_in=v_in.reshape(dims[0],dims[1],dims[2]*dims[3])
            #filled_array=filled_array.reshape(dims[0],dims[1],dims[2]*dims[3])
            #messy_out=Parallel(n_jobs=njobs,verbose=verbose)(delayed(fill_vertical)(tmp_v_in[:,:,position],position) for position in range(dims[2]*dims[3]))
            #for mess in messy_out:
            #  filled_array[:,:,mess[0]]=mess[1]
            #filled_array=filled_array.reshape(dims)
            #filled_array=MW(filled_array>1e10,filled_array)
            #Afilled_array=zeros(filled_array.shape)
            #flatten to optimise optimise parallelisation for horizontal NN inteprolation
            #filled_array=filled_array.reshape(dims[0]*dims[1],dims[2],dims[3])

            #print ('filling land points: %i to go'%(dims[0]*dims[1]))
            #for time in range(dims[0]):
            #    messy_out=Parallel(n_jobs=njobs,verbose=verbose)(delayed(horizontal_NN_interpolation)(lon,lat,filled_array[time,layer],layer) for layer in range(dims[1]))
            #    for mess in messy_out:
            #        Afilled_array[time,mess[0]]=mess[1]
            #Afilled_array=Afilled_array.reshape(dims)
        else: #serial
            for time in range(timelen):
                if mod(time,10)==0: print('filling timestep #%3i out of %3i'%(time,timelen))
                v_in=ncout.variables[varname][time,:,:,:]
                v_in=MW(isnan(v_in),v_in)
                var_out=fill_vertical(v_in)
                dims=v_in.shape
                var_out=MW(var_out>1e10,var_out)
                for z in range(dims[0]):
                    dummy,var_out[z]=horizontal_NN_interpolation(lon,lat,var_out[z],z)
                #mask on the regional model Tmask
                var_out=MW(tmask.reshape(1,Ylen,Xlen).repeat(var_out.shape[0],axis=0)==0,var_out)
                ncout.variables[varname][time,:,:,:]=var_out
    elif ndims==3: #2D variable
        v_in=ncout.variables[varname][:]
        var_out=zeros(v_in.shape)
        if njobs>0:
            messy_out=Parallel(n_jobs=njobs,verbose=verbose)(delayed(horizontal_NN_interpolation)(lon,lat,v_in[time],time) for time in range(timelen))
            for mess in messy_out:
                var_out[mess[0]]=mess[1]
        else:
            for time in range(timelen):
                var_out[t]=horizontal_NN_interpolation(lon,lat,v_in[time],time)
        #mask on the regional model Tmask
        var_out=MW(tmask.reshape(1,Ylen,Xlen).repeat(timelen,axis=0)==0,var_out)
        ncout.variables[varname][:]=var_out

    
    #ncout.variables[varname][:]=Afilled_array
    ncout.close()
    return

def create3Dmask(filename,varname,reg_lat,reg_lon):
    # this function creates a 3D mask of the ESM domain on the regional model resolution (reg_lat,reg_lon)
    
    # read the mask from the variable
    ncin=DS(filename)
    mask_in=ncin.variables[varname][0].mask
    lat_in=ncin.variables['nav_lat'][:]
    lon_in=ncin.variables['nav_lon'][:]
    ncin.close()
    # if longitude is not centered on 0, then correct longitude
    if (lon_in>180).any(): lon_in=where(lon_in>180,lon_in-360,lon_in)
    # if coordinates are single coordinates on regular grid, then creates 2D arrys of coordinates for all gridpoints
    if len(lat_in.shape)==1: lon_in,lat_in=meshgrid(lon_in,lat_in)
    if len(reg_lat.shape)==1:
        reg_x,reg_y=meshgrid(reg_lon,reg_lat)
    else:
        reg_y=reg_lat.copy()
        reg_x=reg_lon.copy()
    Ylen=reg_y.shape[0]
    Xlen=reg_x.shape[1]
    if len(mask_in.shape)==3: #file is 3D
        mask_out=zeros((mask_in.shape[0],Ylen,Xlen))
        for z in range(mask_in.shape[0]):
            xin=lon_in.ravel()
            yin=lat_in.ravel()
            vin=mask_in[z].ravel()
            mask_out[z]=griddata((xin,yin),vin,(reg_x.ravel(),reg_y.ravel()),method='nearest').reshape(Ylen,Xlen)
    else:
        xin=lon_in.ravel()
        yin=lat_in.ravel()
        vin=mask_in.ravel()
        mask_out=griddata((xin,yin),vin,(reg_x.ravel(),reg_y.ravel()),method='nearest').reshape(Ylen,Xlen)
    return mask_out

def plot_check(output_filename,input_filename,varname):
    ncin=DS(input_filename)
    ncout=DS(output_filename)
    timelen=len(ncin.dimensions['time_counter'])
    if 'lev' in ncin.dimensions.keys():
        deplen=len(ncin.dimensions['lev'])
        depth=True
    elif 'depth' in ncin.dimensions.keys():
        deplen=len(ncin.dimensions['deptht'])
        depth=True
    else: depth=False
    if depth:
        depths=(0,int(deplen*.5))
        vin_1=ncin.variables[varname][0,depths[0]]
        vin_2=ncin.variables[varname][0,depths[1]]
        vout_1=ncout.variables[varname][0,depths[0]]
        vout_2=ncout.variables[varname][0,depths[1]]
        t1='original, time=0, depth_lev=0'
        t2='original, time=0, depth_lev=%i'%depths[1]
        t3='interp., time=0, depth_lev=0'
        t4='interp., time=0, depth_lev=%i'%depths[1]
        figname=varname+'_3D.png'
    else:
        vin_1=ncin.variables[varname][0]
        vin_2=ncin.variables[varname][-1]
        vout_1=ncout.variables[varname][0]
        vout_2=ncout.variables[varname][-1]
        t1='original, time=0'
        t2='original, time=last'
        t3='interp., time=0'
        t4='interp., time=last'
        figname=varname+'_surface.png'
    vmin=min(vin_1.min(),vin_2.min(),vout_1.min(),vout_2.min())
    vmax=max(vin_1.max(),vin_2.max(),vout_1.max(),vout_2.max())
    fig=figure(figsize=(12,12))
    fig.add_subplot(2,2,1)
    pcolormesh(vin_1)
    clim(vmin,vmax)
    colorbar()
    title(t1)
    fig.add_subplot(2,2,2)
    pcolormesh(vin_2)
    clim(vmin,vmax)
    colorbar()
    title(t2)
    fig.add_subplot(2,2,3)
    pcolormesh(vout_1)
    clim(vmin,vmax)
    colorbar()
    title(t3)
    fig.add_subplot(2,2,4)
    pcolormesh(vout_2)
    clim(vmin,vmax)
    colorbar()
    title(t4)
    suptitle('Remember that interpolated data below seafloor are not masked\nbottom water values are extended blow seafloor')
    fig.savefig(figname)
    close(fig)
    


if __name__=='__main__':
    if len(sys.argv)>1:
        conf_filename=sys.argv[1]
    else:
        conf_filename='1_interp.yaml'
    configuration_file=open(conf_filename)
    Yconfiguration=Y.load(configuration_file)
    configuration_file.close()
    varname=Yconfiguration['varname']
    mesh_filename=Yconfiguration['meshfilename']
    input_filename=Yconfiguration['inputfilename']
    output_filename=Yconfiguration['outputfilename']
    try:
        method=Yconfiguration['method']
    except:
        method='cubic'
    try: 
        njobs=Yconfiguration['njobs']
    except:
        njobs=-1
    if FORCE_serial:
        njobs=-1
    try:
        verbose=Yconfiguration['verbose']
    except:
        verbose=0
    try:
        plot=Yconfiguration['plotcheck']
    except:
        plot=False
        
    #read lat, lon and T mask of the regional model
    mesh_file=DS(mesh_filename)
    reg_lat=mesh_file.variables['nav_lat'][:]
    reg_lon=mesh_file.variables['nav_lon'][:]
    tmask=mesh_file.variables['tmask'][0,0]
    mesh_file.close()

    print ('extracting ESM mask on regional grid resolution')
    mask3D=create3Dmask(input_filename,varname,reg_lat,reg_lon)

    print ('creating output file')
    create_file(input_filename,output_filename,varname,reg_lat,reg_lon)
    print ('interpolating ESM on regional grid')
    interpolate_on_regional_grid(input_filename,output_filename,varname,reg_lat,reg_lon,mask3D,njobs,method)
    print ('filling masked values to guarantee all domain has reasonable values')
    fill_masked_values(output_filename,varname,tmask)
    if plot:
        print ('plotting picture for control')
        plot_check(output_filename,input_filename,varname)
