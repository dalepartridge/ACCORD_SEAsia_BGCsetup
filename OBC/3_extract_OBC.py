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

from netCDF4 import Dataset as DS
from numpy import zeros, ones, array, mod
from scipy.interpolate import interp1d
import yaml as Y
import sys

def define_bdy(coordsfile,meshfile,grid):
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

def extract_values(BDY,varfile,varname,scaling=1.,bias=0.,fill_mask=False):
    # this function extracte the variable "varname" from the "varfile" netcdf file uding the BDY metadata in the BDY dictionnary
    # if the surface monthly values need to be extracted then the logical switch needs to be set to true and the name of the netcdf file with the surface monthly value need to be read.
    # it also extract the original depth information
    
    # open the files and check what is the name of the depth variable
    fvar=DS(varfile)
    if 'lev' in fvar.variables.keys():
        zvar='lev'
    if 'nav_lev' in fvar.variables.keys():
        zvar='nav_lev'
        zdim='z'
    elif 'depth' in fvar.variables.keys():
        zvar='depth'
    else: 
        print ('No recognised depth array in the file')
        zvar=''
    
    #extract variable
    # if the length of the time dimension is big (i.e. monthly), then extract data one timstep at time to save memory
    print ('extracting values for '+varname)
    timelen=len(fvar.dimensions['time_counter'])
    if zvar!='':
        nlev=len(fvar.dimensions[zdim])
        var=zeros((timelen,nlev,BDY['nbj'].size))
    else:
        var=zeros((timelen,BDY['nbj'].size))
    print ('given the big size the extraction is one timestep at time')
    print ('this might take few minutes')
    for time in range(timelen):
        tmp_var=fvar.variables[varname][time]
        if zvar!='':
            var[time]=tmp_var[:,BDY['nbj'],BDY['nbi']]*scaling-bias
        else:
            var[time]=tmp_var[BDY['nbj'],BDY['nbi']]*scaling-bias
        if mod(time,50)==0:
            print ('%4i out of %4i'%(time,timelen))
            
    # extract depth information
    if zvar!='':
        original_depth=fvar.variables[zvar][:]
    else:
        original_depth=0.
    fvar.close()
    return var,original_depth
    
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

if __name__=='__main__':
    if len(sys.argv)>1:
        conf_filename=sys.argv[1]
    else:
        conf_filename='3_extract_OBC.yaml'
    conf_file=open(conf_filename)
    print ('reading configuration from '+conf_filename)
    Yconfiguration=Y.load(conf_file)
    ystart=Yconfiguration['y0']
    yend=Yconfiguration['y1']
    
    meshfile=Yconfiguration['meshfile']
    coordsfile=Yconfiguration['coordsfile']
    grid=Yconfiguration['grid']
    variables=Yconfiguration['variables']
    
    infile_dict={'PO4':'PO4', 'O2':'OXY', 'NO3':'NO3', 'Si': 'SIL'}  
    out_dict={'TA': 'TA', 'NH4':'ammonium','NO3':'nitrate','PO4':'phosphate','Si':'silicate','temp':'votemper','sal':'vosaline','O2':'oxygen','DIC':'DIC','bioalk':'bioalk','uo':'vozocrtx','vo':'vomecrty','ssh':'sossheig','small_pon':'small_pon','large_poc':'large_poc','calcite_c':'calcite_c','R3_c':'R3_c','nanophytoplankton_Chl':'nanophytoplankton_Chl','small_pop':'small_pop','R2_c':'R2_c','picophytoplankton_n':'picophytoplankton_n','microphytoplankton_p':'microphytoplankton_p','picophytoplankton_c':'picophytoplankton_c','microphytoplankton_c':'microphytoplankton_c','medium_pos': 'medium_pos','microphytoplankton_n':'microphytoplankton_n','picophytoplankton_p':'picophytoplankton_p','nanophytoplankton_p':'nanophytoplankton_p','large_pop':'large_pop','large_pos':'large_pos','diatoms_p':'diatoms_p','diatoms_s':'diatoms_s','diatoms_n':'diatoms_n','microphytoplankton_Chl':'microphytoplankton_Chl','large_pon':'large_pon','small_poc':'small_poc','nanophytoplankton_n':'nanophytoplankton_n','medium_poc':'medium_poc','medium_pon':'medium_pon','nanophytoplankton_c':'nanophytoplankton_c','picophytoplankton_Chl':'picophytoplankton_Chl','diatoms_c':'diatoms_c','diatoms_Chl':'diatoms_Chl','medium_pop': 'medium_pop'}
    
    #here the information about themetadata are read and sotred tin the BDY dictionnary
    BDY=define_bdy(coordsfile,meshfile,grid)

    for varname in variables.keys():
        var_keys=Yconfiguration['variables'][varname]
        print ('Processing %s'%varname)
        
        if 'fixed_value' in var_keys.keys():
            outfile=var_keys['output_file']
            outvar=out_dict[varname]
            value=var_keys['fixed_value']
            print ('Applying fixed value of %s'%str(value))
            
            tmp_file=DS(outfile+'_y%4im%02i.nc'%(ystart,1))
            tmp_dims=tmp_file.variables[outvar][:].shape
            tmp_file.close()
            
            output_val=value*ones(((yend-ystart+1)*12,tmp_dims[1],tmp_dims[3]))
            write_bdy(output_val,outfile,outvar,y0=ystart,y1=yend)
            continue

        elif 'input_file' in var_keys.keys():
            # extract the BGC variable from the reference dataset
            print ('Extracting field from %s'%var_keys['input_file'])
            output_val,input_depth=extract_values(BDY,var_keys['input_file'],infile_dict[varname])
            
            #save the boundary in the appropriate files
            outfile=var_keys['output_file']
            outvar=out_dict[varname]
            write_bdy(output_val,outfile,outvar,y0=ystart,y1=yend)
            continue
