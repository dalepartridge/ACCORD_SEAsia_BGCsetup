from netCDF4 import Dataset as DS
from netCDF4 import default_fillvals as FV
import yaml as Y
from numpy import array,zeros
from datetime import datetime
import sys
import calendar


def create_file (filename,config,coords,gridchar):
    #this function create an empty boundary file to be used in the extract_BDY script
    # it requires as input:
    # filename: filename to create
    # config: the dictionnary containing all configuration options associated to that file from the yaml file
    # coords: a dictionnary containing the coordinates of the bdy points both in i,j,r and lon,lat
    nr=coords['nbr'+gridchar].max()
    dims=coords['nbi'+gridchar].shape
    print ('creating '+filename)
    ncfile=DS(filename,'w')
    ncfile.createDimension('time_counter',None)
    timedim=1
    ncfile.createDimension('xb',dims[1])
    ncfile.createDimension('yb',dims[0])
    var_layers=[]
    for varname in coords.keys():
        #create and copy coordinates variables
        if varname[0]=='n':
            ncfile.createVariable(varname,'int',('yb','xb'))
            ncfile.variables[varname][:]=coords[varname]
        if varname[0]=='g':
            ncfile.createVariable(varname,'f4',('yb','xb'))
            ncfile.variables[varname][:]=coords[varname]
    for varname in config['variables'].keys():
        var_layers+=[config['variables'][varname]['depth'],]
    max_depth=array(var_layers).max()
    check=(var_layers==max_depth)|(var_layers==zeros(len(var_layers)))
    if not(check.all()):
            print ("variables in %s have different numbers of z layers"%filename)
            print ("be sure they are either 0 or the same number")
            sys.exit()
    if max_depth>0:
        ncfile.createDimension('zb',max_depth)
    ncfile.creation_date=datetime.now().isoformat()
    ncfile.rim_width='%iLL'%nr
    for at in config.keys():
        if at=='variables':
            continue
        ncfile.setncattr(at,config[at])
    var_config=config['variables']
    for varname in var_config.keys():
        if var_config[varname]['depth']>0:
            ncfile.createVariable(varname,'f4',('time_counter','zb','yb','xb'),zlib=True,complevel=9,fill_value=FV['f4'])
            ncfile.variables[varname][:]=zeros((timedim,max_depth,dims[0],dims[1]))+FV['f4']
        else:
            ncfile.createVariable(varname,'f4',('time_counter','yb','xb'),zlib=True,complevel=9,fill_value=FV['f4'])
            ncfile.variables[varname][:]=zeros((timedim,dims[0],dims[1]))+FV['f4']
        ncfile.variables[varname].grid='bdy'+gridchar.upper()
        for at in var_config[varname].keys():
            if at=='depth': 
                continue
            ncfile.variables[varname].setncattr(at,var_config[varname][at])
    ncfile.close()
    return

if __name__=='__main__':

    if len(sys.argv)>1:
        conf_filename=sys.argv[1]
    else:
        conf_filename='2_create_OBC_file.yaml'
    conf_file=open(conf_filename)
    print ('reading configuration from '+conf_filename)
    Yconfiguration=Y.load(conf_file)
    ystart=Yconfiguration['ystart']
    yend=Yconfiguration['yend']
    coords_file=Yconfiguration['coords_file']
    out_folder=Yconfiguration['output_folder']
    try:
        monthly=Yconfiguration['monthly']
    except:
        monthly=True
    try:
        daily=Yconfiguration['daily']
    except:
        daily=False
    
    if daily&monthly:
        print ('The script has been configured to be both daily and monthly')
        print ('this is not currently allowed, if you need to generate both types, create two configuration files and run the script twice')
        print ('now only the MONTHLY files will be generated')
        daily=False

    for Yfile in Yconfiguration['files'].keys():
        fcoords=DS(coords_file)
        grid=Yconfiguration['files'][Yfile]['grid']
        gridchar=grid.lower()
        coords={}
        coords['nbi'+gridchar]=fcoords.variables['nbi'+gridchar][:]
        coords['nbj'+gridchar]=fcoords.variables['nbj'+gridchar][:]
        coords['nbr'+gridchar]=fcoords.variables['nbr'+gridchar][:]
        coords['gphi'+gridchar]=fcoords.variables['gphi'+gridchar][:]
        coords['glam'+gridchar]=fcoords.variables['glam'+gridchar][:]
        try:
            rim=min(coords['nbr'+gridchar].max(),Yconfiguration['files'][Yfile]['rim'])
        except:
            rim=coords['nbr'+gridchar].max()
        if rim!=coords['nbr'+gridchar].max():
            dims=coords['nbi'+gridchar].shape
            last=(coords['nbr'+gridchar]<=rim).sum()
            coords['nbi'+gridchar]=coords['nbi'+gridchar][:,:last]
            coords['nbj'+gridchar]= coords['nbj'+gridchar][:,:last]
            coords['nbr'+gridchar]= coords['nbr'+gridchar][:,:last]
            coords['gphi'+gridchar]= coords['gphi'+gridchar][:,:last]
            coords['glam'+gridchar]= coords['glam'+gridchar][:,:last]
        
        fcoords.close()
        for y in range(ystart,yend):
            if monthly:
                for m in range(12):
                    filename=out_folder+'/'+Yfile+'_y%4im%02i.nc'%(y,m+1)
                    create_file(filename,Yconfiguration['files'][Yfile],coords,gridchar)
            elif daily:
                if calendar.isleap(y):
                    days_per_month=(31,29,31,30,31,30,31,31,30,31,30,31)
                else:
                    days_per_month=(31,28,31,30,31,30,31,31,30,31,30,31)
                for m in range(12):
                    for d in range(days_per_month[m]):
                        filename=out_folder+'/'+Yfile+'_y%4im%02id%02i.nc'%(y,m+1,d+1)
                        create_file(filename,Yconfiguration['files'][Yfile],coords,gridchar)
