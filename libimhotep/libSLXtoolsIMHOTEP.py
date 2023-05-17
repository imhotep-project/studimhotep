"""My IMHOTEP tools
This is a collection of  tools i often use when analysing the IMHOTEP project data. 
"""


## standart libraries

import os,sys
import numpy as np

# xarray
import xarray as xr

# plot
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.colors import Colormap

import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.cm as cm
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
from matplotlib.colors import from_levels_and_colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from matplotlib import cm 
from matplotlib.colors import ListedColormap,LinearSegmentedColormap

import cartopy.feature as cfeature

import cmocean

print(f"Name: {__name__}")
print(f"Package: {__package__}")

print("This is a collection of  tools i often use when analysing the IMHOTEP project data.")


#def __init__():
#    print('Init')

def main():
    print('Create example map (global) without any data.')
    #set default grid parameters before plotting
    gridparam = Fpltgridparamdefaults(reg='GLO')
    print('- gridparam')
    print(gridparam)
    # plot
    outo = FpltGLO(0,pltgridparam=gridparam,ty='DIFF',saveo=True,pltshow=False)
    
def Floadmultimb(N,nexp,prefix,varnasss,fo,fity,y1,y2,xselect=[0,20],yselect=[0,20],verbose=False):  
    """ Load data from N ensemble members on a selected subregion
            
        Parameters:
        - N (int): Number of members
        - nexp (str): name of ensemble experiments
        - prefix (str): prexix of experiment
        - varnasss(str): name of variable to read in file
        - fo (str): frequency '1m', '1y'
        - fity (str): type of gridpoint
        - xselect (array of 2): xmin, xmax indices
        - yselect (array of 2): ymin, ymax indices
        
        Returns:
        xarrays: alldat concatenating the data in the subregion points for all the members (e axis)
        """
    # loop on members
    for ie in range(1,N+1):
        # load 1 member
        diridat = Ffindinputdata(nexp,prefix=prefix,fo=fo,onemember=ie)
        files   = Ffindinputdfiles(nexp,diridat,prefix=prefix,fo=fo,onemember=ie)
        if verbose:
            print("read ensemble mb "+str(ie))
            print(diridat)
        data = xr.open_mfdataset(files,decode_times=True)[varnasss].sel(time_counter=slice(y1,y2)).isel(x=slice(xselect[0],xselect[1]),y=slice(yselect[0],yselect[1]))
        # concat with previously loaded members
        if (ie == 1):
            alldat = data
        else:
            alldat = xr.concat([alldat, data], "e")
                                                                                                           
    return alldat

def FloadmultiGP(N,nexp,prefix,varnasss,fo,fity,y1,y2,gpx=[0],gpy=[0],verbose=False):
    """ Load data at a list of gridpoints defined by their x,y indices in N ensemble members
            
        Parameters:
        - N (int): Number of members
        - nexp (str): name of ensemble experiments
        - prefix (str): prexix of experiment
        - varnasss(str): name of variable to read in file
        - fo (str): frequency '1m', '1y'
        - fity (str): type of gridpoint
        - gpx (list of integers): x indices of requested points
        - gpy (list of integers): y indices of requested points.
        
        Returns:
        xarrays: allGP concatenating the data at the requested points for all the members (e axis)
        """
    NGP=len(gpx)
    if len(gpx)!= len(gpy) :
        sys.exit("ERROR: gpx should be of same size as gpy")
    else:
        # loop on gridpoints
        for igp in range(0,NGP):
            GP = Floadmultimb(N,nexp,prefix,varnasss,fo,fity,y1,y2,xselect=[gpx[igp],gpx[igp]+1],yselect=[gpy[igp],gpy[igp]+1],verbose=verbose)      
            # concat with previously loaded members
            if (igp == 0):
                allGP = GP
            else:
                allGP = xr.concat([allGP, GP], "gp")
    return allGP

def FloadcoordGP(dirigrid,gpx=[0],gpy=[0]):
    """ Load lat lon coordinates a list of gridpoints defined by their x,y indices
            
        Parameters:
        - dirigrid (str): grid files where to read nav_lon nav_lat,
        - gpx (list of integers): x indices of requested points
        - gpy (list of integers): y indices of requested points.
        
        Returns:
        xarrays: GPlats,GPlons with latitudes and longitudes of the requested points
        """
    NGP=len(gpx)
    if len(gpx)!= len(gpy) :
        sys.exit("ERROR: gpx should be of same size as gpy")
    else:
        # loop on gridpoints
        for igp in range(0,NGP):
            lat = xr.open_dataset(dirigrid+'mesh_hgr.nc')['nav_lat'].isel(x=slice(gpx[igp],gpx[igp]+1),y=slice(gpy[igp],gpy[igp]+1))
            lon = xr.open_dataset(dirigrid+'mesh_hgr.nc')['nav_lon'].isel(x=slice(gpx[igp],gpx[igp]+1),y=slice(gpy[igp],gpy[igp]+1))     
            # concat with previously loaded members
            if (igp == 0):
                GPlats = lat
                GPlons = lon
            else:
                GPlats = xr.concat([GPlats, lat], "gp")
                GPlons = xr.concat([GPlons, lon], "gp") 
    return GPlats,GPlons

def Ffindij(dirigrid,latval,lonval):
    """ Finds i,j coordinates nearest to latitude and longitude values
            
        Parameters:
        - dirigrid (str): grid files where to read nav_lon nav_lat,
        - latval (float): latitude value requested,
        - lonval (float): longitude value requested.
        
        Returns:
        str: string of 3 characters (integer padded with zeroes).
        """
    lat = xr.open_dataset(dirigrid+'mesh_hgr.nc')['nav_lat']
    lon = xr.open_dataset(dirigrid+'mesh_hgr.nc')['nav_lon']     

    # Find the index of the grid point nearest a specific lat/lon.   
    abslat = np.abs(lat-latval)
    abslon = np.abs(lon-lonval)
    c = np.maximum(abslon, abslat)

    ([j], [i]) = np.where(c == np.min(c))
    
    return i,j
    
    
def Fstrmb(ie):
    """ Reads an integer between 1-100 and returns a string with 2 or 1 zeroes before int so that there are 3 characters in total.
            
        Parameters:
        - ie (int): integer between 1-100
        
        Returns:
        str: name with 2 zeroes before integer.
        """
    if ie<10:
        ena = "00"+str(ie)
    elif ((ie>9)&(ie<100)):
        ena= "0"+str(ie)
    else:
        ena= str(ie)
    return ena

def Ffindinputdata(nexp,prefix="eORCA025.L75-IMHOTEP",fo='1y',allmb=False,onemember='no'):
    """ Find data path depending on which experiment
            
        Parameters:
        - nexp (str): experiment name
        - prefix (str): prefix file name
        - fo (str): input frequency '1y' '1m'
        - allmb (bool): set to True if nexp is an ensemble experiment and you want to read all files from all members
        - onemember (int): set to the  number of the ensemble member (start at 1) you want to load (index starts at 01)
        
        Returns:
        str: path to directory
        """
    
    #  input directory on store
    if nexp=='SC':
        diri="/gpfsstore/rech/cli/commun/IMHOTEP/eORCA025.L75-IMHOTEP.SC-S/"
    else:
        if nexp=='02':
            diri="/gpfsstore/rech/cli/commun/IMHOTEP/eORCA025.L75-IMHOTEP.02/"+fo+"/"
        else:
            if (((nexp=='ES')|(nexp=='EGAI'))|(nexp=='EAI')):
                if allmb:
                    diri="/gpfsstore/rech/cli/rcli002/eORCA025.L75/"+prefix+"."+nexp+".???-S/"
                elif (onemember!='no'):
                    diri="/gpfsstore/rech/cli/rcli002/eORCA025.L75/"+prefix+"."+nexp+"."+Fstrmb(onemember)+"-S/"
                else:
                    diri="/gpfsstore/rech/cli/commun/IMHOTEP/ENSTATS_"+fo+"/"+nexp+"/"
            else:
                diri="/gpfsstore/rech/cli/rcli002/eORCA025.L75/"+prefix+"."+nexp+"-S/"
                
    return diri

def Ffindinputdfiles(nexp,diri,fo,prefix="eORCA025.L75-IMHOTEP",fity='gridTsurf',allmb=False,onemember='no'):
    """ Find data path depending on which experiment
            
        Parameters:
        - nexp (str): experiment name
        - diri (str): path
        - fo (str): data frequency '1m', '1y' '1d'
        - prefix (str): prefix file name
        - fity (str): suffix name of files (file type). Default is gridTsurf.
        
        Returns:
        str: path to files
        """
    
    # data file paths
    if nexp=='02':
        sssfiles = diri+"????/"+prefix+nexp+"*"+fity+".nc"
    else:
        if (((nexp=='ES')|(nexp=='EGAI'))|(nexp=='EAI')): 
            if allmb:
                sssfiles = diri+fo+"/????/"+prefix+"."+nexp+".???*"+fity+".nc"
            elif (onemember!='no'):
                sssfiles = diri+fo+"/????/"+prefix+"."+nexp+"."+Fstrmb(onemember)+"*"+fity+".nc"
            else:
                sssfiles = diri+"ESTATS_"+prefix+"."+nexp+"*"+fity+".nc"
                if (fity=="flxT"):
                    diridat = Ffindinputdata(nexp+".001",prefix=prefix,fo=fo)
                    sssfiles = diridat+fo+"/????/"+prefix+"."+nexp+"*"+fity+".nc"
        else:
            sssfiles = diri+fo+"/????/"+prefix+"."+nexp+"*"+fity+".nc"
    return sssfiles



def Ftrpolyfit2(xrar,obs=False):
    """ Compute linear regression(y = a*t + b) from input data.
    
    Based on  http://atedstone.github.io/rate-of-change-maps/

    Parameters:
    - xar (xarray): input data to compute the linear regression from.
    - obs (str): True if the input data is from obs dataset. Default is False (model dataset). 
        Will change the name of the time coordinate to read (time_counter for model data, and time for obs).
        The obs option is not working yet.
    
    Returns:
    xrtrends:  trend coefficient (xarray). This is the a value from the linear regression in y = a*t + b.
    xrorigins: origin value (xarray). This is the b value from the linear regression in y = a*t + b.
    years: time index in years (xarray). Can be decimal years if input data was monthly.

   """
   
    vals = xrar.values 
    
    if obs==False:
            if hasattr(xrar.time_counter.to_index(), 'calendar'):
                if (xrar.time_counter.to_index().calendar=='noleap'):
                    timeJulian = xrar.time_counter.to_index().to_datetimeindex().to_julian_date()
            else:
                timeJulian = xrar.time_counter.to_index().to_julian_date()

    else:
            timeJulian = xrar.time.to_index().to_julian_date()
        
    # Reshape to an array with as many rows as years and as many columns as there are pixels
    vals2 = vals.reshape(len(timeJulian), -1)

    # Do a first-degree polyfit
    regressions = np.polyfit(timeJulian, vals2, 1)

    # Get the coefficients back
    trends     = regressions[0,:].reshape(vals.shape[1], vals.shape[2])
    origins    = regressions[1,:].reshape(vals.shape[1], vals.shape[2])
    
    # if input  data is from model, the coordinates are x,y
    coordsname = ['x','y']
    
    # if input  data is from ESA-CCI obs, the coordinates are changed to lon,lat
    if obs==True:
        coordsname = ['lon','lat']
        
     # make it back to xarray format
    xrtrends = xr.DataArray(
            data=trends,
            dims=[coordsname[1], coordsname[0]],
            attrs=dict(
                description="trend coeff",
                units="(g/kg)/day"
            ),
    )
    
    xrorigins = xr.DataArray(
            data=origins,
            dims=[coordsname[1], coordsname[0]],
            attrs=dict(
                description="b value",
                units="(g/kg)"
            ),
    )
    return xrtrends,xrorigins,timeJulian


def Ftrseries2(xrar,obs=False):
    """Compute linear trend timeseries
    
    This function calls the Ftrpolyfit function (that computes linear coefficient and origin value, from initial data)
    and uses it to build the trend timeseries for the time values (in Julian days years). 
    
    Parameters:
    - obs (str): True if the input data is from obs dataset. Default is False (model dataset). If obs is True, 
      the time coordinate is 'time', if dat is False, the time coordinate is 'time_counter' (model dataset).
      
    Usage:
       ```
       tr,xrtrends,xrorigins = Ftrseries2(data)
       data_dt = data - tr
       data.isel(x=1050,y=1050).plot()
       tr.isel(x=1050,y=1050).plot()
       ```

    Returns:
    xarray: tr   : trend timeseries over same time period as input data
    xarray: xrtrends  : trend coefficient in unit/year  
    xarray: xrorigins :origin coefficient in same unit as input data.

    """
    
    # calls Ftrpolyfit2 to compute the linear coeffs a and b   (y=ax+b)
    xrtrends,xrorigins,timeJulian = Ftrpolyfit2(xrar,obs=False)

    ic=-1
    # loop on timesteps from input data
    for itj in timeJulian:
        ic = ic+1
        lintr  = itj*xrtrends
        lintro = xrorigins

        # concatenate values in same array
        if obs==False:
            tmp = lintro.expand_dims(dim='time_counter',axis=0) + lintr.expand_dims(dim='time_counter',axis=0)
            if (ic!= 0):
                tr = xr.concat([tr, tmp], dim="time_counter")
            else:
                tr = tmp
        else:
            tmp = lintro.expand_dims(dim='time',axis=0) + lintr.expand_dims(dim='time',axis=0)
            if (ic!= 0):
                tr = xr.concat([tr, tmp], dim="time")
            else:
                tr = tmp
    if obs==False:
        tr = tr.assign_coords({"time_counter": xrar.time_counter})
    else:
        tr = tr.assign_coords({"time": xrar.time})
        
    # convert trend values from unit/day to unit/year    
    xrtrends = xrtrends*365.25
    xrtrends['units'] = "(g/kg)/year"
        
    return tr,xrtrends,xrorigins


def Ftrpolyfit(xrar,fty='1y',dat=False):
    """ Compute linear regression(y = a*t + b) from input data.
        --- NOW REPLACED BY Ftrpolyfit2 ---
    
    Based on  http://atedstone.github.io/rate-of-change-maps/
    Usage: you need to convert to np array the input data before applying this polyfit function.

    Parameters:
    - xar (numpy array): input data to compute the linear regression from.
    - fty (str): frequency of the input data ('1y','1m'). Doesn't work yet for bi-monthtly obs 'bimo'.
    - dat (str): True if the input data is from obs dataset. Default is False (model dataset). Will change the name of the time coordinate to read (time_counter for model data, and time for obs)
    
    Returns:
    xrtrends:  trend coefficient (xarray). This is the a value from the linear regression in y = a*t + b.
    xrorigins: origin value (xarray). This is the b value from the linear regression in y = a*t + b.
    years: time index in years (xarray). Can be decimal years if input data was monthly.

   """
   
    vals = xrar.values 
    if fty=='1y':
        if dat==False:
            years = xrar.time_counter.to_index().year
        else:
            years = xrar.time.to_index().year
    if fty=='1m':
        if dat==False:
            ys = xrar.time_counter.to_index().year.values   # !!! need to check if .values should be removed
            months = xrar.time_counter.to_index().month.values # !!! need to check if .values should be removed
            years = ys+(1./12.)*months
        else:
            ys = xrar.time.to_index().year.values  # !!! need to check if .values should be removed
            months = xrar.time.to_index().month.values  # !!! need to check if .values should be removed
            years = ys+(1./12.)*months

    if fty=='bimo':
        if dat==False:
            sys.exit("Sorry. The bi-monthly frequency is not implemented yet. Only yearly and monthly ('1y' and '1m')")
        else:
            ys = xrar.time.to_index().year.values
            months = xrar.time.to_index().month.values
            days = xr.DataArray(xrar.time.to_index().day.values)
            days0 = days*0.
            days15=days0.where(days==1,0.5)            
            years = ys+(1./12.)*months+days15

    if (((fty!='1m')&(fty!='1y'))&(fty!='bimo')): 
        sys.exit("Sorry. Only yearly, monthly and bi-monthly are implemented yet ('1y' and '1m' and 'bimo')")    

    # Reshape to an array with as many rows as years and as many columns as there are pixels
    vals2 = vals.reshape(len(years), -1)

    # Do a first-degree polyfit
    regressions = np.polyfit(years, vals2, 1)

    # Get the coefficients back
    trends     = regressions[0,:].reshape(vals.shape[1], vals.shape[2])
    origins    = regressions[1,:].reshape(vals.shape[1], vals.shape[2])

    # if input  data is from model, the coordinates are x,y
    coordsname = ['x','y']
    # if input  data is from ESA-CCI obs, the coordinates are changed to lon,lat
    if dat==True:
        coordsname = ['lon','lat']
        
     # make it back to xarray format
    xrtrends = xr.DataArray(
            data=trends,
            dims=[coordsname[1], coordsname[0]],
            attrs=dict(
                description="trend coeff",
                units="(g/kg)/year"
            ),
    )
    
    xrorigins = xr.DataArray(
            data=origins,
            dims=[coordsname[1], coordsname[0]],
            attrs=dict(
                description="b value",
                units="(g/kg)"
            ),
    )
    return xrtrends,xrorigins,years


def Ftrseries(xrtrends,xrorigins,ty,tc,dat=False):
    """Compute linear trend timeseries from linear regression coeffs.
        --- NOW REPLACED BY Ftrseries2 ---
    
    This function takes the result from the Ftrpolyfit function (linear coefficient, origin value, and time in years)
    and uses it to build the trend timeseries for the ty time values (dcimal years). 
    In the output array, the time coordinate is replaced by time values tc so that they are the same as in initial dataset.
    
    Parameters:
    - xrtrends:  trend coefficient (xarray). This is the a value from the linear regression in y = a*t + b.
    - xrorigins: origin value (xarray). This is the b value from the linear regression in y = a*t + b.
    - ty: time index in years (xarray). Can be decimal years if input data was monthly.
    - tc: time index to use in the output detrended data (xarray).
    - dat (str): True if the input data is from obs dataset. Default is False (model dataset). If dat is True, 
      the time coordinate is 'time', if dat is False, the time coordinate is 'time_counter' (model dataset).

    Returns:
    int: Returning value

    """
    ic=-1
    for iy in ty:
        ic = ic+1
        lintr  = iy*xrtrends
        lintro = xrorigins

        if dat==False:
            tmp = lintro.expand_dims(dim='time_counter',axis=0) + lintr.expand_dims(dim='time_counter',axis=0)
            if (ic!= 0):
                tr = xr.concat([tr, tmp], dim="time_counter")
            else:
                tr = tmp
        else:
            tmp = lintro.expand_dims(dim='time',axis=0) + lintr.expand_dims(dim='time',axis=0)
            if (ic!= 0):
                tr = xr.concat([tr, tmp], dim="time")
            else:
                tr = tmp
    if dat==False:
        tr = tr.assign_coords({"time_counter": tc})
    else:
        tr = tr.assign_coords({"time": tc})
    return tr


def Fpltaddfeatures(ax,incrgridlon,incrgridlat,onecohrml='#2E2E2E',alphaland=1,reg='GLO',landedgeco='none'):
    """Add my favourite physical features on map and remove map edges.
    
    Add rivers with 50m scale, and land with color fill.

    Parameters:
    ax (ax): Ax properties from matplotlib.pyplot.axes.
    reg (str): region to plot (default is 'GLO' for global)
    incrgridlon (int): increment of longitude in degrees
    incrgridlat (int): increment of latitude in degrees
    onecohrml (str): htlm color code for land.
    landedgeco (str): html colorcode for land edges.
    
    
    Returns:
    ax: ax (axes properties)
    gl: gl (gridlines properties)

   """
    
    rivers = cartopy.feature.NaturalEarthFeature(
        category='physical', name='rivers_lake_centerlines',
        scale='50m',facecolor='none',edgecolor='b')
    
    # Not used
    #lands = cartopy.feature.NaturalEarthFeature(
    #    category='physical', name='coastline',
    #    scale='50m',facecolor='none',edgecolor='k')

    cl2 = ax.add_feature(cfeature.LAND.with_scale('50m'),facecolor=onecohrml,edgecolor=landedgeco ,linewidth=0.2,alpha=alphaland,zorder=2)
    clr = ax.add_feature(cartopy.feature.RIVERS,alpha=0.5,facecolor='none',edgecolor='#A9E2F3',zorder=6)
    clr2 = ax.add_feature(rivers,alpha=0.5,facecolor='none',edgecolor='#A9E2F3',zorder=6)  ##CEE3F6

    #========= GRIDLINES
    gl =ax.gridlines(xlocs=range(-180,181,incrgridlon), ylocs=range(-90,91,incrgridlat),draw_labels=True,linewidth=1, color='#585858', alpha=0.3, linestyle='--',zorder=8)
    label_style = {'size': 12, 'color': '#BDBDBD', 'weight': 'normal'}
    gl.xlabel_style = label_style
    #gl.xlabels_bottom = False
    gl.ylabel_style = label_style
    #gl.ylabels_right = False

    return ax,gl



def Fshortenvarna(nami):
    """
    Defines short name of the variable to use in file names and to add to plot as string.

    Parameters:
    - nami (str): input var name

    Returns:
    str:  short name

   """
    if nami in ('sossheig','sosaline','sornf','sosstsst','stdev_sossheig','stdev_sosaline','stdev_sosstsst'):
        if nami=='sornf':
            shortvar='RNF'
        if nami=='sosaline':
            shortvar='SSS'
        if nami=='sossheig':
            shortvar='SSH'
        if nami=='sosstsst':
            shortvar='SST'
        if nami=='stdev_sossheig':
            shortvar='STD SSH'
        if nami=='stdev_sosstsst':
            shortvar='STD SST'
        if nami=='stdev_sosaline':
            shortvar='STD SSS'
    else:
        shortvar = nami
    return shortvar


def Fmycolormap(levbounds,cm_base='Spectral_r',cu='w',co='k',istart=0):
    """
    Makes a normalized colormap from min max values and a base colormap.
    
    Also add requested colors for under min and above max values.

    Parameters:
    - levbounds (array): [min, max, incr] min, max and increment of the requested colormap.
    - cm_base (str): name of the base colormap (see colormap list from matplotlib and cmocean).
    - cu (str): color code for under min values.
    - co (str): color code for above max values.
    - istart (int): trick to shift the colormap indices not to start from first color.

    Returns:
    str:  short name

   """
    lmin = levbounds[0]
    lmax = levbounds[1]
    incr = levbounds[2]
    levels = np.arange(lmin,lmax,incr)
    if ( (cm_base=='NCL') | (cm_base=='MJO') | (cm_base=='NCL_NOWI') ):
        nice_cmap = Fmake_SLXcolormap(whichco=cm_base)
    else:
        nice_cmap = plt.get_cmap(cm_base)
    colors = nice_cmap(np.linspace(istart/len(levels),1,len(levels)))[:]
    cmap, norm = from_levels_and_colors(levels, colors, extend='max')
    cmap.set_under(cu)
    cmap.set_over(co)
    return cmap,norm,levels


def Fmake_cmap(colors, position=None, bit=False):
    ''' Produce a colormap with equally spaced colors from a list of RGB colors
    
    This function  takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors. Arrange your tuples so that the first color 
    is the lowest value for the colorbar and the last is the highest.

    Parameters:
    - colors (array): Colors to use in the colormap. The RGB values may either be 
    in 8-bit [0 to 255] (in which bit must be set to True when called) or arithmetic [0 to 1] (default).
    - position (array):  Contains values from 0 to 1 to dictate the location of each color.
    - bit (bool): Set to True if RGB colors are given as 8-bit [0 to 255]. Default is False when RGB 
    colors are given as [0 to 1].
    
    Returns:
    dict:  colormap
    
    '''
    
    import matplotlib as mpl
    import numpy as np
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap



def Fmake_SLXcolormap(reverse=False,whichco='MJO',r=0,g=0,b=0):
    """Define a custom cmap 

    Parameters:
    - Reverse (bool): If true, will create the reverse colormap. Default is False.
    - whichco (str): Some predined colormaps ('MJO', 'NCL', 'NCL_NOWI' bluyello' 'oneco' 'combi'). Default is 'MJO'.
    - r (int): R value from RGB color for a one-color colormap.
    - g (int): G value from RGB color for a one-color colormap.
    - b (int): B value from RGB color for a one-color colormap.

    Returns:
    dict: colormap

   """

    ### colors to include in my custom colormap
    if whichco=='MJO':
        colors_NCLbipo=[(176,17,3,1),(255,56,8,1),(255,196,1,1),(255,255,255,1),(255,255,255,1),(13,176,255,1),(2,88,255,1),(0,10,174,1)]

    if whichco=='NCL':
        colors_NCLbipo=[(11,76,95),(0,97,128),(0,161,191),(0,191,224),(0,250,250),(102,252,252),(153,250,250),(255,255,255),(255,255,255),(252,224,0),(252,191,0),(252,128,0),(252,64,0),(252,33,0),(128,0,0),(0,0,0)]

    if whichco=='NCL_NOWI':
        colors_NCLbipo=[(11,76,95),(0,97,128),(0,161,191),(0,191,224),(0,250,250),(102,252,252),(153,250,250),(255,255,255),(252,224,0),(252,191,0),(252,128,0),(252,64,0),(252,33,0),(128,0,0),(0,0,0)]
        
    if whichco=='bluyello':
        colors_NCLbipo=[(11,76,95),(0,97,128),(0,161,191),(0,191,224),(0,250,250),(102,252,252),(153,250,250),(255,255,255),(252,224,0),(252,191,0),(252,128,0),(252,64,0),(252,33,0),(128,0,0),(0,0,0)]
        
    if whichco=='oneco':
        colors_NCLbipo=[(r,g,b),(0,0,0)]
        
    if whichco=='combi':
        colors_NCLbipo=[(0, 0, 0, 255),(24, 24, 34, 255),(49, 48, 68, 255),(74, 74, 103, 255),(98, 105, 130, 255),(124, 140, 156, 255),(148, 173, 180, 255),(178, 206, 205, 255),(216, 230, 230, 255),(255, 255, 255, 255),(247, 209, 60, 255), (251, 155, 6, 255),(237, 104, 37, 255),(207, 68, 70, 255),(164, 44, 96, 255),(120, 28, 109, 255),(74, 11, 106, 255),(26, 11, 64, 255),(0, 0, 3, 255)]


    ### Call the function make_cmap which returns my colormap
    my_cmap_NCLbipo = Fmake_cmap(colors_NCLbipo[:], bit=True)
    my_cmap_NCLbipo_r = Fmake_cmap(colors_NCLbipo[::-1], bit=True)
    
    if reverse==True:
        my_cmap_NCLbipo = Fmy_cmap_NCLbipo_r

    return(my_cmap_NCLbipo)


def Faddcolorbar(fig,cs,ax,levbounds,levincr=1,tformat="%.2f",tlabel='',shrink=0.45,facmul=1.,orientation='vertical',tc='k',loc='lower right',wth="15%",bbta=(0.08, -0.1,0.9,0.2)):
    """Function to add a customized colorbar on a figure

    Parameters:
    - fig: figure properties
    - cs: plotted field properties
    - ax: ax properties
    - levbounds (array): [min, max, incr] min, max and increment of the requested colormap.
    - levincr (int): increment at which to display colorbar labels.
    - tformat (str): default is "%.2f"
    - tlabel (str): label to display near the colorbar.
    - shrink (float): Proportion in [0,1] for the relative size of the colorbar compared to axis.
    - facmul (float): Multiplication factor if you want to display the colorbar label in a different unit from original unit in the data.     
    - orientation (str): 'horizontal' or vertical'
    - tc (str): color of text labels and ticks
    - loc (str): position of the colorbar (default is 'lower right')
    - wth (str): percentage setting the height of a horizontal colorbar
    - bbta (tuple): bounding box defining the location of the colorbar Default is (0.08, -0.1,0.9,0.2).

    Returns:
    dict:colorbar properties (cb).
   """
    
    lmin = levbounds[0]
    lmax = levbounds[1]
    incr = levbounds[2]
    levels = np.arange(lmin,lmax,incr)
    cblev = levels[::levincr]
    
    if orientation =='horizontal':
        axins1 = inset_axes(ax,
                        height=wth,  # height : 5%
                            width="50%",
                        bbox_to_anchor=bbta,
                        bbox_transform=ax.transAxes,
                        borderpad=0)

    if orientation =='vertical':
        axins1 = inset_axes(ax,
                        height="50%",  # height : 5%
                            width="2%",
                        loc='center left',
                       borderpad=2)

    cb = fig.colorbar(cs,cax=axins1,
                                    extend='both',                   
                                    ticks=cblev,
                                    spacing='uniform',
                                    orientation=orientation,
                                    )
    
    new_tickslabels = [tformat % i for i in cblev*facmul]
    cb.set_ticklabels(new_tickslabels)
    cb.ax.set_xticklabels(new_tickslabels, rotation=70,size=10,color=tc)
    cb.ax.tick_params(labelsize=10,color=tc) 
    cb.set_label(tlabel,size=14,color=tc)
    
    return cb
    

def Fdefnamoplt(imhov,reg,ty,suffix="",imisc=0):
    """
    Defines output name of plot for the given parameters

    Parameters:
    - imhov (imhov instance)
    - reg (str):
    - ty (str): type of plot (DIFF, STD, TR, TM) for difference, Std, Trend, Time-Mean
    - suffix (str): More string to add at the end of the name if needed. Default is none.

    Returns:
    str:  output name of plot

   """
    addna=""
    if imhov.lev==-1:
        pass
    else:
        #suffix="LEV-"+str(int(np.floor(imhov.levm)))+"m_"+suffix
        suffix="LEV-"+suffix
    
    if (((ty=='DIFF')|(ty=='TRDIFF'))|(ty=='TMDIFF')):
        addna="-"+imhov.nexp2
    namo = "JZmap_"+imhov.var+"_"+reg+"_"+ty+"_"+imhov.y1+"-"+imhov.y2+"_"+imhov.nexp+addna+"_"+imhov.fo
     
    #if (ty=="MISC"):
    #    namo = "JZmap_"+imhov.var+"_"+reg+"_"+ty+"_"+imhov.y1+"-"+imhov.y2+"_"+imhov.nexp[4]+"_"+imhov.fo
    if ((ty=="EM")):
        namo = "JZmap_"+imhov.var+"_"+reg+"_"+imhov.nexp+"_"+imhov.fo+"_"+format(imisc, '04d')
    if ((ty=="ESPR")|(ty=="EM")):
        namo = "JZmap_"+ty+"_"+imhov.var+"_"+reg+"_"+imhov.nexp+"_"+imhov.fo+"_"+format(imisc, '04d')
    if (ty=="EMDIFF"):
        namo = "JZmap_"+ty+"_"+imhov.var+"_"+reg+"_"+imhov.nexp+"-"+imhov.nexp2+"_"+imhov.fo+"_"+format(imisc, '04d')
    if suffix=="":
        pass
    else:
        namo = namo +"_"+suffix
    return namo


def Fsaveplt(fig,diro,namo,dpifig=300):
    """Save plot

    Parameters:
    fig: Fig properties to save.
    diro (str): output directory
    namo (str): name of output plot
    dpifig (int): resolution (dpi) of saved plot.
   """
    fig.savefig(diro+namo+".png", facecolor=fig.get_facecolor(),
                edgecolor='none',dpi=dpifig,bbox_inches='tight', pad_inches=0)
    plt.close(fig) 
    
    
def Fpltsetcolorsdefaults(ty= "defaults",
    levbounds=[-1,1.05,0.05],
    cbincr=10,
    cm_base='Spectral_r',
    cu='k',
    co='r' ,   
    onecohrml='#F2F2F2',
    alphaland=1,
    landedgeco='#6E6E6E',
    # rgb color for model land grid points
    rgb=[242,242,242],
    cbformat="%.2f"):
    
    """Set my favourite default color parameters for plotting .

    Print it to remember what are the parameters.
    pltcolparam = Fpltsetcolorsdefaults
    print(pltcolparam)
    
    Parameters:
    - levbounds (array): [min, max, incr] min, max and increment of the requested colormap.
    - cbincr (int): increment at which to display colorbar labels.
    - cm_base='Spectral_r',
    - cu='k',
    - co='r' ,   
    - onecohrml (str):html colorbar
    - landedgeco (str):html colorbar
    - rgb (array): rgb color for model land grid points
    - cbformat (str): Format to print the colorbar labels (how manyDefault is "%.2f".
    
    Returns: 
    dict:  color parameters of the plot

   """
    pltcolparam = {'type':ty,'levbounds':levbounds,'cbincr':cbincr,'cbformat':cbformat,'cm_base':cm_base,'cu':cu,'co':co,'onecohrml':onecohrml,'alphaland':alphaland,'landedgeco':landedgeco,'rgb':rgb}
    return pltcolparam

def Fsetcirclesparamdefaults(
    alpha=1,
    mulmag=1000,
    linewidths=0.3,
    linewidthstr=0.6,
    marker='o',
    facecolors = 'none',
    edgecolors='#08088A',
    edgecolorstrp='#FA5858',
    edgecolorstrn='#01DFD7'):
    
    """Set my favourite default  parameters for plotting circles at river mouths.

    Print it to remember what are the parameters.
    pltcirclesparam = Fsetcirclesparamdefaults
    print(pltcirclesparam)
    
    Parameters:
    alpha (float): transparency [0-1]
    mulmag (int): relative size of circles
    linewidths (float):
    linewidthstr (float):
    marker (str):
    facecolors (str): html color code
    edgecolors (str): html color code
    edgecolorstrp (str): html color code for first group of circles (positive values)
    edgecolorstrn (str): html color code for second group of circles (negative values)
    """
    pltcirclesparam={'alpha':alpha,'mulmag':mulmag,'linewidths':linewidths,'linewidthstr':linewidthstr,'marker':marker,'facecolors':facecolors,'edgecolors':edgecolors,'edgecolorstrp':edgecolorstrp,'edgecolorstrn':edgecolorstrn}
    return pltcirclesparam

def Faddmarkers(trdata,pltmarkersparam):
    """
    Add scatter plot of markers on top of map.

    """ 
    m=plt.scatter(x=pltmarkersparam['x'],
    y=pltmarkersparam['y'],
    alpha=pltmarkersparam['alpha'], #1,
    s=pltmarkersparam['mulmag'],#700
    linewidths=pltmarkersparam['linewidths'], #0.3,
    marker=pltmarkersparam['marker'], #,
    facecolors=color_list_big,#pltmarkersparam['facecolors']color_list,#, #none, 
    edgecolors=color_list_big, #pltmarkersparam['edgecolors'], #'#08088A',
    transform=trdata,
    cmap=color_list,
    zorder=20) 
            
def Fsetmarkersparamdefaults(
    x=0,
    y=5,
    alpha=1,
    mulmag=100,
    linewidths=2,
    marker='x',
    facecolors = '#e74c3c',
    edgecolors='#e74c3c'):
    
    """Set my favourite default  parameters for plotting markers on map

    Print it to remember what are the parameters.
    
    Parameters:
    x (array): x positions,
    y (array): y positions,
    alpha (float): transparency [0-1]
    mulmag (int): relative size of circles
    linewidths (float):
    marker (str):
    facecolors (str): html color code
    edgecolors (str): html color code
    """
    pltmarkersparam={'x':x,'y':y,'alpha':alpha,'mulmag':mulmag,'linewidths':linewidths,'marker':marker,'facecolors':facecolors,'edgecolors':edgecolors}
    return pltmarkersparam

def Fpltgridparamdefaults(reg='GLO',gridl=False,incrgridlon=45,incrgridlat=30,sath=35785831,minlat=-90,maxlat=90,minlon= -180,maxlon=180,loncentr=0,latcentr=10):
    """Set my favourite default  parameters for plotting grid and project map

    Parameters:
    argument1 (int): Description of arg1
    reg (str): region to plot (some regions are predefined: 'GLO' 'gro' 'atl' 'tropatl' 'asia' 'ind' 'bof')
    gridl (bool): not sure if this is still used?
    incrgridlon (int): increment to plot lon labels on map
    incrgridlat (int): increment to plot lon labels on map
    sath (float): height of satellite in specific projection (name?)
    minlat (float):
    maxlat (float):
    minlon (float):
    maxlon (float):
    loncentr (float): lon center of the map in degrees
    latcentr (float): lat center of the map in degrees

    Returns:
    dict : parameters

   """
    # default gridlines parameters
    axextent = [minlon, maxlon, minlat, maxlat]    #

    if reg=='gro':
        loncentr=-35
        latcentr=75
        sath=2085831
        axextent = [minlon, maxlon, minlat, maxlat]

    if reg=='atl':
        loncentr=-35
        latcentr=10
        sath=35785831
        minlat=-20.0
        maxlat=35.0
        minlon= -100
        maxlon=-20
        axextent = [minlon, maxlon, minlat, maxlat]
        
    if reg=='tropatl':     
        loncentr=-50
        latcentr=10
        sath=30000000 
        axextent = [-80, 20, -10, 10]
        minlat=axextent[2]
        maxlat=axextent[3]

    if reg=='asia':       
        sath=30000000
        latcentr=5
        loncentr=110

    if reg=='bob':
        loncentr=90
        latcentr=0
        sath=35785831
        axextent = [77, 98, 8, 25]
        minlat=axextent[2]
        maxlat=axextent[3]
        
    if reg=='bob2':
        loncentr=90
        latcentr=0
        sath=35785831
        axextent = [77, 98, 12, 25]
        minlat=axextent[2]
        maxlat=axextent[3]
        
    if reg=='ind':
        loncentr=90
        latcentr=0
        sath=35785831
        axextent = [30, 130, -21, 35]
        minlat=axextent[2]
        maxlat=axextent[3]

    if reg=='GLO':
        minlat=-70.0
        maxlat=84.0
        minlon= -180
        maxlon=179
        loncentr=0
        incrgridlon=60
        incrgridlat=30 
        axextent = [minlon, maxlon, minlat, maxlat]
        
    pltparam = {'reg':reg,'gridl':gridl,'incrgridlon':incrgridlon,'incrgridlat':incrgridlat,'sath':sath,'minlat':minlat,'maxlat':maxlat,'minlon':minlon,'maxlon':maxlon,'loncentr':loncentr,'latcentr':latcentr,'axextent':axextent}
        
    return pltparam



def Fpltsetcolorsdefaults(ty):
    """
    Set my favourite default color parameters for plot given type of plot.

    Parameters:
    - ty (str): type of plot (DIFF, STD, TR, TM) for difference, Std, Trend, Time-Mean
    
    Returns: (dict) 

   """

    if ty=='default':
        # defaults
        levbounds=[-1,1.05,0.05]
        cbincr=10
        cm_base='Spectral_r'
        cu='k'
        co='r'    

    if ty=='DIFF':
        levbounds=[-0.75,0.77,0.02]
        cbincr=5
        cm_base=Fmake_SLXcolormap(reverse=False,whichco='NCL_NOWI')
        cu='#080449'
        co='#5b2123'
        
    if ty=='TM':
        levbounds=[29,39,0.1]
        cm_base=cmocean.cm.haline
        cu='k'
        co='#FFFF00'
        cbincr=10
        
    if ty=='TMDIFF':
        levbounds=[-1.,1.0,0.05]
        cm_base=Fmake_SLXcolormap(reverse=False,whichco='NCL_NOWI')
        cu='#080449'
        co='#5b2123'
        cbincr=2
    
    if ty=='STD':
        levbounds=[0.05,1.5,0.05]
        cbincr=2
        cm_base='inferno_r'
        cu='w'
        co='k'
        
    if ty=='TR':
        levbounds=[-0.1,0.105,0.005]
        cbincr=2
        cm_base='RdYlBu_r'
        cu='#080449'
        co='#5b2123'
        
    if ty=='TRDIFF':
        levbounds=[-0.1,0.105,0.005]
        cbincr=2
        cm_base='RdYlBu_r'
        cu='#080449'
        co='#5b2123'
        

    onecohrml='#F2F2F2'
    alphaland=1
    landedgeco='#6E6E6E'
    # rgb color for model land grid points
    rgb=[242,242,242]
    cbformat="%.2f"
    
    pltcolparam = {'type':ty,'levbounds':levbounds,'cbincr':cbincr,'cbformat':cbformat,'cm_base':cm_base,'cu':cu,'co':co,'onecohrml':onecohrml,'alphaland':alphaland,'landedgeco':landedgeco,'rgb':rgb}
    
    return pltcolparam


def FselectRNF(sortrnf):  
    """Select points of an xarray when value at this point is above/below a given threshold,
    and output as a dictionnary that contains a  list with lat, lon and value of the selected points. 

    Parameters:
    sortrnf (dict): containing 'out_stack': a 1-d vector with the values to browse, 'nav_lat_stack' and 'nav_lon_stack' the corresponding lat,lon position (1-D arrays), 'threshold' the threshold used to select values above it, 'ty' is a string defining how to select the values compared to threshold value. If 'ty' is 'tr' 2 threshold values are read to select larger values than positive threshold and lower values than negative threshold.  Otherwise only 1 threshhold value is considered and values are selected when above this threshlod.

    Returns:
    dict:

    """
    xrtrends_stack=sortrnf['out_stack']
    nav_lat_stack=sortrnf['nav_lat_stack']
    nav_lon_stack=sortrnf['nav_lon_stack']  
    threshold = sortrnf['threshold'] 
    
    if sortrnf['ty']=='tr':
        rntrposiselect = xrtrends_stack.where(((xrtrends_stack>threshold[1])&((xrtrends_stack>0.))),drop=True)
        latselecttrp = nav_lat_stack.where(((xrtrends_stack>threshold[1])&((xrtrends_stack>0.))),drop=True)
        lonselecttrp = nav_lon_stack.where(((xrtrends_stack>threshold[1])&((xrtrends_stack>0.))),drop=True)

        rntrnegaselect = xrtrends_stack.where(((xrtrends_stack<threshold[0])&((xrtrends_stack<0.))),drop=True)
        latselecttrn = nav_lat_stack.where(((xrtrends_stack<threshold[0])&((xrtrends_stack<0.))),drop=True)
        lonselecttrn = nav_lon_stack.where(((xrtrends_stack<threshold[0])&((xrtrends_stack<0.))),drop=True)

        # for plot purposes, shift lon lat so that they are in the center of the gridpoint instead of on left side
        latselecttrp = latselecttrp+0.125
        lonselecttrp = lonselecttrp+0.125
        latselecttrn = latselecttrn+0.125
        lonselecttrn = lonselecttrn+0.125

        # set magnitude for latter plot (magnitude relative to max runoff)
        magnitudetrp = (rntrposiselect.values)/(rntrposiselect.max().values)
        magnitudetrn = (rntrnegaselect.values)/(-1.*rntrposiselect.max().values)    

        selrnf = {'ty':sortrnf['ty'],'p':sortrnf['p'],'latselecttrn':latselecttrn,'lonselecttrn':lonselecttrn,'latselecttrp':latselecttrp,'lonselecttrp':lonselecttrp,'magnitudetrp':magnitudetrp,'magnitudetrn':magnitudetrn}

    #if ((sortrnf['ty']=='tm')|(sortrnf['ty']=='std')):
    else:
        rnfselect = xrtrends_stack.where(xrtrends_stack>threshold,drop=True)
        latselect = nav_lat_stack.where(xrtrends_stack>threshold,drop=True)
        lonselect = nav_lon_stack.where(xrtrends_stack>threshold,drop=True)

        # for plot purposes, shift lon lat so that they are in the center of the gridpoint instead of on left side
        latselect = latselect+0.125
        lonselect = lonselect+0.125

        # set magnitude for latter plot (magnitude relative to max runoff)
        magnitude = (rnfselect.values)/(rnfselect.max().values)

        selrnf = {'ty':sortrnf['ty'],'p':sortrnf['p'],'latselect':latselect,'lonselect':lonselect,'magnitude':magnitude}
        
    return selrnf
    
    
def Fdeftlabel(indat,ty,strunits,fo='1d',imisc=0):
    """
    Define colorbar textlabel of given plot

    Parameters:
    - indat: imhov class instance
    - strunits (str): variable unit
    
    Return: (str)
    """
    
    un=""
    le=""

    if strunits=="":
        pass
    else:
        un=" ("+strunits+")"

    if indat.lev==-1:
        pass
    else:
        le=" ("+str(int(np.floor(indat.levm)))+" m"+")"

    tlabel=ty+" "+indat.var+le+un+" over "+indat.y1+"-"+indat.y2
    
    if (((ty=="ESPR")|(ty=="EM"))|(ty=="EMDIFF")):
        tyy = indat.data.time_counter.isel(time_counter=imisc).to_dict()['data'].year
        tmm = indat.data.time_counter.isel(time_counter=imisc).to_dict()['data'].month
        tdd = indat.data.time_counter.isel(time_counter=imisc).to_dict()['data'].day
        tlabel="Ensemble "+indat.var+le+un+" on "+str(tyy)+"-"+str(tmm)+"-"+str(tdd)
        if (ty=="EMDIFF"):
            if (fo=='1y'):
                tlabel="DIFF ens. "+indat.var+le+un+" in "+str(tyy)
            else:
                tlabel="DIFF ens. "+indat.var+le+un+" on "+str(tyy)+"-"+str(tmm)+"-"+str(tdd)
            
    
    return tlabel


def Fnospines(ax):
    """
    Remove spines from around the plot
    
    Return: (ax)
    """
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.axis('off')

    
def Ftesttime(s1,s2,tdn='time_counter'):
    """ Test if timeseries have the same time coordinates
    Return (bool)
    """
    if (tdn=='time_counter'):
        testtime = s1.time_counter-s2.time_counter
    else:
        testtime = s1.time-s2.time

    if (testtime.sum()==0):
        is_same_time=True
    else:
        is_same_time=False
    return is_same_time


def FcomputeGM(data,e1,e2):
        """ Compute the global mean of a field at each time step.
            
            Parameters:
            - data (float): xarray of dimension x,y,t 
            - e1 (float): xarray of dimension x,y that is the size of the data mesh in x direction
            - e2 (float): xarray of dimension x,y that is the size of the data mesh in y direction
    
            Return: gm (float) xarray of dimention t containing the global mean value
        """
        
        weighted_data =  data * e1 
        weighted_data =  weighted_data * e2 
        stacked_weighted_data = weighted_data.stack(z=("x", "y"))
        weights = e1 * e2
        stacked_weights = weights.stack(z=("x", "y"))
        summed_weights = stacked_weights.sum(dim='z')
        gm = stacked_weighted_data.sum(dim="z") / summed_weights
        
        return gm
    
def FrmGM(data,gm):
        """ Remove the global mean from a field at each time step.
            
            Parameters:
            - data (float): xarray of dimension x,y,t 
            - gm (float): xarray of dimention t containing the global mean value
    
            Return: data_rmGM (float) xarray of dimention x,y,t containing the data field - gm at each time step
        """
        
        data_rmgm = data - gm.broadcast_like(data)
        
        return data_rmgm

## TODO in imhov class??? 
def Fdata2plot(indat,ty='STD',imisc=0,imhov=True):
    """ Indicate which data to plot as a function of plot type
            
            Parameters:
            - ty (str): plot type
            - indat (imhov class instance)
            - imisc (int): index at which to read data if several timeseteps
            - imhov (bool): True if indat is of class imhov. False if indat is just a simple xarray type.
    
            Return: 
            xarray: 2-D xarray to plot
    """
    if imhov:
        if ty== 'DIFF': 
            data2plot  = indat.std_diff
            print('Plot: imhov.std_diff')
        elif ty== 'STD': 
            data2plot  = indat.std_dt
            print('Plot: imhov.std_dt')
        elif ty== 'TM': 
            data2plot  = indat.tm
            print('Plot: imhov.tm')
        elif ty== 'TMDIFF': 
            data2plot  = indat.tm_diff
            print('Plot: imhov.tm_diff')
        elif ty== 'TR': 
            data2plot  = indat.atr
            print('Plot: imhov.atr')
        elif ty== 'TRDIFF': 
            data2plot  = indat.tr_diff
            print('Plot: imhov.tr_diff')
        elif ((ty== 'ESPR')|(ty== 'EM')): 
            data2plot  = indat.data.isel(time_counter=imisc) 
            print('Plot: imhov.data.isel(time_counter=imisc) ')
        elif (ty== 'EMDIFF'): 
            data2plot  = indat.data_diff.isel(time_counter=imisc)
            print('Plot: imhov.data_diff.isel(time_counter=imisc)')
        else:
            data2plot  = indat.data.mean(dim='time_counter')
    else:
        data2plot  = indat.isel(time_counter=imisc) 
        
    return data2plot
    

def Fpltcircles(selrnf,trdata,pltcirclesparam='no'):
    """
    Add scatter plot of circles on top of map.

   """ 

    # Plot type: time mean or STD: One group of  circles is plotted. All the  same color.  
    if ((selrnf['ty']=='tm')|(selrnf['ty']=='std')): 
            s1=plt.scatter(x=selrnf['lonselect'].values,
            y=selrnf['latselect'].values,
            alpha=pltcirclesparam['alpha'], #1,
            s=selrnf['magnitude']*pltcirclesparam['mulmag'],#700
            linewidths=pltcirclesparam['linewidths'], #0.3,
            marker=pltcirclesparam['marker'], #,
            facecolors=pltcirclesparam['facecolors'], #none, 
            edgecolors=pltcirclesparam['edgecolors'], #'#08088A',
            transform=trdata,
            cmap='inferno_r',
            zorder=20) 
            
    # Plot type: trend plot: Two groups of  circles are plotted. (negative and positive values)
    if (selrnf['ty']=='tr'): 

            s2=plt.scatter(x=selrnf['lonselecttrp'].values,
            y=selrnf['latselecttrp'].values,
            alpha=pltcirclesparam['alpha'], #1,
            s=selrnf['magnitudetrp']*pltcirclesparam['mulmag'],#700
            linewidths=pltcirclesparam['linewidthstr'], #0.3,
            marker=pltcirclesparam['marker'], #,
            facecolors=pltcirclesparam['facecolors'], #none, 
            edgecolors=pltcirclesparam['edgecolorstrp'], #'#08088A',
            transform=trdata,
            cmap='inferno_r',
            zorder=20) 

            s3=plt.scatter(x=selrnf['lonselecttrn'].values,
            y=selrnf['latselecttrn'].values,
            alpha=pltcirclesparam['alpha'], #1,
            s=selrnf['magnitudetrn']*pltcirclesparam['mulmag'],#700
            linewidths=pltcirclesparam['linewidthstr'], #0.3,
            marker=pltcirclesparam['marker'], #,
            facecolors=pltcirclesparam['facecolors'], #none, 
            edgecolors=pltcirclesparam['edgecolorstrn'], #'#08088A',
            transform=trdata,
            cmap='inferno_r',
            zorder=20) 
    
    
    
def FpltGLO(indat,pltgridparam,ty,reg='GLO',strunits="",pltcolparam='no',diro="./",siplt=False,saveo=False,pltshow=True,cbon=False,suffix="",selrnf="no",pltcircles=False,pltcirclesparam='no',fo='1d',imisc=0,imhov=True,cl=True,pltmarkersparam='no',**kwargs):
    """
    Plot Global map

    Parameters:
    - strunits (str): units of variable to be plotted (as to display on plot label)
    - pltgridparam (dict): grid parameters from function pltgridparam.
    - ty (str): type of plot (DIFF, STD, TR, TM) for difference, Std, Trend, Time-Mean
    - pltcolparam (dict): color parameters: Keys:
        - levbounds (vector of size 3) : Min, Max, increment of the colorscale.
        - cm_base : colormap
        - cu (str): html color code for under-values
        - co (str): html color code for over-values
        - onecohrml (str): html color code for lands
        - landedgeco (str): color code for land edges
        - rgb (vector of size 3): r,g,b colors for land cells in model
        - cbincr (int): incr of labels on colorbar
        - cbformat (str): String format of colorbar labels (determines the number of digits). Default is "%.2f"
    - siplt (bool): plot sea ice extent as hatches or not.
    - saveo (bool): save figure as png (default is False)
    - cbon (bool): add colorbar and labels (default is False).
    - suffix (str): string to optionnaly add to plot name.

    Returns: (str) path of output plot


   """

    #---  choose data to plot given the type of plot 'ty'
    if (indat.all()!=0):
        data2plot = Fdata2plot(indat,ty=ty,imisc=imisc,imhov=imhov)
        #---  take care of masking where no data
        # if indat is obs
        if indat.obs:
            maskobs = data2plot.isnull()
            if ty== 'STD':
                data2plot = data2plot.where(data2plot>0) 
        else:    
        # if indat is model data
            data2plot  = data2plot.where(indat.mask!=0).squeeze()  

            m2plt      = indat.mask.where(indat.mask==0)

            if siplt:
                maskedzone = data2plot.where(indat.siextentTM>0)
    else:
        # case where the globe is plotted without data
       print('Plotting map without data.')

    #--- set default color param
    if pltcolparam=='no':
        pltcolparam = Fpltsetcolorsdefaults(ty)


    #--- colormap
    cmap,norm,levels = Fmycolormap(pltcolparam['levbounds'],cm_base=pltcolparam['cm_base'],cu=pltcolparam['cu'],co=pltcolparam['co']) ##080449
    # rgb color for model land grid points
    r=pltcolparam['rgb'][0];g=pltcolparam['rgb'][1];b=pltcolparam['rgb'][2]
    # color for continents from data based (hi-res)
    onecohrml=pltcolparam['onecohrml']
    alphaland=pltcolparam['alphaland']


    #========= CREATE FIGURE
    fig3 = plt.figure(figsize=([18,10]),facecolor='white')

    #========= PLOT DATA

    # Data system proj (if coords are in lat lon, use PlateCarre here)
    trdata  = ccrs.PlateCarree() 

    if (reg=='GLO'):
        ax = plt.axes(projection= ccrs.Mercator(central_longitude=pltgridparam['loncentr'],min_latitude=pltgridparam['minlat'], max_latitude=pltgridparam['maxlat'], globe=None))
    else:
        ax = plt.axes(projection= ccrs.PlateCarree(central_longitude=pltgridparam['loncentr']))
    if (indat!=0):
        cs   = plt.pcolormesh(indat.nav_lon.squeeze(), indat.nav_lat.squeeze(), data2plot,shading='flat',cmap=cmap,transform=trdata,norm=norm)
        if indat.obs:
            csdots  = plt.contourf(indat.nav_lon.squeeze(), indat.nav_lat.squeeze(),maskobs.where(maskobs>=1),cmap='gray',hatches=['...'],transform=trdata,extend='both',alpha=0.)
        else:
            # add grey shading  where ocean mask is 0 (land gridpoints in the model)
            # rgb color for model land grid points
            csland  = plt.pcolormesh(indat.nav_lon.squeeze(), indat.nav_lat.squeeze(), m2plt, shading='flat',cmap=Fmake_SLXcolormap(reverse=False,whichco='oneco',r=r,g=g,b=b),transform=trdata)

        #.where(indat.siextentTM>0)

        if siplt:
            cs2  = plt.contourf(indat.nav_lon.squeeze(),indat.nav_lat.squeeze(),data2plot,cmap='gray',hatches=['//'],transform=trdata,extend='both',alpha=0.)
            cs3  = plt.pcolormesh(indat.nav_lon.squeeze(), indat.nav_lat.squeeze(), data2plot.where(maskedzone.isnull()),shading='flat',cmap=cmap,transform=trdata,norm=norm)
            cs4  = plt.pcolormesh(indat.nav_lon.squeeze(), indat.nav_lat.squeeze(), m2plt.where(maskedzone.isnull()), shading='flat',cmap=Fmake_SLXcolormap(reverse=False,whichco='oneco',r=r,g=g,b=b),transform=trdata)
    else:
        # case where the globe is plotted without data
        pass

    #========= GEOGRAPHICAL FEATURES
    # make plot nice with rivers, continents, grids:
    ax,gl = Fpltaddfeatures(ax,pltgridparam['incrgridlon'],pltgridparam['incrgridlat'],onecohrml=pltcolparam['onecohrml'],alphaland=pltcolparam['alphaland'],reg=pltgridparam['reg'],landedgeco=pltcolparam['landedgeco'])  ##585858#BDBDBD

    if (reg=='GLO'):
        Fnospines(ax)
    else:
        # geographical limits
        ax.set_extent([pltgridparam['axextent'][0],pltgridparam['axextent'][1],pltgridparam['axextent'][2],pltgridparam['axextent'][3]])
        # remove spines from around plot
        Fnospines(ax)
        
    #========= Add SCATTERPLOT OF RUNOFF VALUES   
    if pltcircles:
        # set plotting parameters if not done yet
        if pltcirclesparam=='no':
            pltcirclesparam=Fsetcirclesparamdefaults()
        # plot (scatterplot on top of the current map)
        ci = Fpltcircles(selrnf,trdata,pltcirclesparam=pltcirclesparam)
        # add suffix to name of output plot
        suffix="Circ-"+selrnf['ty']+"-"+str(int(100*selrnf['p']))+suffix
        
    #========= Add SCATTERPLOT OF SOME MARKERS 
    if pltmarkersparam=='no':
        pass
    else:
        mk = Faddmarkers(trdata,pltmarkersparam)

    #========= ADD COLORBAR
    if cbon:
        if (((ty=="ESPR")|(ty=="EM"))|(ty=="EMDIFF")):
            tlabel=Fdeftlabel(indat,ty,strunits,fo=fo,imisc=imisc)
        else:
            tlabel=Fdeftlabel(indat,ty,strunits,fo=fo)
        cb = Faddcolorbar(fig3,cs,ax,pltcolparam['levbounds'],levincr=pltcolparam['cbincr'],tformat=pltcolparam['cbformat'],
                         tlabel=tlabel,facmul=1,orientation='horizontal',tc='k',
                        bbta=(-0.18,-0.25,0.9,0.2))  
        
    #========= ZOOM EXTENT in LAT LON
    if (reg=='GLO'):
        pass
    else:
        #set limits of region
        ax.set_extent([pltgridparam['axextent'][0],pltgridparam['axextent'][1],pltgridparam['axextent'][2],pltgridparam['axextent'][3]])
        # remove spines from around plot
        Fnospines(ax)
        
    #========= PLT SHOW AND SAVE
    if pltshow:
        plt.show()

    if saveo:
        if (indat!=0):
            if (((ty=="ESPR")|(ty=="EM"))|(ty=="EMDIFF")):
                namo = Fdefnamoplt(indat,reg,ty,suffix=suffix,imisc=imisc)
            else:
                namo = Fdefnamoplt(indat,reg,ty,suffix=suffix)
        else:
            namo='examplemap'
            
        Fsaveplt(fig3,diro,namo,dpifig=300)
        outo=diro+namo+'.png'
        print('Plot saved in '+outo)
    else:
        outo="(saveo=False)"
        print('Plot not saved. '+outo)
        
    if cl:
        plt.close(fig3)

    return outo,fig3


#-----------------------------------------    
#-----------------------------------------    
class imhov:
    """The imhov class is a class to process imhotep data
    """

    def __init__(self, origin1, varna,nexp,fo,y1,y2,dirigrid='./',sif='no',diff=False,estatsdiff=False,mbdiff=False,origin2='no',nexp2='no',obs=False,lev=-1):
        """Initialize the imhov class

        Parameters:
        origin1 (str): path of the files to look at 
        varna (str): name of variable to look at
        nexp (str): name of experiment (GAI, S, GA, GI, AI, EGAI, ES, EAI, ESA (for obs dataset) )
        fo (str): frequency at which to look at the data ('1y' or '1m')
        y1 (int): year at which to start considering the data
        y2 (int): year at which to stop considering the data
        dirigrid (str): path of the file containing the grid info
        sif (str): sea ice file (if needed, default is 'no')
        diff (bool): set to True if the aim is to study the difference between 2 experiments 
        estatsdiff (bool):set to True if the aim is to read and study ensemble stat files (precomputed)
        mbdiff (bool): set to True if the aim is to study the difference between two members of a same ensemble experiment.
        origin2 (str): path for the second dataset/experiment if needed
        nexp2(str): name of second experiment (if needed)
        obs (bool): set to True if you aim to deal with obs dataset instead of model dataset
        lev (int): if variable to read is 3-D you should specify here which vertical level to consider

        Returns:
        self : a new instance of the imhov class.

        """
        self.varna = varna
        self.var = Fshortenvarna(self.varna)
        self.origin = origin1
        self.nexp = nexp
        self.y1 = y1
        self.y2 = y2
        self.fo = fo
        self.diff = diff
        self.estatsdiff=estatsdiff
        self.mbdiff = mbdiff
        self.obs = obs
        
        # if ESA-CCI obs dataset for salinity
        if nexp=="ESA":
            self.obs = True
            
        self.sif=sif
        self.lev=int(lev)
        self.levm=0
        self.dirigrid = dirigrid
        
        # if difference is to be computed
        if self.diff:
            self.origin2 = origin2
            self.nexp2 = nexp2 
            
            
        
    def process(self,typ='t'): 
        """Process the imhov instance.

        More precisely, it will apply the following operations:
        - self.loaddata() # load data
        - selfdetrend(fo=self.fo) # detrend datasets
        - self.tmcompute() # compute time-mean of datasets
        - self.stdcompute() # compute time std of datasets
        - self.loadgridinfo(type='t')  # load grid info for later plotting needs
        - if self.diff:
            # compute difference between datasets
            self.diffcompute()

        Parameters:
        self (instance of imhov class)

        Returns:
        Modified instance of the imhov class.

        """
        self.loaddata()
        self.detrend()
        self.tmcompute()
        self.stdcompute()
        self.loadgridinfo(type=typ)  # todo: improve type t later on    
        if self.diff:
                self.diffcompute()
        return self
    
    def processestats(self,typ='t'): 
        """Process the imhov instance in the case estats are loaded.

        More precisely, it will apply the following operations:
        - self.loaddata() # load data
        - self.loadgridinfo(type='t')  # load grid info for later plotting needs
        - self.estatsdiffcompute # Compute the difference between  the ensemble stats of two experiments

        Parameters:
        self (instance of imhov class)

        Returns:
        Modified instance of the imhov class.

        """
        self.loaddata()
        self.loadgridinfo(type=typ)  # todo: improve type t later on  
        self.estatsdiffcompute()
        self.tmcompute()
        self.diffTMcompute()
        
        return self

            
    def loaddata(self):
        """Load data of the imhov instance.

        Parameters:
        self (instance of imhov class)

        Returns:
        Modified instance of the imhov class.

        """
        # if obs dataset to load
        if self.obs:
            d1=str(self.y1)+'-01-01T00:00:00.000000000'
            d2=str(self.y2)+'-12-31T00:00:00.000000000'
            self.data = xr.open_mfdataset(self.origin,decode_times=True)[self.varna].sel(time=slice(d1,d2))
        else:
            
            self.timecod = xr.open_mfdataset(self.origin,decode_times=False)[self.varna]
            # if difference to be compûted, 2 datasets have to be loaded (self.data1 and self.data2)
            if self.diff:
                if self.lev==-1:
                    self.data1 = xr.open_mfdataset(self.origin,decode_times=True)[self.varna].sel(time_counter=slice(self.y1,self.y2))
                    self.data2 = xr.open_mfdataset(self.origin2,decode_times=True)[self.varna].sel(time_counter=slice(self.y1,self.y2))
                else :
                    self.data1 = xr.open_mfdataset(self.origin,decode_times=True)[self.varna].sel(time_counter=slice(self.y1,self.y2)).isel(deptht=self.lev)
                    self.data2 = xr.open_mfdataset(self.origin2,decode_times=True)[self.varna].sel(time_counter=slice(self.y1,self.y2)).isel(deptht=self.lev)
                    self.levm=self.data1.deptht.values
                    self.data1 = self.data1.squeeze()
                    self.data2 = self.data2.squeeze()
                self.data = self.data1
            else:
            #if only 1 dataset has to be loaded (self.data)
                # variable to load is 2d
                if self.lev==-1:
                    self.data = xr.open_mfdataset(self.origin,decode_times=True)[self.varna].sel(time_counter=slice(self.y1,self.y2))
                else:
                # variable to load is 3d and a single vertical level will be read
                    self.data = xr.open_mfdataset(self.origin,decode_times=True)[self.varna].sel(time_counter=slice(self.y1,self.y2)).isel(deptht=self.lev)
                    self.levm=self.data.deptht.values
                    self.data = self.data.squeeze()

            # if sea ice extent to be loaded.
            if self.sif=='no':
                pass
            else:
                ice_1m = xr.open_mfdataset(self.sif,decode_times=True)['simsk15'].sel(time_counter=slice(self.y1,self.y2))
                self.siextentTM = ice_1m.mean(dim='time_counter')

        return self



    def detrend_old(self,fo='1y'):  
        """Detrend data of the imhov instance.
      --- REPLACED BY NEW detrend ---

        Parameters:
        - self (instance of imhov class)
        - fo (str): frequency of data to detrend ('1m' or '1y') 

        Returns:
        Modified instance of the imhov class. self.data_dt is added as well as self.data2_dt if diff set to True.

        """
        # compute linear trend
        self.atr,xrorigins,years = Ftrpolyfit(self.data,fo,dat=self.obs)
        # retrieve trend  timeseries
        if self.obs:
            ts_tr = Ftrseries(self.atr,xrorigins,years,self.data.time,dat=self.obs)
        else:
            ts_tr = Ftrseries(self.atr,xrorigins,years,self.data.time_counter)
        # retrieve detrended timeseries
        self.data_dt = self.data - ts_tr
        
        if self.diff:
            # compute linear trend
            self.atr2,xrorigins2,years = Ftrpolyfit(self.data2,fo,dat=self.obs)
            # retrieve trend  timeseries
            ts_tr2 = Ftrseries(self.atr2,xrorigins2,years,self.data2.time_counter)
            # retrieve detrended timeseries
            self.data2_dt = self.data2 - ts_tr2
   
        return self

    def detrend(self):  
        """Detrend data of the imhov instance.

        Parameters:
        - self (instance of imhov class)

        Returns:
        Modified instance of the imhov class. self.data_dt is added as well as self.data2_dt if diff set to True.

        """
        
        # compute linear trend
        ts_tr,self.atr,xrorigins = Ftrseries2(self.data,obs=self.obs)
                
        # retrieve detrended timeseries
        self.data_dt = self.data - ts_tr
        
        if self.diff:
            # compute linear trend
            ts_tr2,self.atr2,xrorigins2 = Ftrseries2(self.data2,obs=self.obs)

            # retrieve detrended timeseries
            self.data2_dt = self.data2 - ts_tr2
   
        return self
    
    
    def stdcompute(self):
        """Compute time-STD of detrended datasets.

        You need to apply self.detrend() first

        Parameters:
        - self (instance of imhov class)

        Returns:
        Modified instance of the imhov class. self.data_std_dt is added as well as self.std_dt2 if diff set to True.

        """
        if self.obs:
            self.std_dt = self.data_dt.std(dim='time')     
        else:
            self.std_dt = self.data_dt.std(dim='time_counter') 
            if self.diff:
                self.std_dt2 = self.data2_dt.std(dim='time_counter')
        return self
    
    
    def diffcompute(self):
        """Compute the difference between the time-STD, time mean and trend of the 2 datasets.

        diff has to be set to True, otherwise only 1 dataset has been loadded.
        You need to apply self.detrend() and  self.stdcompute() first.

        Parameters:
        - self (instance of imhov class)

        Returns:
        Modified instance of the imhov class. self.data_std_dt is added as well as self.std_dt2 if diff set to True.
        """
        if self.diff:
            self.std_diff = self.std_dt - self.std_dt2
            self.tr_diff = self.atr - self.atr2
            self.tm_diff = self.tm - self.tm2
        else:
            sys.exit('error: data2 not loaded')
        return self
 
    def diffTMcompute(self):
        """Compute the difference between the time-mean of the 2 datasets.

        diff has to be set to True, otherwise only 1 dataset has been loadded.
        

        Parameters:
        - self (instance of imhov class)

        Returns:
        Modified instance of the imhov class. self.tm_diff is added
        """
        if self.diff:
            self.tm_diff = self.tm - self.tm2
        else:
            sys.exit('error: data2 not loaded')
        return self
    
    
    def estatsdiffcompute(self,rmgm=False):
        """ Compute the difference between  the ensemble stats of two experiments 
    
        The difference is computed of eman or estd at each time step. 
        A test checks first if the time series is the same in each experiment.
        Option rmgm set to True will remove the global mean of each field at each tme step before computing the diff.
        Need to apply self.loadgridinfo first as e1 and e2 (horizontal metric might be needed)
            
        Parameters:
        - self (instance of imhov class)
        - rmgm (bool): if you need to remove Global mean prior to computing the mean

        Returns:
        Modified instance of the imhov class. Adds self.data_diff and data_diff_rmGM set to true if global mean was removerd.
        
        """
        if self.estatsdiff:
            if (Ftesttime(self.data1 ,self.data2,tdn='time_counter' )):
            # test if time series are consistent between data1 and data2
                if rmgm:
                # test if need to remove Global Mean first
                    gm1 = FcomputeGM(self.data1,self.e1,self.e2)
                    gm2 = FcomputeGM(self.data2,self.e1,self.e2)
                    self.data_diff = FrmGM(self.data1,gm1) - FrmGM(self.data2,gm2)
                    self.data_diff_rmGM = rmgm
                else:
                    self.data_diff = self.data1 - self.data2
                    self.data_diff_rmGM = rmgm
            else:
                print('ERROR: time is different in data1 and data2')
                pass    
        else:
            print('ERROR: data2 not loaded')
            pass
        return self
    
    def mbdiffcompute(self):
        """ Compute the difference between  2 members of an ensemble
    
        The difference is computed of the 2 members  at each time step. 
        A test checks first if the time series is the same in each experiment.
    
            
        Parameters:
        - self (instance of imhov class)

        Returns:
        Modified instance of the imhov class. Adds self.data_diff
        """
        
        if self.mbdiff:
            if (Ftesttime(self.data1 ,self.data2,tdn='time_counter' )):
                self.data_diff = self.data1 - self.data2
            else:
                print('ERROR: time is different in data1 and data2')
                pass    
        else:
            print('ERROR: data2 not loaded')
            pass
        return self
    
    def tmcompute(self):
        """ Compute the time mean of the datasets
            
        Parameters:
        - self (instance of imhov class)

        Returns:
        Modified instance of the imhov class. Adds self.tm and self.tm2
        """
        if self.obs:
            self.tm = self.data.mean(dim='time') 
        else:
            self.tm = self.data.mean(dim='time_counter') 
            if self.diff:
                    self.tm2 = self.data2.mean(dim='time_counter')    
        return self  
                
    
    def loadgridinfo(self,type='t'):
        """ Load grid info
            
        Parameters:
        - self (instance of imhov class)
        - type (str): grid type ('t','v','u','f')

        Returns:
        Modified instance of the imhov class. Adds self.nav_lon self.nav_lat self.mask self.e1 self.e2
        """
        if self.obs:
            d1=str(self.y1)+'-01-01T00:00:00.000000000'
            d2=str(self.y2)+'-12-31T00:00:00.000000000'
            self.nav_lon = xr.open_mfdataset(self.origin,concat_dim='time',decode_times=True)['lon']
            self.nav_lat = xr.open_mfdataset(self.origin,concat_dim='time',decode_times=True)['lat']
            self.mask = xr.open_mfdataset(self.origin,concat_dim='time',decode_times=True)[self.varna].sel(time=slice(d1,d2)).isnull()
        else:
            diri=self.dirigrid
            self.mask = xr.open_dataset(diri+'mesh_hgr.nc')[type+'mask'][0,0,:,:]
            self.e1    = xr.open_dataset(diri+'mesh_hgr.nc')['e1'+type][0,:,:]
            self.e2    = xr.open_dataset(diri+'mesh_hgr.nc')['e2'+type][0,:,:]
            self.nav_lon = xr.open_dataset(diri+'mesh_hgr.nc')['nav_lon']
            self.nav_lat = xr.open_dataset(diri+'mesh_hgr.nc')['nav_lat']    
        return self

    def convertRNF(self):
        """ Convert RNF from NEMO unit to m3/s
            
        Parameters:
        - self (instance of imhov class)
        
        Returns:
        Modified instance of the imhov class. Adds self.rnfconv
        """
        if self.varna=='sornf':
            data_c = self.data
            data_c = data_c *self.e1*self.e2 
            data_c = data_c/1e3
            data_c.attrs['unit'] = "m3/s" 
            self.data = data_c
            self.rnfconv=True
        return self
   

    def sortoutRNF(self,ty='tm',p=0.95):
        """ Sort out runoffs and get percentile value
            
        Parameters:
        - self (instance of imhov class)
        - ty (str): type of analysis (if 'tm' or 'std' will look for one percentile value (1 side of the pdf), 
        if 'tr' will look for 2 (2 sides of the pdfs))
        
        Returns:
        Modified instance of the imhov class. Return a dict (sortrnf).
        """

        nav_lat_stack = self.nav_lat.stack(z=("x", "y"))
        nav_lon_stack = self.nav_lon.stack(z=("x", "y"))
        
        
        if ty=='tr':
            xrtrends_stack= self.atr.stack(z=("x", "y"))
            xtrend_stack_abs= abs(xrtrends_stack)
            
            #quantile pth high values
            perpos = xrtrends_stack.quantile(p,dim='z')
            perneg = xrtrends_stack.quantile(1-p,dim='z')
            perall = xtrend_stack_abs.quantile(p,dim='z')
            per=np.array([perneg,perpos,perall])
            print('Count occurences above pos. threshold '+str(p)+'(val='+str(perpos.values)+')')
            print(xrtrends_stack.where(xrtrends_stack>perpos,drop=True).count().values)
            print('Count occurences below neg. threshold '+str(1-p)+'(val='+str(perneg.values)+')')
            print(xrtrends_stack.where(xrtrends_stack<perneg,drop=True).count().values)
            out_stack = xrtrends_stack
            
        if ty=='tm':
            tm_stack = self.tm.stack(z=("x", "y"))
            rnf_nonzero = tm_stack.where(tm_stack>0.0,drop=True)
            per = rnf_nonzero.quantile(p,dim='z')
            print('Count occurences above threshold '+str(p)+'(val='+str(per.values)+')')
            print(tm_stack.where(tm_stack>per,drop=True).count().values)
            out_stack = tm_stack
            
        if ty=='std':
            std_stack = self.std_dt.stack(z=("x", "y"))
            rnf_nonzero = std_stack.where(std_stack>0.0,drop=True)
            per = rnf_nonzero.quantile(p,dim='z')
            print('Count occurences above threshold '+str(p)+'(val='+str(per.values)+')')
            print(std_stack.where(std_stack>per,drop=True).count().values)
            out_stack = std_stack
        
        sortrnf = {'ty':ty,'threshold':per,'nav_lat_stack':nav_lat_stack,'nav_lon_stack':nav_lon_stack,'out_stack':out_stack,'p':p}
        return sortrnf 
    

    def saveTR(self,trno=2,diro='./'):
        """ Save Linear Trend in nc file
            
        Parameters:
        - self (instance of imhov class)
        - trno: save only var 1 or also var2
        
        Returns:
        
        """
        fina = "TR_"+self.nexp+"_"+self.y1+"-"+self.y2+".nc"
        self.atr.name = 'TR'
        self.atr.attrs['exp'] =  self.nexp
        self.atr.to_netcdf(diro+fina,mode='w',format='NETCDF4')

        if trno==2:
            fina2 = "TR_"+self.nexp2+"_"+self.y1+"-"+self.y2+".nc"
            self.atr2.name = 'TR'
            self.atr2.attrs['exp'] =  self.nexp2
            self.atr2.to_netcdf(diro+fina2,mode='w',format='NETCDF4')

            
        return 

if __name__ == "__main__":
    main()
