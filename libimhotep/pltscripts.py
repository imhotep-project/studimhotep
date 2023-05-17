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

from libimhotep import libSLXtoolsIMHOTEP as li

print(f"Name: {__name__}")
print(f"Package: {__package__}")

print("This is a collection of  basic scripts to plot the different GLOBAL plots with my favourite settings.")


def __init__():
    print('Init')

def main():
    print('test')
    
def pltDIFF(DIFF,diro,pltshow=False):
    #set default grid parameters before plotting
    gridparam = li.Fpltgridparamdefaults(reg='GLO')
    print('- gridparam')
    print(gridparam)

    # set default color parameters before plotting
    pltcolparam = li.Fpltsetcolorsdefaults('DIFF')
    print('- pltcolparam')
    print(pltcolparam)

    # set default parameters for circles (runoffs)
    pltparam =  li.Fsetcirclesparamdefaults()
    print('- pltcirclesparam')
    print(pltparam)

    # plot
    outo = li.FpltGLO(DIFF,pltgridparam=gridparam,ty='DIFF',diro=diro,siplt=False,saveo=True,pltshow=pltshow,cbon=True,strunits="g/kg",pltcolparam=pltcolparam,selrnf='no',pltcircles=False,pltcirclesparam=pltparam)
    
    
def pltSTD(DIFF,diro,pltshow=False):
    #set default grid parameters before plotting
    gridparam = li.Fpltgridparamdefaults(reg='GLO')
    print('- gridparam')
    print(gridparam)

    # set default parameters for circles (runoffs)
    pltparam =  li.Fsetcirclesparamdefaults()
    print('- pltcirclesparam')
    print(pltparam)
    
    # set default color parameters before plotting
    pltcolparam = li.Fpltsetcolorsdefaults('STD')
    
    outo = li.FpltGLO(DIFF,pltgridparam=gridparam,ty='STD',siplt=False,saveo=True,pltshow=False,diro=diro,cbon=True,strunits="g/kg",pltcolparam=pltcolparam,selrnf='no',pltcircles=False,pltcirclesparam=pltparam)
    
def pltTR(DIFF,diro,selrnf='no',pltcircles=False,pltshow=False):
    #set default grid parameters before plotting
    gridparam = li.Fpltgridparamdefaults(reg='GLO')
    print('- gridparam')
    print(gridparam)
    
    # set default color parameters before plotting
    pltcolparam = li.Fpltsetcolorsdefaults('TR')
    pltcolparam

    # set default parameters for circles (runoffs)
    pltcolparam = li.Fpltsetcolorsdefaults('TR')
    pltcolparam['levbounds']=[-0.03,0.031,0.001]
    pltcolparam['cbincr']=5
    pltcolparam['cbformat'] = '%.3f'
    pltparam = li.Fsetcirclesparamdefaults()
    pltparam['linewidthstr']=1
    pltparam['mulmag']=600
    pltparam

    #selrnf=selrnftr
    outo = li.FpltGLO(DIFF,pltgridparam=gridparam,ty='TR',siplt=False,saveo=True,pltshow=pltshow,diro=diro,cbon=True,strunits="g/kg/yr",pltcolparam=pltcolparam,selrnf=selrnf,pltcircles=pltcircles,pltcirclesparam=pltparam)

    
def pltTRDIFF(DIFF,diro,selrnf='no',pltcircles=False,pltshow=False):
    #set default grid parameters before plotting
    gridparam = li.Fpltgridparamdefaults(reg='GLO')
    print('- gridparam')
    print(gridparam)
    
    # set default color parameters before plotting
    pltcolparam = li.Fpltsetcolorsdefaults('TRDIFF')
    pltcolparam['levbounds']=[-0.01,0.011,0.001]
    #pltcolparam['levbounds']=[-0.03,0.031,0.001]
    pltcolparam['cbincr']=5
    pltcolparam['cbformat'] = '%.3f'
    
    # set default parameters for circles (runoffs)
    pltparam = li.Fsetcirclesparamdefaults()
    pltparam['linewidthstr']=1
    pltparam['mulmag']=600

    pltparam
    outo = li.FpltGLO(DIFF,pltgridparam=gridparam,ty='TRDIFF',siplt=False,saveo=True,pltshow=pltshow,diro=diro,cbon=True,strunits="g/kg/yr",pltcolparam=pltcolparam,selrnf=selrnf,pltcircles=pltcircles,suffix="colorbar",pltcirclesparam=pltparam)

    
def pltTM(DIFF,diro,reg='GLO',lev='def',cbincr='def'):
    #set default grid parameters before plotting
    gridparam = li.Fpltgridparamdefaults(reg=reg)
    print('- gridparam')
    print(gridparam)

    # set default color parameters before plotting
    pltcolparam = li.Fpltsetcolorsdefaults('TM')
    if (lev=='def'):
        pass
    else:
        pltcolparam['levbounds']= lev
    if (cbincr=='def'):
        pass
    else:
        pltcolparam['cbincr']= cbincr

    # set default parameters for circles (runoffs)
    pltparam = li.Fsetcirclesparamdefaults()
    pltparam['linewidths']=0.6
    pltparam['edgecolors']='#FA5858'

    outo = li.FpltGLO(DIFF,pltgridparam=gridparam,ty='TM',siplt=False,saveo=True,cbon=True,pltshow=pltshow,strunits="g/kg",pltcolparam=pltcolparam,selrnf='no',pltcircles=False,pltcirclesparam=pltparam)
    
    
def pltTMDIFF(DIFF,diro,pltshow=False):
    #set default grid parameters before plotting
    gridparam = li.Fpltgridparamdefaults(reg='GLO')
    print('- gridparam')
    print(gridparam)

    # set default color parameters before plotting
    pltcolparam = li.Fpltsetcolorsdefaults('TMDIFF')
    pltcolparam['levbounds']=[-0.5,0.51,0.01]
    pltcolparam['cbincr']=5

    # set default parameters for circles (runoffs)
    pltparam = li.Fsetcirclesparamdefaults()
    pltparam['linewidths']=0.3
    #pltparam['edgecolors']='#FA5858'

    outo =li.FpltGLO(DIFF,pltgridparam=gridparam,ty='TMDIFF',siplt=False,saveo=True,pltshow=pltshow,cbon=True,strunits="g/kg",pltcolparam=pltcolparam,selrnf='no',pltcircles=False,pltcirclesparam=pltparam)
    
    

if __name__ == "__main__":
    main()
