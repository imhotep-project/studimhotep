"""Script to plot the basic analyses from a given experiment (either 1 year or 1 month frequency data). 

It calls functions from  my library of tools in libimhotep/libSLXclassIMHOTEP.py
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

import warnings
warnings.filterwarnings('ignore')

# custom tools for plotting

from libimhotep import libSLXtoolsIMHOTEP as li
from libimhotep import pltscripts as pltGLO

import cmocean
plt.rcParams.update({'hatch.color': '#086A87'})




# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# NAMELIST
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Which plots?
PDIFF   = False  # diff of time-std of the 2 exp (after detrending)
PSTD    = False  # time-std of exp1 (imhov.data1)
PTR     = True  # linear trend of exp1 (imhov.data1)
PTRDIFF = True  # diff of linear trend of the 2 exp 
PTM     = False  # time-mean of exp1 (imhov.data1)
PTMDIFF = False  # diff of time-mean of the 2 exp

prefix = "eORCA025.L75-IMHOTEP"
nexp = "EAI"
nexpREF = "ES"

#load runoffs?
RNFL = True

varnarnf = 'sornf'
varnasss = 'sosaline'


# years to read data from:
y1='1980'
y2='2018'

fo="1y" # output frequency

# path for data (experiments)
diridat = li.Ffindinputdata(nexp,prefix=prefix,fo=fo)
diridref = li.Ffindinputdata(nexpREF,prefix=prefix,fo=fo)

# input directory on work for grid info
diri="/gpfswork/rech/cli/rcli002/eORCA025.L75/eORCA025.L75-I/"

# plot directory
diro="/gpfswork/rech/cli/regi915/PLT/dec2022/"+fo+"/"

# data output directory
dirdat="/gpfswork/rech/cli/regi915/DAT/"

# full file names
sssfiles = li.Ffindinputdfiles(nexp,diridat,fo,prefix=prefix,fity="gridTsurf")
sssfilesREF = li.Ffindinputdfiles(nexpREF,diridref,fo,prefix=prefix,fity="gridTsurf")

#  Runoffs file list
rnffiles = li.Ffindinputdfiles(nexp,diridat,fo,prefix=prefix,fity="flxT")
print(rnffiles)

# ICE
if nexp=='02':
    sifiles = diridat+"????/"+prefix+nexp+"*icemod.nc"
else:
    sifiles = diridat+"1m/????/"+prefix+"."+nexp+"*icemod.nc"


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# LOAD and PROCESS DATA
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
print("======")
print('--- LOAD and PROCESS DATA over time period:')
# SSS DIFF
DIFF = li.imhov(sssfiles, varnasss,nexp,fo,y1,y2,dirigrid=diri,diff=True,sif='no',origin2=sssfilesREF,nexp2=nexpREF)

DIFF.process()

print(DIFF.data.time_counter.isel(time_counter=0).values)
print(DIFF.data.time_counter.isel(time_counter=-1).values)
print(DIFF.data.time_counter.size)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# LOAD and PROCESS runoffs
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if RNFL:
    print("======")
    print('--- LOAD and PROCESS RUNOFFS:')

    RNF1 = li.imhov(rnffiles, varnarnf,nexp,'1y',y1,y2,diff=False,dirigrid=diri)
    RNF1.process()

    # sort out runoffs
    selrnf = li.FselectRNF(RNF1.sortoutRNF(ty='tm',p=0.95))
    selrnftr = li.FselectRNF(RNF1.sortoutRNF(ty='tr',p=0.975))
    selrnfstd = li.FselectRNF(RNF1.sortoutRNF(ty='std',p=0.95))
        
        
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PLOT DIFF of time-STD
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

if PDIFF:
    print("======")
    print('--- PLOT STD DIFF map over the period:')
    pltGLO.pltDIFF(DIFF,diro)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PLOT DIFF of time-STD
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

if PSTD:
    print("======")
    print('--- PLOT STD of first EXP over the period:')
    pltGLO.pltSTD(DIFF,diro)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PLOT Global trend map (TR)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


if PTR:
    print("======")
    print('--- PLOT TR of first EXP over the period:')
    pltGLO.pltTR(DIFF,diro,selrnf=selrnftr,pltcircles=True)


#------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PLOT Global trend DIFF map (TRDIFF)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


if PTRDIFF:
    print("======")
    print('--- PLOT TR DIFF over the period:')
    pltGLO.pltTRDIFF(DIFF,diro,selrnf=selrnftr,pltcircles=True)




#------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PLOT Global Mean (TM)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


if PTM:
    print("======")
    print('--- PLOT TM of first exp over the period:')
    pltGLO.pltTM(DIFF,diro)




#------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PLOT Global Mean DIFF (TMDIFF)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

if PTMDIFF:
    print("======")
    print('--- PLOT TM DIFF over the period:')
    pltGLO.pltTMDIFF(DIFF,diro)
