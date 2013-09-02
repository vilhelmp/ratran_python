

### Imports
# python builtins
import os as _os
import sys as _sys
import subprocess as _subprocess

# python extra modules
import scipy as _scipy
from matplotlib import pyplot as _pl
from matplotlib.transforms import blended_transform_factory as _btf

# adapy internal imports
#~ from ...libs import cgsconst as _cgs
from .. import moldata as _moldata
#~ from ...tools import get_colors as _gc

#~ from ...views import set_rc
#~ _pl.rcParams['text.usetex'] = True
#~ _pl.rcParams['savefig.dpi'] = 300
#~ _pl.rcParams['figure.dpi'] = 150
_pl.rcParams['figure.facecolor'] = 'w'
_pl.rcParams['font.size'] = 10
#~ set_rc()

from .. import helpers

import ratout as ratout
import ratin as ratin
# read in main functions for
# reading and making Ratran runs
from ratout import Read
from ratin import Make
