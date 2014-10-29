import os as _os
from .. import cgsconst as _cgs
from ..helpers import *
import scipy as _sp

def read_ratraninput(modelfile = "transphere.mdl"):
    """
    Read the input model to Ratran and plots all parameters 
    agains ra/r1 column.
    
        
    #Example model file, can be more lines before @
    
    # Ratran input file based on Transphere results
    rmax=1.757e+15
    ncell=51
    tcmb=2.728
    columns=id,ra,rb,nh,nm,tk,td,db,vr
    gas:dust=75.0
    @
       1 0.000e+00 3.501e+12 0.000e+00 0.000e+00 2.259e+02 2.259e+02 2.400e+00 0.000e+00
     ...
    """
    # imports...
    from scipy import array, where
    
    class Mdl: pass
    
    #with open(modelfile, 'r') as f:
    #    lines = f.read().split('\n')
    #table_start_indicator = where(array(lines) == '@')[0][0]
    ## get the header/preamble with the singel values
    #header = lines[:table_start_indicator]
    ## first get the comments in the header
    #Mdl.comments = [i for i in header if i[0] == "#"]
    ## second get the keyword - value pairs
    #values = [i for i in header if i[0] != "#"]
    #key_val = array([i.split('=') for i in values])
    ## add a checker function for input?
    ## put into class attributes
    #for key, val in zip(key_val[:,0], key_val[:,1]):
    #    if key in ['columns']:
    #        Mdl.__dict__[key] = val.split(',')
    #    elif key in ['ncell']:
    #        Mdl.__dict__[key.replace(':','2')] = int(val)
    #    else:
    #        Mdl.__dict__[key.replace(':','2')] = float(val)
    
    with open(modelfile, 'r') as f:
        line = f.readline()
        Mdl.comments = []
        while not line.startswith('@'):
            if line.startswith('#'):
                Mdl.comments.append(line)
                line = f.readline()
                pass
            else:
                keyval = line.strip('\n').split('=')
                try:
                    setattr(Mdl, keyval[0].replace(':','_'), float(keyval[1]))
                except(ValueError):
                    setattr(Mdl, keyval[0].replace(':','_'), keyval[1])
                line = f.readline()
        Mdl.columns = Mdl.columns.split(',')
        lines = f.readlines()
    
    # now get the whole data table
    #table = lines[table_start_indicator + 1:]
    table = lines
    Mdl.table = array([i.split() for i in table]).astype('float')
    Mdl.table = Mdl.table.transpose()
    # put each column in to its own attribute of the class
    for col,vals in zip(Mdl.columns, Mdl.table):
        Mdl.__dict__[col] = vals
    # convert to AUs
    Mdl.ra_au = Mdl.ra * 100 / _cgs.AU # input in m, so times 100 and devide by cm/AU
    Mdl.rb_au = Mdl.rb * 100 / _cgs.AU
    Mdl.r_au = (Mdl.ra_au + Mdl.rb_au) / 2.
    
    # convert to relative abundances
    Mdl.nh[0] = 1.0
    try:
        Mdl.ne[0] = 1.0
    except (AttributeError):
        Mdl.ne = 0
    Mdl.nm_rel = Mdl.nm / (Mdl.nh + Mdl.ne)
    Mdl.nh[0] = 0.0
    try:
        Mdl.ne[0] = 0.0
    except (AttributeError, TypeError):
        pass
    Mdl.nm_rel[0] = 0.0
    return Mdl

def plot_ratraninput(directory = '', modelfile = "transphere.mdl"):
    import matplotlib.pyplot as pl
    from scipy import arange, array
    pl.ion()
    #~ from adavis import set_rc
    from matplotlib import rc
    rc('axes', linewidth=1)
    rc('patch', linewidth=1.5)
    rc('lines', linewidth=1.5, markeredgewidth=1)
    # read in the model file
    Ratran_mdl = read_ratraninput(_os.path.join(directory, modelfile))
    
    # -2 because we dont want the first two cols, id and ra
    N = len(Ratran_mdl.columns) - 2
    # now how many subplots do we need?
    if N % 3:
        np = N / 3 + 1
    else:
        np = N / 3
    pl.close()
    fig = pl.figure(1)#num = 1, figsize = ())
    # create dictionary to send/return the axes back
    plots = dict()
    # loop over the columns and plot
    for (i, dat) in zip(arange(N), Ratran_mdl.table[2:]):
        ax = 'ax{0}'.format(i)
        ylbl = Ratran_mdl.columns[i + 2]

        pln =  i +1
        if pln >1:
            plots[ax] = fig.add_subplot(np, 3, pln, sharex=plots['ax0'])
        else:
            plots[ax] = fig.add_subplot(np, 3, pln)
        x = (Ratran_mdl.ra[:] + Ratran_mdl.rb[:])/2*100/_cgs.AU
      
        if ylbl in ['db', 'vr']:
            plots[ax].semilogx(x, dat, '.')
            plots[ax].semilogx(x, dat, '-')
        elif ylbl in ['rb']:
            plots[ax].loglog(x, dat*100/_cgs.AU, '.')
            plots[ax].loglog(x, dat*100/_cgs.AU, '-')
        else:
            plots[ax].loglog(x, dat, '.')
            plots[ax].loglog(x, dat, '-')
        # set the right ylabel
        if ylbl in ['nh', 'nm', 'ne']:
            ylbl += ' (cm-3)'
        elif ylbl in ['tk', 'td', 'te']:
            ylbl += ' (K)'
        elif ylbl in ['db', 'vr']:
            ylbl += ' (km s-1)'
        elif ylbl in ['rb']:
            ylbl += ' (AU)'
        plots[ax].set_ylabel(ylbl)
        plots[ax].grid()
    [plots['ax{0}'.format(i)].set_xlabel('ra (AU)') for i in arange(N-3, N)]
    #~ [plots['ax{0}'.format(i)].xaxis.set_visible(0) for i in arange(N-3)]
    #~ [plots['ax{0}'.format(i)].xaxis.set_ticklabels([]) for i in arange(N-3)]


    #~ plots['ax{0}'.format(N-1)].set_xlabel('ra (AU)')
    fig.subplots_adjust(left=0.11, right= 0.97, bottom=0.11, top=0.96, wspace=0.43, hspace=0.15)
    return plots, fig

def create_molecular_abundance(temperature, 
                                abund_type = 'jump', 
                                Tjump = 100, 
                                Xs = [1E-4, 1E-9],
                                smooth = 0):
    """
    IMPORTANT:
    - assumes there is only one point where T > 100 K
    - and that this point has the highest cell number of all cells with
     T > 100 K
    """
    # calculates the molecular abundance from a predefined type and 
    # options
    # 'jump' abundance creates a jump from Xout to Xin 
    # where (temperature > Tjump)
    #
    # smooth : create jump that spans several 'smooth' numer of cells
    # -> should be an even number...
    #
    from scipy import where, ones, diff, sign, array, arange, exp
    [Xin, Xout] = Xs
    # jump abundance
    if 'jump' in abund_type:
        i100k = where(temperature >= Tjump)[0]
        mol_abundance = ones(len(temperature)) * Xout
        mol_abundance[i100k] = Xin
        if not smooth:
            print('discontinuity')
            # if it should be a discontinuity
            abund_param = False
        if smooth:
            # first make sure it is an even number of cellls to 
            # smooth over
            #~ if smooth % 2 != 0:
                #~ smooth += 1
                #~ print('Argument \'smooth\' needs to be an even number',
                    #~ 'adding one.')
            # use sigmoid function to connect the discontinuities        
            # assumes that the abundance is constant before 
            # and after jump
            ijump = max(i100k) - _sp.ceil(smooth/2.)
            Xdiff = _sp.log10(Xin-Xout) # log difference
            dX = Xdiff/float(smooth)    # increase in abundance for each cell
            cells_x = [ijump, ijump + smooth]
            #~ Xs_diff = diff(Xs)                  # outer - inner
            #~ d = direction = sign(Xs_diff)[0]    # direction of function
            #~ B = height = abs(Xs_diff)[0]        # height of the sigmoid
            #~ A = constant = max(Xs)              # curve start?-> min value
            #~ x0 = center = ijump+0.5
            #~ a = width = abs(diff(cells_x)[0])/10. # the width needs
                                #~ # to be very small, should investigate
                                #~ # and express more analytically
            #~ sigmoid = lambda x: A + d * B / (1 + exp(-(x - x0) / a))
            #~ y_interp = splineinterpolate1d(x, y, k=2)
            #~ x = arange(cells_x[0],cells_x[1],1);
            #~ y = 1.0 / (1 + exp(-x / 0.1))
            
            #~ for i  in arange(x[0], x[1], 1):
                #~ mol_abundance[i] = y_interp(i)
            #~ mol_abundance[:x[0]] = Xin
            j = 1
            for i in arange(cells_x[0], cells_x[1],1):
                #~ mol_abundance[i] = sigmoid(i)
                mol_abundance[i] = mol_abundance[cells_x[0]] - 10**(dX * j)
                #~ print mol_abundance[i]
                j += 1
            #~ abund_param = dict(constant = A, direction = d, height = B, center = x0, width = a)
    # add more abundance types later on
    else:
        raise Exception('No abundance type given.')
    # send back the calculated abundance
    
    return mol_abundance #, abund_param


def ratran_environment_check():
    # check RATRAN path
    # get it from the environment variable
    # "RATRANRATRAN", that needs to be set 
    try:
        ratran_path = _os.environ['RATRAN']
    except (KeyError):
        ratran_path = False
    if ratran_path:
        ratran_bin = _os.path.join(ratran_path, 'bin')
        global RUN_AMC
        RUN_AMC = _os.path.join(ratran_bin, 'amc')
        global RUN_SKY
        RUN_SKY = _os.path.join(ratran_bin, 'sky')
        return True
    else:
        print('Create an environment variable called RATRAN, '
        'otherwise the binaries cannot be found.')
        return False

class Make(object):
    def __init__(self, **kwargs):
        """
        r = 0.0, 
        rho_dust = 0.0, 
        temp = 0.0, 
        db = 0.0, 
        abund = 0.0, 
        vr = 0.0, 
        tdust = 0.0, 
        dustonly = 0, 
        modelfile = 'transphere.mdl', 
        dpc = 0.0, 
        imsize = 129, 
        pixel = 0.5, 
        trans = '220.0e9', 
        writeonly = 0, 
        skyonly = 0, 
        molfile='', 
        ncell = 50, 
        outputfile="ratranResult", 
        snr=20, 
        fixset=1e-6, 
        minpop=1e-4, 
        unit='Jypx', 
        opstates=0, 
        gas2dust=100, 
        nphot=1000, 
        temp_limit=10.0, 
        rho_limit=1E-4, 
        pxl_radius = 32, 
        los = 2
        
        
        Only supports one collision partner. If both para- and ortho-H2
        is tabulated in \'molfile\', please check the validity and if
        the gas-to-dust ratio has to be changed.

        Input :

        """

        from scipy import array
        # imports
        #~ import cgsconst as _cgs
        from scipy import zeros, array, logspace, linspace, log10
        from scipy import where, pi, exp, zeros, diff
        import scipy.interpolate
        import sys
        import os
        import numpy as np
        from time import time

        # Checking input parameters
        #
        # 'abundance',     'jump, 100, 1E-4, 1E-9', 'type, tjump, in, out', 'str'
        params = [
        'r',                            0,          'cm',    'array',   # Radial points
        'rin',                          0,          'AU',    'float',   # Inner radius
        'rout',                         0,          'AU',    'float',   # Outer radius (where to cut off)
        'rhodust',                      0,       'g/cm3',    'array',   # Dust density
        'molfile',      'ph2-18o-ph2.dat',            '',      'str',   # Name of moldata file
        'modelfile',     'transphere.mdl',            '',      'str',   # Name of model input file
        'outputfile',      'ratranresult',            '',      'str',   # Name of output file
        'db',                         0.0,        'km/s',    'float',   # 1/e half-width of line profile (Doppler b-parameter) (km s-1)
        'abundtype',               'jump',            '',      'str',   # Molecular abundance type ['jump', '?']
        'tjump',                    100.0,           'K',    'float',   # If 'jump' profile, at what T should the jump be
        'smoothjump',                   0,         'pxl',      'int',
        #~ 'collapse_radius',         1000.0,          'AU',    'float',   # radii where the collapse has proceeded a*t where a is the sound speed and t time since collapse start, or it is where the knee is in the shu model
        'xs',                [1E-4, 1E-9],    'relative',     'list',   # If 'jump' profile, what is [inner, outer] relative abundance
        'rrefs',                      [0],          'AU',     'list',   # at what radii the change in ncell occur
        #~ 'npsref',                     [0],            '',     'list',   # how many points to create in each rrefs interval
        'spacing',               ['log'],            '',     'list',    # what type spacing for the grid
        'vr',                         0.0,        'km/s',    'array',   # Radial velocity
        'velocitydirection',       'None',            '',      'str',   # Velocity direction if vr given 'infall', 'outflow'
        #'velocityfield',     lambda,              ,      'str',   # Velocity model 'shu_infall', 'db'
        'tdust',                      0.0,           'K',    'array',   # Dust temperature profile
        'outformat',               'fits',            '',      'str',   # output format, for SKY
        'dustonly',                 False,            '',     'bool',   # 
        'skyonly',                  False,            '',     'bool',   # 
        'writeonly',                False,            '',     'bool',   # Only write the input files (obsolete)
        'kappa',           'jena,thin,e6',            '',      'str',   # powerlaw,NU0,KAPPA0,BETA  OR  jena,(bare|thin|thick),(no|e5|e6|e7|e8)
        'temp',                         0,           'K',    'array',   # A power law emissivity model, kappa=KAPPA0*(nu/NU0)^BETA, where NU0 is in Hz, KAPPA0 in cm2/g_dust, and BETA is the power law index.
        'templim',                    8.0,           'K',    'float',   # Daniel uses 8 K
        'nh2lim',                     1E4,        'cm-3',    'float',   #
        'Tconstouter',              False,            '',     'bool',   # Constant temperature (=templim) where temp < templim
        'trans',                      '7',            '',      'str',   # Transition number(s) as string. If the input 'molfile' is defined, trans contains the transition numbers to be calculated. These are the numbers at the start of lines (10+NLEV) to (10+NLEV+NLIN) in the molecular data file.
        'dpc',                        0.0,          'pc',    'float',   # Distance to source
        'imsize',                     129,      'pixels',      'int',   # Number of pixels in the output image
        'pixel',                      0.5,    'asec/pxl',    'float',   # Pixel size in arcseconds
        'pxlradius',                   32,            '',      'int',   # Region (in numbers of pixels radius w.r.t. image center) over which to use multiple lines of sight (los)
        'los',                          2,            '',      'int',   # Number of lines of sight
        'opr',                         -1,            '',    'float',   # ortho-to-para ratio, or in RATRAN language : ne/nh ratio i.e. 2nd coll partner numerator
        'chans',                       50,            '',      'int',   # number of velocity channels
        'chwidth',                    0.2,            '',    'float',   # Channel width, in km/s
        'unit',                    'Jypx',            '',      'str',   # Output units ['Jypx', 'K', 'Wm2Hzsr']
        'snr',                       10.0,       'ratio',    'float',   # Requested minimum signal-to-noise
        'fixset',                    1E-6,            '',    'float',   # Convergence requirement for first stage
        'minpop',                    1E-4,            '',    'float',   # Minimum population to include in S/N calculation
        'nphot',                     1000,            '',      'int',   # Number of photons
        'opstates',                 False,            '',     'bool',   # ?
        'ncell',                     [20],            '',     'list',   # Number of grid cells in each section of rref
        'gas2dust',                 100.0,       'ratio',    'float',   # Gas to dust ratio to be used in the run
        'directory',       'ratr_model_1',    'dir name',      'str']   # Directory to work in
        # if loadfile is input, drop everythin and just load the file 
        # and, if it exists, the output
        
        if ratran_environment_check():
            pass
        elif not ratran_environment_check():
            _sys.exit('Path error, need to define variables RATRAN in '
                        'your .bashrc or similar.')
            
        print ('Model created with the following parameters:')
        
        # check input parameters
        class Input: pass
        param_zip = zip(params[0::4],params[1::4], params[2::4], params[3::4])
        for par, stdval, unit, typ in param_zip:
            if par in kwargs: # if input was given, save it
                value = kwargs[par]
                #~ printvalue = kwargs[par]
                #~ if par in ['r', 'rho_dust', 'temp']:
                    #~ printvalue = printvalue[0]
                #~ # print the value that was input for that given parameter
                #~ print '   {0:9} : {1} {2}'.format(par, str(printvalue).strip('\"'), unit)
                
            elif par not in kwargs: # if not input was given, use default
                if stdval != None:
                    # print the default value for the parameter
                    #~ print '   {0:9} : {1:15} {2:8} \t(default)'.format(par, stdval, unit)
                    value = stdval
                else:
                    raise Exception('Wrong parameter input/handling'
                                    'Check input and/or code...')
            # Check what type it should be and store it 
            #string
            if typ == 'str':
                #~ kwargs[par] = '{0}'.format(kwargs[par])
                Input.__dict__[par] = str(value)
                self.__dict__[par] = str(value)
            #integer (not sure about 'itypemw' perhaps boolean)
            elif typ == 'int':
                Input.__dict__[par] = int(value)
                self.__dict__[par] = int(value)
            #float
            elif typ == 'float':
                Input.__dict__[par] = float(value)
                self.__dict__[par] = float(value)
            #bool
            elif typ == 'bool':
                Input.__dict__[par] = bool(value)
                self.__dict__[par] = bool(value)
            # array
            elif typ == 'array':
                Input.__dict__[par] = array(value)
                self.__dict__[par] = array(value)
            elif typ == 'list':
                Input.__dict__[par] = list(value)
                self.__dict__[par] = list(value)
            elif typ == 'mixed':
                #~ try:
                Input.__dict__[par] = value
                self.__dict__[par] = value
                #~ except ValueError:
                    #~ try:
                        #~ Input.__dict__[par] = int(value)
                    #~ except ValueError:
                        #~ Input.__dict__[par] = str(value)
            else:
                Input.__dict__[par] = value
                self.__dict__[par] = value
        #
        Input.trans = Input.trans.replace(' ', '') # remove whitespaces!
        
        # input parameters contains all the input needed to 
        # create this class again
        self.Input = Input
        #
        # if pixel size is smaller (in AU) than rin, print warning
        #~ self.pixel
                
        # copy important parameters to the main class
        #~ self.r = self.Input.r
        #~ self.rhodust = self.Input.rhodust
        #~ self.temp = self.Input.temp
        #~ self.tdust = self.Input.tdust        
        #~ self.directory = self.Input.directory
        #~ self.temp = self.Input.temp
        #~ self.db = self.Input.db
        #~ self.vr = self.Input.vr
        
        #~ self.collapse_radius *= _cgs.AU # convert to cm
        # come up with a nicer solution than this...
        # check that none of the required arrays are not empty, or
        # to short
        # doesnt work for integers, which is default input...
        #
        # check input for arrays and such
        if len(self.r)<5:
            raise Exception('neeed array as r, rho_dust, abund and temp')
        if len(self.rhodust)<5:
            raise Exception('neeed array as r, rho_dust, abund and temp')
        if len(self.temp)<5:
            raise Exception('neeed array as r, rho_dust, abund and temp')
        
        if self.Tconstouter:
            i = where(self.temp < self.templim)
            self.temp[i] = self.templim
        
        # If no dust temperature given, assume it is in equilibrium 
        # with the gas
        if self.tdust == 0.0:
            self.tdust = self.temp
        
        if self.rin:
            self.rin *= _cgs.AU # convert rin to cm
            _index = max(where(self.r<self.rin)[0])
            #~ print _index
            self.r = self.r[_index:]
            #~ print self.r
            self.temp = self.temp[_index:]
            self.tdust = self.tdust[_index:]
            self.rhodust = self.rhodust[_index:]
            if type(self.db) == type(1.0):
                pass
            else:                
                self.db = self.db[_index:]
            self.vr = self.vr[_index:]
        elif not self.rin:
            self.rin = self.r[0]
        
        if self.rout:
            self.rout *= _cgs.AU # convert rout to cm
            _index = min(where(self.r>self.rout)[0])
            #~ print _index
            self.r = self.r[:_index]
            #~ print self.r
            self.temp = self.temp[:_index]
            self.tdust = self.tdust[:_index]
            self.rhodust = self.rhodust[:_index]
            
            if type(self.db) == type(1.0):
                pass
            else:                
                self.db = self.db[:_index]
            self.vr = self.vr[:_index]
        elif not self.rout:
            self.rout = self.r[-1]
    
        # calculate the radial dependence of the molecular
        # abundance depends on what type of abundance type is choosen
        #~ self.abund, self.abund_param =  create_molecular_abundance(self.temp, 
        self.abund =  create_molecular_abundance(self.temp, 
                                abund_type = self.Input.abundtype, 
                                Tjump = self.Input.tjump, 
                                Xs = self.Input.xs,
                                smooth = self.Input.smoothjump)
        #~ return None
        #~ self.abund = self.Input.abund
        #
        # CHECK if directory exists
        input_dir_path = os.path.join(os.getcwd(), self.directory)
        # if the directory exists
        if not make_dirs(input_dir_path): # will raise error if things go south (i.e., permissions not correct [I think...])
            print('Directory exists, continuing.')
        
        save_ratran(self)
        
        # rewrite this part when changing to object oriented
        # now it is very hack-ish and non pythonic
        ################################################################
        ################################################################
        # MOLECULAR H2 NUMBER DENSITY
        # 
        # TODO : why time gas2dust here? <- changed to times 100
        # what is correct, 100 for the standard assumption
        # rhodust is g/cm3
        # from dust density to number density (cm-3)
        #   cm-3 =  g/cm3 * 100 / (muH2 * g)
        # 100 is gas:dust, but the input gas2dust is only for the 
        # run itself, i.e. only related to the molecules
        self.nh2 = self.rhodust * 100  / (_cgs.MUH2 * _cgs.MP)
        # nh2 is number density of H2
        # rhodust * 100 = rhogas (g/cm3) (all gas H2+He+Metals)
        # MUH2 * MP =  molecular mass (H2+He+Metals) in g
        # so
        # rhogas / molecular mass = number density of H2 in cm-3
        
        ################################################################
        ################################################################
        # Velocity grid
        #~ 
        #~ self.vr = -1 * sqrt(2 * cgs.GG * self.mstar / (ratr.r)) * 1E-5 # to go over to km/s
        #~ self.vr[(self.r / cgs.AU) > self.rref] = 0.0  # r_inf in Crimier+2010
        #~ if self.collapse_radius: # if we have a maximum collapse radius
            #~ self.vr_int[(self.rr > self.collapse_radius)] = 0.0
        # so that we can run log10(self.vr) (these values are rounded off to 0.0 before writing to input file)
        self.vr[self.vr == 0.0] = 1E-20
        #~ vr_negative = where(self.vr < 0.0)[0]
        #~ self.vr[vr_negative] *= -1
        ################################################################
        ################################################################
        # ENVELOPE CUT OFF LIMIT
        #Find the envelope cut off for T and n
        #
        # if T goes below tepmlim (10K def) somewhere in the model
        try:
           ind_T = where(self.temp < self.templim)[0].min()
        #~ ind = where(r<(1.2E4*cgs.AU))[0].max()
        except (ValueError):
            ind_T = False
        # if n goes below rholim (1E4 def /cm3) somewhere in the model
        try:
            ind_n = where((self.nh2) < self.nh2lim)[0].min()
        except (ValueError):
            ind_n = False
        # T or n strongest constraints on radius
        # Neither T nor n constrain the radius
        # thus just use the last element
        if ind_n == False and ind_T == False:
            self.r_constraint = None
            ind = len(self.r)-1
        # Both constraint, which comes first
        elif ind_n != False and ind_T != False:
            # ind_n comes first
            ind = min((ind_n, int_T))
            # what if both have the same...
            # it will pick T, ok
            self.r_constraint = ['n', 'T'][ind_n < ind_T]
        elif ind_n != False:
            ind = ind_n
            self.r_constraint = 'n'
        elif ind_T != False:
            ind = ind_T
            self.r_constraint = 'T'
        
        # get values at cut off
        self.r_10k = self.r[ind]
        self.rhodust_10k = self.rhodust[ind]
        self.nh2_10k = self.rhodust_10k * 100 / _cgs.MUH2 / _cgs.MP
        self.temp_10k = self.temp[ind]
        #~ self.Y = self.Input.r.max() / self.Input.r.min()
        #~ print self.r_10k, self.r.min()
        self.Y = self.r_10k / self.r.min()
        self.ind = ind
        #
        # get values at r = 1000 AU
        ind_r1000 = where(self.r > 1000 * _cgs.AU)[0].min()
        self.rhodust_r1000 = self.rhodust[ind_r1000]
        self.nh2_r1000 =  self.nh2[ind_r1000]
        self.temp_r1000 = self.temp[ind_r1000]
        #
        # cut off where T<templim OR nh2<nh2lim (8-10 K and 1E4 cm-3)
        # first we have to remove all cells where T<templim K
        # RATRAN does not work well with them
        # after this you use the self.parameter.
        # TODO : perhaps not the best tactics..?
        # even if we set Tconstouter, it could still cut off due
        # to the density, which is good (?).
        self.r = self.r[:ind]           
        self.rhodust = self.rhodust[:ind]
        self.temp = self.temp[:ind]
        self.tdust = self.tdust[:ind]
        self.abund = self.abund[:ind]
        self.nh2 = self.nh2[:ind]
        self.vr = self.vr[:ind]
        ################################################################
        ################################################################
        # Refinement, for the refinement, easiest way is to redefine rx!
        # isn't it weird to first create a grid in transphere, and then 
        # another here?
        # need to be able to create a refinement grid, so perhaps just 
        # a handfull of cells outside of Tjump, and alot inside of it
        #
        # TODO new grid making method
        # TODO  refactor code!
        """
        What needs to be done here:
        - refinement around the Tjump, if smoothjump is True
        - refinement inside of 100 K
        -> several regions with different cell-density
        """
        # 'rrefs',                      [0],          'AU',     'list',   # what intervals to boost the number of points
        # 'npsref',                     [0],            '',     'list',   # how many points to create in each rrefs interval
        # 'refspace',               ['log'],            '',     'list',   # what type spacing for the reference grid
        # 'ncell',                       20,            '',      'int',   # Number of grid cells
        #
        # how to input
        # if ncell = [20], then its just like below, the whole range of radii
        # if ncell =[10,20], then 'rrefs' has to be input with one value
        # so between rin and rrefs[0] you get ncell[0] cells.
        # so having
        # 'ncell' = [10,20]
        # 'rrefs' = [50]
        # as input would mean a 10 point grid from 'rin' to 50 AU
        # and 20 point grid from 50 AU to 'rout'
        #
        # ONE grid
        # if its 0 = only one grid
        if not self.Input.rrefs[0]:
            self.rx = logspace( log10( self.r[0] ),
                        log10( self.r[-1] ),
                        num = self.ncell[0] + 1,
                        endpoint = True
                        )
        # SEVERAL grids
        # if its not 0, then we have n_grids > 1
        elif self.Input.rrefs[0]:
            #~ from scipy import linspace
            self.rrefs = [i * _cgs.AU for i in self.Input.rrefs] # AU to cm
            # check if spacing is long enough for the number spaces
            # if rrefs = [10], we need two spacing, e.g. spacing=['log','log']
            if len(self.Input.spacing)-1 != len(self.Input.rrefs):
                print('Warning, to few refspace supplied, '
                        'not as many as rrefs, assuming the first/default '
                        'is the same for all.')
                self.spacing = [self.Input.spacing[0] for i in range(len(self.Input.rrefs)+1)]

            # pick out the start and stop intervals
            r_grids = [] # list of radius grids to concatenate later

            # if rrefs is [10,30,100]
            # create array that is [10,10,30,30,100,100]
            self.rrefs = _sp.vstack([self.rrefs, self.rrefs])
            self.rrefs = self.rrefs.transpose().flatten()
            # then add rin and rout
            # so that it would be [rin, 10, 10, 30, 30, 100, 100, rout]
            self.rrefs = _sp.concatenate(([[self.rin], self.rrefs, [self.rout]]))
            # create start and stop values
            starts, stops = self.rrefs[0::2], self.rrefs[1::2]
            # start the construction of the grids
            print (stylify('Refinement with', fg='r'))
            for start, stop, npoints, spacing in zip(starts, stops, self.ncell, self.spacing):
                # if its the last point, include the endpoint
                endpoint_bool =  (stop == stops[-1])
                print('{0} cells : {1} - {2:.1} AU'.format(npoints, start/_cgs.AU, stop/_cgs.AU))
                # is the grid linear or log spaced?
                
                if spacing.lower() in ['lin', 'linear', 'linspace']:
                    # create that part of the grid and change rx accordingly
                    r_grids.extend(linspace(start, stop, 
                                        num = npoints, 
                                        endpoint = endpoint_bool)
                                )
                elif spacing.lower() in ['log', 'logarithm', 'logarithmic', 'logspace']:
                    # create that part of the grid and change rx accordingly
                    r_grids.extend(logspace(log10(start), log10(stop), 
                                        num = npoints, 
                                        endpoint = endpoint_bool)
                                    )
            self.rx = array(r_grids)
            #~ self.rx = self.rx
            #~ # units should be in cm
            #~ self.rrefs_cm = [i*_cgs.AU for i in self.Input.rrefs]
        else:
            print('no refinement!')
        #~ return None
        self.rx = np.insert(self.rx, 0, 0)
        self.r1 = self.rx[0:-1]
        self.r2 = self.rx[1:]
        
        #~ from scipy import dstack
        #~ self.rr = dstack((self.r1, self.r2)).ravel()
        
        #~ r1=np.insert(r[0:-1],0,0)
        #~ r2=np.array(r)
        self.rr = zeros(len(self.rx)-1, float)
        
        ################################################################
        ################################################################
        # INTERPOLATION of values
        # grid point distances allways have to be logarithmically spaces
        # does linear even make sense?
        # TODO : rr needs to be the the averaged radius in that cell!!
        # i.e. center of mass!!!
        # create rr array which is just the mean radius of each cell
        # assumes all spacings are logarithmic!!!
        self.rr[1:] = 10**( (log10(self.r1[1:]) + log10(self.r2[1:])) / 2.0 )
        
        
        ##### what if part needs to be interpolated linearly??
        # this whole implementation is more complex than what it needs to be!
        # need to rewrite this, so it is more my own code...
        
        # and the first cell has the same value as the first,
        # so the interpolated values are just copied and identical 
        # in cell 0 and 1 (except nh2int)
        self.rr[0] = self.rr[1]
        # Interpolate the values to 'ncell' cells
        self.nh2f = scipy.interpolate.interp1d(log10(self.r), log10(self.nh2))
        self.tkf = scipy.interpolate.interp1d(log10(self.r), log10(self.temp))
        self.tdf = scipy.interpolate.interp1d(log10(self.r), log10(self.tdust))
        self.abund_f = scipy.interpolate.interp1d(log10(self.r), log10(self.abund))
        self.vr_f = scipy.interpolate.interp1d(log10(self.r), log10(self.vr))
        #
        # Convert logarithms to floats
        self.nh2int = 10**self.nh2f(log10(self.rr))
        self.tkint = 10**self.tkf(log10(self.rr))
        self.tdint = 10**self.tdf(log10(self.rr))
        self.abund_int = 10**self.abund_f(np.log10(self.rr))
        self.vr_int = 10**self.vr_f(np.log10(self.rr))
        # if it is infall, multiply by -1 
        # log10 does not work all that well
        if self.velocitydirection in ['infall']:
            self.vr_int *= -1
        # round off the array so that 1E20 is 0E15
        self.vr_int = array([round(i, 15) for i in self.vr_int])
        ############################
        #~ nhint_p = nhint*2/4.
        #~ nhint_p[0] = 0.0
        #~ nhint_o = nhint*2/4.
        #~ nhint_o[0] = 0.0
        #~ teint = 10**tkf(log10(rr))
        ############################
        # nh2 needs to start at 0
        self.nh2int[0] = 0.0
        
        self.teint = self.tkint
        ################################################################
        ################################################################
        # O/P ratio of H2
        # define teint, ortho, para
        if self.Input.opr >= 0:
            self.opr = self.Input.opr
            print ('opr = {0}'.format(self.opr))
            para = 1. / (1. + self.opr)
            ortho = self.opr * para
            self.para = _sp.ones_like(self.teint) * para
            self.ortho = _sp.ones_like(self.teint) * ortho
        elif self.Input.opr == -1:
            print (stylify('Temperature dependent ortho/para', fg='r'))
            self.opr = 9.0 * exp(-170.6 / self.teint)
            self.opr = np.clip(self.opr, 1.0E-3, 3.0)
            self.para = 1.0 / (1 + self.opr)
            self.ortho = 1 - self.para
        else:
            print('Could not understand the \'op\' parameter, please correct.')
        ################################################################
        ################################################################
        # mass 
        # mass of it all
        #~ vol=[]
        #~ mass=[]
        # V = 4*pi*r**3/3
        # r in cm (?)
        #~ V = 4 * pi * (self.r2**3 - self.r1**3) / 3     # cm3
        #~ V = 4 * pi * ((self.r2**3  - self.r1**3 )) / 3     # cm3
        #~ self.M = V * self.nh2int * _cgs.MUH2 * _cgs.MP # g = cm3 * g/cm3
        # proper calculation of mass
        # TODO : the density is slightly wrong though, but shouldn't matter
        # too much for these simple models and lines.
        # need to correct this in the future though...
        rho = self.nh2int * _cgs.MUH2 * _cgs.MP
        dr = self.r2 - self.r1
        r = (self.r1 + self.r2)/2.
        self.M = 4 * _sp.pi * _sp.sum(r**2 * rho * dr)
        self.M /= _cgs.MSUN                            # Msun
        
        # to get the column density, integrate over radius r1 to r_10k
        #r_10k * 2 nh2_10k
    
        ################################################################
        #  print info. -> move to __str__ method
        print ('M_10K   : {0:<7.2f} Msun\n'
                'R_10K   : {1:<7.0f} AU\n'
                'nH2_10K : {2:<7.1e} cm-3\n'
                'Y       : {3:<7.0f}\n'
                'T       : {4:<7.1f} K\n'.format(self.M.sum(),
                                            self.r_10k/_cgs.AU,
                                            self.nh2_10k,
                                            self.Y,
                                            self.temp_10k))
        print 'Constraining the envelope : ', self.r_constraint
        print ('nH2_r1000   : {0:<7.1e} cm-3\n'
                'T_r1000     : {1:7.1f} K\n'.format(self.nh2_r1000,
                                             self.temp_r1000))
        
        print('printing input files')
        
        """
        id    : shell number
        ra,rb : inner & outer radius (m)
        za,zb : lower & upper height (m) (2D only)
        nh    : density (cm-3) of main collision partner (usually H2)
        nm    : density (cm-3) of molecule
        ne    : density (cm-3) of second collision partner (e.g. electrons)
        tk    : kinetic temperature (K) 
        td    : dust temperature (K)
        te    : electron/second coll. partner temperature (K)
        db    : 1/e half-width of line profile (Doppler b-parameter) (km s-1)
        vr    : radial velocity (km s-1)
        """
        with open(_os.path.join(self.directory, self.modelfile),'w') as f:
            f.write('# Ratran input file based on Transphere results'+'\n')
            if self.skyonly: 
                f.write('# ... intended for (SKY) continuum calculations only.'+'\n')
            f.write("rmax={0:.5E}\n".format( self.r2[-1] / 100 ))       # rmax in METERS (convert from cm i.e. / 100)
            f.write("ncell={0:}\n".format(len(self.r2)))
            f.write("tcmb=2.735\n")
            f.write("columns=id,ra,rb,nh,nm,ne,tk,td,te,db,vr\n")
            f.write("gas:dust={0}\n".format(self.gas2dust))
            if self.skyonly: 
                f.write("kappa={0}\n".format(self.kappa))
            f.write('@\n')
            # r1/r2 in meter (convert from cm)
            for ii in range(0, len(self.r1)):
                test = ("{0:4} "                         #  1 id : shell number
                        "{1:12.5E} "                     #  2 ra : inner radius (m)  
                        "{2:12.5E} "                     #  3 rb : outer radius (m)
                        "{3:12.5E} "                     #  4 nh : density (cm-3) of main coll. partner (usually H2)
                        "{4:12.5E} "                     #  5 nm : density (cm-3) of molecule
                        "{5:12.5E} "                     #  6 ne : density (cm-3) of second coll. partner (e.g. electrons)
                        "{6:12.5E} "                     #  7 tk : kinetic temperature (K) 
                        "{7:12.5E} "                     #  8 td : dust temperature (K)
                        "{8:12.5E} "                     #  9 te : second coll. partner temperature (K)
                        "{9:12.5E} "                     # 10 db : 1/e half-width of line profile (Doppler b-parameter) (km s-1)
                        "{10:12.5E}\n")                  # 11 vr : radial velocity (km s-1)
                # now print the whole shebang
                f.write(test.format(ii + 1,               #  1 id : shell number
                    self.r1[ii] / 100.0,                  #  2 ra : inner radius (m)  
                    self.r2[ii] / 100.0,                  #  3 rb : outer radius (m)
                    self.nh2int[ii] * self.para[ii],      #  4 nh : density (cm-3) of main coll. partner (usually p-H2)
                    self.nh2int[ii] * self.abund_int[ii], #  5 nm : density (cm-3) of molecule
                    self.nh2int[ii] * self.ortho[ii],     #  6 ne : density (cm-3) of second coll. partner (e.g, e^-, o-H2)
                    self.tkint[ii],                       #  7 tk : kinetic temperature (K) 
                    self.tdint[ii],                       #  8 td : dust temperature (K)
                    self.teint[ii],                       #  9 te : second coll. partner temperature (K)
                    self.db,                              # 10 db : 1/e half-width of line profile (Doppler b-parameter) (km s-1)
                    round(self.vr_int[ii], 15))           # 11 vr : radial velocity (km s-1)
                        )          
                                    
                                    
                                    
                #~ f.write("%4i %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E" % (ii+1, self.r1[ii]/100.0, self.r2[ii]/100.0, self.nh2int[ii], self.nh2int[ii]*self.abund_int[ii], self.tkint[ii], self.tdint[ii], self.db, self.vr[ii])+'\n')
            
        if not self.skyonly:
            if self.molfile == '':
                sys.exit('Error: for AMC calculations the molecular datafile (molfile) needs to be set.')
            with open(_os.path.join(self.directory, "amc.inp"),'w') as f:
                f.write("source={0}\n".format(self.modelfile))
                f.write("outfile=populations.pop\n")
                f.write("molfile={0}\n".format(self.molfile))
                f.write("snr={0}\n".format(self.snr))
                f.write("velo=grid\n")
                # velo=grid if velocity vector is given in input model?
                f.write("nphot={0}\n".format(self.nphot))
                f.write("kappa={0}\n".format(self.kappa))
                f.write("minpop={0:3.2E}\n".format(self.minpop))
                #~ f.write("seed=1971\n")
                f.write("fixset={0:3.2E}\n".format(self.fixset))
                f.write("trace=on\n")
                f.write("go\n")
                f.write("q\n")
                f.write("\n")


        with open(_os.path.join(self.directory, "sky.inp"),'w') as f:
            if self.skyonly:
                f.write("source={0}\n".format(self.modelfile))
            else:
                f.write("source=populations.pop\n")                     # just use the AMC output file (always set to populations.pop above)
            f.write("format={0}\n".format(self.outformat))
            f.write("outfile="+self.outputfile+"\n")
            f.write("trans={0}\n".format(self.trans))
            f.write("pix={0},{1:f},{2},{3}\n".format(self.imsize, self.pixel, self.pxlradius, self.los))
            if self.skyonly:
                f.write("chan=1,1.0\n")
            else:
                f.write("chan={0},{1:f}\n".format(self.chans, self.chwidth))
            f.write("distance={0}\n".format(self.dpc))
            f.write("units={0}\n".format(self.unit))
            f.write("go\n")
            f.write("q\n")
            f.write("\n")
            
    def __str__(self):
        
        print ('M_10K   : {0:<7.2f} Msun\n'
        'R_10K   : {1:<7.0f} AU\n'
        'nH2_10K : {2:<7.1e} cm-3\n'
        'Y       : {3:<7.0f}\n'
        'T       : {4:<7.1f} K\n'.format(self.M.sum(),
                                    self.r_10k/_cgs.AU,
                                    self.nh2_10k,
                                    self.Y,
                                    self.temp_10k))
        print 'Constraining the envelope : ', self.r_constraint
        print ('nH2_r1000   : {0:<7.1e} cm-3\n'
                'T_r1000     : {1:7.1f} K\n'.format(self.nh2_r1000,
                                             self.temp_r1000))
        return 'Info. about the model'

    def save(self, filename='ratran_input.pick'):
        import pickle
        with open(os.path.join(), 'w') as f:
            pickle.dump(vars(self.Input), f)

    def run(self):
        # run ratran with the setup in the directory
        #
        # catch : ### WARNING
        # AMC: minimum S/N  |  converged  |     photons  |  increase to
        # AMC: -------------|-------------|--------------|-------------
        # and the output after this..
        # e.g.
        # AMC: 1.38711E+00  |      78.43% |       89000  |      116000
        # import modules, subprocess replaces os
        import subprocess
        from time import time, sleep
        import sys
        # 
        if not self.writeonly:
            if not self.skyonly:
                #~ print "Starting AMC calculation..."
                with ChangeDirectory(self.directory):
                    f = open('amc.log', 'w')
                    f.close()
                    t1 = time()
                    proc = subprocess.Popen([RUN_AMC, 'amc.inp'],
                                    stdout = subprocess.PIPE, 
                                    stderr = subprocess.STDOUT)
                    #~ sys.stdout.write('Iteration no : ')
                    #~ sys.stdout.flush() # flush output so we can write again
                    #~ amc_out = []
                    self.amc_output = []
                    step = 0
                    sys.stdout.write('AMC Fixset converged : ')
                    sys.stdout.flush()
                    while True:
                        # first : if process is done, break the loop
                        if proc.poll() != None: 
                            break
                        nextline = proc.stdout.readline()
                        #~ amc_out.append(nextline)
                        self.amc_output.append(nextline)
                        #~ sys.stdout.write(nextline+'\n')
                        #~ sys.stdout.flush()
                        if "% converged" in nextline:
                            conv_perc = nextline[-18:-10]
                            sys.stdout.write('{0:7}{1}'.format(conv_perc, '\b'*8)) # print out a point
                            sys.stdout.flush()
                        if "FIXSET convergence reached..." in nextline:
                            step += 1
                            sys.stdout.write('{0:7} {1}\n'.format('100.00%', 'FIXSET done, starting RANDOM'))
                            sys.stdout.flush()
                            sys.stdout.write('{0}'.format('Percent converged : '))
                            sys.stdout.flush()
                        if "% | " in nextline and step == 1:
                            conv_perc = nextline[24:31]
                            sys.stdout.write('{0:8}{1}'.format(conv_perc, '\b'*8)) # print out a point
                            sys.stdout.flush()
                        if "Warning" in nextline or "Fatal Error" in nextline:
                            sys.stdout.write('{0}'.format(nextline)) # print out a point
                            sys.stdout.flush()
                        open('amc.log', 'a').write('{0}'.format(nextline))
                            #~ # grab the iteration number
                            #~ iter_no = int(nextline[10:])
                            #~ sys.stdout.write('{0:3} \b\b\b\b'.format(iter_no)) # print out a point
                            #~ sys.stdout.flush()
                            
                        #~ if "Error" in nextline:
                            # if it is the error of the first iteration
                            # we grab the error of it
                            #~ if iter_no == 1:
                                #~ start_error = float(nextline[13:])
                                #~ first = 0 # now we are not at the first any longer
                            # if it is not the first error, calculate per cent done
                            #~ else:
                                #~ current_error = float(nextline[13:])
                                #~ diff1 = start_error - self.convcrit
                                #~ diff2 = start_error - current_error
                                #~ p_done = (diff2 / diff1 * 100)
                                #~ sys.stdout.write(' {0} : {1}\r'.format(iter_no, p_done)) # print out a point
                                #~ sys.stdout.flush()    # flush output so we can write again
                        #~ sleep(0.5)            # wait for 0.5 second
                        #~ sys.stdout.write('Downloading File FooFile.txt [%d%%]\r'%i)
                        #~ sys.stdout.flush()
                    print('\nAMC took {0:2.1f} seconds'.format((time()-t1)))
                    f.close()
                    #~ self.amc_output = ''.join(amc_out)
            # now run SKY
            with ChangeDirectory(self.directory):
                f = open('sky.log', 'w')
                f.close()
                t1 = time()
                proc = subprocess.Popen([RUN_SKY, 'sky.inp'],
                                    stdout = subprocess.PIPE, 
                                    stderr = subprocess.STDOUT)
                #~ sys.stdout.write('Iteration no : ')
                #~ sys.stdout.flush() # flush output so we can write again
                #~ sky_out = []
                self.sky_output = []
                while True:
                    # first : if process is done, break the loop
                    if proc.poll() != None: 
                        break
                    nextline = proc.stdout.readline()
                    #~ print(nextline)
                    #~ sky_out.append(nextline)
                    self.sky_output.append(nextline)
                    #~ sys.stdout.write(nextline+'\n'); sys.stdout.flush()
                    open('sky.log', 'a').write('{0}'.format(nextline))
                    
                    #~ if "Iteration" in nextline:
                        #~ # grab the iteration number
                        #~ iter_no = int(nextline[10:])
                        #~ sys.stdout.write('{0:3} \b\b\b\b'.format(iter_no)) # print out a point
                        #~ sys.stdout.flush()
                        
                    #~ if "Error" in nextline:
                        # if it is the error of the first iteration
                        # we grab the error of it
                        #~ if iter_no == 1:
                            #~ start_error = float(nextline[13:])
                            #~ first = 0 # now we are not at the first any longer
                        # if it is not the first error, calculate per cent done
                        #~ else:
                            #~ current_error = float(nextline[13:])
                            #~ diff1 = start_error - self.convcrit
                            #~ diff2 = start_error - current_error
                            #~ p_done = (diff2 / diff1 * 100)
                            #~ sys.stdout.write(' {0} : {1}\r'.format(iter_no, p_done)) # print out a point
                            #~ sys.stdout.flush()    # flush output so we can write again
                    #~ sleep(0.5)            # wait for 0.5 second
                    
                    #~ sys.stdout.write('Downloading File FooFile.txt [%d%%]\r'%i)
                    #~ ys.stdout.flush()
                print('\nSKY took {0:2.1f} seconds'.format((time()-t1)))
                f.close()
                #~ self.sky_output = ''.join(sky_out)
        #~ # now get the output as a string, so we can check stuff if
        #~ # we want to
        #~ self.transphere_output = ''.join(trans_out)
        #~ # read in the output-files, for checks, plots etc
        #~ self.Envstruct, self.Convhist = read_transphereoutput(self)




class MakeSky(object):
    
    def __init__():
        params = [
        'r',                            0,          'cm',    'array',   # Radial points
        'rout',                         0,          'AU',    'float',   # Outer radius (where to cut off)
        'rhodust',                      0,       'g/cm3',    'array',   # Dust density
        'modelfile',     'transphere.mdl',            '',      'str',   # Name of model input file
        'outputfile',      'ratranresult',            '',      'str',   # Name of output file
        'tdust',                      0.0,           'K',    'array',   # Dust temperature profile
        'outformat',               'fits',            '',      'str',   # output format, for SKY
        'kappa',           'jena,thin,e6',            '',      'str',   # powerlaw,NU0,KAPPA0,BETA  OR  jena,(bare|thin|thick),(no|e5|e6|e7|e8)
        'temp',                         0,           'K',    'array',   # A power law emissivity model, kappa=KAPPA0*(nu/NU0)^BETA, where NU0 is in Hz, KAPPA0 in cm2/g_dust, and BETA is the power law index.
        'templim',                    8.0,           'K',    'float',   # Daniel uses 8 K
        'nh2lim',                     1E4,        'cm-3',    'float',   #
        'Tconstouter',              False,            '',     'bool',   # Constant temperature (=templim) where temp < templim
        'frequency',              203.4E9,          'Hz',    'float',   # Transition number(s) as string. If the input 'molfile' is defined, trans contains the transition numbers to be calculated. These are the numbers at the start of lines (10+NLEV) to (10+NLEV+NLIN) in the molecular data file.
        'dpc',                        0.0,          'pc',    'float',   # Distance to source
        'imsize',                     129,      'pixels',      'int',   # Number of pixels in the output image
        'pixel',                      0.5,    'asec/pxl',    'float',   # Pixel size in arcseconds
        'pxlradius',                   32,            '',      'int',   # Region (in numbers of pixels radius w.r.t. image center) over which to use multiple lines of sight (los)
        'los',                          2,            '',      'int',   # Number of lines of sight
        'unit',                    'Jypx',            '',      'str',   # Output units ['Jypx', 'K', 'Wm2Hzsr']
        'snr',                       10.0,       'ratio',    'float',   # Requested minimum signal-to-noise
        'nphot',                     1000,            '',      'int',   # Number of photons
        'gas2dust',                 100.0,       'ratio',    'float',   # Gas to dust ratio to be used in the run
        'directory',              'sky_1',    'dir name',      'str']   # Directory to work in


        if ratran_environment_check():
            pass
        elif not ratran_environment_check():
            _sys.exit('Path error, need to define variables RATRAN in '
                        'your .bashrc or similar.')
      
        if Tconstouter:
            i = where(temp < templim)
            temp[i] = templim
        
        # If no dust temperature given, assume it is in equilibrium 
        # with the gas
        if tdust == 0.0:
            tdust = temp
        
        # calculate the radial dependence of the molecular
        # abundance depends on what type of abundance type is choosen
        #~ self.abund, self.abund_param =  create_molecular_abundance(self.temp, 
        self.abund =  create_molecular_abundance(self.temp, 
                                abund_type = self.Input.abundtype, 
                                Tjump = self.Input.tjump, 
                                Xs = self.Input.xs,
                                smooth = self.Input.smoothjump)
        #~ return None
        #~ self.abund = self.Input.abund
        #
        # CHECK if directory exists
        input_dir_path = os.path.join(os.getcwd(), self.directory)
        # if the directory exists
        if not make_dirs(input_dir_path): # will raise error if things go south (i.e., permissions not correct [I think...])
            print('Directory exists, continuing.')
        
        save_ratran(self)
        
        # rewrite this part when changing to object oriented
        # now it is very hack-ish and non pythonic
        ################################################################
        ################################################################
        # MOLECULAR H2 NUMBER DENSITY
        # 
        # TODO : why time gas2dust here? <- changed to times 100
        # what is correct, 100 for the standard assumption
        # rhodust is g/cm3
        # from dust density to number density (cm-3)
        #   cm-3 =  g/cm3 * 100 / (muH2 * g)
        # 100 is gas:dust, but the input gas2dust is only for the 
        # run itself, i.e. only related to the molecules
        self.nh2 = self.rhodust * 100  / (_cgs.MUH2 * _cgs.MP)
        # nh2 is number density of H2
        # rhodust * 100 = rhogas (g/cm3) (all gas H2+He+Metals)
        # MUH2 * MP =  molecular mass (H2+He+Metals) in g
        # so
        # rhogas / molecular mass = number density of H2 in cm-3
        
        ################################################################
        ################################################################
        # Velocity grid
        #~ 
        #~ self.vr = -1 * sqrt(2 * cgs.GG * self.mstar / (ratr.r)) * 1E-5 # to go over to km/s
        #~ self.vr[(self.r / cgs.AU) > self.rref] = 0.0  # r_inf in Crimier+2010
        #~ if self.collapse_radius: # if we have a maximum collapse radius
            #~ self.vr_int[(self.rr > self.collapse_radius)] = 0.0
        # so that we can run log10(self.vr) (these values are rounded off to 0.0 before writing to input file)
        self.vr[self.vr == 0.0] = 1E-20
        #~ vr_negative = where(self.vr < 0.0)[0]
        #~ self.vr[vr_negative] *= -1
        ################################################################
        ################################################################
        # ENVELOPE CUT OFF LIMIT
        #Find the envelope cut off for T and n
        #
        # if T goes below tepmlim (10K def) somewhere in the model
        try:
           ind_T = where(self.temp < self.templim)[0].min()
        #~ ind = where(r<(1.2E4*cgs.AU))[0].max()
        except (ValueError):
            ind_T = False
        # if n goes below rholim (1E4 def /cm3) somewhere in the model
        try:
            ind_n = where((self.nh2) < self.nh2lim)[0].min()
        except (ValueError):
            ind_n = False
        # T or n strongest constraints on radius
        # Neither T nor n constrain the radius
        # thus just use the last element
        if ind_n == False and ind_T == False:
            self.r_constraint = None
            ind = len(self.r)-1
        # Both constraint, which comes first
        elif ind_n != False and ind_T != False:
            # ind_n comes first
            ind = min((ind_n, int_T))
            # what if both have the same...
            # it will pick T, ok
            self.r_constraint = ['n', 'T'][ind_n < ind_T]
        elif ind_n != False:
            ind = ind_n
            self.r_constraint = 'n'
        elif ind_T != False:
            ind = ind_T
            self.r_constraint = 'T'
        
        # get values at cut off
        self.r_10k = self.r[ind]
        self.rhodust_10k = self.rhodust[ind]
        self.nh2_10k = self.rhodust_10k * 100 / _cgs.MUH2 / _cgs.MP
        self.temp_10k = self.temp[ind]
        #~ self.Y = self.Input.r.max() / self.Input.r.min()
        #~ print self.r_10k, self.r.min()
        self.Y = self.r_10k / self.r.min()
        self.ind = ind
        #
        # get values at r = 1000 AU
        ind_r1000 = where(self.r > 1000 * _cgs.AU)[0].min()
        self.rhodust_r1000 = self.rhodust[ind_r1000]
        self.nh2_r1000 =  self.nh2[ind_r1000]
        self.temp_r1000 = self.temp[ind_r1000]
        #
        # cut off where T<templim OR nh2<nh2lim (8-10 K and 1E4 cm-3)
        # first we have to remove all cells where T<templim K
        # RATRAN does not work well with them
        # after this you use the self.parameter.
        # TODO : perhaps not the best tactics..?
        # even if we set Tconstouter, it could still cut off due
        # to the density, which is good (?).
        self.r = self.r[:ind]           
        self.rhodust = self.rhodust[:ind]
        self.temp = self.temp[:ind]
        self.tdust = self.tdust[:ind]
        self.abund = self.abund[:ind]
        self.nh2 = self.nh2[:ind]
        self.vr = self.vr[:ind]
        ################################################################
        ################################################################
        # Refinement, for the refinement, easiest way is to redefine rx!
        # isn't it weird to first create a grid in transphere, and then 
        # another here?
        # need to be able to create a refinement grid, so perhaps just 
        # a handfull of cells outside of Tjump, and alot inside of it
        #
        # TODO new grid making method
        # TODO  refactor code!
        """
        What needs to be done here:
        - refinement around the Tjump, if smoothjump is True
        - refinement inside of 100 K
        -> several regions with different cell-density
        """
        # 'rrefs',                      [0],          'AU',     'list',   # what intervals to boost the number of points
        # 'npsref',                     [0],            '',     'list',   # how many points to create in each rrefs interval
        # 'refspace',               ['log'],            '',     'list',   # what type spacing for the reference grid
        # 'ncell',                       20,            '',      'int',   # Number of grid cells
        #
        # how to input
        # if ncell = [20], then its just like below, the whole range of radii
        # if ncell =[10,20], then 'rrefs' has to be input with one value
        # so between rin and rrefs[0] you get ncell[0] cells.
        # so having
        # 'ncell' = [10,20]
        # 'rrefs' = [50]
        # as input would mean a 10 point grid from 'rin' to 50 AU
        # and 20 point grid from 50 AU to 'rout'
        #
        # ONE grid
        # if its 0 = only one grid
        if not self.Input.rrefs[0]:
            self.rx = logspace( log10( self.r[0] ),
                        log10( self.r[-1] ),
                        num = self.ncell[0] + 1,
                        endpoint = True
                        )
        # SEVERAL grids
        # if its not 0, then we have n_grids > 1
        elif self.Input.rrefs[0]:
            #~ from scipy import linspace
            self.rrefs = [i * _cgs.AU for i in self.Input.rrefs] # AU to cm
            # check if spacing is long enough for the number spaces
            # if rrefs = [10], we need two spacing, e.g. spacing=['log','log']
            if len(self.Input.spacing)-1 != len(self.Input.rrefs):
                print('Warning, to few refspace supplied, '
                        'not as many as rrefs, assuming the first/default '
                        'is the same for all.')
                self.spacing = [self.Input.spacing[0] for i in range(len(self.Input.rrefs)+1)]

            # pick out the start and stop intervals
            r_grids = [] # list of radius grids to concatenate later

            # if rrefs is [10,30,100]
            # create array that is [10,10,30,30,100,100]
            self.rrefs = _sp.vstack([self.rrefs, self.rrefs])
            self.rrefs = self.rrefs.transpose().flatten()
            # then add rin and rout
            # so that it would be [rin, 10, 10, 30, 30, 100, 100, rout]
            self.rrefs = _sp.concatenate(([[self.rin], self.rrefs, [self.rout]]))
            # create start and stop values
            starts, stops = self.rrefs[0::2], self.rrefs[1::2]
            # start the construction of the grids
            print (stylify('Refinement with', fg='r'))
            for start, stop, npoints, spacing in zip(starts, stops, self.ncell, self.spacing):
                # if its the last point, include the endpoint
                endpoint_bool =  (stop == stops[-1])
                print('{0} cells : {1} - {2:.1} AU'.format(npoints, start/_cgs.AU, stop/_cgs.AU))
                # is the grid linear or log spaced?
                
                if spacing.lower() in ['lin', 'linear', 'linspace']:
                    # create that part of the grid and change rx accordingly
                    r_grids.extend(linspace(start, stop, 
                                        num = npoints, 
                                        endpoint = endpoint_bool)
                                )
                elif spacing.lower() in ['log', 'logarithm', 'logarithmic', 'logspace']:
                    # create that part of the grid and change rx accordingly
                    r_grids.extend(logspace(log10(start), log10(stop), 
                                        num = npoints, 
                                        endpoint = endpoint_bool)
                                    )
            self.rx = array(r_grids)
            #~ self.rx = self.rx
            #~ # units should be in cm
            #~ self.rrefs_cm = [i*_cgs.AU for i in self.Input.rrefs]
        else:
            print('no refinement!')
        #~ return None
        self.rx = np.insert(self.rx, 0, 0)
        self.r1 = self.rx[0:-1]
        self.r2 = self.rx[1:]
        
        #~ from scipy import dstack
        #~ self.rr = dstack((self.r1, self.r2)).ravel()
        
        #~ r1=np.insert(r[0:-1],0,0)
        #~ r2=np.array(r)
        self.rr = zeros(len(self.rx)-1, float)
        
        ################################################################
        ################################################################
        # INTERPOLATION of values
        # grid point distances allways have to be logarithmically spaces
        # does linear even make sense?
        # TODO : rr needs to be the the averaged radius in that cell!!
        # i.e. center of mass!!!
        # create rr array which is just the mean radius of each cell
        # assumes all spacings are logarithmic!!!
        self.rr[1:] = 10**( (log10(self.r1[1:]) + log10(self.r2[1:])) / 2.0 )
        
        
        ##### what if part needs to be interpolated linearly??
        # this whole implementation is more complex than what it needs to be!
        # need to rewrite this, so it is more my own code...
        
        # and the first cell has the same value as the first,
        # so the interpolated values are just copied and identical 
        # in cell 0 and 1 (except nh2int)
        self.rr[0] = self.rr[1]
        # Interpolate the values to 'ncell' cells
        self.nh2f = scipy.interpolate.interp1d(log10(self.r), log10(self.nh2))
        self.tkf = scipy.interpolate.interp1d(log10(self.r), log10(self.temp))
        self.tdf = scipy.interpolate.interp1d(log10(self.r), log10(self.tdust))
        self.abund_f = scipy.interpolate.interp1d(log10(self.r), log10(self.abund))
        self.vr_f = scipy.interpolate.interp1d(log10(self.r), log10(self.vr))
        #
        # Convert logarithms to floats
        self.nh2int = 10**self.nh2f(log10(self.rr))
        self.tkint = 10**self.tkf(log10(self.rr))
        self.tdint = 10**self.tdf(log10(self.rr))
        self.abund_int = 10**self.abund_f(np.log10(self.rr))
        self.vr_int = 10**self.vr_f(np.log10(self.rr))
        # if it is infall, multiply by -1 
        # log10 does not work all that well
        if self.velocitydirection in ['infall']:
            self.vr_int *= -1
        # round off the array so that 1E20 is 0E15
        self.vr_int = array([round(i, 15) for i in self.vr_int])
        ############################
        #~ nhint_p = nhint*2/4.
        #~ nhint_p[0] = 0.0
        #~ nhint_o = nhint*2/4.
        #~ nhint_o[0] = 0.0
        #~ teint = 10**tkf(log10(rr))
        ############################
        # nh2 needs to start at 0
        self.nh2int[0] = 0.0
        
        self.teint = self.tkint
        ################################################################
        ################################################################
        # O/P ratio of H2
        # define teint, ortho, para
        if self.Input.opr >= 0:
            self.opr = self.Input.opr
            print ('opr = {0}'.format(self.opr))
            para = 1. / (1. + self.opr)
            ortho = self.opr * para
            self.para = _sp.ones_like(self.teint) * para
            self.ortho = _sp.ones_like(self.teint) * ortho
        elif self.Input.opr == -1:
            print (stylify('Temperature dependent ortho/para', fg='r'))
            self.opr = 9.0 * exp(-170.6 / self.teint)
            self.opr = np.clip(self.opr, 1.0E-3, 3.0)
            self.para = 1.0 / (1 + self.opr)
            self.ortho = 1 - self.para
        else:
            print('Could not understand the \'op\' parameter, please correct.')
        ################################################################
        ################################################################
        # mass 
        # mass of it all
        #~ vol=[]
        #~ mass=[]
        # V = 4*pi*r**3/3
        # r in cm (?)
        #~ V = 4 * pi * (self.r2**3 - self.r1**3) / 3     # cm3
        #~ V = 4 * pi * ((self.r2**3  - self.r1**3 )) / 3     # cm3
        #~ self.M = V * self.nh2int * _cgs.MUH2 * _cgs.MP # g = cm3 * g/cm3
        # proper calculation of mass
        # TODO : the density is slightly wrong though, but shouldn't matter
        # too much for these simple models and lines.
        # need to correct this in the future though...
        rho = self.nh2int * _cgs.MUH2 * _cgs.MP
        dr = self.r2 - self.r1
        r = (self.r1 + self.r2)/2.
        self.M = 4 * _sp.pi * _sp.sum(r**2 * rho * dr)
        self.M /= _cgs.MSUN                            # Msun
        
        # to get the column density, integrate over radius r1 to r_10k
        #r_10k * 2 nh2_10k
    
        ################################################################
        #  print info. -> move to __str__ method
        print ('M_10K   : {0:<7.2f} Msun\n'
                'R_10K   : {1:<7.0f} AU\n'
                'nH2_10K : {2:<7.1e} cm-3\n'
                'Y       : {3:<7.0f}\n'
                'T       : {4:<7.1f} K\n'.format(self.M.sum(),
                                            self.r_10k/_cgs.AU,
                                            self.nh2_10k,
                                            self.Y,
                                            self.temp_10k))
        print 'Constraining the envelope : ', self.r_constraint
        print ('nH2_r1000   : {0:<7.1e} cm-3\n'
                'T_r1000     : {1:7.1f} K\n'.format(self.nh2_r1000,
                                             self.temp_r1000))
        
        print('printing input files')
        
        """
        id    : shell number
        ra,rb : inner & outer radius (m)
        za,zb : lower & upper height (m) (2D only)
        nh    : density (cm-3) of main collision partner (usually H2)
        nm    : density (cm-3) of molecule
        ne    : density (cm-3) of second collision partner (e.g. electrons)
        tk    : kinetic temperature (K) 
        td    : dust temperature (K)
        te    : electron/second coll. partner temperature (K)
        db    : 1/e half-width of line profile (Doppler b-parameter) (km s-1)
        vr    : radial velocity (km s-1)
        """
        with open(_os.path.join(self.directory, self.modelfile),'w') as f:
            f.write('# Ratran input file based on Transphere results'+'\n')
            if self.skyonly: 
                f.write('# ... intended for (SKY) continuum calculations only.'+'\n')
            f.write("rmax={0:.5E}\n".format( self.r2[-1] / 100 ))       # rmax in METERS (convert from cm i.e. / 100)
            f.write("ncell={0:}\n".format(len(self.r2)))
            f.write("tcmb=2.735\n")
            f.write("columns=id,ra,rb,nh,nm,ne,tk,td,te,db,vr\n")
            f.write("gas:dust={0}\n".format(self.gas2dust))
            if self.skyonly: 
                f.write("kappa={0}\n".format(self.kappa))
            f.write('@\n')
            # r1/r2 in meter (convert from cm)
            for ii in range(0, len(self.r1)):
                test = ("{0:4} "                         #  1 id : shell number
                        "{1:12.5E} "                     #  2 ra : inner radius (m)  
                        "{2:12.5E} "                     #  3 rb : outer radius (m)
                        "{3:12.5E} "                     #  4 nh : density (cm-3) of main coll. partner (usually H2)
                        "{4:12.5E} "                     #  5 nm : density (cm-3) of molecule
                        "{5:12.5E} "                     #  6 ne : density (cm-3) of second coll. partner (e.g. electrons)
                        "{6:12.5E} "                     #  7 tk : kinetic temperature (K) 
                        "{7:12.5E} "                     #  8 td : dust temperature (K)
                        "{8:12.5E} "                     #  9 te : second coll. partner temperature (K)
                        "{9:12.5E} "                     # 10 db : 1/e half-width of line profile (Doppler b-parameter) (km s-1)
                        "{10:12.5E}\n")                  # 11 vr : radial velocity (km s-1)
                # now print the whole shebang
                f.write(test.format(ii + 1,               #  1 id : shell number
                    self.r1[ii] / 100.0,                  #  2 ra : inner radius (m)  
                    self.r2[ii] / 100.0,                  #  3 rb : outer radius (m)
                    self.nh2int[ii] * self.para[ii],      #  4 nh : density (cm-3) of main coll. partner (usually p-H2)
                    self.nh2int[ii] * self.abund_int[ii], #  5 nm : density (cm-3) of molecule
                    self.nh2int[ii] * self.ortho[ii],     #  6 ne : density (cm-3) of second coll. partner (e.g, e^-, o-H2)
                    self.tkint[ii],                       #  7 tk : kinetic temperature (K) 
                    self.tdint[ii],                       #  8 td : dust temperature (K)
                    self.teint[ii],                       #  9 te : second coll. partner temperature (K)
                    self.db,                              # 10 db : 1/e half-width of line profile (Doppler b-parameter) (km s-1)
                    round(self.vr_int[ii], 15))           # 11 vr : radial velocity (km s-1)
                        )          
                                    
        with open(_os.path.join(self.directory, "sky.inp"),'w') as f:
            if self.skyonly:
                f.write("source={0}\n".format(self.modelfile))
            else:
                f.write("source=populations.pop\n")                     # just use the AMC output file (always set to populations.pop above)
            f.write("format={0}\n".format(self.outformat))
            f.write("outfile="+self.outputfile+"\n")
            f.write("trans={0}\n".format(self.trans))
            f.write("pix={0},{1:f},{2},{3}\n".format(self.imsize, self.pixel, self.pxlradius, self.los))
            if self.skyonly:
                f.write("chan=1,1.0\n")
            else:
                f.write("chan={0},{1:f}\n".format(self.chans, self.chwidth))
            f.write("distance={0}\n".format(self.dpc))
            f.write("units={0}\n".format(self.unit))
            f.write("go\n")
            f.write("q\n")
            f.write("\n")
            


    def run(self):
        # run ratran with the setup in the directory
        import subprocess
        from time import time, sleep
        import sys
        # Run SKY
        with ChangeDirectory(self.directory):
            f = open('sky.log', 'w')
            f.close()
            t1 = time()
            proc = subprocess.Popen([RUN_SKY, 'sky.inp'],
                                stdout = subprocess.PIPE, 
                                stderr = subprocess.STDOUT)
            self.sky_output = []
            while True:
                # first : if process is done, break the loop
                if proc.poll() != None: 
                    break
                nextline = proc.stdout.readline()
                self.sky_output.append(nextline)
                open('sky.log', 'a').write('{0}'.format(nextline))
            print('\nSKY took {0:2.1f} seconds'.format((time()-t1)))
            f.close()






