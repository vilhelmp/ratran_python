### Imports
# python builtins
import os as _os
import sys as _sys
import subprocess as _subprocess

# python extra modules
import scipy as _sp
from matplotlib import pyplot as _pl; _pl.ion()
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
from .. import cgsconst as _cgs
from ..helpers import get_colors as _gc
from ..helpers import *


class Read(object):
    def __init__(self, popfile = '', directory = '', molfile = '', kappa = '', skylogfile = ''):
        """
        This class should take a directory and ratran output file
        and just read in whatever it can. Should have reasonable
        assumption.
        
        popfile : just name of popfile
        directory : directory where the Ratran output is located
        molfile : path to the moldata file (if different from that 
                    of popfile(s))
        """
        
        ### Parse directory input
        self._check_directory(directory)
        ##### The popfile is mandatory output of RATRAN so it is 
        ##### kind of needed for this class to be of any use...
        ### Parse popfile input
        self._check_popfile(popfile)
        ### create the full path to the popfile
        self.path_popfile = _os.path.join(self.directory, self.popfile)
        ### Read the population output from Ratran
        self._read_population_output()
        ### Read the AMC input file
        self.path_amc_input_file  = _os.path.join(self.directory, 'amc.inp')
        if _os.path.isfile(self.path_amc_input_file):
            self._read_amc_input()
        ### Read the SKY input file
        self.path_sky_input_file  = _os.path.join(self.directory, 'sky.inp')
        if _os.path.isfile(self.path_sky_input_file):
            self._read_sky_input()
        ### check sky log file
        if self._check_sky_log(skylogfile):
            self._read_sky_logfile()
        ### Check fits output
        if self._check_fits:
            self._read_fits()
        ##### If there are no history files present, don't load anything.
        ### Check if history files present in directory
        if self._check_history_files(): # if so, read them
            self._read_population_history()
        else:
            pass
        ### Check molfile,
        if self._check_molfile(molfile):
            self.Moldata  = _moldata.Read(self.molfile)
        ### Dust opacity
        if self._check_kapfile(kappa):
            if 'powerlaw' in self.Kappa.kappath:
                # if we need a powerlaw
                self._kappa_powerlaw()
            else:
                # else it is jena
                self._read_kappa()

    #### checking functions
    def _check_directory(self, directory):
        if directory:   # if directory input
            self.directory = directory
        else:
            self.directory = _os.getcwd()
        
    def _check_popfile(self, popfile):
        if popfile:
            self.popfile = popfile
        else:
            try:
                self.popfile = [i for i in _os.listdir(self.directory) if i[-3:] == 'pop'][0]
            except (IndexError):
                raise Exception ('No popfile given or found, or bad directory input.')
    
    def _check_molfile(self, molfile):
        if molfile: 
            # if moldata file input
            if not _os.path.isfile(molfile):
                print('Molfile does not exists/wrong path.')
                return False
            else:
                self.molfile = molfile
                return True
        elif self.Pop.__dict__.has_key('molfile'): 
            # moldata file in popfile
            print self.Pop.molfile
            if not _os.path.isfile(self.Pop.molfile):
                print('Molfile does not exists/wrong path.')
                return False
            self.molfile = self.Pop.molfile
            return True
        else: # moldata not input, nor in Pop
            # NO moldata read, no Moldata nested class exists
            # raise warning?
            print('Molfile neither input nor in population output.')
            return False

    def _check_sky_log(self, skylogfile):
        if not skylogfile: # if no name is supplied, try with 'sky.log'
            skylogfile = 'sky.log'
            # now, check if it exists
        if not _os.path.isfile(skylogfile):
            print('Sky log file does not exists/wrong path.')
            return False
        else:
            self.skylogfile = skylogfile            
            return True
     
    def _check_kapfile(self, kappa):
        if kappa: # if kapfile was input
            if not _os.path.isfile(kappa):
                print('Kappa file does not exists/wrong path.')
                return False
            else: # it is a path to a file
                self.kappath = kappa
                return True
        elif self.Pop.__dict__.has_key('kappa'):
            # no file was input, get it from RATRAN
            if 'powerlaw' in self.Pop.kappa: # powerlaw model of opacity
                class Kappa:
                    pass
                Kappa.kappath = 'powerlaw'
                self.Kappa = Kappa
                return True
            elif self.Pop.__dict__.has_key('ratran_path') and 'jena' in self.Pop.kappa: # else its jena..
                kapfile = self.Pop.kappa.replace(',', '_')
                kapfile = kapfile.__add__('.tab')           # file format
                kapdir = _os.path.join(self.Pop.ratran_path, 'kappa')
                kappath = _os.path.join(kapdir, kapfile)
                print kappath
                if not _os.path.isfile(kappath):
                    print('Opacity (kappa) file does not exists/wrong path.')
                    return False
                else:
                    class Kappa:
                        pass
                    Kappa.kapdir = kapdir
                    Kappa.kappath = kappath
                    self.Kappa = Kappa
                    return True
        else:
            print('No dust opacity found (kappa)')
            return False

    def _check_history_files(self):
        hisfiles = [i for i in _os.listdir(self.directory) if i[-3:] == 'his']
        if not hisfiles:
            return False
        else:
            return True
    
    def _check_fits(self):
        if not hasattr(self, 'Sky'):
            return False
        else:
            return self.Sky.output

    #### reading functions    
    def _read_population_output(self):
        from scipy import arange, array
        #~ if not directory:   # if no directory was given
            #~ directory = _os.getcwd()
        class Pop(object): pass
        with open(self.path_popfile) as f:
            line = f.readline()
            Pop.comments = []
            while not line.startswith('@'):
                if line.startswith('#'):
                    Pop.comments.append(line)
                    line = f.readline()
                    pass
                else:
                    keyval = line.strip('\n').split('=')
                    try:
                        setattr(Pop, keyval[0].replace(':','_'), float(keyval[1]))
                    except(ValueError):
                        setattr(Pop, keyval[0].replace(':','_'), keyval[1])
                    line = f.readline()
            Pop.columns = Pop.columns.split(',')
            lines = f.readlines()
        # try to get the RATRAN path
        if Pop.__dict__.has_key('molfile'):
            i = Pop.molfile.split('/').index('molec')
            path = Pop.molfile.split('/')[:i]
            path.append('')
            Pop.ratran_path = '/'.join(path)
        lines = array([i.strip().split() for i in lines], dtype='float')
        lines = lines.transpose()
        for colname, i in zip(Pop.columns, arange(len(Pop.columns))):
            if colname == 'lp':
                # the lp are all the columns left,
                # the number of columns amount to the numbers of 
                # levels of that particulat molecule (see moldata file)
                setattr(Pop, colname, lines[i:]) 
            else:
                setattr(Pop, colname, lines[i])
        self.r = (Pop.rb + Pop.ra) / 2.
        self.r_au = (Pop.rb + Pop.ra) / 2. * 100 / _cgs.AU
        #[setattr(self, col, i) in zip(self.columns,arange(len(self.columns)-1))]
        self.Pop = Pop
        
    def _read_population_history(self):
        hisfiles = [i for i in _os.listdir(self.directory) if i[-3:] == 'his']
        self.hisfiles = sorted(hisfiles)
        class Pop_his(object): pass
        Pop_his.hisfiles = self.hisfiles
        # copy pasted below, check and refine
        from scipy import loadtxt, zeros, array
        tables = []
        onetime = False
        for f in self.hisfiles:    
            with open(_os.path.join(self.directory, f)) as text:
                lines = text.readlines()
                fname = text.name
            #~ print('Reading file : {0}'.format(f))
            rawhdr = [c for c in lines if len(c)<100 and not c.startswith('@')]
            preamble = [i.strip().split('=') for i in rawhdr if not i.startswith('#')]
            while not onetime:
                for keyval in preamble:
                    try:
                        setattr(Pop_his, keyval[0].replace(':','_'), float(keyval[1]))
                    except(ValueError):
                        setattr(Pop_his, keyval[0].replace(':','_'), keyval[1])
                onetime = True
                Pop_his.preamble = preamble
                #~ self.molfile = rawhdr[]
            dat = [c for c in lines if len(c)>100 and not c.startswith('@')]
            # 'raw' table data
            coldata = loadtxt(dat).transpose()
            popfile = dict()
            popfile.__setitem__('filename', str(fname))
            # get the individual settings written in the header/preamble
            [popfile.__setitem__(i[0], i[1]) for i in preamble]
            # from the settings get the columns present
            # the script assumes that the last column is level populations
            columns = [i for i in popfile['columns'].split(',')]
            [popfile.__setitem__(i, j) for (i,j) in zip(columns[:-1], coldata[:len(columns)-1])]
            # convergence percentage
            convergence = [i for i in rawhdr if i.startswith('#') and 'converged' in i]
            convergence = convergence[0].strip().strip('#AMC:').strip('converged').strip('\% ')
            popfile.__setitem__('convergence', float(convergence))
            # get the level populations table
            popfile.__setitem__('lp', coldata[len(columns)-1:])
            # write the raw data and header to the dictionary as well
            popfile.update(dict(header = rawhdr, rawdata = coldata))   
            # append to list of population file info
            tables.append(popfile)
        print('History files read.')
        Pop_his.pop_tables = tables
        self.Pop_his = Pop_his

    def _read_amc_input(self):
        class Amc:
            pass
        with open(self.path_amc_input_file) as f:
            line = f.readline()
            Amc.comments = []
            while line != 'go\n':
                #~ print line
                if line.startswith('#'):
                    Amc.comments.append(line)
                    line = f.readline()
                    pass
                else:
                    keyval = line.strip('\n').split('=')
                    try:
                        setattr(Amc, keyval[0].replace(':','_'), float(keyval[1]))
                    except(ValueError):
                        setattr(Amc, keyval[0].replace(':','_'), keyval[1])
                    line = f.readline()
            if hasattr(Amc,'kappa'):
                Amc.kappa = Amc.kappa.split(',')
        self.Amc = Amc

    def _read_sky_input(self):
        ## kind of duplicate of _read_amc_input, share and shorten?
        class Sky:
            pass
        with open(self.path_sky_input_file) as f:
            line = f.readline()
            Sky.comments = []
            while line != 'go\n':
                #~ print line
                if line.startswith('#'):
                    Sky.comments.append(line)
                    line = f.readline()
                    pass
                else:
                    keyval = line.strip('\n').split('=')
                    try:
                        setattr(Sky, keyval[0].replace(':','_'), float(keyval[1]))
                    except(ValueError):
                        setattr(Sky, keyval[0].replace(':','_'), keyval[1])
                    line = f.readline()
            # direct approach, ok for these three
            if hasattr(Sky,'pix'):
                Sky.pix = Sky.pix.split(',')
            if hasattr(Sky,'chan'):
                Sky.chan = Sky.chan.split(',')
                Sky.chan = [float(i) for i in Sky.chan]
            if hasattr(Sky,'trans'):
                try:
                    Sky.trans = Sky.trans.split(',')
                    Sky.trans = [int(i) for i in Sky.trans]
                except:
                    pass
        self.Sky = Sky

    def _read_sky_logfile(self):
        #TODO : expand to read errors, msgs etc
        # read in the whole sky log file, shouldn't be big
        f = open(self.skylogfile)
        lines = f.readlines()
        f.close()
        dust = [line.split()[1:] for line in lines if line.startswith('dtau_dust')]
        line = [line.split()[1:] for line in lines if line.startswith('dtau_line')]
        dust = _sp.array(dust, dtype='float')
        line = _sp.array(line, dtype='float')
        transitions = _sp.unique(dust[:,0])
        shells = _sp.unique(dust[:,1])
        dtau_dust = dict()
        dtau_line = dict()
        dtau_tot = dict()
        for t in transitions:
            d = []
            l = []
            for s in shells:
                d.append( _sp.mean([i[2] for i in dust if ((i[0]==t) * (i[1]==s))]) )
                l.append( _sp.mean([i[2] for i in line if ((i[0]==t) * (i[1]==s))]) )
            dtau_dust[t] = _sp.copy(d)
            dtau_line[t] = _sp.copy(l)
            dtau_tot[t] = _sp.array(d) + _sp.array(l)
        # create object to store in main class
        class Tau(object):pass
        Tau.dtau_dust = dtau_dust
        Tau.dtau_line = dtau_line
        Tau.dtau_tot = dtau_tot
        Tau.transitions = transitions
        Tau.shells = shells
        self.Tau = Tau
       
    def _read_fits(self):
        print('Fits reading not implemented.')

    def kappa(self, nu, inunit='freq'):
        if inunit == 'freq':
            cm = _cgs.CC / nu
        elif inunit == 'mum':
            cm = nu * 1e-4
        
        from scipy import interpolate
        value = interpolate.splev(_sp.log10(cm), self.Kappa.tck)
        return 10**(value)

    def _read_kappa(self):
        self.Kappa.table, self.Kappa.tck = helpers.get_kappa_jena(self.Kappa.kappath)

    ### other functions
    def _kappa_powerlaw(self):
        # redefine kappa to be a powerlaw function
         self.Kappa.tck = helpers.kappa_powerlaw( nu0 = 1e12, kappa0 = 15, beta = 1.8 )

    #### plotting functions
    ## for final populations
    def plot_structure(self):
        pass
        
    #~ def plot_tau(self, trans = [12, 10], width = 1.0E4):
        #~ """
        #~ Plot opacity/optical depth (tau) and the optical depth in LTe
        #~ of each cell
        #~ """
#~ 
        #~ radii = self.r_au
        #~ 
        #~ # ra, rb in m, convert to cm
        #~ # width in cm/s
        #~ itup, itdown = trans[0] - 1, trans[1] - 1
#~ 
        #~ transIndex = helpers.get_transindex(self, trans)
        #~ 
        #~ # First, get the opacity from the model outpupt
        #~ class Tau: pass
        #~ self.Tau = Tau
        #~ # the density of the upper energy level
        #~ # from the model data
        #~ self.Tau.ds = (self.Pop.rb - self.Pop.ra) * 100 # in cm, not m
        #~ self.Tau.Nu = Nu = self.Tau.ds * self.Pop.nm * self.Pop.lp[itup]
        #~ self.Tau.Nu = nu = self.Pop.nm * self.Pop.lp[itup]
        #~ self.Tau.Nl = Nl = self.Tau.ds * self.Pop.nm * self.Pop.lp[itdown]
        #~ self.Tau.Nl = nl = self.Pop.nm * self.Pop.lp[itdown]
        #~ self.Tau.Au = Aul = self.Moldata.radtrans[transIndex]['aul']
        #~ self.Tau.freq = freq = self.Moldata.radtrans[transIndex]['freq']
#~ 
        #~ self.Tau.Eu = Eu = self.Moldata.radtrans[transIndex]['eu']
        #~ self.Tau.gu = gu = self.Moldata.elev[itup]['weight']
        #~ self.Tau.gl = gl = self.Moldata.elev[itdown]['weight']
#~ 
        #~ self.Tau.T = T = self.Pop.tk
        #~ # now calculate the opacity
        #~ # tau = helpers.calc_dtau(Nu, Nl, Aul, freq, width, gu, gl, T)
            #calc_alpha_line()
#~ 
        #~ # self.nh2 = self.rhodust * 100  / (_cgs.MUH2 * _cgs.MP)
        #~ 
        #~ 
        #~ alpha_line = calc_alpha_line(nu, nl, Aul, freq, gu, gl, self.Pop.db)
        #~ alpha_dust = self.kappa(freq) * self.Pop.nh * _cgs.MUH2 * _cgs.MP/ self.Pop.gas_dust
        #~ # tau = calc_tau(Nu, Aul, freq, width, Eu, T)
        #~ # inull = _sp.where(tau == 0)
        #~ # tau[inull] = 1E-15
        #~ tau = (alpha_line + alpha_dust) * self.Tau.ds
        #~ self.Tau.tau = tau
        #~ self.Tau.alpha_line = alpha_line
        #~ self.Tau.alpha_dust = alpha_dust
        #~ 
        #~ 
        #~ # plotting
        #~ _pl.ion()
        #~ # fig = _pl.figure(num=1, figsize=(3.5,3))
        #~ fig = _pl.figure(num=1)
        #~ ax = fig.add_subplot(111)
        #~ ax.loglog(radii , alpha_line*self.Tau.ds, label=r'line', color='#559922', lw=2, marker='o', ms=3, mew=0)
        #~ ax.loglog(radii , alpha_dust*self.Tau.ds, label=r'dust', color='#992255', lw=2, marker='o', ms=3, mew=0)
        #~ ax.loglog(radii , tau, label=r' - '.join([self.Moldata.get_lvl(i, tex = 1) for i in trans]), color=_gc(1).next(), lw=2, marker='o', ms=3, mew=0)
#~ 
        #~ 
        #~ # Second, calculate the opacity if it is pure LTE conditions
        #~ 
        #~ try:
            #~ # first initialise the get_partition method
            #~ # to get the partition function values
            #~ # at all temperatures
            #~ self.Moldata.get_partition()
            #~ plot_lte = False
        #~ except:
            #~ print ('cannot get the partition function values.')
            #~ plot_lte = False
            #~ 
        #~ if plot_lte:
            #~ # Nu  = nm * (rb - ra) * gu / Qrot(T) * exp(-Eu/T)
            #~ nm = self.Pop.nm
            #~ ra = self.Pop.ra * 100 # in cm
            #~ rb = self.Pop.rb * 100 # in cm
            #~ gu = self.Moldata.elev[itup]['weight']
            #~ qrot = self.Moldata.qrot(T)
            #~ self.Tau.Nu_lte = Nu_lte = nm * (rb - ra) * gu / qrot * _sp.exp(-Eu/T)
            #~ tau_lte = helpers.calc_tau_lte(Nu_lte,  Aul, freq, width, Eu, T)
            #~ # inull = _sp.where(tau_lte == 0)
            #~ # tau_lte[inull] = 1E-15
            #~ self.Tau.tau_lte = tau_lte
            #~ # plot it
            #~ _pl.loglog(radii, tau_lte, label=r' LTE', color='#008822', lw=2, marker='o', ms=3, mew=0)
#~ 
        #~ # plot end of model
        #~ trans1 = _btf(ax.transData, ax.transAxes)
        #~ linesetting = dict(color=_gc(1).next(), transform=trans1, lw=1,
                           #~ ls='dashed')
        #~ ax.semilogx([radii[-1], radii[-1]], [0, 1], **linesetting)
        #~ ax.semilogx([radii[0], radii[0]], [0, 1],**linesetting )
        #~ # labels, legend and grid
        #~ ax.set_xlabel('Radius [AU]')
        #~ ax.set_ylabel(r'$d\tau$')
        #~ ax.legend()
        #~ ax.grid()

    def plot_populations(self, levels = [], runjump = 10, leveljump = 10):
        """ 
        Function to plot the populations level for 
        every 'runjump' iteration and every 'leveljump' population 
        level
        """
        from scipy import linspace, log10, logspace, array
        #~ import matplotlib.pyplot as pl; 
        _pl.ion()
        from matplotlib import cm
        #~ from adapy.libs import cgsconst as cgs

        radii = self.r_au
               
        lenlvls = len(self.Pop.lp)
        if not levels:
            levels = _sp.arange(1, lenlvls+1, leveljump)
        #
        if hasattr(self, 'Pop_his'):
            tables_select = self.Pop_his.pop_tables[::10]
            lp = array([[i['lp'][j-1] for j in levels] for i in self.Pop_his.pop_tables[::int(runjump)]])
            # plot all iterations
            [[_pl.loglog(radii , i, color=str(c), lw=1, ls=':', marker='o', ms=3, mew=0) for i in j] for (j,c) in zip(lp, linspace(0.7, 0.2, len(lp)))]
        #
        # should plot the resulting 
        # a bit to complicated list comprehension... reformat to better
        # programming style...
        [_pl.loglog(radii , self.Pop.lp[i], label=self.Moldata.get_lvl(i+1, tex = 1), color=c, lw=1.5, marker='o', ms=3, mew=0) for (i,c) in zip(array(levels)-1, _gc(len(levels)))]
        _pl.legend(loc=3, numpoints=1, ncol=len(levels)/8+1)
        _pl.grid()
        
    def plot_tex(self, trans = [12, 10], runjump = 10, history = True):    # specify which transition
        """
        Plots the excitation temperature for every 'runjump' iteration 
        and every 'leveljump' population level
        """
        from scipy import linspace, array
        #~ import matplotlib.pyplot as pl; 
        _pl.ion()
        #~ from matplotlib import cm
        #~ from adapy.libs import cgsconst as cgs
        
        radii = (self.Pop.ra + self.Pop.rb)/2. * 100 / _cgs.AU
        if hasattr(self, 'Pop_his') and history:
            ### get the populations for the levels from the different runs 
            runs = array([i for i in self.Pop_his.pop_tables[::int(runjump)]]) # every 10th run
            lp = [[j['lp'][i-1] for i in trans] for j in runs]
        
        ### 
        gweights = [self.Moldata.elev[trans[0]-1]['weight'], self.Moldata.elev[trans[1]-1]['weight']] 
        print('   Molecular weights : upper {0[0]}  lower {0[1]}'.format(gweights))
        #~ nu = _cgs.CC*abs(self.Moldata.elev[trans[0]-1]['energies'] - self.Moldata.elev[trans[1]-1]['energies'])
        transIndex = get_transindex(self, trans)
        eu =  self.Moldata.radtrans[transIndex]['eu']
        nu = self.Moldata.radtrans[transIndex]['freq']
        print('   Frequency : {0:3.3f} E09 GHz'.format(nu*1e-9))
        trans_str = [self.Moldata.get_lvl(i) for i in trans]
        #trans_str = [self.elev[trans[0] ,.elev[trans[1]-1]['j']]
        print('   Transition : {0[0]} - {0[1]}'.format(trans_str))
        if hasattr(self, 'Pop_his') and history:
            tex = []
            for levels in lp:
                t = temp_pop(levels[1], levels[0], gweights[1], gweights[0], eu)
                tex.append(t)
            [_pl.loglog(radii, j, color=str(c), lw=1, marker='o', ms=4, mew=0)
                for (j, c) in zip(tex, linspace(0.7, 0, len(tex)))]
        #
        ### get the final tex curve
        tex_final = temp_pop(self.Pop.lp[trans[1]-1], self.Pop.lp[trans[0]-1],
                             gweights[1], gweights[0], eu)
        str_transitions = [self.Moldata.get_lvl(i, tex=1) for i in trans]
        _pl.loglog(radii, tex_final, label=r' - '.join(str_transitions),
                   color=_gc(1).next(), lw=2, marker='o', ms=5, mew=0)
        #
        _pl.legend()
        _pl.grid()

        print(' Return the calculated Tex')
        return tex_final
        
        # tex for all transitions? why do that...
        # perhaps be able to supply several pairs of
        # transitions but not more
        #~ tex = [[temp_pop(ldown, lup, gweights[1], gweights[0], nu) for
            #~ (ldown, lup) in zip(run[1], run[0])] for (run) in lp]
        #~ return tex, lp
        #~ for runlp in tex:
            #~ for

    def plot_radiation(self, trans=[12,10]):
        """
            Plot the radiation field in each cell
            and the LTE radiation field
        """
        nu = self.Pop.lp[trans[0]-1]
        nl = self.Pop.lp[trans[1]-1]
        gu = self.Moldata.elev[trans[0]-1]['weight']
        gl = self.Moldata.elev[trans[1]-1]['weight']
        
        itup, itdown = trans[0] - 1, trans[1] - 1
        try:
            transIndex = [i['trans'] for i in self.Moldata.radtrans if i['up'] == trans[0] and i['down'] == trans[1]][0]
        except (IndexError):
            print('No such transition for molecule : {0}'.format(self.Moldata.molecule))
            print('Check your moldata file at : {0}'.format(self.Moldata.molfile))
            return False
        transIndex -= 1
        
        freq = self.Moldata.radtrans[transIndex]['freq']
        
        
        snu = calc_radiation(nu, nl, freq, gu, gl)
        radii = self.r_au
        _pl.ion()
        #~ fig = _pl.figure(num=1, figsize=(3.5,3))
        fig = _pl.figure(num=1)
        ax = fig.add_subplot(111)
        _pl.loglog(radii, snu, color='#008822', lw=2, marker='o', ms=3, mew=0)

        # plot end of model
        trans1 = _btf(ax.transData, ax.transAxes)
        linesetting = dict(color=_gc(1).next(), transform=trans1, lw=1,
                           ls='dashed')
        ax.loglog([radii[-1], radii[-1]], [0, 1], **linesetting)
        ax.loglog([radii[0], radii[0]], [0, 1],**linesetting )
        # labels, legend and grid
        ax.set_xlabel('Radius [AU]')
        ax.set_ylabel(r'S$_\nu$')
        #~ ax.legend()
        ax.grid()

    def plot_tau(self):
        N = len(self.Tau.transitions)
        if N % 3:
            np = N / 3 + 1
        else:
            np = N / 3
        ncol = 3
        if N <3:
            ncol = N
        for i in self.Tau.transitions:
            _pl.subplot(np, ncol, i)
            _pl.step(self.r_au[1:], self.Tau.dtau_tot[i], 'k', label='Total', lw=1.5)
            _pl.step(self.r_au[1:], self.Tau.dtau_line[i], 'b', label='Line', lw=1.5)
            _pl.step(self.r_au[1:], self.Tau.dtau_dust[i], 'g', label='Dust', lw=1.5)
            _pl.title(str(int(i)))
            _pl.xlabel('r (AU)')
            _pl.ylabel(r'$\tau$')
            _pl.xscale('log')
            _pl.yscale('log')
            _pl.legend()

def _read_populations_simple(popfile):
    """
        Function to read in the populations file from RATRAN.
        Outputs the comments and tables as separate lists of strings.
    """
    with open(popfile) as f:
        line = f.readline()
        comments = []
        # get the comments
        while not line.startswith('@'):
            comments.append(line)
            line = f.readline()
        table = f.readlines()
    return comments, table


def _run_cell(args):
    folder = args['folder']
    cells = args['cells']
    comments = args['comments'][:]
    table = args['table'][:]
    skyinp = args['skyinp']
    popfile = args['popfile']
    sky_text = args['sky_text'][:]
    # peel of cells
    table = table[:cells]
    # change rmax
    irmax = [i for i in xrange(len(comments)) if comments[i][:4] == 'rmax'][0]
    comments[irmax] = 'rmax={0}\n'.format(table[-1].split()[2])
    # change ncell
    incell = [i for i in xrange(len(comments)) if comments[i][:5] == 'ncell'][0]
    comments[incell] = 'ncell={0}\n'.format(len(table))
    # create the popfile text to write
    poptext = comments[:]
    poptext.append('@\n')
    temp = [poptext.append(i) for i in table]
    # write sky.inp and populations.pop
    with helpers.ChangeDirectory(folder):
        # write sky.inp
        with open(skyinp, 'w') as f:
            f.writelines(sky_text)
        with open(popfile, 'w') as f:
            f.writelines(poptext)
        # call sky on the setup
        #~ _os.system('sky {0}'.format(skyinp))
        import subprocess
        RUN_SKY = '/home/magnusp/applications/Ratran/bin/sky'
        #~ proc = subprocess.Popen([RUN_SKY, skyinp])
        #~ _os.system(RUN_SKY+' '+skyinp)
        f = open('sky.log', 'w')
        f.close()
        print('Running sky in {0}'.format(folder))
        proc = subprocess.Popen([RUN_SKY, skyinp],
                            stdout = subprocess.PIPE, 
                            stderr = subprocess.STDOUT)
        while True:
            # first : if process is done, break the loop
            if proc.poll() != None: 
                break
            nextline = proc.stdout.readline()
            #~ print nextline
            open('sky.log', 'a').write('{0}'.format(nextline))


def _get_intensity_simple(fits_file):
    from pyfits import getdata
    data = getdata(fits_file)
    contsub = data - data[0]
    intensity = contsub.sum()
    return intensity

    
def peel_onion( skyinp='sky.inp', popfile='populations.pop', fits_file='_007.fits' ):

    #~ from scipy import arange
    
    comments, table  = _read_populations_simple(popfile)

    t = table[:]
    t = [i.strip('\n') for i in t]
    t = [i.split() for i in t]
    from scipy import array
    ra = array(t, dtype='float').transpose()[1] * 100/_cgs.AU
    rb = array(t, dtype='float').transpose()[2] * 100/_cgs.AU
    r = (ra+rb)/2
    with open(skyinp) as f:
        sky_text = f.read()
    # make the folders
    no_cells = range(1,len(table) + 1)
    folders = ['onion/{0}'.format(i) for i in no_cells]
    try:
        temp = [_os.makedirs('{0}'.format(i)) for i in folders]
    except(OSError):
        intensities = []
        for folder in folders:
            file_to_open = _os.path.join(folder, fits_file)
            print file_to_open
            intensity  = _get_intensity_simple(file_to_open)
            intensities.append(intensity)
        return intensities, table, comments, r
        
    #~ [_os.makedirs(i) for i in folders]


    
    #~ # run in parallell!
    #~ try:
        #~ import pprocess
        #~ import multiprocessing
        #~ multi = True
    #~ except(ImportError):
        #~ multi = False
    #~ if multi:
#~ 
        #~ print('running in parallel (pprocess)')
        #~ from time import time
        #~ t1 = time()
        #~ list_of_args = [dict(folder=i,
                #~ cells=j,
                #~ comments=comments,
                #~ table=table,
                #~ skyinp=skyinp,
                #~ popfile=popfile,
                #~ sky_text=sky_text) for (i,j) in zip(folders, no_cells)]
        #~ nproc = pprocess.get_number_of_cores() - 1 
        #~ results = pprocess.Map(limit=nproc, reuse=1)
        #~ parallel_function = results.manage(pprocess.MakeReusable(_run_cell))
        #~ [parallel_function(args) for args in list_of_args]
        #~ parallel_results = results
        #~ return parallel_results
        #~ 
    #~ elif not multi:
    print('not running in parallel (not implemented)')
    from time import time
    t1 = time()
    for folder, cells in zip( folders, no_cells):
        #~ print folder, cells, comments, table, skyinp, popfile , sky_text
        
        to_send = dict(folder=folder,
            cells=cells,
            comments=comments,
            table=table,
            skyinp=skyinp,
            popfile=popfile,
            sky_text=sky_text)
        # I've got a peeling!
        _run_cell(to_send)
        
    print 'This took: '
    print (time()-t1)
    print 'seconds'

    intensities = []
    for folder in folders:
        file_to_open = _os.path.join(folder, fits_file)
        print file_to_open
        intensity  = _get_intensity_simple(file_to_open)
        intensities.append(intensity)
   
    return intensities, table, comments, r


def calc_alpha_line(nu, nl, Aul, freq, gu, gl, db):
    """
    Calculate the absorption (alpha) for given parameters
    from Rybicki & Lightman.
    """
    part1 = _cgs.CC**3 * Aul / ( 8 * _sp.pi * freq**4 * db)
    part2 = nl * gu / gl - nu
    alpha = part1 * part2
    
    #~ alpha = _cgs.HH * nu / (4 * _scipiy.pi)
    #~ dtau = _cgs.CC**2 / (8 * _sp.pi * freq**2) * Aul * (Nl * gu / gl - Nu) 
    #~ dtau = alpha * ds
    
    
    #~ part1 = _cgs.CC**3 * Aul / (8 * _sp.pi * freq**3 * width) 
    #~ part2 = (_sp.exp(Eu/T) - 1)
    #~ part2 = (Nl * gu / float(gl) - Nu)
    #~ return part1 * part2
    return alpha


def plot_onion(inp='_007.fits'):

    intens,tab,com,r = peel_onion(fits_file=inp)

    fig = _pl.figure()
    fig.clf()
    ax = fig.add_subplot(111)

    inten_percent = intens/max(intens)*100
    ax.step(r[1:], inten_percent[1:], where='post', lw=4, color='0.6')
    ax.step(r[1:], inten_percent[1:], where='post', lw=2, color='#33AA33')
    ax.set_xscale('log')
    ax.set_ylim((-5, 105))
    ax.grid(which='minor', axis='x')

    ax2 = ax.twinx()
    inten_add = [intens[i]-intens[i-1] for i in xrange(1,len(intens))]
    #inten_add = [intens[i+1]-intens[i] for i in xrange(0,len(intens))]
    ax2.step(r[1:], inten_add, where='post', lw=2, color='0.6')
    ax2.step(r[1:], inten_add, where='post', lw=1, color='#3333AA')

    ymin,ymax = min(inten_add), max(inten_add)
    ymin, ymax = ymin-0.15*abs(ymin), ymax+0.35*abs(ymax)
    ax2.set_ylim((ymin, ymax))

    ax.set_xlim((r[1]*0.95, r[-1]*1.05))

    ax.set_xlabel('AU')
    ax.set_ylabel('% of max intensity', weight='bold', size=12, color='#33AA33')
    ax2.set_ylabel('Added intensity', weight='bold', size=12,color='#3333AA')
    return ax, ax2




