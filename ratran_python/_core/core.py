

class Input(object):pass

class Output(object):
    def __init__(self, popfile = '', directory = '', molfile = ''):
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
            print('Reading file : {0}'.format(f))
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
        Pop_his.pop_tables = tables
        self.Pop_his = Pop_his

    def _read_amc_input(self):
        class Amc:
            pass
        with open(self.path_amc_input_file) as f:
            line = f.readline()
            Amc.comments = []
            while line != 'go\n':
                print line
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
                print line
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
                Sky.trans = Sky.trans.split(',')
                Sky.trans = [int(i) for i in Sky.trans]
        self.Sky = Sky

    def _read_fits(self):
        print('Tjo')
    
    #### plotting functions
    ## for final populations
    def plot_structure(self):
        pass
        
    def plot_tau(self, trans = [12, 10], width = 1.0E4):
        """
        Plot opacity/optical depth (tau) and the optical depth in LTe
        of each cell
        """

        radii = self.r_au
        
        # ra, rb in m, convert to cm
        # width in cm/s
        itup, itdown = trans[0] - 1, trans[1] - 1

        transIndex = get_transindex(self, trans)
        
        # First, get the opacity from the model outpupt
        class Tau: pass
        self.Tau = Tau
        # the density of the upper energy level
        # from the model data
        self.Tau.Nu = Nu = (self.Pop.rb - self.Pop.ra) * 100.0 * self.Pop.nm * self.Pop.lp[itup]
        self.Tau.Au = Aul = self.Moldata.radtrans[transIndex]['aul']
        self.Tau.freq = freq = self.Moldata.radtrans[transIndex]['freq']
        self.Tau.Eu = Eu = self.Moldata.radtrans[transIndex]['eu']
        self.Tau.T = T = self.Pop.tk
        # now calculate the opacity
        tau = calc_tau(Nu, Aul, freq, width, Eu, T)
        #~ inull = _scipy.where(tau == 0)
        #~ tau[inull] = 1E-15
        self.Tau.tau = tau
        
        
        
        # plotting
        _pl.ion()
        #~ fig = _pl.figure(num=1, figsize=(3.5,3))
        fig = _pl.figure(num=1)
        ax = fig.add_subplot(111)
        ax.loglog(radii , tau, label=r' - '.join([self.Moldata.get_lvl(i, tex = 1) for i in trans]), color=_gc(1).next(), lw=2, marker='o', ms=3, mew=0)

        
        # Second, the calculate the opacity if it is pure LTE conditions
        
        try:
            # first initialise the get_partition method
            # to get the partition function values
            # at all temperatures
            self.Moldata.get_partition()
            plot_lte = True
        except:
            print ('cannot get the partition function values.')
            plot_lte = False
            
        if plot_lte:
            # Nu  = nm * (rb - ra) * gu / Qrot(T) * exp(-Eu/T)
            nm = self.Pop.nm
            ra = self.Pop.ra * 100 # in cm
            rb = self.Pop.rb * 100 # in cm
            gu = self.Moldata.elev[itup]['weight']
            qrot = self.Moldata.qrot(T)
            self.Tau.Nu_lte = Nu_lte = nm * (rb - ra) * gu / qrot * _scipy.exp(-Eu/T)
            tau_lte = calc_tau(Nu_lte, Aul, freq, width, Eu, T)
            #~ inull = _scipy.where(tau_lte == 0)
            #~ tau_lte[inull] = 1E-15
            self.Tau.tau_lte = tau_lte
            # plot it
            _pl.loglog(radii, tau_lte, label=r' LTE', color='#008822', lw=2, marker='o', ms=3, mew=0)

        # plot end of model
        trans1 = _btf(ax.transData, ax.transAxes)
        linesetting = dict(color=_gc(1).next(), transform=trans1, lw=1,
                           ls='dashed')
        ax.loglog([radii[-1], radii[-1]], [0, 1], **linesetting)
        ax.loglog([radii[0], radii[0]], [0, 1],**linesetting )
        # labels, legend and grid
        ax.set_xlabel('Radius [AU]')
        ax.set_ylabel(r'$\tau$')
        ax.legend()
        ax.grid()
    
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
            levels = _scipy.arange(1, lenlvls+1, leveljump)
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
