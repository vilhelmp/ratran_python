
from .. import cgsconst as _cgs

import scipy as _scipy
# perhaps put help functions in separate file?
# help.py, extra.py,



### Help functions
def bplanck(nu, T):
    """
        Returns the Spectral Radiance (Planck curve)
    """
    from scipy import exp, constants
    try:
        x = _cgs.HH * nu / (_cgs.KK * T)
        bpl = (2.0 * _cgs.HH * nu**3 / _cgs.CC**2) / (exp(x)-1.0)
        return bpl
    except (ZeroDivisionError):
        return 0

def find_intensity(fitsfile, interval = [], nsig = 3):
    from scipy import array, arange, where, log, pi, meshgrid
    import matplotlib.pyplot as pl
    from pyfits import getdata, getheader
    from adapy.adacore import gaussfit2d
    pl.ion()
    class SkyModel: pass

    data = getdata(fitsfile)
    SkyModel.continuum = data[0]
    #~ data -= data[0] # remove continuum
    SkyModel.data = data - SkyModel.continuum
    header = getheader(fitsfile)
    
    # get some of the stuff from the header
    SkyModel.bunit = header['BUNIT']
    SkyModel.restfreq = header['RESTFREQ']
    
    # create the velocity array, just for fun
    v_cdelt = header['CDELT3']*1e-3     # in km/s
    v_crpix = header['CRPIX3']
    v_crval = header['CRVAL3']
    v_naxis = header['NAXIS3']
    
    v_array = arange(v_naxis) - v_crpix
    v_array *= v_cdelt
    SkyModel.v_array = v_array + v_crval
    SkyModel.v_cdelt = v_cdelt         # in km/s
    
    
    SkyModel.ra_cdelt = header['CDELT1']*3600
    SkyModel.dec_cdelt = header['CDELT2']*3600
    SkyModel.ra_array = ((arange(header['NAXIS1']) - header['CRPIX1']) * header['CDELT1']*3600) + header['CRVAL1']
    SkyModel.dec_array = ((arange(header['NAXIS2']) - header['CRPIX2']) * header['CDELT2']*3600) + header['CRVAL2']
    
    # assume model peak in center
    z, y ,x = SkyModel.data.shape
    SkyModel.spectrum = SkyModel.data[:, y/2, x/2]
    
    if len(SkyModel.data.shape) < 3: # needs to be a cube for analysis to work
        print("Wrong data shape of input fits file")
        return 0
    if interval == []: # if no interval given, need to do it interactively
        from adapy.adacore import fit_gauss1d as gaussfit
        #~ from pylab import ginput
        #~ from matplotlib.widgets import Cursor
        #~ fig = pl.figure()
        #~ ax = fig.add_subplot(111) 
        #~ ax.plot(SkyModel.v_array, SkyModel.spectrum)
        #~ cursor = Cursor(ax, color='red', linewidth=2 )
        #~ print('Click on lower limit')
        #~ x_low  = ginput()[0][0]
        #~ print('Click on upper limit')
        #~ x_high  = ginput()[0][0]
        #~ fig.close()
        #~ SkyModel.mom0 = 
        # simple guesses/assumptions, 
        # perhaps extend to calculate moment0/1 for the position
        # and width of the distribution?
        datamax = SkyModel.spectrum.max()
        # where is the max?
        ii = where(SkyModel.spectrum.max() == SkyModel.spectrum)[0] 
        width_estimate = (SkyModel.v_array.max() - SkyModel.v_array.min()) * 0.2
        # fit a 1D Gaussian
        results_1d = gaussfit((SkyModel.v_array, SkyModel.spectrum), 
                            params=(
                                    datamax,
                                    SkyModel.v_array[ii], 
                                    width_estimate
                                    ),
                            verbose=0
                            )[0]
        SkyModel.results_1d = results_1d
        amplitude_1d = results_1d[2]
        position_1d = results_1d[1]
        fwhm_1d = results_1d[2]
        sigmafromfwhm = 1 / (2 * (2 * log(2))**.5)
        sigma_1d = fwhm_1d * sigmafromfwhm
        
        interval = position_1d + array([-1, 1]) * nsig * sigma_1d
        print("Integration interval : 1D Gaussian fit"
                " (+/- {0} sigma)".format(nsig))
        SkyModel.interval = interval
    else:
        SkyModel.interval = interval
        print("Integration interval : input")
        
    indices = where(
                        (SkyModel.v_array >= interval[0]) * 
                        (SkyModel.v_array <= interval[1])
                        )
    
    SkyModel.zero = SkyModel.data[indices].sum(axis=0) * abs(ModelData.v_cdelt)
    X, Y = meshgrid(arange(header['NAXIS1']),arange(header['NAXIS2']))
    #~ results_2d = gaussfit2d((X, Y, ModelData.zero), params=(0.0, (0.1, 64, 64, 2, 2, 0)))[0]
    results_2d = gaussfit2d((X, Y, SkyModel.zero), fitheight=0)[0]
    # Volume (integral) of the Gaussian
    # V = 2 pi Amp sigma1 sigma2
    SkyModel.amplitude_2d = results_2d[0]
    SkyModel.sigma_x = results_2d[3]
    SkyModel.sigma_y = results_2d[4]
    SkyModel.results_2d = results_2d
    
    SkyModel.intensity = pi * 2 * SkyModel.amplitude_2d * SkyModel.sigma_x * SkyModel.sigma_y
    print('Integrated intensity : {0:.2f} Jy'.format(SkyModel.intensity))
    return SkyModel

def _tex_transition(n1, n2, g1, g2, nu):
    """ 
    Calculate the excitation temperature of transition from 2 to 1
    """
    
    numer = _cgs.HH * nu
    denom = _cgs.KK * _scipy.log(n1 * g2 / (n2 * g1))
    return numer / denom

def _pdfcheck(pdf):
        if not pdf: 
            _plt.ion()
        elif pdf:
            _plt.ioff()

def _pdfsave(pdf, pdfname, **kwargs):
    if pdf:
            _plt.savefig('{0}.pdf'.format(str(pdfname)), kwargs)

def read_molfile(filepath, output='dictionary'):
    return True


def convolve():
    pass

def get_colors(color):
    import colorsys as _cs
    for hue in range(color):
        hue = 1. * hue / color
        col = [int(x) for x in _cs.hsv_to_rgb(hue, 1.0, 230)]
        yield "#{0:02x}{1:02x}{2:02x}".format(*col)

class ChangeDirectory:
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = _os.getcwd()
        _os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        _os.chdir(self.savedPath)
    

### Imports
# python builtins
import os as _os
import sys as _sys
import subprocess as _subprocess

def test_thick(abund, flux, removelast = 1):
    """
        Test the linearity of some values.
    """
    from scipy import polyfit
    from scipy import arange
    print('removelast={0} points in fit'.format(removelast))
    p = polyfit(abund[:-1*int(removelast)], flux[:-1*int(removelast)], deg=1)
    line = lambda x: p[0]*x + p[1]    
    import matplotlib.pyplot as pl
    pl.ion()


    off = (line(abund[-1]) - flux[-1])/line(abund[-1])*100
    print('last value off by : {0:.2f}% from fit'.format(off))
    
    pl.plot(abund, flux, 'or', label='data')
    X = arange(abund[0]*0.9, abund[-1]*1.1, abund[0]*0.1)
    
    pl.plot(X, line(X),'-g', label='fit')
    pl.legend(loc=4)


########################################################################
# OLD HELP FUNCTIONS

def read_transphereoutput(self, ext = 0):
    from scipy import array
    if ext == 0: ext=''
    filetoread = 'envstruct' + ext + '.dat'
    path_to_file = _os.path.join(self.directory, filetoread)
    with open(path_to_file, 'r') as f:
        lines = f.read().split('\n')
    nr = int(lines[0])
    dat_envstruct = array([i.split() for i in lines[2:nr + 2]], dtype='float')


    filetoread = 'spectrum.dat'
    path_to_file = _os.path.join(self.directory, filetoread)
    with open(path_to_file, 'r') as f:
        lines = f.read().split('\n')
    nr = int(lines[0])
    dat_spectrum = array([i.split() for i in lines[2:nr+2]], dtype='float')

    #~ nr = int(f.readline().strip())
    #~ dat = np.zeros((3,nr),float)
    #~ for ii in range(0,3):
        #~ for jj in range(0,nr):
            #~ dum = f.readline().strip()
            #~ dat[ii,jj] = dum
    #~ f.close()

    class Envstruct:
        r = dat_envstruct[:,0]
        rho_dust = dat_envstruct[:,1]
        temp = dat_envstruct[:,2]

    #~ class Spectrum:
    Envstruct.frequency = dat_spectrum[:,0]
    Envstruct.intensity = dat_spectrum[:,1]
    #~ Envstruct.Spectrum = Spectrum
    #~ self.Envstruct = Envstruct

    #~ import numpy as np
    filetoread = 'convhist.info'
    path_to_file = _os.path.join(self.directory, filetoread)
    f = open(path_to_file, 'r')
    nn = int(f.readline().strip().split()[0])
    f.close()

    # Convergence history
    filetoread = 'convhist.dat'
    path_to_file = _os.path.join(self.directory, filetoread)
    with open(path_to_file, 'r') as f:
        lines = f.read().split('\n')
    nr = int(lines[0].strip())
    if nr == 0: raise Exception('Nothing run, no convergence history.')
    x1 = nr+1

    #These need to depend on value of nr
    dat1 = array([i.split() for i in lines[1:x1]], dtype='float')
    dat2 = array([i.split() for i in lines[x1+1:x1*2]], dtype='float')
    dat3 = array([i.split() for i in lines[x1*2+1:x1*3]], dtype='float')

    dat = array([dat1,dat2,dat3])

    #~ f = open('convhist.dat','r')
    #~ nr = int(f.readline().strip())
    #~ dat = np.zeros((9,nn,nr),float)
    #~ for jj in range(0,nn):
        #~ for kk in range(0,nr):
            #~ dum = f.readline().strip().split()
            #~ if dum == []: dum=f.readline().strip().split()
            #~ dat[0:9,jj,kk]=np.array(dum,dtype=float)
    #~ f.close()

#    if nn gt 1 then idx=[1,2,0] else idx=[1,0]. Note transpose commands not executed...
    class Convhist:
        temp=dat[:,:,0]
        jjme=dat[:,:,1]
        hhme=dat[:,:,2]
        jj=  dat[:,:,3]
        hh=  dat[:,:,4]
        kapt=dat[:,:,5]
        kapj=dat[:,:,6]
        kaph=dat[:,:,7]
        fj=  dat[:,:,8]
    #~ self.Convhist = Convhist

    #~ f = open('envstruct.inp')
    #~ nr = int(f.readline().strip())
    #~ dat = np.zeros((3,nr),float)
    #~ for ii in range(0,nr):
        #~ dum=f.readline().strip().split()
        #~ if dum == []: dum=f.readline().strip().split()
        #~ dat[0:3,ii]=np.array(dum,dtype=float)
    #~ r=dat[0,:]
    #~ f.close()

    #~ convhist={'r': r, 'temp': temp, 'jjme': jjme, 'hhme': hhme, 'jj': jj, 'hh': hh, 'kapt': kapt, 'kapj': kapj, 'kaph': kaph, 'fj': fj}
    #~ self.Envstruct = envstruct
    #~ self.convhist = convhist
    return Envstruct, Convhist

def create_grid(r_in, r_out, nshell, space = 'powerlaw1', end = True):
    # function to create grid
    if space == 'log10':
        from scipy import log10, logspace
        # get the exponent of the start- and
        # stop-radius in input units
        start = [log10(r_in), 0][r_in == 0]
        stop = log10(r_out)
        radii = logspace(start, stop, num=nshell, endpoint=end)
    elif space == "powerlaw1":
        from scipy import arange
        radii = r_in * (r_out/r_in)**(arange(nshell)/(nshell - 1.0))
    elif space == 'linear':
        from scipy import linspace
        # linearly spaced grid
        radii = linspace(r_in, r_out, num=nshell, endpoint=end)
    elif space == 'powerlaw2':
        from scipy import linspace
        # first check if coefficients to the power-law was given
        #~ if 'exp' in kwargs:
            #~ p_exp = kwargs['exp']
        #~ else: # if not, set it to 2, i.e. r^2
            #~ p_exp = 2
        radii = r_in + (r_out - r_in)*(linspace(r_in, r_out, num=nshell, endpoint=end)/(r_out))**2
        #pr_int('Not implemented yet.')
        #raise ParError(spaced)
    else:
        raise Exception(space)
    return radii

# FIXME, does not work tries to import old adavis module
def plot_spectrum(freq, intensity, dpc = 0, jy = 0, pstyle = '', xlog = 1, ylog = 1):
    import sys
    import matplotlib.pyplot as pl
    pl.ion()
    #~ from ..views import set_rc
    #~ set_rc
    xcoord = 1.0e4 * _cgs.CC / freq

    if dpc == 0: sys.exit('Error: distance needs to be set when plotting flux')

    distfact = 1.e0/ (dpc**2)

    if jy != 0:
        lumfact = 1e+23
    else:
        lumfact = freq

    pl.plot(xcoord, distfact * lumfact * intensity, pstyle)
    pl.xlabel(r'$\lambda\, [\mu \mathrm{m}]$')

    if jy != 0:
        pl.ylabel(r'$F_\nu$\, [Jy]')
    else:
        pl.ylabel(r'$\nu F_\nu \, [\mathrm{erg cm}^{-2}\, \mathrm{s}^{-1}]$')

    if xlog == 1: pl.xscale('log')
    if ylog == 1: pl.yscale('log')

def write_ratraninput(self):
    # ugly hack, need to restructure the whole radiative transfer module
    from scipy import log, log10, pi, array, arange
    
    
    input_dir_path = _os.path.join(_os.getcwd(), self.directory)
    # if the directory exists
    if not make_dirs(input_dir_path): # will raise error if things go south (i.e., permissions not correct [I think...])
        print('Directory exists, continuing.')
    
    V = 4 * pi * ((self.ra**3  - self.ra**3 )) / 3     # cm3
    self.M = V * (self.nh + self.ne ) * _cgs.MUH2 * _cgs.MP # g = cm3 * g/cm3
    self.M /= _cgs.MSUN                            # Msun
    
    # to get the column density, integrate over radius r1 to r_10k
    #r_10k * 2 nh2_10k
    
    ################################################################
    #  print info. -> move to __str__ method
    #~ print ('M_10K   : {0:<7.2f} Msun\n'
            #~ 'R_10K   : {1:<7.0f} AU\n'
            #~ 'nH2_10K : {2:<7.1e} cm-3\n'
            #~ 'Y       : {3:<7.0f}\n'
            #~ 'T       : {4:<7.1f} K\n'.format(self.M.sum(),
                                        #~ self.r_10k/_cgs.AU,
                                        #~ self.nh2_10k,
                                        #~ self.Y,
                                        #~ self.temp_10k))
    #~ print 'Constraining the envelope : ', self.r_constraint
    #~ print ('nH2_r1000   : {0:<7.1e} cm-3\n'
            #~ 'T_r1000     : {1:7.1f} K\n'.format(self.nh2_r1000,
                                         #~ self.temp_r1000))
    
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
    # the model file!
    with open(_os.path.join(self.directory, self.modelfile),'w') as f:
        f.write('# Ratran input file based on Transphere results'+'\n')
        if self.skyonly: 
            f.write('# ... intended for (SKY) continuum calculations only.'+'\n')
        f.write("rmax={0:.5E}\n".format( self.rb[-1] ))       # rmax in METERS (convert from cm i.e. / 100)
        f.write("ncell={0:}\n".format(len(self.rb)))
        f.write("tcmb=2.735\n")
        f.write("columns=id,ra,rb,nh,nm,ne,tk,td,te,db,vr\n")
        f.write("gas:dust={0}\n".format(self.gas2dust))
        if self.skyonly: 
            f.write("kappa={0}\n".format(self.kappa))
        f.write('@\n')
        # r1/r2 in meter (convert from cm)
        for ii in range(0, len(self.rb)):
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
            f.write(test.format(ii + 1,              #  1 id : shell number
                self.ra[ii],                         #  2 ra : inner radius (m)  
                self.rb[ii],                         #  3 rb : outer radius (m)
                self.nh[ii],                         #  4 nh : density (cm-3) of main coll. partner (usually p-H2)
                self.nm[ii],                         #  5 nm : density (cm-3) of molecule
                self.ne[ii],                         #  6 ne : density (cm-3) of second coll. partner (e.g, e^-, o-H2)
                self.tk[ii],                         #  7 tk : kinetic temperature (K) 
                self.td[ii],                         #  8 td : dust temperature (K)
                self.te[ii],                         #  9 te : second coll. partner temperature (K)
                self.db[ii],                         # 10 db : 1/e half-width of line profile (Doppler b-parameter) (km s-1)
                round(self.vr[ii], 15))              # 11 vr : radial velocity (km s-1)
                    )          
                                
                                
                                
            #~ f.write("%4i %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E" % (ii+1, self.r1[ii]/100.0, self.r2[ii]/100.0, self.nh2int[ii], self.nh2int[ii]*self.abund_int[ii], self.tkint[ii], self.tdint[ii], self.db, self.vr[ii])+'\n')
    # the AMC.inp file
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
            #~ f.write("seed=1971\n")       # NOT GOOD, avoid!
            f.write("fixset={0:3.2E}\n".format(self.fixset))
            f.write("go\n")
            f.write("q\n")
            f.write("\n")
    
    # the SKY.inp file
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

def check_input(input_dictionary, input_defaults):
    return 0

def make_dirs(path):
    import os
    import errno
    # the makedirs function will raise a EEXIST error 
    # if the directory already exists, if so it returns False
    # if some other error is raise, it will raise that 
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        else:
            return False
    return True

class ChangeDirectory:
    """
    # Now you can enter the directory like this:
    #~ import subprocess
    #~ with cd("~/Library"):
        #~ # we are in ~/Library
        #~ run some code
        #~ subprocess.call("ls")
    """
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = _os.getcwd()
        _os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        _os.chdir(self.savedPath)

# temporary function
# needs to be more modular
def plot_envstruct(self, mol_abundance = '', mark100k = True, **kawargs):
    if not hasattr(self, 'Envstruct'):
        raise Exception('you havent read in the transphere output')
    import matplotlib.pyplot as pl
    from matplotlib.ticker import ScalarFormatter, LogFormatter
    pl.ion()
    pl.close()
    fig = pl.figure(1, figsize=(8,6))
    ax1 = fig.add_subplot(111)
    pl.grid()
    ax2 = ax1.twinx()
    # Density
    p1 = ax1.loglog(self.Envstruct.r/_cgs.AU, self.n_h2, label='n_H2', **kawargs)
    ax1.set_xlabel('Radius (AU)')
    #~ ax1.set_xscale('log')
    ax1.set_ylabel('Number Density (cm-3)')
    # Temperature
    p2 = ax2.loglog(self.Envstruct.r/_cgs.AU, self.Envstruct.temp, color='r', label='Temp', **NiceLineSettings)

    ax2.yaxis.set_major_formatter(ScalarFormatter())
    
    ax2.set_ylabel('Temp (K)', color='r')
    if mol_abundance != '':
        def make_patch_spines_invisible(ax):
            ax.set_frame_on(True)
            ax.patch.set_visible(False)
            for sp in ax.spines.itervalues():
                sp.set_visible(False)
        ax3 = ax1.twinx()
        #~ ylims = ax3.get_ylim()
        #ax3.set_ylim(-0.05E-7, 1.85E-7)

        #~ p3 = ax3.loglog(self.Envstruct.r/_cgs.AU, mol_abundance, 'g')
        p3 = ax3.semilogx(self.Envstruct.r/_cgs.AU, mol_abundance, color='g', label='Mol Abund', **NiceLineSettings)
        ax3.spines["right"].set_position(("axes", 1.2))
        make_patch_spines_invisible(ax3)
        ax3.spines["right"].set_visible(True)
        ax3.set_ylabel('Rel. Abundance', color='g')
        #ax1.legend([p1, p2, p3], ['Density', 'Temp', 'Rel. Abundance'])
        #~ ax3.yaxis.set_major_formatter(())
        #~ ax3.xticks(['1E-9','1E-8','1E-7','1E-6','1E-5','1E-4'],[1E-9,1E-8,1E-7,1E-6,1E-5,1E-4])
        #~ ax3.set_yticks([1E-9,1E-8,1E-7,1E-6,1E-5,1E-4], minor=True)
        #~ ax3.tick_params(axis='y', direction='in')
        fig.subplots_adjust(right = 0.75)
    
    if mark100k:
        from scipy import where
        # where is the value closest to 100 K?
        i_100k = where(abs(100 - self.Envstruct.temp).round(2) == round(min(abs(100 - self.Envstruct.temp)), 2))[0][0]
        r_100k = self.Envstruct.r[i_100k]/_cgs.AU
        t_100k = self.Envstruct.temp[i_100k]
        ax2.annotate('T = {0:.1f} K\nR = {1:.1f} AU'.format(round(t_100k,2), r_100k),
                xy=(r_100k, t_100k), xycoords='data',
                xytext=(-30, -100), textcoords='offset points', fontsize=12,
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=-.2"))
        #~ ax2.plot([r_100k], [t_100k] , 'o',color='r', ms=4, mew=0)
        #~ pl.legend('n_H2', 'Temp', 'Mol Abund')
    #~ else:
    ax1.xaxis.set_major_formatter(ScalarFormatter())
    if mol_abundance == '':
        #Create custom artists
        simArtist = pl.Line2D((0,1),(0,0), color='b')
        anyArtist = pl.Line2D((0,1),(0,0), color='r')
        
        #Create legend from custom artist/label lists
        ax1.legend([simArtist,anyArtist],
                  ['Density', 'Temperature'])
    elif mol_abundance != '':
        #Create custom artists
        simArtist = pl.Line2D((0,1),(0,0), color='b')
        anyArtist = pl.Line2D((0,1),(0,0), color='r')
        molArtist = pl.Line2D((0,1),(0,0), color='g')
        
        #Create legend from custom artist/label lists
        ax1.legend([simArtist, anyArtist, molArtist],
                  ['Density', 'Temperature', 'Mol. abundance'])
    
# test functions
def cleanup_ratran():
    import os
    filelist = ['populations.pop', 'amc.inp', 'sky.inp', 'transphere.mdl']
    os.system('rm -Rf _007')
    os.system('rm -Rf image.conv')
    for f in filelist:
        os.system('rm {0}'.format(f))

def cleanup_transphere():
    import os
    filelist = ['convhist.info', 'external_meanint.inp' , 'spectrum.dat', 'transphere.dat', 'dustopac_1.inp',  'envstruct.dat', 'starinfo.inp', 'transphere.inp', 'convhist.dat',  'dustopac.inp', 'envstruct.inp', 'starspectrum.inp']
    for f in filelist:
        os.system('rm {0}'.format(f))

# ratrarun obsolete. saved to be able to check later
"""
def run_ratran(r = 0.0, rho_dust = 0.0, temp = 0.0, db = 0.0, abund = 0.0, vr = 0.0, tdust = 0.0, dustonly = 0, mdl_file = 'transphere.mdl', dpc = 0.0, imsize = 129, pixel = 0.5, trans = '220.0e9', writeonly = 0, skyonly = 0, molfile='', ncell = 50, outputfile="ratranResult", snr=20, fixset=1e-6, minpop=1e-4, unit='Jypx', opstates=0, gas2dust=100, nphot=1000, temp_limit=10.0, rho_limit=1E-4, pxl_radius = 32, los = 2):

    #~ TODO!!!!: create a class, that saves everything,
    #~ envelope structure, model, etc before running RATRAN
    #~ (type : line or continuum)
    #~ TODO : validate the interpolation
    #~ TODO : check the input files for common errors(?)
    #~ TODO : calculate the column density


    from scipy import zeros, array, logspace, log10, pi, where
    import scipy.interpolate
    import sys
    import os
    import numpy as np
    from time import time

    # rewrite this part when changing to object oriented
    # now it is very hack-ish and non pythonic
        
    # Find the envelope cut off for T and n
    #
    # see if T goes below 10K (def) somewhere in the model
    try:
       ind_T = where(temp<temp_limit)[0].min()
    #~ ind = where(r<(1.2E4*cgs.AU))[0].max()
    except (ValueError): 
        ind_T = False
    # see if n goes below 1E-4 (def) somewhere in the model
    try:
        n_h2 = rho_dust * gas2dust  / _cgs.MUH2 / _cgs.MP
        ind_n = where((n_h2) < rho_limit)[0].min()
    except (ValueError):
        ind_n = False
    # T or n strongest constraints on radius
    # Neither T nor n constrain the radius
    # thus just use the last element
    if ind_n == False and ind_T == False:
        r_constraint = None
        ind = len(r)-1
    # Both constraint, which comes first
    elif ind_n != False and ind_T != False:
        # ind_n comes first
        ind = min((ind_n, int_T))
        # what if both have the same...
        # it will pick T, ok
        r_constraint = ['n', 'T'][ind_n < ind_T]
    elif ind_n != False:
        ind = ind_n
        r_constraint = 'n'
    elif ind_T != False:
        ind = ind_T
        r_constraint = 'T'
    
    r_10k = r[ind]
    print ind
    rho_dust_10k = rho_dust[ind]
    nh2_10k = rho_dust_10k * 100 / _cgs.MUH2 / _cgs.MP
    temp_10k = temp[ind]
    Y = r.max() / r.min()

    ind_r1000 = where(r > 1000 * _cgs.AU)[0].min()
    rho_dust_r1000 = rho_dust[ind_r1000]
    nh2_r1000 = rho_dust_r1000 * 100 / _cgs.MUH2 / _cgs.MP
    temp_r1000 = temp[ind_r1000]

    ###### cut off where T<10 K
    # first we have to remove all cells where T<10 K
    # RATRAN does not work well with them
    r = r[:ind]
    rho_dust = rho_dust[:ind]
    temp = temp[:ind]
    abund = abund[:ind]

    if tdust == 0.0:
        tdust = temp
    # from dust density to number density
    #    g/cm-3 * 100 / g
    # ! rho is dust density !
    nh = rho_dust * 100.0 / _cgs.MUH2 / _cgs.MP

    print 'Writing model in {0}'.format(mdl_file)
    #
    # for the refinement, easiest way is to redefine rx!
    #
    rx = logspace(log10(r[0]), log10(r[-1]), num=ncell+1, endpoint=True)
    rx = np.insert(rx, 0, 0)
    r1 = rx[0:-1]
    r2 = rx[1:]
    #~ r1=np.insert(r[0:-1],0,0)
    #~ r2=np.array(r)

    rr = zeros(ncell+1,float)
    rr[1:] = 10**( (np.log10(r1[1:]) + np.log10(r2[1:])) / 2.0 )
    rr[0] = rr[1]
    # Interpolate the values to 'ncell' cells
    nhf = scipy.interpolate.interp1d(log10(r), log10(nh))
    tkf = scipy.interpolate.interp1d(log10(r), log10(temp))
    tdf = scipy.interpolate.interp1d(log10(r), log10(tdust))
    abund_f = scipy.interpolate.interp1d(log10(r), log10(abund))
    #
    # Convert logarithms to floats
    nhint = 10**nhf(log10(rr))
    tkint = 10**tkf(log10(rr))
    tdint = 10**tdf(log10(rr))
    abund_int = 10**abund_f(np.log10(rr))
    ############################
    #~ nhint_p = nhint*2/4.
    #~ nhint_p[0] = 0.0
    #~ nhint_o = nhint*2/4.
    #~ nhint_o[0] = 0.0
    #~ teint = 10**tkf(log10(rr))
    ############################
    nhint[0] = 0.0

    # mass of it all
    #~ vol=[]
    #~ mass=[]
    # V = 4*pi*r**3/3
    # r in cm (?)
    V = 4 * pi * (r2**3-r1**3) / 3         # cm3
    M = V * nhint * _cgs.MUH2 * _cgs.MP # cm3 * g/cm3 ?
    M /= _cgs.MSUN
    
    # to get the column density, integrate over radius r1 to r_10k
    #r_10k * 2 nh2_10k

    print ('M_10K   : {0:<7.2f} Msun\n'
            'R_10K   : {1:<7.0f} AU\n'
            'nH2_10K : {2:<7.1e} cm-3\n'
            'Y       : {3:<7.0f}\n'
            'T       : {4:<7.1f} K\n'.format(M.sum(),
                                        r_10k/_cgs.AU,
                                        nh2_10k,
                                        Y,
                                        temp_10k))
    print 'Constraining the envelope : ', r_constraint
    print ('nH2_r1000   : {0:<7.1e} cm-3\n'
            'T_r1000     : {1:7.1f} K\n'.format(nh2_r1000,
                                         temp_r1000))
    #~ raise Exception('Test')
    #~ return r2,r1
    #
    # INPUT MODEL
    #
    f = open(mdl_file,'w')
    f.write('# Ratran input file based on Transphere results'+'\n')
    if skyonly == 1:
        f.write('# intended for (SKY) continuum calculations only.'+'\n')
    f.write('rmax={0:-5.3e}\n'.format( max(r2) / 1.0e2 ))
    f.write('ncell=%i'  % (len(r2))+'\n')
    f.write('tcmb=2.728\n')
    #~ f.write('columns=id,ra,rb,nh,nm,tk,td,db,vr\n')
    f.write('columns=id,ra,rb,nh,nm,tk,td,db,vr\n')
    f.write('gas:dust={0}\n'.format(gas2dust))
    if skyonly == 1:
        f.write('kappa=jena,thin,e6\n')
    f.write('@\n')
    for ii in range(0,len(r1)):
        f.write("%4i %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e" % (ii+1, r1[ii]/100.0, r2[ii]/100.0, nhint[ii], nhint[ii]*abund_int[ii], tkint[ii], tdint[ii], db, vr)+'\n')
    #~ for ii in range(0,len(r1)):
        #~ f.write("%4i %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E" % (ii+1, r1[ii]/100.0, r2[ii]/100.0, nhint_p[ii], nhint[ii]*abund_int[ii], nhint_o[ii], tkint[ii], tdint[ii], teint[ii], db, vr)+'\n')
    f.close()
    #
    # AMC input file
    #
    if skyonly == 0:
        if molfile == '': sys.exit('Error: for AMC calculations the molecular datafile (molfile) needs to be set.')
        f = open("amc.inp",'w')
        f.write("source="+mdl_file+'\n')
        f.write("outfile=populations.pop\n")
        f.write("molfile="+molfile+'\n')
        f.write("snr={0}\n".format(snr))
        f.write("nphot={0}\n".format(nphot))
        f.write("kappa=jena,thin,e6\n")
        f.write("minpop={0}\n".format(minpop))
        f.write("seed=1971\n")
        f.write("fixset={0}\n".format(fixset))
        f.write("go\n")
        f.write("q\n")
        f.write(" \n")
        f.close()
    #
    # SKY input file
    #
    f = open("sky.inp",'w')
    if skyonly == 0:
        f.write("source=populations.pop\n")
    else:
        f.write("source="+mdl_file+"\n")
    f.write("format=miriad\n")
    f.write("outfile="+outputfile+"\n")
    f.write("trans="+trans+"\n")
    #~ f.write("pix="+str(imsize)+","+str(pixel)+",32,8\n")
    f.write("pix={0},{1},{2},{3}\n".format(imsize, pixel, pxl_radius, los))
    if skyonly == 0:
        f.write("chan=100,0.2\n")
    else:
        f.write("chan=1,1.0\n")
    f.write("distance="+str(dpc)+"\n")
    f.write("units={0}\n".format(unit))
    f.write("go\n")
    f.write("q\n")
    f.close()
    #
    #
    #
    if writeonly == 0:
        if skyonly == 0:
            #~ print "Starting AMC calculation..."
            t1 = time()
            os.system('amc amc.inp')
            print 'AMC took :  {0:2.2f} hours'.format((time()-t1)/3600.0)
            os.system('alert \"AMC has finished running.\"')
        #~ print "Starting SKY calculation..."
        t1 = time()
        os.system("sky sky.inp")
        print ' SKY took : {0:2.2f} seconds'.format(time()-t1)
        os.system('alert \"SKY has finished running.\"')
"""

def save_ratran(Obj, filename = 'ratranmodel.pickle'):
    # take input object
    # and save it as a dictionary to filename in directory
    import pickle
    # take all the original input parameters and put into a dictionary
    inp = vars(Obj.Input)
    # now a dictionary is a non dynamic structure
    # as opposed to a dynamically created object attribute
    # e.g., with the __dict__ method (-> doesn't work with pickle)
    with ChangeDirectory(Obj.directory):
        with open(filename, 'w') as f:
            pickle.dump(inp, f)

def load_ratran(directory = '', filename = 'ratranmodel.pickle'):
    # load input object from filename in directory
    # and create the ratran object
    import pickle
    with open(_os.path.join(directory, filename), 'r') as f:
        inputdict = pickle.load(f)
    Obj = Ratran(**inputdict)
    # IDEA : add so that it loads the output(?) as well?
    return Obj

def save_transphere(Obj, filename = 'transpheremodel.pickle'):
    # take input object
    # and save it as a dictionary to filename in directory
    import pickle
    # take all the original input parameters and put into a dictionary
    inp = vars(Obj.Input)
    # now, a dictionary is a non dynamic structure
    # as opposed to a dynamically created object attribute
    # e.g., with the __dict__ method (-> doesn't work with pickle)
    with ChangeDirectory(Obj.directory):
        with open(filename, 'w') as f:
            pickle.dump(inp, f)

def load_transphere(directory = '', filename = 'transpheremodel.pickle'):
    # load input object from filename in directory
    # and create the transphere object
    import pickle
    with open(_os.path.join(directory, filename), 'r') as f:
        inputdict = pickle.load(f)
    Obj = Transphere(**inputdict)
    # IDEA : add so that it loads the output(?) as well?
    return Obj

def temp_pop(nl, nu, gl, gu, eu):
    """ 
    Calculate the excitation temperature of transition from 2 to 1
    """
    
    numer = eu
    denom = _scipy.log(nl * gu / (nu * gl))
    return numer / denom

def calc_nu(N, Q, gu, Eu, T):
    """
    Calculate the column density of the upper energy level
    assuming LTE 
    """
    part1 = N / Q * gu
    part1 = _scipy.exp(-Eu / (_cgs.KK * T))
    return part1 * part2

def kappa_powerlaw( nu0, kappa0, beta ):
    #~ kappa=powerlaw,NU0,KAPPA0,BETA
    #~ where NU0 is in Hz, KAPPA0 in cm2/g_dust, and BETA is freq.index.
    #~ ...in m2/kg_dust: (0.1 converts cm2/g to m2/kg)
    #~ kappa=0.1* kappa0*(nu/nu0)**beta
    #~ table = _scipy.loadtxt(filepath).transpose()
    kappa = lambda nu : kappa0 * (nu/nu0)**beta
    from scipy import arange
    lamtab = _scipy.log10( arange(1E-4,10,1e-3) ) # to cgs unit (cm)
    kaptab = _scipy.log10( kappa(_cgs.CC / arange(1E-4,10,1e-3)) ) # input Hz
    
    from scipy import interpolate
    
    tck = interpolate.splrep(lamtab, kaptab, k=1, s=0)
        
    return tck

def get_kappa_jena(filepath):

    table = _scipy.loadtxt(filepath).transpose()
    
    lamtab = _scipy.log10( table[0] * 1.0E-4 ) # to cgs unit (cm)
    kaptab = _scipy.log10( table[1] )
    
    from scipy import interpolate
    
    tck = interpolate.splrep(lamtab, kaptab, k=1, s=0)
    
    #~ def function(nu):
        #~ return interpolate.splev(nu, tck)
    
    #~ tck = scipy.interpolate.splrep(X, Y, k=1, s=0)
    #~ scipy.interpolate.splev(6, tck)
    
    # function should check if within file table limit
    
    return _scipy.array([lamtab, kaptab]), tck

def calc_dtau(Nu, Nl, Aul, freq, gu, gl):
    """
    Calculate the optical depth (tau) for given parameters
    from Goldsmith and Langer (1999) or not?
    """
    
    #~ alpha = _cgs.HH * nu / (4 * _scipiy.pi)
    dtau = _cgs.CC**2 / (8 * _scipy.pi * freq**2) * Aul * (Nl * gu / gl - Nu) 
    #~ dtau = alpha * ds
    
    
    #~ part1 = _cgs.CC**3 * Aul / (8 * _scipy.pi * freq**3 * width) 
    #~ part2 = (_scipy.exp(Eu/T) - 1)
    #~ part2 = (Nl * gu / float(gl) - Nu)
    #~ return part1 * part2
    return dtau

def calc_dtau_lte(Nu, Aul, freq, width, Eu, T):
    """
    Calculate the optical depth (tau) for given parameters
    from Goldsmith and Langer (1999) or not?
    """
    part1 = Nu * _cgs.CC**3 * Aul / (8 * _scipy.pi * freq**3 * width) 
    part2 = (_scipy.exp(Eu/T) - 1)
    return part1 * part2

def calc_radiation(nu, nl, freq, gu, gl):

    part1 = 2 * _cgs.HH * freq**3 / _cgs.CC**2
    part2 = nl * gu / (nu * gl)
    sv = part1 * 1 / (part2 - 1)
    return sv

def get_transindex(Obj, trans):
    try:
        transIndex = [i['trans'] for i in Obj.Moldata.radtrans if i['up'] == trans[0] and i['down'] == trans[1]][0]
    except (IndexError):
        print('No such transition for molecule : {0}'.format(self.Moldata.molecule))
        print('Check your moldata file at : {0}'.format(self.Moldata.molfile))
        return False
    transIndex -= 1
    return transIndex

# END OLD HELP FUNCTIONS
########################################################################
