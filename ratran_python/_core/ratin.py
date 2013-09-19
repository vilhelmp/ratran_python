import os as _os
from .. import cgsconst as _cgs

class Make(object):pass

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
    Mdl.ne[0] = 1.0
    Mdl.nm_rel = Mdl.nm / (Mdl.nh + Mdl.ne)
    Mdl.nh[0] = 0.0
    Mdl.ne[0] = 0.0
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
    pl.close()
    fig = pl.figure(1)#num = 1, figsize = ())
    plots = dict()
    for (i, dat) in zip(arange(N), Ratran_mdl.table[2:]):
        ax = 'ax{0}'.format(i)
        ylbl = Ratran_mdl.columns[i + 2]
        # now how many subplots do we need?
        if N % 3:
            np = N / 3 + 1
        else:
            np = N / 3
        pln =  i +1
        if pln >1:
            plots[ax] = fig.add_subplot(np, 3, pln, sharex=plots['ax0'])
        else:
            plots[ax] = fig.add_subplot(np, 3, pln)
        x = Ratran_mdl.ra*100/_cgs.AU
        if ylbl in ['db', 'vr']:
            plots[ax].semilogx(x, dat, '.')
            plots[ax].semilogx(x, dat, '-')
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
            ylbl += ' (m)'
        plots[ax].set_ylabel(ylbl)
        plots[ax].grid()
    [plots['ax{0}'.format(i)].set_xlabel('ra (AU)') for i in arange(N-3, N)]
    #~ [plots['ax{0}'.format(i)].xaxis.set_visible(0) for i in arange(N-3)]
    #~ [plots['ax{0}'.format(i)].xaxis.set_ticklabels([]) for i in arange(N-3)]


    #~ plots['ax{0}'.format(N-1)].set_xlabel('ra (AU)')
    fig.subplots_adjust(left=0.11, right= 0.97, bottom=0.11, top=0.96, wspace=0.43, hspace=0.15)
    return plots, fig
