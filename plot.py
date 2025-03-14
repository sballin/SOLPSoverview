#!/bin/env python3
'''
Create a PDF summarizing a run.

Usage: in run directory, do `plot.py` to generate a file called plots.pdf.
'''
import os
import copy
import subprocess
import re
import datetime
from heapq import nsmallest
import traceback
import glob
import shutil
import textwrap
import gzip
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erfc
from scipy.io import netcdf
import scipy.interpolate as interpolate
from labellines import labelLines
import netCDF4 as nc
import f90nml
import warnings
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import collections as mc
from matplotlib.backends.backend_pdf import PdfPages
import aurora # suppress cryptography deprecation warning with pip install cryptography==2.5
from aurora import coords
import quixote
import solpspy

np.seterr(divide='ignore', invalid='ignore') # suppress RuntimeWarning: invalid value encountered in true_divide
# suppress f90nml warning
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
warnings.filterwarnings("ignore", message=".*Value .* is not assigned to any variable and has been removed.*")
parser = f90nml.Parser() # It works, but it doesn't like things like userfluxparm(1,1) = bunch of numbers. It wants userfluxparm(:,1) = bunch of numbers


def get_cmd_output(arr):
    proc = subprocess.Popen(arr, stdout=subprocess.PIPE)
    return proc.communicate()[0].decode('utf-8')


def argnear(array, value):
    return np.argmin(np.abs(np.array(array)-value))


def parseDat(filename):
    out = {}
    with open(filename, 'r') as f:
        txt = f.read()
    lines = txt.split('\n')
    for l in lines:
        try:
            if l[0] == "'":
                l = l.split('#')[0]
                key, val = l.split() 
                key = key.replace("'", "")
                val = float(val.replace("'", ""))
                out[key] = val
        except Exception as e:
            pass
    return out


def showPoloidalLoc(ix):
    """
    Get polygons describing B2 grid at poloidal location ix.
    Adapted from aurora method get_b2_patches.
    """
    xx = so.data("crx").transpose(2, 1, 0)
    yy = so.data("cry").transpose(2, 1, 0)
    NY = int(so.data("ny"))
    NX = int(so.data("nx"))

    if so.form == "files":
        # eliminate boundary cells
        xx = xx[1:-1, 1:-1, :]
        yy = yy[1:-1, 1:-1, :]

    patches = []
    for iy in np.arange(0, NY):
        rr = np.atleast_2d(xx[ix, iy, [0, 1, 3, 2]]).T
        zz = np.atleast_2d(yy[ix, iy, [0, 1, 3, 2]]).T
        patches.append(matplotlib.patches.Polygon(np.hstack((rr, zz)), True, linewidth=3))

    # collect all patches
    p = matplotlib.collections.PatchCollection(patches, False, fc="w", edgecolor="red", linewidth=0.5)
    plt.gca().add_collection(p)


def rmid(so, geqdsk):
    '''Modified from Aurora.'''
    # # find rhop along midplane grid chords
    # rhop_chord_HFS = coords.get_rhop_RZ(so.data("cr")[:, JXI], so.data("cz")[:, JXI], so.geqdsk)
    rhop_chord_LFS = coords.get_rhop_RZ(so.data("cr")[:, JXA], so.data("cz")[:, JXA], geqdsk)
    r = coords.rad_coord_transform(rhop_chord_LFS, 'rhop', 'Rmid', geqdsk) # Without deepcopy you get error: 'solps_case' object has no attribute 'geqsk'
    rsep = (r[so.ny//2-1]+r[so.ny//2])/2
    return r-rsep


def plot_wall_geometry(so, ax, mark=None):
    """Method to plot vessel wall segment geometry from wall_geometry field in fort.44 file"""

    out = so.load_fort44()
    wall_geometry = out["wall_geometry"]

    Wall_Seg = []
    RR = wall_geometry[0::2]
    ZZ = wall_geometry[1::2]
    NLIM = out["NLIM"]

    for i in range(0, NLIM):
        line = [(RR[2 * i], ZZ[2 * i]), (RR[2 * i + 1], ZZ[2 * i + 1])]
        Wall_Seg.append(line)
    Wall_Collection = mc.LineCollection(Wall_Seg, colors="k", linewidth=0.5, zorder=0)
    ax.add_collection(Wall_Collection)
    
    # Show a user-specified element
    if mark != None:
        Wall_Seg2 = []
        for i in [mark]:
            line = [(RR[2 * i], ZZ[2 * i]), (RR[2 * i + 1], ZZ[2 * i + 1])]
            Wall_Seg2.append(line)
        Wall_Collection2 = mc.LineCollection(Wall_Seg2, colors="r", linewidth=0.5, zorder=1)
        ax.add_collection(Wall_Collection2)
    
    ax.set_xlim(RR.min() - 0.05, RR.max() + 0.05)
    ax.set_ylim(ZZ.min() - 0.05, ZZ.max() + 0.05)
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.set_aspect("equal")

    so.WS = Wall_Seg
    so.WC = Wall_Collection
    
    
def plotmesh():
    plt.triplot(so.xnodes, so.ynodes, so.triangles, '-', color='gray', lw=0.2, zorder=0)
    p = so.get_b2_patches()
    p.set_edgecolor('C0')
    p.set_facecolor('white')
    p.set_zorder(1)
    plt.gca().add_collection(p)


def getTransportCoeffs():
    '''Source: https://github.com/ORNL-Fusion/SOLPSxport/blob/87107473c9ecb7c002082e4dfc8b199412af6572/SOLPSutils.py#L880'''
    dn = np.zeros((sp.ny,sp.ns))
    dp = np.zeros((sp.ny,sp.ns))
    chii = np.zeros((sp.ny,sp.ns))
    chie = np.zeros(sp.ny)
    vlax = np.zeros((sp.ny,sp.ns))
    vlay = np.zeros((sp.ny,sp.ns))
    vsa = np.zeros((sp.ny,sp.ns))
    sig = np.zeros(sp.ny)
    alf = np.zeros(sp.ny)    

    nml = parser.read('b2.transport.parameters')
    this = nml['transport']['parm_dna']
    for ispec in range(len(this)):
        dn[:,ispec] = this[ispec]

    dsa = np.loadtxt('dsa')
    nml = parser.read('b2.transport.inputfile')
    tdata = nml['transport']['tdata']
    nS = len(tdata)
    for ispec in range(nS):
        nKinds = len(tdata[ispec])
        for jkind in range(nKinds):
            this = tdata[ispec][jkind]

            # Check if this kind was filled with none by f90nml (not defined)
            test = [i[1] for i in this]
            if all(val is None for val in test):
                continue

            # Read and interpolate back on dsa
            xRaw = [i[0] for i in this] 
            yRaw = [i[1] for i in this]
            yInterp = np.interp(dsa,xRaw,yRaw)

            if jkind+1 == 1:
                dn[:,ispec] = yInterp
            elif jkind+1 == 2:
                dp[:,ispec] = yInterp
            elif jkind+1 == 3:
                chii[:,ispec] = yInterp
            elif jkind+1 == 4:
                chie[:] = yInterp
            elif jkind+1 == 5:
                vlax[:,ispec] = yInterp                        
            elif jkind+1 == 6:
                vlay[:,ispec] = yInterp                        
            elif jkind+1 == 7:
                vsa[:,ispec] = yInterp                        
            elif jkind+1 == 8:
                sig[:] = yInterp                        
            elif jkind+1 == 9:
                alf[:] = yInterp                        
    return {'dn': dn, 'dp': dp, 'chii': chii, 'chie': chie, 'vlax': vlax, 'vlay': vlay, 'vsa': vsa, 'sig': sig, 'alf': alf}


def plotTransportCoeffs(ax):
    dsa = np.loadtxt('dsa')
    coeffs = getTransportCoeffs()
    ax.plot(dsa*1000, coeffs['dn'][:,0], marker='.', label='$D$')
    ax.plot(dsa*1000, coeffs['chii'][:,0], marker='.', label='$\chi_i$')
    ax.plot(dsa*1000, coeffs['chie'], marker='.', label='$\chi_e$')
    ax.set_xlabel(yyclabel)
    ax.set_ylabel('[m$^2$/s]')
    ax.set_yscale('log')
    ax.legend()


def getvals(file, key, valtype):
    '''e.g. 
    getvals('b2.boundary.parameters', 'bccon(0,1)', int) 
    returns
    [2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8] in a case with D+Ne
    '''
    with open(file, 'r') as f:
        for l in f:
            if key in l:
                valstring = l.split('=')[1].replace('\n','').replace(' ','').replace("'",'')
                return [valtype(i) for i in valstring.split(',') if i != '']


def get_b2tallies(var):
    '''
    returns (times, vals)
    fixed up from solpspy: https://gitlab.com/ipar/solpspy/-/blob/master/solpspy/classes/rundir.py#L1684'''
    exeb2tallies = os.path.join('b2mn.exe.dir/b2tallies.nc')
    b2tallies = os.path.join('b2tallies.nc')

    if os.access(exeb2tallies, os.R_OK):
        ftallies=netcdf.netcdf_file(exeb2tallies,'r')
    elif os.access(b2tallies, os.R_OK):
        ftallies=netcdf.netcdf_file(b2tallies,'r')

    times = ftallies.variables['times'][:].copy()
    data =  ftallies.variables[var][:].copy()
    ftallies.close()
    return times, data


# divertor pressure
# Get density right in front of semi-transparent surface
def get_pumped_fluxes():
    i = 0

    D_pumped_pump = []
    i_D_pumped_pump = []
    T_pumped_pump = []
    i_T_pumped_pump = []

    D_pumped_core = []
    i_D_pumped_core = []
    T_pumped_core = []
    i_T_pumped_core = []

    influx = []
    i_influx = []

    ipumped = 0
    iinflux = 0

    Dpumpiter = 0
    Tpumpiter = 0
    Dcoreiter = 0
    Tcoreiter = 0
    
    if os.path.isfile('run.log.gz'):
        f = gzip.open('run.log.gz','rt') # if in working directory right after run
    elif os.path.isfile('run.log'):
        f = open('run.log') # if opening [date].tar.gz save
    reading_pumped = False
    for l in f:
        # Empty line signals end of flux readout
        if reading_pumped and len(l.split()) == 0:
            reading_pumped = False
            D_pumped_pump.append(Dpumpiter)
            i_D_pumped_pump.append(ipumped)
            T_pumped_pump.append(Tpumpiter)
            i_T_pumped_pump.append(ipumped)

            D_pumped_core.append(Dcoreiter)
            i_D_pumped_core.append(ipumped)
            T_pumped_core.append(Tcoreiter)
            i_T_pumped_core.append(ipumped)
            ipumped += 1
            continue
        
        if reading_pumped:
            surf, species, flux = l.split()
            if '10' in surf and 'D' in species:
                if '2' in species:
                    Dpumpiter += 2*float(flux)/1.602e-19
                else:
                    Dpumpiter += float(flux)/1.602e-19
            if '10' in surf and 'T' in species:
                if '2' in species:
                    Tpumpiter += 2*float(flux)/1.602e-19
                else:
                    Tpumpiter += float(flux)/1.602e-19
            if '-1' in surf and 'D' in species:
                if '2' in species:
                    Dcoreiter += 2*float(flux)/1.602e-19
                else:
                    Dcoreiter += float(flux)/1.602e-19
            if '-1' in surf and 'T' in species:
                if '2' in species:
                    Tcoreiter += 2*float(flux)/1.602e-19
                else:
                    Tcoreiter += float(flux)/1.602e-19
            continue
                    
        if 'NO. SPECIES  PUMPED FLUX' in l:
            reading_pumped = True
            Dpumpiter = 0
            Tpumpiter = 0
            Dcoreiter = 0
            Tcoreiter = 0
            continue
        
        if 'INFLUX=' in l:
            _, flux = l.split()
            influx.append(float(flux)/1.602e-19)
            i_influx.append(iinflux)
            iinflux += 1
            continue
    f.close()
            
    D_pumped_pump = np.array(D_pumped_pump)
    i_D_pumped_pump = np.array(i_D_pumped_pump)
    T_pumped_pump = np.array(T_pumped_pump)
    i_T_pumped_pump = np.array(i_T_pumped_pump)
    D_pumped_core = np.array(D_pumped_core)
    i_D_pumped_core = np.array(i_D_pumped_core)
    T_pumped_core = np.array(T_pumped_core)
    i_T_pumped_core = np.array(i_T_pumped_core)
    return np.mean(D_pumped_pump), np.mean(T_pumped_pump)


def get_triangle_centers():
    '''so.triangles contains indices of nodes of each triangle'''
    xc = []
    yc = []
    for t in so.triangles:
        xc.append(np.mean([so.xnodes[int(t[i])] for i in [0,1,2]]))
        yc.append(np.mean([so.ynodes[int(t[i])] for i in [0,1,2]]))
    xc = np.array(xc)
    yc = np.array(yc)
    return xc, yc


def getConfigText():
    txt = r'$\bf{Path}$ ' + os.getcwd()
    # Date created
    date = datetime.datetime.now().strftime('%I:%M %p %a %d %b %Y')
    txt += '\n' + r'$\bf{Plots\ created}$ ' + date
    solpstop = os.environ['SOLPSTOP']
    txt += '\n' + r'$\bf{SOLPSTOP}$ ' + solpstop
    # SOLPS version
    proc = subprocess.Popen(['git','describe'], stdout=subprocess.PIPE)
    out = proc.communicate()[0].decode('utf-8').replace('\n', '')
    txt += '\n' + r'$\bf{SOLPS\ version}$ ' + out

    b2mn = parseDat('b2mn.dat')    
    bdry = parser.read('b2.boundary.parameters')['boundary']
    trans = parser.read('b2.transport.parameters')['transport']

    # Timesteps
    timesteps = int(b2mn["b2mndr_ntim"])
    resids = getResiduals()
    iterations = np.shape(resids)[0]
    try:
        txt += '\n' + r'$\bf{Timesteps}$ ' + f'{iterations}/{timesteps} at {b2mn["b2mndr_dtim"]:.2e} s'
    except Exception as e:
        print(traceback.format_exc()) 
        txt += '\n' + r'$\bf{Timesteps}$ ' + f"{type(e).__name__}: {e}"

    # Species
    try:
        txt += '\n' + r'$\bf{Species}$ ' + ', '.join(species_nice) + f' ({len(species_nice)})'
    except Exception as e:
        print(traceback.format_exc()) 
        txt += '\n' + r'$\bf{Species}$ ' + f"{type(e).__name__}: {e}"
        
    # Core power
    try:
        txt += '\n' + r'$\bf{Core\ power}$ ' + f'$P_e$={bdry["enepar"][0][0]/1e6:.3g} MW, $P_i$={bdry["enipar"][0][0]/1e6:.3g} MW'
    except Exception as e:
        print(traceback.format_exc()) 
        txt += '\n' + r'$\bf{Core\ power}$ ' + f"{type(e).__name__}: {e}"

    # Core density
    try:
        bccons = getvals('b2.boundary.parameters', 'bccon(0,1)', int)
        bccondict = {1: 'density', 2: 'gradient', 8: 'flux', 27: 'feedback flux'}
        unitdict = {1: 'm$^{-3}$', 2: 'm$^{-4}$', 8: '/s', 27: '/s'}
        conpars = getvals('b2.boundary.parameters', 'conpar(0,1,1)', float)
        bcs = []
        for i, n in enumerate(nuclei):
            if n == 'H':
                spec = 'H'
                if sp.am[i] == 2:
                    spec = 'D'
                elif sp.am[i] == 3:
                    spec = 'T'
                if sp.za[i] == 1:
                    spec += '+'
                if bccons[i] in bccondict.keys():
                    bctype = bccondict[bccons[i]]
                else:
                    bctype = f'bccon={bccons[i]}'
                if bccons[i] in bccondict.keys():
                    units = unitdict[bccons[i]]
                else:
                    units = ''
                bcs.append(f'{spec} {bctype} {conpars[i]} {units}')
        txt += '\n' + r'$\bf{Core\ BCs}$ ' + ', '.join(bcs)
    except Exception as e:
        print(traceback.format_exc()) 
        txt += '\n' + r'$\bf{Core\ BCs}$ ' + f"{type(e).__name__}: {e}"

    # Puffs
    try:
        strat_specs = getvals('b2.neutrals.parameters', 'species_start', int)
        crcstra = getvals('b2.neutrals.parameters', 'crcstra', str)
        fluxes = getvals('b2.neutrals.parameters', 'userfluxparm(1,1)', float)
        puffs = []
        for i in range(len(strat_specs)):
            if crcstra[i] == 'C': # puff
                puffs.append(f'{nuclei_nice[strat_specs[i]]} {fluxes[i]} /s')
        txt += '\n' + r'$\bf{Puffs}$ ' + ', '.join(puffs)
    except Exception as e:
        print(traceback.format_exc()) 
        txt += '\n' + r'$\bf{Puffs}$ ' + f"{type(e).__name__}: {e}"

    # Diffusivity
    try:
        txt += '\n' + r'$\bf{Anomalous\ ion\ diffusivity}$ ' + f'{trans["parm_dna"][1]:.3g} ' + r'm$^2$/s' # 1 is index of ions
    except Exception as e:
        print(traceback.format_exc()) 
        txt += '\n' + r'$\bf{Anomalous\ ion\ diffusivity}$ ' + f'{trans["parm_dna"][1]:.3g} ' + r'm$^2$/s' # 1 is index of ions

    # Pinch velocity
    try:
        txt += '\n' + r'$\bf{Anomalous\ ion\ pinch\ velocity}$ ' + f'{trans["parm_vla"][1]:.3g} ' + r'm/s' # 1 is index of ions
    except Exception as e:
        print(traceback.format_exc()) 
        txt += '\n' + r'$\bf{Anomalous\ ion\ pinch\ velocity}$ ' + f"{type(e).__name__}: {e}"

    # Thermal diffusivity
    try:
        txt += '\n' + r'$\bf{Anomalous\ thermal\ diffusivity}$ ' + f'$\chi_e$={trans["parm_hce"]:.3g} m$^2$/s, $\chi_i$={trans["parm_hci"][1]:.3g} m$^2$/s' # 1 is index of ions
    except Exception as e:
        print(traceback.format_exc()) 
        txt += '\n' + r'$\bf{Anomalous\ thermal\ diffusivity}$ ' + f"{type(e).__name__}: {e}"

    # Viscosity
    try:
        txt += '\n' + r'$\bf{Anomalous\ viscosity}$ ' + f'{trans["parm_vsa"][1]:.3g} ' # 1 is index of ions
    except Exception as e:
        print(traceback.format_exc()) 
        txt += '\n' + r'$\bf{Anomalous\ viscosity}$ ' + f"{type(e).__name__}: {e}"
    txt += '\n'
    
    try:
        try:
            txt += '\n' + r'$\bf{Impurity\ fixed\ fraction}$ ' + str(b2mn["b2sqel_artificial_radiation"])
        except Exception as e: 
            txt += '\n' + r'$\bf{Impurity\ fixed\ fraction}$ ' + str(b2mn["art_rad"])
    except Exception as e:
        print(traceback.format_exc()) 

    try:
        nlim = so.fort44['NLIM']
        pumpabs = np.max(so.fort44['wlabsrp(A)'][0,:nlim])
        pumpa = np.sum(so.fort44['wlpump(A)'][0,:nlim])
        pumpm = np.sum(so.fort44['wlpump(M)'][0,:nlim])
        txt += '\n' + r'$\bf{Pump\ surface\ absorption}$ ' + f'{pumpabs*100:.3g}%'
        txt += '\n' + r'$\bf{Pumping\ rate}$ ' + f'{pumpa:.3g} atoms/s, {pumpm:.3g} molecules/s'
    except Exception as e:
        print(traceback.format_exc())
        txt += '\n' + r'$\bf{Pumping\ rate}$ ' + f"{type(e).__name__}: {e}"

    if has_impurities:
        try:
            # Pumping S not really valid for full impurity model cases because of pdena[ntriang,natm=2 (D+Imp)]
            D_pumped_flux, T_pumped_flux = get_pumped_fluxes()
            xc, yc = get_triangle_centers()
            mask = (1.78<xc)&(xc<1.835)&(-1.405<yc)&(yc<-1.385)
            vals = np.sum(sp.pdena[:, :], axis=1)+2*np.sum(sp.pdenm[:, :], axis=1)
            ngdiv = np.mean(vals.flatten()[mask])
            pump_surface_area = 0.719 # m^-2
            S_pump = (T_pumped_flux+D_pumped_flux)/pump_surface_area/ngdiv
            txt += '\n' + f'$\\bf{{S_{{pump}}}}$ {S_pump:.3g}'

            # Compression for full impurity model cases assuming D+Imp only
            pumpMask = (1.78<xc)&(xc<1.9)&(-1.37<yc)&(yc<-1.33)
            nH0 = sp.pdena[:, 0]+2*sp.pdenm[:,0]
            nImp0 = sp.pdena[:, 1]
            nH0Pump = np.mean(nH0.flatten()[pumpMask])
            nImp0Pump = np.mean(nImp0.flatten()[pumpMask])
            nHSep = float(interpolate.interp1d(yyc, sp.na[JXA+1,1:-1,1], kind='linear')(0))
            nImpSep = float(interpolate.interp1d(yyc, sp.na[:,:,3:].sum(axis=2)[JXA+1,1:-1], kind='linear')(0))
            txt += '\n' + f'$\\bf{{H\\ compression}}$ {nH0Pump/nHSep:.3g} (n0_div={nH0Pump:.3g}/ni_sep,omp={nHSep:.3g})'
            txt += '\n' + f'$\\bf{{{n}\\ compression}}$ {nImp0Pump/nImpSep:.3g} (n0_div={nImp0Pump:.3g}/ni_sep,omp={nImpSep:.3g})'
        except Exception as e:
            print('full impurity model compression calculation error:', traceback.format_exc())
            ngdiv = np.nan
            S_pump = np.nan
            D_pumped_flux = np.nan
            T_pumped_flux = np.nan
        
    try:
        t, f = get_b2tallies('fnayreg')
            
        for n in nucset:
            n_ion_indices = np.where((np.array(nuclei)==n)*(sp.za >= 1))[0]
            end_flux_core = np.sum(f[-1,n_ion_indices,2])
            start_flux_core = np.sum(f[0,n_ion_indices,2])
            end_flux_sep = np.sum(f[-1,n_ion_indices,4])
            start_flux_sep = np.sum(f[0,n_ion_indices,4])
            end_flux_mc = np.sum(f[-1,n_ion_indices,6])
            start_flux_mc = np.sum(f[0,n_ion_indices,6])
            txt += '\n' + f'$\\bf{{{n}\\ ion\\ flux}}$ core {end_flux_core:.3e} /s ({100*(end_flux_core/start_flux_core-1):+.2g}%), sep {end_flux_sep:.3e} /s ({100*(end_flux_sep/start_flux_sep-1):+.2g}%), mc {end_flux_mc:.3e} /s ({100*(end_flux_mc/start_flux_mc-1):+.2g}%)'
    except Exception as e:
        print(traceback.format_exc())      
        txt += '\n' + r'$\bf{Particle\ fluxes}$ ' + f"{type(e).__name__}: {e}"


    if has_impurities:
        try:
            imp_ion_indices = np.where((np.array(nuclei)==impurity)*(sp.za >= 1))[0]
            imp_dens = np.sum(sp.na[:,:,imp_ion_indices], axis=2)
            ecount_sol = np.sum((sp.ne*sp.vol)[sp.masks.sol])
            impcount_sol = np.sum((imp_dens*sp.vol)[sp.masks.sol])
            ecount_core = np.sum((sp.ne*sp.vol)[sp.masks.core])
            impcount_core = np.sum((imp_dens*sp.vol)[sp.masks.core])
            imp_percent_core = 100*impcount_core/ecount_core
            imp_percent_sol = 100*impcount_sol/ecount_sol
            txt += '\n' + f'$\\bf{{{impurity}\\ fraction}}$ core {imp_percent_core:.3g}% , SOL {imp_percent_sol:.3g}%'
        except Exception as e:
            print(traceback.format_exc())      
            txt += '\n' + r'$\bf{Impurity\ fraction}$ ' + f"{type(e).__name__}: {e}"

    try:
        Psep = np.sum((sp.fhey+sp.fhiy)[sp.ixp:sp.oxp+1,sp.sep])
        Prad = sp.line_radiation.sum(axis=2)[sp.masks.core].sum()
        txt += '\n' + r'$\bf{Power\ crossing\ separatrix}$ ' + f'{Psep/1e6:.3g} MW, {Prad/1e6:.3g} MW radiated in core'
    except Exception as e:
        print(traceback.format_exc()) 
        txt += '\n' + r'$\bf{Power\ crossing\ separatrix}$ ' + f"{type(e).__name__}: {e}"

    try:
        Pin = np.sum((sp.fhey+sp.fhiy)[sp.ixp:sp.oxp+1,1])
        fntx = sp.fnt[:,:,0]
        fnty = sp.fnt[:,:,1]
        myfhtx = sp.fhex + sp.fhix + np.sum(sp.fhpx, 2) + np.sum(sp.fhmx, 2) + sp.fhjx #+ fntx
        myfhty = sp.fhey + sp.fhiy + np.sum(sp.fhpy, 2) + np.sum(sp.fhmy, 2) + sp.fhjy #+ fnty
        Pinner = np.sum(-myfhtx[1,1:-1])
        Pouter = np.sum(myfhtx[-1,1:-1])
        Pwalls = np.sum(myfhty[1:-1,-1]) + -np.sum(myfhty[sp.oxp+1:,1]) + -np.sum(myfhty[:sp.ixp,1])
        Prad = sp.line_radiation.sum() 
        Pout = Pinner + Pouter + Pwalls + Prad

        txt += '\n' + r'$\bf{Power\ balance}$ ' + f'{Pin/1e6:.3g} MW in, {Pout/1e6:.3g} MW out\n({Pinner/1e6:.3g} MW to inner target, {Pouter/1e6:.3g} MW to outer target, {Pwalls/1e6:.3g} MW to walls, {Prad/1e6:.3g} MW radiation)'
    except Exception as e:
        print(traceback.format_exc())
        txt += '\n' + r'$\bf{Power\ balance}$ ' + f"{type(e).__name__}: {e}"
            
    try:
        jobid = os.environ['SLURM_JOB_ID']
        #with open('current_run.txt', 'r') as f:
        #    current_run = f.read()
        #    jobid = int(current_run.split()[1])
        jobinfo = '\n'.join(get_cmd_output(['sacct', '-j', jobid, '--format=JobID,JobName,Partition,State,Elapsed,NodeList']).split('\n')[:3])
        txt += '\n\n' + jobinfo
        line = jobinfo.split('\n')[2]
        time_elapsed = line.split()[4]
        hours, minutes, seconds = map(int, time_elapsed.split(':'))
        seconds_elapsed = hours * 3600 + minutes * 60 + seconds
        txt += '\n\n' + f'Seconds per iteration: {seconds_elapsed/iterations:.2g}'
    except Exception as e:
        print(traceback.format_exc()) 
        txt += '\n\n' + r'$\bf{Job\ info}$ ' + f"{type(e).__name__}: {e}"
        
    try:
        txt += '\n\n' + get_cmd_output([f'{os.path.dirname(os.path.abspath(__file__))}/checkup'])
    except Exception as e:
        print(traceback.format_exc())   
        txt += '\n\n' + r'$\bf{Checkup}$ ' + f"{type(e).__name__}: {e}"
    
    return txt


def wrapText(text, width=100):
    newLines = []
    for l in text.split('\n'):
        if l == '':
            # Keep intentional newlines in original text
            newLines.append('')
        else:
            newLines.extend(textwrap.wrap(l, width=width))
    return '\n'.join(newLines)


def page0():
    fig = plt.figure(figsize=(8.5, 11), dpi=dpi)
    plt.axes(frameon=False)
    plt.gca().get_xaxis().set_ticks([])
    plt.gca().get_yaxis().set_ticks([])
    txt = wrapText(getConfigText())
    fig.text(0.1, 0.95, txt, transform=fig.transFigure, size=11, horizontalalignment='left', verticalalignment='top', fontfamily='monospace', fontsize=9)

    
def page1():
    # Plot D and chi using https://github.com/Reksoatr/SOLPS_Scripts/blob/cb54a040abeccd3fde71a8354f8927df3ab7669d/SOLPS_Scripts/B2TransportParser.py
    
    fig = plt.figure(figsize=(8.5, 11), dpi=dpi)
    plt.suptitle('Equilibrium: $%d\\times%d$, $B_T=%+.2g$ T @ $R=%.2g$ m, $I_P=%+.2g$ MA' % (so.geqdsk['NW'], so.geqdsk['NH'], so.geqdsk['BCENTR'], so.geqdsk['RCENTR'], so.geqdsk['CURRENT']/1e6))
    gs = fig.add_gridspec(2,2)
    ax = fig.add_subplot(gs[:, 0])
    plotmesh()
    plot_wall_geometry(so, ax)
    plt.scatter(so.geqdsk['RMAXIS'], so.geqdsk['ZMAXIS'], marker='+', c='salmon', s=200, linewidths=0.5, zorder=10)
    plt.axhline(so.geqdsk['ZMAXIS'], lw=0.5, ls=':', c='black')
    # showPoloidalLoc(JXA)
    # showPoloidalLoc(JXI)
    # showPoloidalLoc(sp.ixp)
    # showPoloidalLoc(sp.oxp)

    def plot2D(what, ax, levels=None, Z_in=None, **kw):
        if levels is None:
            if what in ['PHIRZ_NORM', 'RHOpRZ', 'RHORZ', 'PSIRZ_NORM']:
                levels = np.r_[0.1:10:0.1]
                label_levels = levels[:9]
            else:
                levels = np.linspace(np.nanmin(so.geqdsk['AuxQuantities'][what]), np.nanmax(so.geqdsk['AuxQuantities'][what]), 20)
                label_levels = levels
        else:
            label_levels = levels

        label = kw.pop('label', None)  # Take this out so the legend doesn't get spammed by repeated labels

        # use this to set up the plot key word args, get the next line color, and move the color cycler along
        (l,) = ax.plot(so.geqdsk['AuxQuantities']['R'], so.geqdsk['AuxQuantities']['R'] * np.nan, **kw)
        # contours
        cs = ax.contour(
            so.geqdsk['AuxQuantities']['R'],
            so.geqdsk['AuxQuantities']['Z'],
            so.geqdsk['AuxQuantities'][what],
            levels,
            colors=[l.get_color()] * len(levels),
            linewidths=l.get_linewidth(),
            alpha=l.get_alpha(),
            linestyles=l.get_linestyle(),
        )

        # get the color
        kw1 = copy.copy(kw)
        kw1['linewidth'] = kw['linewidth'] + 1
        kw1.setdefault('color', ax.lines[-1].get_color())

        # boundary
        #ax.plot(self['RBBBS'], self['ZBBBS'], label=label, **kw1)
        
    plot2D('RHORZ', plt.gca(), color='salmon', ls='--', linewidth=0.5)
    
    # Profiles of pressure and safety factor
    # Adapted from https://omfit.io/_modules/omfit_classes/omfit_eqdsk.html#OMFITgeqdsk.plot
    kw = {}
    xlabel_in_legend = False
    
    xName = '$\\rho$'
    if 'RHOVN' in so.geqdsk and np.sum(so.geqdsk['RHOVN']):
        x = so.geqdsk['RHOVN']
    else:
        x = so.geqdsk['AuxQuantities']['RHO']

    if 'label' not in kw:
        kw['label'] = (' '.join([a.strip() for a in so.geqdsk['CASE'][3:]])).strip()
        if not len(kw['label']):
            kw['label'] = (' '.join([a.strip() for a in so.geqdsk['CASE']])).strip()
            if not len(kw['label']):
                kw['label'] = os.path.split(so.geqdsk.filename)[1]
    if xlabel_in_legend:
        kw['label'] += ' vs ' + xName

    ax = fig.add_subplot(gs[0, 1])
    ax.plot(x, so.geqdsk['PRES'], **kw)
    kw.setdefault('color', ax.lines[-1].get_color())
    ax.set_title(r'$\,$ Pressure')
    ax.ticklabel_format(style='sci', scilimits=(-1, 2), axis='y')
    ax.set_xlabel((not xlabel_in_legend) * xName)

    ax = fig.add_subplot(gs[1, 1], sharex=ax)
    ax.plot(x, so.geqdsk['QPSI'], **kw)
    ax.set_title('Safety Factor')
    ax.set_xlabel((not xlabel_in_legend) * xName)
    
    
def getResiduals():
    file_in = "b2mn.exe.dir/b2ftrace" 
    if not os.path.isfile(file_in):
        file_in = "b2ftrace"
        if not os.path.isfile(file_in):
            print("Error: Could not find b2ftrace")

    with open(file_in) as f:
        lines = f.readlines()

    mydatalist = []
    counter = 0
    read_data = False
    for line in lines:
        if "data" in line:
            counter += 1
            read_data = True
            continue
        if read_data:
            line = line.split()
            part_list = [float(i) for i in line]
            mydatalist.extend(part_list)
        
    mydata = np.array(mydatalist)
    mydata = mydata.reshape(counter,int(len(mydatalist)/counter))
    return mydata


def page2():
    colnamesD = ['conD0', 'conD1', 'momD0', 'momD1', 'totmom', 'ee', 'ei', 'phi']
    colnamesDT = ['conD0', 'conD1', 'conT0', 'conT1', 'momD0', 'momD1', 'momT0', 'momT1', 'totmom', 'ee', 'ei', 'phi']
    colnames = colnamesDT
    
    fig, axs = plt.subplots(3, 2, figsize=(8.5,11), dpi=dpi)
    lw = 1
    
    try:
        ax = axs[0,0]
        ax.set_title('Residuals')

        resids = getResiduals()
        its = np.shape(resids)[0]

        for i in range(len(colnames)):
            col = colnames[i]
            if 'D0' in col or 'T0' in col or 'totmom' in col:
                continue
            ax.plot(resids[:,i+2], label=col, lw=lw)
        ax.set_yscale('log')
        ax.set_xlabel('Iterations')
        labelLines(ax.get_lines(), align=False)
    except Exception as e:
        print(traceback.format_exc())
    
    try:
        ax = axs[1,0]
        # Negative value of b2mndr_stim in b2mn.dat means this file gets appended to
        ds = nc.Dataset('b2time.nc')
        times = np.array(ds['timesa'])
        
        # ax.set_title('Plasma variables at separatrix [keV, $10^{20}$ m$^{-3}$]')
        # ax.plot(times[-its:]*1e3, np.array(ds['tesepm'])[-its:]/1e3, label='$T_{e\ MP}$', ls='-', lw=lw)
        # ax.plot(times[-its:]*1e3, np.array(ds['tisepm'])[-its:]/1e3, label='$T_{i\ MP}$', ls='-', lw=lw)
        # ax.plot(times[-its:]*1e3, np.array(ds['nesepm'])[-its:]/1e20, label='$n_{e\ MP}$', ls='-', lw=lw)
        # ax.plot(times[-its:]*1e3, np.array(ds['tesepa'])[-its:]/1e3, label='$T_{e\ OT}$', ls='--', lw=lw)
        # ax.plot(times[-its:]*1e3, np.array(ds['tisepa'])[-its:]/1e3, label='$T_{i\ OT}$', ls='--', lw=lw)
        # ax.plot(times[-its:]*1e3, np.array(ds['nesepa'])[-its:]/1e20, label='$n_{e\ OT}$', ls='--', lw=lw)
        
        ax.set_title('Change in plasma variables at separatrix')
        ax.plot(np.array(ds['tesepm'])[:]/np.array(ds['tesepm'])[0], label='$T_{e\ MP}$', ls='-', lw=lw)
        ax.plot(np.array(ds['tisepm'])[:]/np.array(ds['tisepm'])[0], label='$T_{i\ MP}$', ls='-', lw=lw)
        ax.plot(np.array(ds['nesepm'])[:]/np.array(ds['nesepm'])[0], label='$n_{e\ MP}$', ls='-', lw=lw)
        ax.plot(np.array(ds['tesepa'])[:]/np.array(ds['tesepa'])[0], label='$T_{e\ OT}$', ls='--', lw=lw)
        ax.plot(np.array(ds['tisepa'])[:]/np.array(ds['tisepa'])[0], label='$T_{i\ OT}$', ls='--', lw=lw)
        ax.plot(np.array(ds['nesepa'])[:]/np.array(ds['nesepa'])[0], label='$n_{e\ OT}$', ls='--', lw=lw)
        
        # ax.plot(sp.inn.te_sep, label='$T_{e\ IT}$', ls=':')
        # ax.plot(sp.inn.ti_sep, label='$T_{i\ IT}$', ls=':')
        # ax.plot(sp.inn.ne_sep/1e20, label='$n_{e\ IT}$', ls=':')
        ax.set_yscale('log')
        ax.set_xlabel('Iterations')
        labelLines(ax.get_lines(), align=False)
    except Exception as e:
        print(traceback.format_exc())
    
    try:
        ax = axs[2,0]
        # Modified from scripts/energy_balance.py
        # Negative value of b2mndr_stim in b2mn.dat means this file gets appended to
        if os.access('b2mn.exe.dir/b2tallies.nc', os.R_OK):
            f=nc.Dataset('b2mn.exe.dir/b2tallies.nc','r')
        else:
            f=nc.Dataset('b2tallies.nc','r')
            vreg=f.dimensions['vreg'].size
            xreg=f.dimensions['xreg'].size
            yreg=f.dimensions['yreg'].size
            ns=f.dimensions['ns'].size
            time=f.dimensions['time'].size
            times=f.variables['times']
            fhexreg=f.variables['fhexreg']
            fheyreg=f.variables['fheyreg']
            fhixreg=f.variables['fhixreg']
            fhiyreg=f.variables['fhiyreg']

        scale = 1e6
        # ax.set_title('Energy fluxes [MW]')
        # ax.plot(times[-its:]*1e3,(-fhexreg[-its:,1]-fhixreg[-its:,1])/scale, label='-W', lw=lw)
        # ax.plot(times[-its:]*1e3,(-fhexreg[-its:,2]-fhixreg[-its:,2])/scale, label='-w', lw=lw)
        # ax.plot(times[-its:]*1e3,( fhexreg[-its:,3]+fhixreg[-its:,3])/scale, label='e', lw=lw)
        # ax.plot(times[-its:]*1e3,( fhexreg[-its:,4]+fhixreg[-its:,4])/scale, label='E', lw=lw)
        # ax.plot(times[-its:]*1e3,( fheyreg[-its:,2]+fhiyreg[-its:,2])/scale, label='core', lw=lw)
        # ax.plot(times[-its:]*1e3,( fheyreg[-its:,4]+fhiyreg[-its:,4])/scale, label='sep', lw=lw)
        # ax.plot(times[-its:]*1e3,(np.sum(fheyreg[-its:,5:8],axis=1)+np.sum(fhiyreg[-its:,5:8],axis=1))/scale, label='mcw', lw=lw)
        # ax.plot(times[-its:]*1e3,(-fheyreg[-its:,1]-fhiyreg[-its:,1]-fheyreg[-its:,3]-fhiyreg[-its:,3])/scale, label='-pfw', lw=lw)
        
        ax.set_title('Change in energy fluxes')
        ax.plot((-fhexreg[:,1]-fhixreg[:,1])/(-fhexreg[0,1]-fhixreg[0,1]), label='-W', lw=lw)
        ax.plot((-fhexreg[:,2]-fhixreg[:,2])/(-fhexreg[0,2]-fhixreg[0,2]), label='-w', lw=lw)
        ax.plot(( fhexreg[:,3]+fhixreg[:,3])/( fhexreg[0,3]+fhixreg[0,3]), label='e', lw=lw)
        ax.plot(( fhexreg[:,4]+fhixreg[:,4])/( fhexreg[0,4]+fhixreg[0,4]), label='E', lw=lw)
        ax.plot(( fheyreg[:,2]+fhiyreg[:,2])/( fheyreg[0,2]+fhiyreg[0,2]), label='core', lw=lw)
        ax.plot(( fheyreg[:,4]+fhiyreg[:,4])/( fheyreg[0,4]+fhiyreg[0,4]), label='sep', lw=lw)
        ax.plot((np.sum(fheyreg[:,5:8],axis=1)+np.sum(fhiyreg[:,5:8],axis=1))/(np.sum(fheyreg[0,5:8])+np.sum(fhiyreg[0,5:8])), label='mcw', lw=lw)
        ax.plot((-fheyreg[:,1]-fhiyreg[:,1]-fheyreg[:,3]-fhiyreg[:,3])/(-fheyreg[0,1]-fhiyreg[0,1]-fheyreg[0,3]-fhiyreg[0,3]), label='-pfw', lw=lw)

        labelLines(ax.get_lines(), align=False)
        ax.set_xlabel('Iterations')
        ax.set_yscale('log')
    except Exception as e:
        print(traceback.format_exc())

    try:
        # Modified from scripts/particle_balance.py
        ax = axs[0,1]
        ax.set_title('Normalized particle error')
        if os.access('b2mn.exe.dir/b2tallies.nc', os.R_OK):
            f = nc.Dataset('b2mn.exe.dir/b2tallies.nc','r')
        else:
            f = nc.Dataset('b2tallies.nc','r')
        vreg = f.dimensions['vreg'].size
        xreg = f.dimensions['xreg'].size
        yreg = f.dimensions['yreg'].size
        ns = f.dimensions['ns'].size
        time = f.dimensions['time'].size
        times = f.variables['times']
        species_names=f.variables['species']
        species=[b''.join(species_names[i,:]).strip().decode('utf-8') for i in range(species_names.shape[0])]
        elements=[re.sub('[-+0-9]','',species[i]) for i in range(len(species))]
        mask=[re.match('[a-zA-Z]+0$',species[i])!=None for i in range(len(species))]
        s=0
        bounds=[]
        for i in range(mask.count(True)-1):
            e=mask.index(True,s+1)
            bounds+=[[s,e]]
            s=e
        bounds+=[[s,len(mask)]]

        fnaxreg = f.variables['fnaxreg']
        fnayreg = f.variables['fnayreg']
        b2stbr_sna_reg = f.variables['b2stbr_sna_reg']
        b2sext_sna_reg = f.variables['b2sext_sna_reg']

        if vreg == 5:
            FULL_X = np.array([0,1,0,0,-1,0,0])
            FULL_Y = np.array([0,1,1,1,0,-1,-1,-1])
        elif vreg ==2:
            FULL_X = np.array([0,1,-1])
            FULL_Y = np.array([0,1,-1])
        else:
            raise ValueError('Value of vreg=%s not currently coded' % vreg)

        for i in range(len(bounds)):
            fnax = (fnaxreg[:,bounds[i][0]:bounds[i][1],:].sum(axis=1)*FULL_X).sum(axis=1)
            fnay = (fnayreg[:,bounds[i][0]:bounds[i][1],:].sum(axis=1)*FULL_Y).sum(axis=1)
            b2stbr = b2stbr_sna_reg[:,bounds[i][0]:bounds[i][1],0].copy()
            b2stbr[:,0] = 0 # remove the neutral
            if b2stbr.shape[1] > 2: b2stbr[:,-1] = 0 # remove the fully stripped species for He and up
            b2stbr = b2stbr.sum(axis=1)

            b2sext = b2sext_sna_reg[:,bounds[i][0]:bounds[i][1],0].copy()
            b2sext = b2sext.sum(axis=1)

            fna_norm = np.max([np.max(np.abs(fnaxreg[:,bounds[i][0]:bounds[i][1],:])),np.max(np.abs(fnayreg[:,bounds[i][0]:bounds[i][1],:]))])
            ax.plot((fnax+fnay+b2stbr+b2sext)/fna_norm, label=elements[bounds[i][0]], lw=lw)
        ax.axhline(0.01, ls='--', lw=1, c='red')
        ax.axhline(-0.01, ls='--', lw=1, c='red')
        ax.set_xlabel('Iterations')
        labelLines(ax.get_lines(), align=False);
    except Exception as e:
        print(traceback.format_exc())
    
    try:
        ax = axs[1,1]
        ax.set_title('$Z_{eff}$ for each region')
        # Modified from scripts/zeffreg.py

        if os.access('b2mn.exe.dir/b2tallies.nc', os.R_OK):
            f=nc.Dataset('b2mn.exe.dir/b2tallies.nc','r')
        else:
            f=nc.Dataset('b2tallies.nc','r')
            vreg=f.dimensions['vreg'].size
            times=f.variables['times']
            nereg=f.variables['nereg']
            ne2reg=f.variables['ne2reg']

        for i in np.arange(0, vreg):
            ax.plot(ne2reg[:,i]/nereg[:,i], label=str(i), lw=lw)
        ax.set_xlabel('Iterations')
        labelLines(ax.get_lines(), align=False)
    except Exception as e:
        print(traceback.format_exc())

    try:
        ax = axs[2,1]
        # ax.set_title('Total particles and energy [$10^{20}$, $10^{20}$ keV]')
        # ax.plot(times[-its:]*1e3, np.array(ds['tmne'])[-its:]/1e20, label='num. particles', ls='-', lw=lw)
        # ax.plot(times[-its:]*1e3, np.array(ds['tmte'])[-its:]/1e20/1e3, label='e energy', ls='-', lw=lw)
        # ax.plot(times[-its:]*1e3, np.array(ds['tmti'])[-its:]/1e20/1e3, label='i energy', ls='-', lw=lw)
        
        ax.set_title('Change in total particles and energy')
        ax.plot(np.array(ds['tmne'])/np.array(ds['tmne'])[0], label='num. electrons', ls='-', lw=lw)
        ax.plot(np.array(ds['tmte'])/np.array(ds['tmte'])[0], label='e energy', ls='-', lw=lw)
        ax.plot(np.array(ds['tmti'])/np.array(ds['tmti'])[0], label='i energy', ls='-', lw=lw)

        ax.set_yscale('log')
        ax.set_xlabel('Iterations')
        labelLines(ax.get_lines(), align=False)
    except Exception as e:
        print(traceback.format_exc())

    plt.tight_layout()
        
    
def page3():
    fig, axs = plt.subplots(3, 3, figsize=(8.5,11), dpi=dpi)
    plotargs = {'marker': '+'}
    c0 = 'C0'
    c1 = 'C1'
    c2 = 'C2'
    lw1 = 1
    lw2 = 3
    
    ax = axs[0,0]
    ax.set_title('Midplane $n$ [$10^{20}$ m$^{-3}$]')
    try:
        nesep = interpolate.interp1d(yyc, so.ne[:,JXA]/1e20, kind='linear')(0)
        ax.plot(yyc, so.ne[:,JXA]/1e20, label=f'$n_e$ ($n_{{e\\ sep}}$={nesep:.3g}e20 m$^{{-3}}$)', **plotargs)
        ndsep = interpolate.interp1d(yyc, so.na[1,:,JXA]/1e20, kind='linear')(0)
        ax.plot(yyc, so.na[1,:,JXA]/1e20, label=f'$n_D$ ($n_{{D\\ sep}}$={ndsep:.3g}e20 m$^{{-3}}$)', **plotargs)
        if has_impurities:
            nimp = so.na[3:,:,JXA].sum(axis=0)
            nimpsep = interpolate.interp1d(yyc, nimp/1e20, kind='linear')(0)
            ax.plot(yyc, nimp/1e20, label=f'$n_{{{impurity}}}$ ($n_{{{impurity}\\ sep}}$={nimpsep*1e20:.3g} m$^{{-3}}$)', **plotargs)
        expfun = lambda x, A, lamda_inv: A*np.exp(-x*lamda_inv) # needs to be in this form for curve_fit to work
        mmstart = -0.2
        mmend = 2
        istart = argnear(yyc, mmstart)
        iend = argnear(yyc, mmend)
        isep = argnear(yyc, 0)
        qofit, _ = curve_fit(expfun, yyc[istart:iend], so.ne[istart:iend,JXA]/1e20, p0=[np.max(so.ne[istart:iend,JXA]/1e20),1], bounds=(0, np.inf))
        lamda = 1/qofit[1] # lamda_q in mm
        ax.plot(yyc[istart:iend], expfun(np.array(yyc[istart:iend]), *qofit), c='k', ls='-', label=f'$\lambda_{{ne}}$={lamda:.3g} mm')
        ax.legend(fontsize=7);
        ax.set_ylim([0,np.max([so.na[:,:,JXA].max(),so.ne[:,JXA].max()])/1e20*1.1])
    except Exception as e:
        print(traceback.format_exc())
    ax.set_xlabel(yyclabel)
    
    ax = axs[0,1]
    ax.set_title('Midplane temperature [eV]')
    ax.plot(yyc, so.te[:,JXA], label='$T_e$', c=c1, lw=lw1, **plotargs)
    ax.plot(yyc, so.ti[:,JXA], label='$T_i$', c=c1, lw=lw2, **plotargs)
    try:
        expfun = lambda x, A, lamda_inv: A*np.exp(-x*lamda_inv) # needs to be in this form for curve_fit to work
        mmstart = -0.2
        mmend = 2
        istart = argnear(yyc, mmstart)
        iend = argnear(yyc, mmend)
        isep = argnear(yyc, 0)
        qofit, _ = curve_fit(expfun, yyc[istart:iend], so.te[istart:iend,JXA], p0=[np.max(so.te[istart:iend,JXA]),1], bounds=(0, np.inf))
        lamda = 1/qofit[1] # lamda_q in mm
        tesep = interpolate.interp1d(yyc, so.te[:,JXA], kind='linear')(0)
        ax.plot(yyc[istart:iend], expfun(np.array(yyc[istart:iend]), *qofit), c='k', ls='-', label=f'$\lambda_T$={lamda:.3g} mm'+'\n$T_{e,sep}$='+f'{tesep:.2g} eV')
        ax.legend(fontsize=7);
    except Exception as e:
        print(traceback.format_exc())
    ax.set_xlabel(yyclabel)
    ax.legend()

    ax = axs[0,2]
    ax.set_title('Transport coefficients')
    try:
        plotTransportCoeffs(ax)
    except Exception as e:
        print(traceback.format_exc())
    
    ax = axs[1,0]
    ax.set_title('Inner target $n$ [$10^{20}$ m$^{-3}$]')
    for i in range(sp.na.shape[2]):
        if sp.za[i] > 0:
            alpha = 1-sp.za[i]/(sp.zn[i]+2)
            ax.plot(yyct, sp.na[0,1:-1,i], label=species_nice[i], c=C0123[nuclei_nice_set.index(nuclei_nice[i])]+(alpha,))
    ax.set_yscale('log')
    labelLines(ax.get_lines(), align=True, xvals=[yyct[2]]*sp.na.shape[2], fontsize=8, outline_width=3);
    ax.set_xlabel(yyclabel)
    
    ax = axs[1,1]
    ax.set_title('Inner target temperature [eV]')
    ax.plot(yyct, so.te[:,0], label='$T_e$', c=c1, lw=lw1, **plotargs)
    ax.plot(yyct, so.ti[:,0], label='$T_i$', c=c1, lw=lw2, **plotargs)
    ax.set_xlabel(yyclabel)
    ax.legend()

    ax = axs[1,2]
    ax.set_title('Inner target $j_{sat}$ [MA/m$^2$]')
    try:
        jsat = 1.6e-19*sp.na[:,:,isp]*np.sqrt((sp.zeff*sp.te*1.6e-19+sp.ti*1.6e-19)/sp.am_kg[isp])
        ax.plot(yyct, jsat[0,1:-1]/1e6, c=c2, **plotargs)
    except Exception as e:
        print(traceback.format_exc())
    ax.set_xlabel(yyclabel)
    
    ax = axs[2,0]
    ax.set_title('Outer target $n$ [$10^{20}$ m$^{-3}$]')
    for i in range(sp.na.shape[2]):
        if sp.za[i] > 0:
            alpha = 1-sp.za[i]/(sp.zn[i]+2)
            ax.plot(yyct, sp.na[-1,1:-1,i], label=species_nice[i], c=C0123[nuclei_nice_set.index(nuclei_nice[i])]+(alpha,))
    ax.set_yscale('log')
    labelLines(ax.get_lines(), align=True, xvals=[yyct[2]]*sp.na.shape[2], fontsize=8, outline_width=3);
    ax.set_xlabel(yyclabel)
    
    ax = axs[2,1]
    ax.set_title('Outer target temperature [eV]')
    ax.plot(yyct, so.te[:,-1], label='$T_e$', c=c1, lw=lw1, **plotargs)
    ax.plot(yyct, so.ti[:,-1], label='$T_i$', c=c1, lw=lw2, **plotargs)
    ax.set_xlabel(yyclabel)
    ax.legend()

    ax = axs[2,2]
    ax.set_title('Outer target $j_{sat}$ [MA/m$^2$]')
    try:
        ax.plot(yyct, jsat[-1,1:-1]/1e6, c=c2, **plotargs)
    except Exception as e:
        print(traceback.format_exc())
    ax.set_xlabel(yyclabel)
    
    plt.tight_layout()
    

def page4():
    fig = plt.figure(figsize=(8.5,11), dpi=dpi)
    gs = fig.add_gridspec(3,2)
    plotargs = {'marker': '+'}

    pstatic = sp.psi+sp.pse #sp.ne*(sp.te+sp.ti)*sp.qe
    pdyn = sp.pdynd #sp.mi[:,:,isp]*sp.ne*sp.ua[:,:,isp]**2 # ua=cs at 0 and -1 so use those instead of reasoning by cell face flux
    ptot = pstatic+pdyn
    ptot = sp.psi+sp.pse+sp.pdynd
    ptot = sp.ptot
    iy = sp.sep-1
    ixp = spq.leftcut[0]+1
    oxp = spq.rightcut[0]+1

    # Resources on heat flux components:
    # https://repo.tok.ipp.cas.cz/kripner/solps-doc/-/blob/master/doc/Katka's_notes_on_SOLPS-ITER.md#solpspy-heat-flux-components
    # https://repo.tok.ipp.cas.cz/kripner/solps-doc/-/blob/master/doc/Comparing_experiment_to_model.md#parallel-heat-flux-q_parallel
    # Appendix A https://iopscience.iop.org/article/10.1088/1741-4326/ac9917/pdf
    # Manual page 471
    Ptot = spq.fhtx
    Ptotei = spq.fhex + spq.fhix
    qpol = Ptot/sp.sx # wrong, sx is different here
    qpar = qpol/sp.pitch
    qsurf = Ptot/sp.sx # not sure if sx is east or west, sx[0,:] is all zeros, and difference between sx[1,:] and sx[2,:] is ~1%
    qsurfei = Ptotei/sp.sx

    ax = fig.add_subplot(gs[0, 0])
    ax.set_title('Inboard static+dynamic pressure [Pa]')
    ax.plot(yyc[iy-1:], pstatic[JXI,iy:-1]+pdyn[JXI,iy:-1], c='C1', lw=1, label='Inner midplane', **plotargs)
    ax.plot(yyc[iy-1:], pstatic[0,iy:-1]+pdyn[0,iy:-1], c='C1', lw=3, label='Inner target', **plotargs)
    ax.set_xlabel(yyclabel)
    ax.set_yscale('log')
    ax.legend()

    ax = fig.add_subplot(gs[0, 1])
    ax.set_title('Outboard static+dynamic pressure [Pa]')
    ax.plot(yyc[iy-1:], pstatic[JXA,iy:-1]+pdyn[JXA,iy:-1], c='C0', lw=1, label='Outer midplane', **plotargs)
    ax.plot(yyc[iy-1:], pstatic[-1,iy:-1]+pdyn[-1,iy:-1], c='C0', lw=3, label='Outer target', **plotargs)
    ax.set_xlabel(yyclabel)
    ax.set_yscale('log')
    ax.legend()

    def fit_exp(yyc, q, color):
        def expfun(x, A, lamda_inv):
            """
            Exponential function for curve fitting.
            
            Parameters:
            x: Array of x values
            A: Amplitude of exponential
            lamda_inv: Inverse decay length (1/lambda) (necessary for curve_fit to work nicely)
            """
            return A*np.exp(-x*lamda_inv)
        
        istart = np.argmax(q)
        mmstart = yyc[istart]
        mmend = mmstart+1 # fit up to 1 mm away from peak
        iend = argnear(yyc, mmend)
        # if heat flux profile is non-monotonic, only fit part before it starts rising again
        for i in range(istart+1, iend+1):
            if q[i] > q[i-1]:
                iend = i-1
                mmend = yyc[iend]
                break
        qofit, _ = curve_fit(expfun, yyc[istart:iend+1], q[istart:iend+1], p0=[np.max(q),1], bounds=(0, np.inf))
        lamda = 1/qofit[1] # lamda_q in mm
        ax.plot(yyc[istart:iend+1], expfun(np.array(yyc[istart:iend+1]), *qofit), c=color, alpha=0.4, label=f'$\lambda_q$={lamda:.3g} mm')
        return lamda
    
    # plot and fit qpar e+i at X-point
    ax = fig.add_subplot(gs[1, 0])
    ax.set_title('$q_{||\ e+i}$ at X-point [GW/m$^2$]')
    ptot = spq.fhex+spq.fhix
    sxpol = sp.sx*sp.qc
    q = ptot/sxpol*wbbl[:,:,3]/np.abs(wbbl[:,:,0])
    ax = plt.gca()
    ax.plot(yyct[sp.sep-1:], q[oxp,sp.sep:-1]/1e9, label=f'Outboard ({np.sum(ptot[oxp,sp.sep:-1])/1e6:.3g} MW)', **plotargs)
    ax.plot(yyct[sp.sep-1:], -q[ixp,sp.sep:-1]/1e9, label=f'Inboard ({np.sum(-ptot[ixp,sp.sep:-1])/1e6:.3g} MW)', **plotargs)
    try:
        lamda_xo = fit_exp(yyct[sp.sep-1:], q[oxp,sp.sep:-1]/1e9, 'C0')
    except Exception as e:
        print(traceback.format_exc())
        lamda_xo = None
    try:
        lamda_xi = fit_exp(yyct[sp.sep-1:], -q[ixp,sp.sep:-1]/1e9, 'C1')
    except Exception as e:
        print(traceback.format_exc())
        lamda_xi = None
    ax.set_yscale('log')
    ax.set_xlabel(yyclabel)
    ax.legend(fontsize=8, loc='upper right')

    # plot and fit qpar tot at X-point
    ax = fig.add_subplot(gs[1, 1])
    ax.set_title('$q_{||\ tot}$ at X-point [GW/m$^2$]')
    ptot = spq.fhtx
    sxpol = sp.sx*sp.qc
    q = ptot/sxpol*wbbl[:,:,3]/np.abs(wbbl[:,:,0])
    ax = plt.gca()
    ax.plot(yyct[sp.sep-1:], q[oxp,sp.sep:-1]/1e9, label=f'Outboard ({np.sum(ptot[oxp,sp.sep:-1])/1e6:.3g} MW)', **plotargs)
    ax.plot(yyct[sp.sep-1:], -q[ixp,sp.sep:-1]/1e9, label=f'Inboard ({np.sum(-ptot[ixp,sp.sep:-1])/1e6:.3g} MW)', **plotargs)
    try:
        lamda_xo = fit_exp(yyct[sp.sep-1:], q[oxp,sp.sep:-1]/1e9, 'C0')
    except Exception as e:
        print(traceback.format_exc())
        lamda_xo = None
    try:
        lamda_xi = fit_exp(yyct[sp.sep-1:], -q[ixp,sp.sep:-1]/1e9, 'C1')
    except Exception as e:
        print(traceback.format_exc())
        lamda_xi = None
    ax.set_yscale('log')
    ax.set_xlabel(yyclabel)
    ax.legend(fontsize=8, loc='upper right')

    # plot and fit qsurf total
    ax = fig.add_subplot(gs[2, :])
    ax.set_title('$q_{surf\ tot}$ [MW/m$^2$]')

    #-radial profile of qpar below entrance to the outer leg
    psurfo = Ptot[-1,:]
    qsurfo = psurfo[1:-1]/sp.sx[-1,1:-1]
    #-radial profile of qpar below entrance to the inner leg
    psurfi = -Ptot[1,:]
    qsurfi = psurfi[1:-1]/sp.sx[1,1:-1]

    def qEich(rho, q0, S, lqi, qbg, rho_0):
            rho = rho - rho_0
            # lqi is inverse lamda_q
            return q0/2*np.exp((S*lqi/2)**2-rho*lqi)*erfc(S*lqi/2-rho/S)+qbg

    if lamda_xo != None:
        lqoGuess = lamda_xo
    else:
        lqoGuess = 1 # mm
    if lamda_xi != None:
        lqiGuess = lamda_xi
    else:
        lqiGuess = 1 # mm
    # lamda_q fits
    yycm = yyct/1000
    yycml = np.linspace(yycm[0], yycm[-1], 500)
    bounds = ([0,0,0,0,yycm[0]], [np.inf,np.inf,np.inf,np.inf,yycm[-1]])
    # Fit outer
    oguess = (np.max(qsurfo)-np.min(qsurfo), lqoGuess/1000/2, 1000/lqoGuess, np.min(qsurfo), 0)
    try:
        qsurfol = interpolate.interp1d(yycm, qsurfo)(yycml)
        qsofit, _ = curve_fit(qEich, yycml, qsurfol, p0=oguess, bounds=bounds)
        lqeo, So = 1000/qsofit[2], qsofit[1]*1000 # lamda_q and S in mm
    except Exception as e:
        print('qsurf outer fit failed:', e)
        qsofit = None
    iguess = (np.max(qsurfi)-np.min(qsurfi), lqiGuess/1000/2, 1000/lqiGuess, np.min(qsurfi), 0)
    # Fit inner
    try:
        qsurfil = interpolate.interp1d(yycm, qsurfi)(yycml)
        qsifit, _ = curve_fit(qEich, yycml, qsurfil, p0=iguess, bounds=bounds)
        lqei, Si = 1000/qsifit[2], qsifit[1]*1000 # lamda_q and S in mm 
    except Exception as e:
        print('qsurf inner fit failed:', e)
        qsifit = None

    ax.plot(yyct, qsurf[-1,1:-1]/1e6, c='C0', label=f'Outer target ({np.sum(Ptot[-1,1:-1])/1e6:.3g} MW)', **plotargs)
    ax.plot(yyct, -qsurf[1,1:-1]/1e6, c='C1', label=f'Inner target ({np.sum(-Ptot[1,1:-1])/1e6:.3g} MW)', **plotargs)

    if np.any(qsofit):
        ax.plot(yycml*1000, qEich(yycml, *qsofit)/1e6, c='C0', alpha=0.4,
                label=r'Eich fit: $\lambda_q$ = %.3g mm, $S$ = %.3g mm' % (lqeo, So))
    if np.any(qsifit):
        ax.plot(yycml*1000, qEich(yycml, *qsifit)/1e6, c='C1', alpha=0.4,
                label=r'Eich fit: $\lambda_q$ = %.3g mm, $S$ = %.3g mm' % (lqei, Si))

    # ax.set_yscale('log')
    ax.legend(fontsize=8, loc='upper right')
    ax.set_xlabel(yyclabel)
    plt.tight_layout()


def plot2d_b2(so, vals, ax=None, scale="log", label="", lb=None, ub=None, cmap='viridis', minratio=1e-4, **kwargs):
    """Modified from aurora.
    Method to plot 2D fields on B2 grids. 
    Colorbars are set to be manually adjustable, allowing variable image saturation.

    Parameters
    ----------
    vals : array (self.data('ny'), self.data('nx'))
        Data array for a variable of interest.
    ax : matplotlib Axes instance
        Axes on which to plot. If left to None, a new figure is created.
    scale : str
        Choice of 'linear','log' and 'symlog' for matplotlib.colors.
    label : str
        Label to set on the colorbar. No label by default.
    lb : float
        Lower bound for colorbar. If left to None, the minimum value in `vals` is used.
    ub : float
        Upper bound for colorbar. If left to None, the maximum value in `vals` is used.
    kwargs
        Additional keyword arguments passed to the `PatchCollection` class.
    """
    if ax is None:
        fig, ax = plt.subplots(1, figsize=(9, 11))

    # get polygons describing B2 grid
    p = so.get_b2_patches()
    p.set_cmap(cmap)
    p.set_edgecolor(None)

    if (
        so.form == "files"
        and np.prod(vals.shape)
        == so.data("crx").shape[1] * so.data("crx").shape[2]
    ):
        vals = vals[1:-1, 1:-1]

    # fill patches with values
    _vals = vals.flatten()

    # Avoid zeros that may derive from low Monte Carlo statistics
    if np.any(_vals == 0):
        nsm = nsmallest(2, np.unique(_vals))
        if len(nsm) > 1:
            _vals[_vals == 0.0] = nsm[1]
        else: # all zeros
            scale = "linear"

    p.set_array(np.array(_vals))

    if lb is None:
        lb = np.min(_vals)
    if ub is None:
        ub = np.max(_vals)

    extend = 'neither'
    if scale == "linear":
        p.set_clim([lb, ub])
    elif scale == "log":
        if lb < 0:
            lb = ub/1e3
        if lb == 0:
            lb = 1e-10
        if lb < ub*minratio:
            lb = ub*minratio
            extend = 'min'
        p.norm = matplotlib.colors.LogNorm(vmin=lb, vmax=ub)
    elif scale == "symlog":
        ub = np.max(np.abs(_vals))
        p.norm = matplotlib.colors.SymLogNorm(
            linthresh=ub / 100.0, base=10, linscale=1, vmin=-ub, vmax=ub
        )
    else:
        raise ValueError("Unrecognized scale parameter")

    ax.add_collection(p)

    cbar = plt.colorbar(
        p,
        ax=ax,
        pad=0.01,
        # ticks=[ub, ub / 10, lb / 10, lb] if scale == "symlog" else None,
        orientation='horizontal',
        extend=extend
    )
    cbar.ax.tick_params(labelsize=8)

    ax.set_title(label)
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.axis("scaled")


def plot2d_eirene(so, vals, ax=None, scale="log", label="", lb=None, ub=None,
    replace_zero=True, lcfs=False, minratio=1e-4, **kwargs):
    """Method to plot 2D fields from EIRENE.

    Parameters
    ----------
    vals : array (self.triangles)
        Data array for an EIRENE variable of interest.
    ax : matplotlib Axes instance
        Axes on which to plot. If left to None, a new figure is created.
    scale : str
        Choice of 'linear','log' and 'symlog' for matplotlib.colors.
    label : str
        Label to set on the colorbar. No label by default.
    lb : float
        Lower bound for colorbar. If left to None, the minimum value in `vals` is used.
    ub : float
        Upper bound for colorbar. If left to None, the maximum value in `vals` is used.
    replace_zero : boolean
        If True (default), replace all zeros in 'vals' with minimum value in 'vals'
    kwargs
        Additional keyword arguments passed to the `tripcolor` function.

    """
    # Avoid zeros that may derive from low Monte Carlo statistics
    # np.nan_to_num(vals,copy=False)

    if replace_zero:
        if np.any(vals == 0):
            vals[vals == 0.0] = nsmallest(2, np.unique(vals.flatten()))[1]

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 11))

    if lcfs and hasattr(so, "geqdsk") and "RBBBS" in so.geqdsk:
        # plot LCFS
        ax.plot(so.geqdsk["RBBBS"], so.geqdsk["ZBBBS"], c="k")

    if lb is None:
        lb = np.nanmin(vals)
    if ub is None:
        ub = np.nanmax(vals)

    # given quantity is on EIRENE triangulation
    extend = 'neither'
    if scale == "linear":
        norm = None
    elif scale == "log":
        if lb < ub*minratio:
            lb = ub*minratio
            extend = 'min'
        norm = matplotlib.colors.LogNorm(vmin=lb, vmax=ub)
    elif scale == "symlog":
        norm = matplotlib.colors.SymLogNorm(
            linthresh=ub / 10.0, base=10, linscale=0.5, vmin=lb, vmax=ub
        )
    else:
        raise ValueError("Unrecognized scale parameter")

    cntr = ax.tripcolor(
        so.xnodes,
        so.ynodes,
        so.triangles,
        facecolors=vals.flatten(),
        norm=norm,
        cmap=cmapdefault,
        **kwargs,
    )
    cbar = plt.colorbar(cntr, ax=ax, pad=0.01, orientation='horizontal', extend=extend)
    cbar.ax.tick_params(labelsize=8)

    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.set_title(label)
    ax.axis("scaled")
    
    
def page5():
    
    def common(ax, wallGeometry=True, cellEdges=False):
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)
        ax.set_facecolor('gray')
        ax.axis('equal') # 'scaled' before left a lot of whitespace
        if wallGeometry:
            plot_wall_geometry(so, ax)
        if not cellEdges:
            try:
                ax.collections[0].set_edgecolor(None)
            except:
                pass
        
    fig, axs = plt.subplots(2,4, figsize=(8.5,11), dpi=dpi)
    
    try:
        ax = axs[0,0]
        plot2d_b2(so, so.data('ne'), ax=ax, label=r'$n_e$ [m$^{-3}$]')
        # plot2d_b2(so, sp.ne[1:-1,1:-1].T, ax=ax, label=r'$n_e$ [m$^{-3}$]')
        common(ax)
    except Exception as e:
        ax.axis('off')
        print(traceback.format_exc())
    
    try:
        ax = axs[0,1]
        # atomic density (sum over species)
        if sp.eirene_names == ['H']:
            # If only hydrogen-type species are being modeled (includes D, T), sum over them
            val = np.sum(sp.pdena[:, :], axis=1)
        else:
            # If non-H species are included too, show the first of the H species
            val = sp.pdena[:, sp.eirene_species('H')]
        plot2d_eirene(so, val, scale="log", label=r"$n_{H\ atoms}$ [$m^{-3}$]", ax=ax)
        common(ax, wallGeometry=True, cellEdges=True)
    except Exception as e:
        ax.axis('off')
        print(traceback.format_exc())

    try:
        ax = axs[0,2]
        # molecular density (sum over species)
        if sp.eirene_names == ['H']:
            # If only hydrogen-type species are being modeled (includes D, T), sum over them
            val = np.sum(sp.pdenm[:, :], axis=1)
        else:
            # If non-H species are included too, show the first of the H species
            val = sp.pdenm[:, sp.eirene_species('H')]
        plot2d_eirene(so, val, scale="log", label=r"$n_{H\ molecules}$ [$m^{-3}$]", ax=ax)
        common(ax, wallGeometry=True, cellEdges=True)
    except Exception as e:
        ax.axis('off')
        print(traceback.format_exc())

    try:
        ax = axs[0,3]
        # https://gitlab.com/ipar/solpspy/-/blob/master/public/apidoc/rundir_cdn.html#L2804
        scale = 1/float(sp.read_switch('b2mndr_rescale_neutrals_sources'))
        plot2d_b2(so, (sp.rsana[:,:,0]*sp.dab2[:,:,0]/(sp.na[:,:,0]*scale)/sp.vol)[1:-1,1:-1].T, ax=ax, label=r'Ionization rate [s$^{-1}$ m$^{-3}$]', scale="log")
        # plot2d_b2(so, sp.ne[1:-1,1:-1].T, ax=ax, label=r'$n_e$ [m$^{-3}$]')
        common(ax)
    except Exception as e:
        ax.axis('off')
        print(traceback.format_exc())

    try:
        ax = axs[1,0]
        plot2d_b2(so, so.data('te'), ax=ax, label=r'$T_e$ [eV]')
        common(ax)
    except Exception as e:
        ax.axis('off')
        print(traceback.format_exc())

    try:
        ax = axs[1,1]
        plot2d_b2(so, so.data('ti'), ax=ax, label=r'$T_i$ [eV]')
        common(ax)
    except Exception as e:
        ax.axis('off')
        print(traceback.format_exc())

    try:
        ax = axs[1,2]
        plot2d_b2(so, (sp.po-sp.po.min()+1)[1:-1,1:-1].T, ax=ax, label=r'$\phi$ [V]')
        common(ax)
    except Exception as e:
        ax.axis('off')
        print(traceback.format_exc())

    try:
        ax = axs[1,3]
        cs = np.sqrt((sp.te+sp.ti)*1.6e-19/(2*1.66e-27))
        mach = sp.ua[:,:,1]/cs
        plot2d_b2(so, mach[1:-1,1:-1].T, ax=ax, label=r'Mach number', scale="symlog", cmap=plt.cm.bwr)
        common(ax)
    except Exception as e:
        ax.axis('off')
        print(traceback.format_exc())

    # ax = axs[0,2]
    # For ionization see "SOLPS variables" note
    # common(ax)
    
    plt.tight_layout()


def page6():
    
    def common(ax, wallGeometry=True, cellEdges=False):
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)
        ax.set_facecolor('gray')
        ax.axis('equal') # 'scaled' before left a lot of whitespace
        if wallGeometry:
            plot_wall_geometry(so, ax)
        if not cellEdges:
            try:
                ax.collections[0].set_edgecolor(None)
            except:
                pass
        
    fig, axs = plt.subplots(2,4, figsize=(8.5,11), dpi=dpi)

    cmap = 'viridis'

    try:
        ax = axs[0,0]
        plot2d_b2(so, (sp.fchx[1:-1,1:-1]/sp.sx[1:-1,1:-1]).T, ax=ax, label=r'Poloidal current [A/m$^2$]', scale='symlog', cmap=plt.cm.bwr)
        common(ax)
    except Exception as e:
        ax.axis('off')
        print(traceback.format_exc())

    try:
        ax = axs[0,1]
        plot2d_b2(so, (sp.rqbrm[1:-1,1:-1,:].sum(2)/sp.vol[1:-1,1:-1]).T/1e6, ax=ax, label=f'Bremsstrahlung [MW/m$^3$]\n({np.sum(sp.rqbrm[1:-1,1:-1,:].sum(2))/1e6:.3g} MW)', cmap=cmapdefault)
        common(ax)
    except Exception as e:
        ax.axis('off')
        print(traceback.format_exc())

    try:
        ax = axs[0,2]
        plot2d_b2(so, (sp.line_radiation.sum(axis=2)/sp.vol)[1:-1,1:-1].T/1e6, ax=ax, label=f'Line radiation [MW/m$^3$]\n({sp.line_radiation.sum()/1e6:.3g} MW)', cmap=cmapdefault)
        common(ax)
    except Exception as e:
        ax.axis('off')
        print(traceback.format_exc())
        
    try:
        ax = axs[0,3]
        plot2d_b2(so, (sp.rqahe[1:-1,1:-1].sum(axis=2)/sp.vol[1:-1,1:-1]).T/1e6, ax=ax, label=f'Electron cooling [MW/m$^3$]\n({sp.rqahe[1:-1,1:-1].sum()/1e6:.3g} MW)', cmap=cmapdefault)
        common(ax)
    except Exception as e:
        ax.axis('off')
        print(traceback.format_exc())

    try:
        ax = axs[1,0]
        if has_impurities:
            iimp = sp.eirene_names.index(impurity)
            plot2d_eirene(so, sp.pdena[:, iimp], scale="log", label=f"{impurity}"+ r"$^0$ dens. [m$^{-3}$]", ax=ax)
            common(ax, wallGeometry=True, cellEdges=True)
        else:
            ax.axis('off')
    except Exception as e:
        ax.axis('off')
        print(traceback.format_exc())

    try:
        ax = axs[1,1]
        if has_impurities:
            b2_impurity_indices = [i for i in range(len(sp.species_names)) if impurity in sp.species_names[i] and sp.za[i] > 0]
            plot2d_b2(so, np.sum(sp.na[:,:,b2_impurity_indices], axis=2)[1:-1,1:-1].T, ax=ax, label=f'{impurity}'+r'$^{+...}$ dens. [m$^{-3}$]', cmap=cmap)
            common(ax)
        else:
            ax.axis('off')
    except Exception as e:
        ax.axis('off')
        print(traceback.format_exc())

    try:
        ax = axs[1,2]
        if has_impurities:
            imp_ion_indices = np.where((np.array(nuclei)==impurity)*(sp.za >= 1))[0]
            # h_ion_indices = np.where((np.array(nuclei)=='H')*(sp.za >= 1))[0]
            # h_dens = np.sum(sp.na[:,:,h_ion_indices], axis=2)
            imp_dens = np.sum(sp.na[:,:,imp_ion_indices], axis=2)
            imp_percent = 100*imp_dens/sp.ne
            plot2d_b2(so, imp_percent[1:-1,1:-1].T, ax=ax, label=f'{impurity} ion %', scale='linear', cmap=cmap)
            common(ax)
        else:
            ax.axis('off')
    except Exception as e:
        ax.axis('off')
        print(traceback.format_exc())

    try:
        ax = axs[1,3]
        if has_impurities:
            plot2d_b2(so, np.abs(sp.zeff[1:-1,1:-1].T), ax=ax, label=r'$Z_{eff}$', scale='linear', cmap=cmap)
            common(ax)
        else:
            ax.axis('off')
    except Exception as e:
        ax.axis('off')
        print(traceback.format_exc())


def page7():
    fig = plt.figure(figsize=(8.5,11), dpi=dpi)
    gs = fig.add_gridspec(2, 2)
    plotargs = {'marker': ''}
    legend_fontsize = 8

    fhtx = spq.fhtx
    fhp = spq.fhp.sum(3) # sum over species
    fhm = spq.fhm.sum(3) # sum over species
    myfht = spq.fhe + spq.fhi + spq.fhj + fhp + fhm + spq.fnt
    fhpx = fhp[:,:,0]
    fhmx = fhm[:,:,0]
    fntx = spq.fnt[:,:,0]
    myfhtx = myfht[:,:,0]

    def common(ix, sign, ax):
        #ax.plot(yyc, sign*(sp.fhtx/sp.sx)[ix,1:-1]/1e6, c='k', label=f'Total: {sign*np.sum(sp.fhtx[ix,1:-1])/1e6:.3g} MW', marker='+')
        ax.plot(yyct, sign*(fhtx/sp.sx)[ix,1:-1]/1e6, c='k', label=f'Total: {sign*np.sum(fhtx[ix,1:-1])/1e6:.3g} MW', **plotargs)
        ax.plot(yyct, sign*(spq.fhex/sp.sx)[ix,1:-1]/1e6, label=f'Electron thermal energy: {sign*np.sum(spq.fhex[ix,1:-1])/1e6:.3g} MW', **plotargs)
        ax.plot(yyct, sign*(spq.fhix/sp.sx)[ix,1:-1]/1e6, label=f'Ion thermal energy: {sign*np.sum(spq.fhix[ix,1:-1])/1e6:.3g} MW', **plotargs)
        ax.plot(yyct, sign*(fhpx/sp.sx)[ix,1:-1]/1e6, label=f'Potential/ionization: {sign*np.sum(fhpx[ix,1:-1])/1e6:.3g} MW', ls='-', **plotargs)
        ax.plot(yyct, sign*(fhmx/sp.sx)[ix,1:-1]/1e6, label=f'Kinetic energy: {sign*np.sum(fhmx[ix,1:-1])/1e6:.3g} MW', ls='--',**plotargs)
        ax.plot(yyct, sign*(spq.fhjx/sp.sx)[ix,1:-1]/1e6, label=f'Electrostatic energy: {sign*np.sum(spq.fhjx[ix,1:-1])/1e6:.3g} MW', ls=':',**plotargs)
        ax.plot(yyct, sign*(fntx/sp.sx)[ix,1:-1]/1e6, label=f'Viscous corrections: {sign*np.sum(fntx[ix,1:-1])/1e6:.3g} MW', ls='-.', **plotargs)
        ax.set_xlabel(yyclabel)
        ax.set_ylabel('$q_{surf}$ [MW/m$^2$]')
        #ax.axhline(0, c='gray', lw=0.5, ls=':')
        ax.legend(fontsize=legend_fontsize, bbox_to_anchor=(0.5, -0.15), loc='upper center')

    # Resources on heat flux components:
    # https://repo.tok.ipp.cas.cz/kripner/solps-doc/-/blob/master/doc/Katka's_notes_on_SOLPS-ITER.md#solpspy-heat-flux-components
    # https://repo.tok.ipp.cas.cz/kripner/solps-doc/-/blob/master/doc/Comparing_experiment_to_model.md#parallel-heat-flux-q_parallel
    # Appendix A https://iopscience.iop.org/article/10.1088/1741-4326/ac9917/pdf
    # Manual page 471
    try:
        Ptot = sp.fhtx #sp.fhex + sp.fhix
        myfhtx = sp.fhex + sp.fhix + np.sum(sp.fhpx, 2) + np.sum(sp.fhmx, 2) + sp.fhjx #+ sp.fnt[:,:,0]
        Ptotei = sp.fhex + sp.fhix
        qpol = Ptot/sp.sx # wrong, sx is different here
        qpar = qpol/sp.pitch
        qsurf = Ptot/sp.sx # not sure if sx is east or west, sx[0,:] is all zeros, and difference between sx[1,:] and sx[2,:] is ~1%
        myqsurf = myfhtx
        qsurfei = Ptotei/sp.sx
        
        ax = fig.add_subplot(gs[0, 0])
        ax.set_title('Inner target')
        common(1, -1, ax)

        ax = fig.add_subplot(gs[0, 1])
        ax.set_title('Outer target')
        common(-1, 1, ax)
    except Exception as e:
        print(traceback.format_exc())

    # wlld plots
    # Descriptions of terms here: https://solps.pages.tok.ipp.cas.cz/solps-doc/Energy_fluxes_deep_dive/#divertor-target-loads-from-wlld
    try:
        if has_wlld:
            ax = fig.add_subplot(gs[1, 0])
            ax.set_title('Inner target')
            sx = sp.sx[1,1:-1]
            ax.plot(yyct, it['Wtot']/1e6, label=f'Total ({np.sum(it["Wtot"]*sx)/1e6:.2g} MW)', c='black')
            ax.plot(yyct, it['Whtpl']/1e6, label=f'Plasma particles impinging ({np.sum(it["Whtpl"]*sx)/1e6:.2g} MW)')
            ax.plot(yyct, it['Wptpl']/1e6, label=f'Plasma recombination ({np.sum(it["Wptpl"]*sx)/1e6:.2g} MW)', ls='--')
            ax.plot(yyct, it['Wneut']/1e6, label=f'Neutral particles ({np.sum(it["Wneut"]*sx)/1e6:.2g} MW)')
            ax.plot(yyct, it['Wrad']/1e6, ls=':', label=f'Radiation ({np.sum(it["Wrad"]*sx)/1e6:.2g} MW)')
            ax.legend(fontsize=legend_fontsize, bbox_to_anchor=(0.5, -0.15), loc='upper center')
            ax.set_ylabel('$q_{surf}$ [MW/m$^2$]')
            ax.set_xlabel(yyclabel);

            ax = fig.add_subplot(gs[1, 1])
            ax.set_title('Outer target')
            sx = sp.sx[-1,1:-1]
            ax.plot(yyct, ot['Wtot']/1e6, label=f'Total ({np.sum(ot["Wtot"]*sx)/1e6:.2g} MW)', c='black')
            ax.plot(yyct, ot['Whtpl']/1e6, label=f'Plasma particles impinging ({np.sum(ot["Whtpl"]*sx)/1e6:.2g} MW)')
            ax.plot(yyct, ot['Wptpl']/1e6, label=f'Plasma recombination ({np.sum(ot["Wptpl"]*sx)/1e6:.2g} MW)', ls='--')
            ax.plot(yyct, ot['Wneut']/1e6, label=f'Neutral particles ({np.sum(ot["Wneut"]*sx)/1e6:.2g} MW)')
            ax.plot(yyct, ot['Wrad']/1e6, ls=':', label=f'Radiation ({np.sum(ot["Wrad"]*sx)/1e6:.2g} MW)')
            ax.legend(fontsize=legend_fontsize, bbox_to_anchor=(0.5, -0.15), loc='upper center')
            ax.set_ylabel('$q_{surf}$ [MW/m$^2$]')
            ax.set_xlabel(yyclabel);
        else:
            pass
    except Exception as e:
        print(traceback.format_exc())

    plt.tight_layout()


def page8():
    fig = plt.figure(figsize=(8.5,11), dpi=dpi)
    gs = fig.add_gridspec(3, 3)
    plotargs = {'marker': ''}
    legend_fontsize = 8

    def common(ax, wallGeometry=True, cellEdges=False):
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)
        ax.set_facecolor('gray')
        ax.axis('equal') # 'scaled' before left a lot of whitespace
        if wallGeometry:
            plot_wall_geometry(so, ax)
        if not cellEdges:
            try:
                ax.collections[0].set_edgecolor(None)
            except:
                pass
        ax.set_xlim([1.290, 1.75])
        ax.set_ylim([-1.385, -0.834])

    try:
        ax = fig.add_subplot(gs[0, 0])
        plt.plot(yyct[sp.sep-1:], ((sp.na_eirene[-1,1:-1,0]+2*sp.dmb2[-1,1:-1,0])/sp.na[JXA+1,1:-1,1])[sp.sep-1:], marker='+', label='H')
        if has_impurities:
            plt.plot(yyct[sp.sep-1:], (sp.na_eirene[-1,1:-1,2]/sp.na[JXA+1,1:-1,3:].sum(axis=1))[sp.sep-1:], marker='+', label=impurity)
        plt.yscale('log')
        plt.legend(fontsize=7)
        plt.ylabel('Compression ($n^0_{ot}/n^{i}_{omp}$)')
        plt.xlabel(yyclabel)
    except Exception as e:
        print(traceback.format_exc())

    try:
        ax = fig.add_subplot(gs[0, 1])
        plt.plot(yyc, (sp.na_eirene[:,:,0]+2*sp.dmb2[:,:,0])[JXA+1,1:-1], label='$n^0_{H\ omp}$', ls='--', c='C0', lw=2)
        plt.plot(yyct, (sp.na_eirene[:,:,0]+2*sp.dmb2[:,:,0])[-1,1:-1], label='$n^0_{H\ ot}$', ls='-', c='C0', lw=2)
        plt.plot(yyc, sp.na[JXA+1,1:-1,1], label='$n^i_{H\ omp}$', ls='--', c='C0', lw=1)
        plt.plot(yyct, sp.na[-1,1:-1,1], label='$n^i_{H\ ot}$', ls='-', c='C0', lw=1)

        if has_impurities:
            imp_ion_indices = np.where((np.array(nuclei)==impurity)*(sp.za >= 1))[0]
            imp_dens = np.sum(sp.na[:,:,imp_ion_indices], axis=2)
            plt.plot(yyc, sp.na_eirene[JXA+1,1:-1,2], label=f'$n^0_{{{impurity}\\ omp}}$', ls='--', c='C1', lw=2)
            plt.plot(yyct, sp.na_eirene[-1,1:-1,2],   label=f'$n^0_{{{impurity}\\ ot}}$', ls='-', c='C1', lw=2)
            plt.plot(yyc, imp_dens[JXA+1,1:-1],       label=f'$n^i_{{{impurity}\\ omp}}$', ls='--', c='C1', lw=1)
            plt.plot(yyct, imp_dens[-1,1:-1],         label=f'$n^i_{{{impurity}\\ ot}}$', ls='-', c='C1', lw=1)
        plt.yscale('log')
        plt.legend(fontsize=7, loc='upper left')
        plt.ylabel('Density [m$^{-3}$]')
        plt.xlabel(yyclabel);
    except Exception as e:
        print(traceback.format_exc())

    try:
        if has_impurities:
            ax = fig.add_subplot(gs[0, 2])
            imp_ion_indices = np.where((np.array(nuclei)==impurity)*(sp.za >= 1))[0]
            imp_dens = np.sum(sp.na[:,:,imp_ion_indices], axis=2)
            imp_percent = 100*imp_dens/sp.ne
            plt.plot(yyc, imp_percent[JXA+1,1:-1], label='Midplane')
            plt.plot(yyct, imp_percent[-1,1:-1], label='Outer target')
            plt.legend(fontsize=7)
            plt.ylabel(f'Impurity concentration ($n^i_{{{impurity}}}/n_e$) [%]')
            plt.xlabel(yyclabel)
    except Exception as e:
        print(traceback.format_exc())
            
    try:
        ax = fig.add_subplot(gs[1, 0])
        plot2d_b2(so, (sp.line_radiation.sum(axis=2)/sp.vol)[1:-1,1:-1].T/1e6, ax=ax, label=f'Line radiation [MW/m$^3$]\n({sp.line_radiation.sum()/1e6:.3g} MW)', cmap='viridis')
        common(ax)
    except Exception as e:
        print(traceback.format_exc())

    try:
        ax = fig.add_subplot(gs[1, 1])
        plot2d_b2(so, (sp.rqahe[1:-1,1:-1].sum(axis=2)/sp.vol[1:-1,1:-1]).T/1e6, ax=ax, label=f'Electron cooling [MW/m$^3$]\n({sp.rqahe[1:-1,1:-1].sum()/1e6:.3g} MW)', cmap='viridis')
        common(ax)
    except Exception as e:
        print(traceback.format_exc())
        
    try:
        ax = fig.add_subplot(gs[1, 2])
        plot2d_b2(so, so.data('te'), ax=ax, label=r'$T_e$ [eV]', cmap='viridis', lb=1, ub=sp.te[sp.masks.sol].max())
        common(ax)
    except Exception as e:
        print(traceback.format_exc())

    try:
        ax = fig.add_subplot(gs[2, 0])
        if has_impurities:
            b2_impurity_indices = [i for i in range(len(sp.species_names)) if impurity in sp.species_names[i] and sp.za[i] > 0]
            plot2d_b2(so, np.sum(sp.na[:,:,b2_impurity_indices], axis=2)[1:-1,1:-1].T, ax=ax, label=f'{impurity}'+r' ion dens. [m$^{-3}$]', cmap='viridis')
            common(ax)
        else:
            ax.axis('off')
    except Exception as e:
        ax.axis('off')
        print(traceback.format_exc())

    try:
        ax = fig.add_subplot(gs[2, 1])
        if has_impurities:
            imp_ion_indices = np.where((np.array(nuclei)==impurity)*(sp.za >= 1))[0]
            imp_dens = np.sum(sp.na[:,:,imp_ion_indices], axis=2)
            imp_percent = 100*imp_dens/sp.ne
            plot2d_b2(so, imp_percent[1:-1,1:-1].T, ax=ax, label=f'{impurity} ion %', scale='linear', cmap='viridis')
            common(ax)
        else:
            ax.axis('off')
    except Exception as e:
        ax.axis('off')
        print(traceback.format_exc())

    try:
        ax = fig.add_subplot(gs[2, 2])
        plot2d_b2(so, so.data('ne'), ax=ax, label=r'$n_e$ [m$^{-3}$]', cmap='viridis')
        common(ax)
    except Exception as e:
        print(traceback.format_exc())

    plt.tight_layout()

def page9():
    fig = plt.figure(figsize=(8.5,11), dpi=dpi)
    gs = fig.add_gridspec(4,2)
    ms = 3
    lw = 1.5
    kw = {'marker':'.', 'markersize': ms, 'lw': lw}

    def legend(ax, **args):
        #ax.legend(fontsize=6, bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        ax.legend(fontsize=6, **args)

    try:
        fhtx = spq.fhtx
        fhp = spq.fhp.sum(3) # sum over species
        fhm = spq.fhm.sum(3) # sum over species
        myfht = spq.fhe + spq.fhi + spq.fhj + fhp + fhm + spq.fnt
        fhpx = fhp[:,:,0]
        fhmx = fhm[:,:,0]
        fntx = spq.fnt[:,:,0]
        myfhtx = myfht[:,:,0]

        ax = fig.add_subplot(gs[0, 0])
        ax.set_title('Inboard X-point')
        ax.set_ylabel('$q_{\parallel}$ [GW/m$^2$]')
        ax.set_xlabel('$R-R_{sep}$ [mm]')
        ix = sp.ixp-1
        iy = sp.sep
        sxpol = spq.sx*spq.qz[:,:,1]
        qfac = -(1/sxpol*spq.bb[:,:,3]/np.abs(spq.bb[:,:,0]))[ix,iy:-1]/1e9
        ax.plot(yyct[iy-1:], qfac*spq.fhtx[ix,iy:-1], label=f'fht ({-np.sum(spq.fhtx[ix,iy:-1])/1e6:.2f} MW)', c='k', **kw)
        ax.plot(yyct[iy-1:], qfac*myfhtx[ix,iy:-1], label=f'e+i+j+p+m+nt ({-np.sum(myfhtx[ix,iy:-1])/1e6:.2f} MW)', c='k', alpha=0.5, **kw)
        ax.plot(yyct[iy-1:], qfac*spq.fhex[ix,iy:-1], label=f'fhe (e thermal, {-np.sum(spq.fhex[ix,iy:-1])/1e6:.2f} MW)', **kw)
        ax.plot(yyct[iy-1:], qfac*spq.fhix[ix,iy:-1], label=f'fhi (i thermal, {-np.sum(spq.fhix[ix,iy:-1])/1e6:.2f} MW)', **kw)
        ax.plot(yyct[iy-1:], qfac*spq.fhjx[ix,iy:-1], label=f'fhj (electrostatic, {-np.sum(spq.fhjx[ix,iy:-1])/1e6:.2f} MW)', **kw)
        ax.plot(yyct[iy-1:], qfac*fhpx[ix,iy:-1], label=f'fhp (ionization, {-np.sum(fhpx[ix,iy:-1])/1e6:.2f} MW)', **kw)
        ax.plot(yyct[iy-1:], qfac*fhmx[ix,iy:-1], label=f'fhm (kinetic, {-np.sum(fhmx[ix,iy:-1])/1e6:.2f} MW)', **kw)
        ax.plot(yyct[iy-1:], qfac*fntx[ix,iy:-1], label=f'fnt (viscous, {-np.sum(fntx[ix,iy:-1])/1e6:.2f} MW)', **kw)
        legend(ax)

        ax = fig.add_subplot(gs[0, 1])
        ax.set_title('Outboard X-point')
        ax.set_ylabel('$q_{\parallel}$ [GW/m$^2$]')
        ax.set_xlabel('$R-R_{sep}$ [mm]')
        ix = sp.oxp
        iy = sp.sep
        sxpol = spq.sx*spq.qz[:,:,1]
        qfac = (1/sxpol*spq.bb[:,:,3]/np.abs(spq.bb[:,:,0]))[ix,iy:-1]/1e9
        ax.plot(yyct[iy-1:], qfac*spq.fhtx[ix,iy:-1], label=f'fht ({np.sum(spq.fhtx[ix,iy:-1])/1e6:.2f} MW)', c='k', **kw)
        ax.plot(yyct[iy-1:], qfac*myfhtx[ix,iy:-1], label=f'e+i+j+p+m+nt ({np.sum(myfhtx[ix,iy:-1])/1e6:.2f} MW)', c='k', alpha=0.5, **kw)
        ax.plot(yyct[iy-1:], qfac*spq.fhex[ix,iy:-1], label=f'fhe (e thermal, {np.sum(spq.fhex[ix,iy:-1])/1e6:.2f} MW)', **kw)
        ax.plot(yyct[iy-1:], qfac*spq.fhix[ix,iy:-1], label=f'fhi (i thermal, {np.sum(spq.fhix[ix,iy:-1])/1e6:.2f} MW)', **kw)
        ax.plot(yyct[iy-1:], qfac*spq.fhjx[ix,iy:-1], label=f'fhj (electrostatic, {np.sum(spq.fhjx[ix,iy:-1])/1e6:.2f} MW)', **kw)
        ax.plot(yyct[iy-1:], qfac*fhpx[ix,iy:-1], label=f'fhp (ionization, {np.sum(fhpx[ix,iy:-1])/1e6:.2f} MW)', **kw)
        ax.plot(yyct[iy-1:], qfac*fhmx[ix,iy:-1], label=f'fhm (kinetic, {np.sum(fhmx[ix,iy:-1])/1e6:.2f} MW)', **kw)
        ax.plot(yyct[iy-1:], qfac*fntx[ix,iy:-1], label=f'fnt (viscous, {np.sum(fntx[ix,iy:-1])/1e6:.2f} MW)', **kw)
        legend(ax)

        ax = fig.add_subplot(gs[1, 0])
        ax.set_title('Inner target')
        ax.set_ylabel('$q_{surf}$ [MW/m$^2$]')
        ax.set_xlabel('$R-R_{sep}$ [mm]')
        ix = 1
        iy = 1
        qfac = -1/spq.sx[ix,iy:-1]/1e6
        ax.plot(yyct[iy-1:], qfac*spq.fhtx[ix,iy:-1], label=f'fht ({-np.sum(spq.fhtx[ix,iy:-1])/1e6:.2f} MW)', c='k', **kw)
        ax.plot(yyct[iy-1:], qfac*myfhtx[ix,iy:-1], label=f'e+i+j+p+m+nt ({-np.sum(myfhtx[ix,iy:-1])/1e6:.2f} MW)', c='k', alpha=0.5, **kw)
        ax.plot(yyct[iy-1:], qfac*spq.fhex[ix,iy:-1], label=f'fhe (e thermal, {-np.sum(spq.fhex[ix,iy:-1])/1e6:.2f} MW)', **kw)
        ax.plot(yyct[iy-1:], qfac*spq.fhix[ix,iy:-1], label=f'fhi (i thermal, {-np.sum(spq.fhix[ix,iy:-1])/1e6:.2f} MW)', **kw)
        ax.plot(yyct[iy-1:], qfac*spq.fhjx[ix,iy:-1], label=f'fhj (electrostatic, {-np.sum(spq.fhjx[ix,iy:-1])/1e6:.2f} MW)', **kw)
        ax.plot(yyct[iy-1:], qfac*fhpx[ix,iy:-1], label=f'fhp (ionization, {-np.sum(fhpx[ix,iy:-1])/1e6:.2f} MW)', **kw)
        ax.plot(yyct[iy-1:], qfac*fhmx[ix,iy:-1], label=f'fhm (kinetic, {-np.sum(fhmx[ix,iy:-1])/1e6:.2f} MW)', **kw)
        ax.plot(yyct[iy-1:], qfac*fntx[ix,iy:-1], label=f'fnt (viscous), {-np.sum(fntx[ix,iy:-1])/1e6:.2f} MW', **kw)
        legend(ax)

        ax = fig.add_subplot(gs[1, 1])
        ax.set_title('Outer target')
        ax.set_ylabel('$q_{surf}$ [MW/m$^2$]')
        ax.set_xlabel('$R-R_{sep}$ [mm]')
        ix = -1
        iy = 1
        qfac = 1/spq.sx[ix,iy:-1]/1e6
        ax.plot(yyct[iy-1:], qfac*spq.fhtx[ix,iy:-1], label=f'fht ({np.sum(spq.fhtx[ix,iy:-1])/1e6:.2f} MW)', c='k', **kw)
        ax.plot(yyct[iy-1:], qfac*myfhtx[ix,iy:-1], label=f'e+i+j+p+m+nt ({np.sum(myfhtx[ix,iy:-1])/1e6:.2f} MW)', c='k', alpha=0.5, **kw)
        ax.plot(yyct[iy-1:], qfac*spq.fhex[ix,iy:-1], label=f'fhe (e thermal, {np.sum(spq.fhex[ix,iy:-1])/1e6:.2f} MW)', **kw)
        ax.plot(yyct[iy-1:], qfac*spq.fhix[ix,iy:-1], label=f'fhi (i thermal, {np.sum(spq.fhix[ix,iy:-1])/1e6:.2f} MW)', **kw)
        ax.plot(yyct[iy-1:], qfac*spq.fhjx[ix,iy:-1], label=f'fhj (electrostatic, {np.sum(spq.fhjx[ix,iy:-1])/1e6:.2f} MW)', **kw)
        ax.plot(yyct[iy-1:], qfac*fhpx[ix,iy:-1], label=f'fhp (ionization, {np.sum(fhpx[ix,iy:-1])/1e6:.2f} MW)', **kw)
        ax.plot(yyct[iy-1:], qfac*fhmx[ix,iy:-1], label=f'fhm (kinetic, {np.sum(fhmx[ix,iy:-1])/1e6:.2f} MW)', **kw)
        ax.plot(yyct[iy-1:], qfac*fntx[ix,iy:-1], label=f'fnt (viscous), {np.sum(fntx[ix,iy:-1])/1e6:.2f} MW', **kw)
        legend(ax)

        ax = fig.add_subplot(gs[2, 0])
        ax.set_title('X-point to inner target deltas');
        ax.set_ylabel('$P_{pol\ targ}-P_{pol\ Xpt}$ [MW]')
        ax.set_xlabel('iy')
        ixs = 1
        ixe = sp.ixp-1
        # Define y-values as variables and add their sums to labels
        x = range(1, sp.ny - 1)
        delta_fhtx = (spq.fhtx[ixe, 1:-1] - spq.fhtx[ixs, 1:-1]) / 1e6
        ax.plot(x, delta_fhtx, label=f'fht ({delta_fhtx.sum():+.2f} MW)', c='k', **kw)
        delta_myfhtx = (myfhtx[ixe, 1:-1] - myfhtx[ixs, 1:-1]) / 1e6
        ax.plot(x, delta_myfhtx, label=f'e+i+j+p+m+nt ({delta_myfhtx.sum():+.2f} MW)', c='k', alpha=0.5, **kw)
        delta_fhex = (spq.fhex[ixe, 1:-1] - spq.fhex[ixs, 1:-1]) / 1e6
        ax.plot(x, delta_fhex, label=f'fhe (e thermal, {delta_fhex.sum():+.2f} MW)', **kw)
        delta_fhix = (spq.fhix[ixe, 1:-1] - spq.fhix[ixs, 1:-1]) / 1e6
        ax.plot(x, delta_fhix, label=f'fhi (i thermal, {delta_fhix.sum():+.2f} MW)', **kw)
        delta_fhjx = (spq.fhjx[ixe, 1:-1] - spq.fhjx[ixs, 1:-1]) / 1e6
        ax.plot(x, delta_fhjx, label=f'fhj (electrostatic, {delta_fhjx.sum():+.2f} MW)', **kw)
        delta_fhpx = (fhpx[ixe, 1:-1] - fhpx[ixs, 1:-1]) / 1e6
        ax.plot(x, delta_fhpx, label=f'fhp (ionization, {delta_fhpx.sum():+.2f} MW)', **kw)
        delta_fhmx = (fhmx[ixe, 1:-1] - fhmx[ixs, 1:-1]) / 1e6
        ax.plot(x, delta_fhmx, label=f'fhm (kinetic, {delta_fhmx.sum():+.2f} MW)', **kw)
        delta_fntx = (fntx[ixe, 1:-1] - fntx[ixs, 1:-1]) / 1e6
        ax.plot(x, delta_fntx, label=f'fnt (viscous, {delta_fntx.sum():+.2f} MW)', **kw)
        ax.axvline(sp.sep, c='k', ls=':', lw=lw)
        ax.axhline(0, c='k', ls=':', lw=lw)
        legend(ax, loc='upper left')

        ax = fig.add_subplot(gs[2, 1])
        ax.set_title('X-point to outer target deltas');
        ax.set_ylabel('$P_{pol\ targ}-P_{pol\ Xpt}$ [MW]')
        ax.set_xlabel('iy')
        ixs = sp.oxp + 1
        ixe = sp.nx-1
        # Define y-values as variables and add their sums to labels
        x = range(1, sp.ny - 1)
        delta_fhtx = (spq.fhtx[ixe, 1:-1] - spq.fhtx[ixs, 1:-1]) / 1e6
        ax.plot(x, delta_fhtx, label=f'fht ({delta_fhtx.sum():+.2f} MW)', c='k', **kw)
        delta_myfhtx = (myfhtx[ixe, 1:-1] - myfhtx[ixs, 1:-1]) / 1e6
        ax.plot(x, delta_myfhtx, label=f'e+i+j+p+m+nt ({delta_myfhtx.sum():+.2f} MW)', c='k', alpha=0.5, **kw)
        delta_fhex = (spq.fhex[ixe, 1:-1] - spq.fhex[ixs, 1:-1]) / 1e6
        ax.plot(x, delta_fhex, label=f'fhe (e thermal, {delta_fhex.sum():+.2f} MW)', **kw)
        delta_fhix = (spq.fhix[ixe, 1:-1] - spq.fhix[ixs, 1:-1]) / 1e6
        ax.plot(x, delta_fhix, label=f'fhi (i thermal, {delta_fhix.sum():+.2f} MW)', **kw)
        delta_fhjx = (spq.fhjx[ixe, 1:-1] - spq.fhjx[ixs, 1:-1]) / 1e6
        ax.plot(x, delta_fhjx, label=f'fhj (electrostatic, {delta_fhjx.sum():+.2f} MW)', **kw)
        delta_fhpx = (fhpx[ixe, 1:-1] - fhpx[ixs, 1:-1]) / 1e6
        ax.plot(x, delta_fhpx, label=f'fhp (ionization, {delta_fhpx.sum():+.2f} MW)', **kw)
        delta_fhmx = (fhmx[ixe, 1:-1] - fhmx[ixs, 1:-1]) / 1e6
        ax.plot(x, delta_fhmx, label=f'fhm (kinetic, {delta_fhmx.sum():+.2f} MW)', **kw)
        delta_fntx = (fntx[ixe, 1:-1] - fntx[ixs, 1:-1]) / 1e6
        ax.plot(x, delta_fntx, label=f'fnt (viscous, {delta_fntx.sum():+.2f} MW)', **kw)
        ax.axvline(sp.sep, c='k', ls=':', lw=lw)
        ax.axhline(0, c='k', ls=':', lw=lw)
        legend(ax, loc='upper left')

        ax = fig.add_subplot(gs[3, 0])
        ax.set_title('Inner divertor leg')
        ax.set_ylabel('$P_{pol\ SOL+PFR}$ [MW]')
        ax.set_xlabel('ix')
        ixs = 1
        ixe = sp.ixp
        iys = 1
        iye = -1
        ax.plot(range(ixs, ixe), -spq.fhtx[ixs:ixe, iys:iye].sum(1)/1e6, c='k', label='fht', **kw)
        ax.plot(range(ixs, ixe), -myfhtx[ixs:ixe, iys:iye].sum(1)/1e6, label='e+i+j+p+m+nt', c='k', alpha=0.5, **kw)
        ax.plot(range(ixs, ixe), -spq.fhex[ixs:ixe, iys:iye].sum(1)/1e6, label='fhe (e thermal)', **kw)
        ax.plot(range(ixs, ixe), -spq.fhix[ixs:ixe, iys:iye].sum(1)/1e6, label='fhi (i thermal)', **kw)
        ax.plot(range(ixs, ixe), -spq.fhjx[ixs:ixe, iys:iye].sum(1)/1e6, label='fhj (electrostatic)', **kw)
        ax.plot(range(ixs, ixe), -fhpx[ixs:ixe, iys:iye].sum(1)/1e6, label='fhp (ionization)', **kw)
        ax.plot(range(ixs, ixe), -fhmx[ixs:ixe, iys:iye].sum(1)/1e6, label='fhm (kinetic)', **kw)
        ax.plot(range(ixs, ixe), -fntx[ixs:ixe, iys:iye].sum(1)/1e6, label='fnt (viscous)', **kw)
        legend(ax, loc='upper center')

        ax = fig.add_subplot(gs[3, 1])
        ax.set_title('Outer divertor leg');
        ax.set_ylabel('$P_{pol\ SOL+PFR}$ [MW]')
        ax.set_xlabel('ix')
        ixs = sp.oxp+1
        ixe = sp.nx
        iys = 1
        iye = -1
        ax.plot(range(ixs, ixe), spq.fhtx[ixs:ixe, iys:iye].sum(1)/1e6, c='k', label='fht', **kw)
        ax.plot(range(ixs, ixe), myfhtx[ixs:ixe, iys:iye].sum(1)/1e6, label='e+i+j+p+m+nt', c='k', alpha=0.5, **kw)
        ax.plot(range(ixs, ixe), spq.fhex[ixs:ixe, iys:iye].sum(1)/1e6, label='fhe (e thermal)', **kw)
        ax.plot(range(ixs, ixe), spq.fhix[ixs:ixe, iys:iye].sum(1)/1e6, label='fhi (i thermal)', **kw)
        ax.plot(range(ixs, ixe), spq.fhjx[ixs:ixe, iys:iye].sum(1)/1e6, label='fhj (electrostatic)', **kw)
        ax.plot(range(ixs, ixe), fhpx[ixs:ixe, iys:iye].sum(1)/1e6, label='fhp (ionization)', **kw)
        ax.plot(range(ixs, ixe), fhmx[ixs:ixe, iys:iye].sum(1)/1e6, label='fhm (kinetic)', **kw)
        ax.plot(range(ixs, ixe), fntx[ixs:ixe, iys:iye].sum(1)/1e6, label='fnt (viscous)', **kw)
        legend(ax, loc='upper center')
    except Exception as e:
        print(e)

    plt.tight_layout()


def page10():
    fig = plt.figure(figsize=(8.5,11), dpi=dpi)
    gs = fig.add_gridspec(3,2)
    plotargs = {'marker': '+'}

    iy = sp.sep-1
    ixp = spq.leftcut[0]+1
    oxp = spq.rightcut[0]+1

    # Resources on heat flux components:
    # https://repo.tok.ipp.cas.cz/kripner/solps-doc/-/blob/master/doc/Katka's_notes_on_SOLPS-ITER.md#solpspy-heat-flux-components
    # https://repo.tok.ipp.cas.cz/kripner/solps-doc/-/blob/master/doc/Comparing_experiment_to_model.md#parallel-heat-flux-q_parallel
    # Appendix A https://iopscience.iop.org/article/10.1088/1741-4326/ac9917/pdf
    # Manual page 471
    Ptot = spq.fhtx
    Ptotei = spq.fhex + spq.fhix
    qpol = Ptot/sp.sx # wrong, sx is different here
    qpar = qpol/sp.pitch
    qsurf = Ptot/sp.sx # not sure if sx is east or west, sx[0,:] is all zeros, and difference between sx[1,:] and sx[2,:] is ~1%
    qsurfei = Ptotei/sp.sx
    

    def fit_exp(yyc, q, color):
        def expfun(x, A, lamda_inv):
            """
            Exponential function for curve fitting.
            
            Parameters:
            x: Array of x values
            A: Amplitude of exponential
            lamda_inv: Inverse decay length (1/lambda) (necessary for curve_fit to work nicely)
            """
            return A*np.exp(-x*lamda_inv)
        
        istart = np.argmax(q)
        mmstart = yyc[istart]
        mmend = mmstart+1 # fit up to 1 mm away from peak
        iend = argnear(yyc, mmend)
        # if heat flux profile is non-monotonic, only fit part before it starts rising again
        for i in range(istart+1, iend+1):
            if q[i] > q[i-1]:
                iend = i-1
                mmend = yyc[iend]
                break
        qofit, _ = curve_fit(expfun, yyc[istart:iend+1], q[istart:iend+1], p0=[np.max(q),1], bounds=(0, np.inf))
        lamda = 1/qofit[1] # lamda_q in mm
        ax.plot(yyc[istart:iend+1], expfun(np.array(yyc[istart:iend+1]), *qofit), c=color, alpha=0.4, label=f'$\lambda_q$={lamda:.3g} mm')
        return lamda
    
    # plot and fit qpar e+i at X-point
    ax = fig.add_subplot(gs[0, 0])
    ax.set_title('$q_{||\ e+i}$ at X-point [GW/m$^2$]')
    ptot = spq.fhex+spq.fhix
    sxpol = sp.sx*sp.qc
    q = ptot/sxpol*wbbl[:,:,3]/np.abs(wbbl[:,:,0])
    ax = plt.gca()
    ax.plot(yyct[sp.sep-1:], q[oxp,sp.sep:-1]/1e9, label=f'Outboard ({np.sum(ptot[oxp,sp.sep:-1])/1e6:.3g} MW)', **plotargs)
    ax.plot(yyct[sp.sep-1:], -q[ixp,sp.sep:-1]/1e9, label=f'Inboard ({np.sum(-ptot[ixp,sp.sep:-1])/1e6:.3g} MW)', **plotargs)
    try:
        lamda_xo = fit_exp(yyct[sp.sep-1:], q[oxp,sp.sep:-1]/1e9, 'C0')
    except Exception as e:
        print(traceback.format_exc())
        lamda_xo = None
    try:
        lamda_xi = fit_exp(yyct[sp.sep-1:], -q[ixp,sp.sep:-1]/1e9, 'C1')
    except Exception as e:
        print(traceback.format_exc())
        lamda_xi = None
    ax.set_yscale('log')
    ax.set_xlabel(yyclabel)
    ax.legend(fontsize=8, loc='upper right')

    # plot and fit qpar tot at X-point
    ax = fig.add_subplot(gs[0, 1])
    ax.set_title('$q_{||\ tot}$ at X-point [GW/m$^2$]')
    ptot = spq.fhtx
    sxpol = sp.sx*sp.qc
    q = ptot/sxpol*wbbl[:,:,3]/np.abs(wbbl[:,:,0])
    ax = plt.gca()
    ax.plot(yyct[sp.sep-1:], q[oxp,sp.sep:-1]/1e9, label=f'Outboard ({np.sum(ptot[oxp,sp.sep:-1])/1e6:.3g} MW)', **plotargs)
    ax.plot(yyct[sp.sep-1:], -q[ixp,sp.sep:-1]/1e9, label=f'Inboard ({np.sum(-ptot[ixp,sp.sep:-1])/1e6:.3g} MW)', **plotargs)
    try:
        lamda_xo = fit_exp(yyct[sp.sep-1:], q[oxp,sp.sep:-1]/1e9, 'C0')
    except Exception as e:
        print(traceback.format_exc())
        lamda_xo = None
    try:
        lamda_xi = fit_exp(yyct[sp.sep-1:], -q[ixp,sp.sep:-1]/1e9, 'C1')
    except Exception as e:
        print(traceback.format_exc())
        lamda_xi = None
    ax.set_yscale('log')
    ax.set_xlabel(yyclabel)
    ax.legend(fontsize=8, loc='upper right')

    # plot and fit qsurf total
    ax = fig.add_subplot(gs[2, 0])
    ax.set_title('$q_{surf\ tot}$ [MW/m$^2$]')

    #-radial profile of qpar below entrance to the outer leg
    psurfo = Ptot[-1,:]
    qsurfo = psurfo[1:-1]/sp.sx[-1,1:-1]
    #-radial profile of qpar below entrance to the inner leg
    psurfi = -Ptot[1,:]
    qsurfi = psurfi[1:-1]/sp.sx[1,1:-1]

    def qEich(rho, q0, S, lqi, qbg, rho_0):
            rho = rho - rho_0
            # lqi is inverse lamda_q
            return q0/2*np.exp((S*lqi/2)**2-rho*lqi)*erfc(S*lqi/2-rho/S)+qbg

    if lamda_xo != None:
        lqoGuess = lamda_xo
    else:
        lqoGuess = 1 # mm
    if lamda_xi != None:
        lqiGuess = lamda_xi
    else:
        lqiGuess = 1 # mm
    # lamda_q fits
    yycm = yyct/1000
    yycml = np.linspace(yycm[0], yycm[-1], 500)
    bounds = ([0,0,0,0,yycm[0]], [np.inf,np.inf,np.inf,np.inf,yycm[-1]])
    # Fit outer
    oguess = (np.max(qsurfo)-np.min(qsurfo), lqoGuess/1000/2, 1000/lqoGuess, np.min(qsurfo), 0)
    try:
        qsurfol = interpolate.interp1d(yycm, qsurfo)(yycml)
        qsofit, _ = curve_fit(qEich, yycml, qsurfol, p0=oguess, bounds=bounds)
        lqeo, So = 1000/qsofit[2], qsofit[1]*1000 # lamda_q and S in mm
    except Exception as e:
        print('qsurf outer fit failed:', e)
        qsofit = None
    iguess = (np.max(qsurfi)-np.min(qsurfi), lqiGuess/1000/2, 1000/lqiGuess, np.min(qsurfi), 0)
    # Fit inner
    try:
        qsurfil = interpolate.interp1d(yycm, qsurfi)(yycml)
        qsifit, _ = curve_fit(qEich, yycml, qsurfil, p0=iguess, bounds=bounds)
        lqei, Si = 1000/qsifit[2], qsifit[1]*1000 # lamda_q and S in mm 
    except Exception as e:
        print('qsurf inner fit failed:', e)
        qsifit = None

    ax.plot(yyct, qsurf[-1,1:-1]/1e6, c='C0', label=f'Outer target ({np.sum(Ptot[-1,1:-1])/1e6:.3g} MW)', **plotargs)

    def qEichGivenLq(rho, q0, S, qbg, rho_0):
        rho = rho - rho_0
        # lqi is inverse lamda_q
        lqi = 1000/lamda_xo
        return q0/2*np.exp((S*lqi/2)**2-rho*lqi)*erfc(S*lqi/2-rho/S)+qbg

    # Fit outer using Lq found at X-point
    oguess2 = (np.max(qsurfo)-np.min(qsurfo), lamda_xo/1000/2, np.min(qsurfo), 0) # try max for rho_0
    bounds2 = ([0,0,0,yycm[0]], [np.inf,np.inf,np.inf,yycm[-1]])
    qsofit2, pcov = curve_fit(qEichGivenLq, yycm, qsurfo, p0=oguess2, bounds=bounds2)
    ax.plot(yycml*1000, qEichGivenLq(yycml, *qsofit2)/1e6, c='#d180ff', alpha=0.4,
                label=r'Fit: $\lambda_q$ = %.3g mm, $S$ = %.3g mm' % (lamda_xo, qsofit2[1]*1000))

    if np.any(qsofit):
        ax.plot(yycml*1000, qEich(yycml, *qsofit)/1e6, c='C0', alpha=0.4,
                label=r'Fit: $\lambda_q$ = %.3g mm, $S$ = %.3g mm' % (lqeo, So))

    # ax.set_yscale('log')
    ax.legend(fontsize=7, loc='upper right')
    ax.set_xlabel(yyclabel)
    ax.set_xlim([-1, yyc[-1]])
    plt.tight_layout()

    # 2D map of S and LQ fitting
    ax = fig.add_subplot(gs[2, 1])
    lqgiven = None
    Sgiven = None
    def qEichGivenLqAndS(rho, q0, qbg, rho_0):
        rho = rho - rho_0
        # lqi is inverse lamda_q
        lqi = 1000/lqgiven
        S = Sgiven/1000
        return q0/2*np.exp((S*lqi/2)**2-rho*lqi)*erfc(S*lqi/2-rho/S)+qbg

    lqarr = np.linspace(lqeo/5, lqeo*5, 30)
    Sarr = np.linspace(So/5, So*5, 30)
    SSEs = np.zeros((len(lqarr), len(Sarr)))
    for ilq, lqgiven in enumerate(lqarr):
        for iS, Sgiven in enumerate(Sarr):
            try:
                qsofit3, pcov = curve_fit(qEichGivenLqAndS, yycm, qsurfo, p0=[np.max(qsurfo), np.min(qsurfo), 0], bounds=([0,0,yycm[0]], [np.inf,np.inf,yycm[-1]]))
                SSEs[ilq, iS] = np.sum((qEichGivenLqAndS(yycm, *qsofit3)-qsurfo)**2)
            except Exception as e:
                print(e)
                SSEs[ilq, iS] = np.nan
    # Normalize
    SSEs /= np.nanmin(SSEs)

    plt.contourf(lqarr, Sarr, SSEs.T, cmap='jet', levels=list(range(1,51)), vmin=1, vmax=50, extend='max') # contour needs .T!
    plt.colorbar(label='SSE/SSE$_{min}$')
    contour10 = plt.contour(lqarr, Sarr, SSEs.T, colors='white', levels=[10]) # contour needs .T!
    plt.clabel(contour10, inline=True, fmt={10: '10'}, colors='white')

    # Show Lq, S values at which SSE was calculated
    Lqm, Sm = np.meshgrid(lqarr, Sarr)
    plt.scatter(Lqm.flatten(), Sm.flatten(), marker='.', s=1, c='k', alpha=0.1)

    plt.scatter(lqeo, So, c='red', marker='.', label='Fit $\lambda_q$ and $S$')
    plt.scatter(lamda_xo, qsofit2[1]*1000, c='red', marker='x', label='Fit using Xpt $\lambda_q$')
    lqminind, Sminind  = np.unravel_index(np.nanargmin(SSEs), SSEs.shape)
    plt.scatter(lqarr[lqminind], Sarr[Sminind], c='lime', marker='*', label='minimum SSE')
    plt.legend(fontsize=7)
    plt.gca().set_facecolor('#890100')

    plt.xlabel('$\lambda_q$ [mm]')
    plt.ylabel('S [mm]')
    plt.tight_layout()


def finishPage(pdf):
    try:
        plt.tight_layout()
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()
    except Exception as e:
        print(traceback.format_exc())
        
dpi = 150
cmapdefault = 'viridis'
scriptdir = '/home/sballinger/solps_utils'

if __name__ == '__main__':
    matplotlib.use('agg')
    
    if os.path.isfile('plots.pdf'):
        os.remove('plots.pdf')
    if os.path.isfile('b2fplasmf') and os.path.getmtime('b2fplasmf') < os.path.getmtime('b2fstate'):
        print('Removing b2fplasmf file older than b2fstate')
        os.remove('b2fplasmf')

    geqdsk = glob.glob('../baserun/*.geq*')[0]
    so = aurora.solps_case(b2fstate_path="./b2fstate", b2fgmtry_path="../baserun/b2fgmtry", geqdsk=geqdsk)
    sp = solpspy.SolpsData('.', mag='lsn')
    spq = quixote.SolpsData('.')

    if os.path.isfile('b2fplasmf'):
        if not np.all(sp.ne[1:-1,1:-1].T/so.ne == 1.0):
            raise Exception('Non-matching b2fstate and b2fplasmf files. Do b2run b2uf or remove b2fplasmf.')

    # Modified from Aurora
    # now obtain also the simple poloidal grid slice near the midplane (LFS and HFS)
    # These are commonly used for SOLPS analysis, using the JXA and JXI indices (which we re-compute here)
    Z_core = so.data("cz")[so.unit_r:2*so.unit_r, so.unit_p:3*so.unit_p]
    R_core = so.data("cr")[so.unit_r:2*so.unit_r, so.unit_p:3*so.unit_p]

    # find indices of poloidal grid nearest to Z=0 in the innermost radial shell
    midplane_LFS_idx = np.argmin(np.abs(Z_core[0, R_core[0, :] > so.geqdsk["RMAXIS"]]))
    midplane_HFS_idx = np.argmin(np.abs(Z_core[0, R_core[0, :] < so.geqdsk["RMAXIS"]]))

    # convert to indices on so.data('cz') and so.data('cr')
    JXI = so.unit_p + np.arange(Z_core.shape[1])[R_core[0, :] < so.geqdsk["RMAXIS"]][midplane_HFS_idx] # HFS_mid_pol_idx
    JXA = so.unit_p + np.arange(Z_core.shape[1])[R_core[0, :] > so.geqdsk["RMAXIS"]][midplane_LFS_idx] # LFS_mid_pol_idx
    print('JXA', JXA)

    nuclei = [re.findall(r'[a-zA-Z]+', s)[0] for s in sp.species_names] # e.g. ['H', 'H', 'H', 'H'] for DT
    nuclei_nice = [] # e.g. ['D', 'D', 'T', 'T'] for DT
    for i in range(len(sp.zn)):
        if sp.zn[i] == 1:
            spec = 'H'
            if sp.am[i] == 2:
                spec = 'D'
            elif sp.am[i] == 3:
                spec = 'T'
        else:
            spec = nuclei[i]
        nuclei_nice.append(spec)
    nucset = list(set(nuclei))
    nuclei_nice_set = list(set(nuclei_nice))
    if nucset == ['H']:
        has_impurities = False
        impurity = None
    else:
        has_impurities = True
        impurity = [i for i in nucset if i != 'H'][0]
    species_nice = [sp.species_names[i].replace(nuclei[i], nuclei_nice[i]) for i in range(len(sp.species_names))] # e.g. ['D', D+', 'T', 'T+'] for DT
    C0123 = [(0.12156863, 0.46666667, 0.70588235), (1.00000000, 0.49803922, 0.00000000), (0.17254902, 0.62745098, 0.17254902), (0.83921569, 0.15294118, 0.15686275)]

    has_wlld = False
    if has_impurities: # multi-species
        # Do wlld analysis and load results
        # Descriptions of terms here: https://solps.pages.tok.ipp.cas.cz/solps-doc/Energy_fluxes_deep_dive/#divertor-target-loads-from-wlld
        try:
            shutil.rmtree('b2pl.exe.dir')
        except Exception as e:
            pass
        try:
            result = subprocess.run('echo "0100 wlld" | b2plot', shell=True, executable='/bin/tcsh', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            cols = ['x','ne','Te','Ti','Wtot','Wpart','Wrad','Wpls','Wneut','Wheat','Wpot','Whtpl','Wptpl','xMP','Rclfc','Zclfc','Wpar_xpt','WWpar','WWtrg','Lcnnt','Lcnnx']
            datai = np.loadtxt('b2pl.exe.dir/ld_tg_i.dat').T
            it = {cols[i]: datai[i] for i in range(len(cols))}
            datao = np.loadtxt('b2pl.exe.dir/ld_tg_o.dat').T
            ot = {cols[i]: datao[i] for i in range(len(cols))}
            has_wlld = True
        except Exception as e:
            print('wlld analysis failed:', e)

    yyc = np.loadtxt('dsa')[1:-1]*1000 #rmid(so, so.geqdsk)*1000
    yyclabel = '$R_{omp}-R_{sep}$ [mm]'
    isp = sp.species_names.index('H+') # ion species index

    i = sp.oxp
    ig = 0
    f = interpolate.interp1d(so.geqdsk.get2D('PSIRZ_NORM', sp.cr[i,1:-1], sp.cz[i,1:-1], interp='cubic'), yyc, fill_value='extrapolate', kind='cubic')
    yyct = f(so.geqdsk.get2D('PSIRZ_NORM', sp.cr[i+1,1:-1], sp.cz[i+1,1:-1], interp='cubic'))
    wbbl = quixote.tools.parsers.b2('../baserun/b2fgmtry', 'wbbl', spq.bb.shape)
    savefile = 'plots.pdf'
    
    def doPage(pageFunction, pdf):
        print(pageFunction.__name__)
        pageFunction()
        finishPage(pdf)
    
    with PdfPages(savefile) as pdf:
        doPage(page0, pdf)
        doPage(page1, pdf)
        doPage(page2, pdf)
        doPage(page3, pdf)
        doPage(page4, pdf)
        doPage(page5, pdf)
        doPage(page6, pdf)
        doPage(page7, pdf)
        doPage(page8, pdf)
        doPage(page9, pdf)
        doPage(page10, pdf)
