import matplotlib
matplotlib.use('Agg')
import numpy as np
import os
import sys
from copy import deepcopy as dcopy
from glob import glob
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from scipy.special import erf

import get_dspec as dspec
import cand_plots as cplt
import pulse_plt as sp_plt



#############################
## PULSE OBJECTS / SORTING ##
#############################

class Pulse(object):
    def __init__(self, beam, line):
        dm, sigma, time, sample, dfact = line.split()
        dfact.rstrip('\n')
        self.beam   = beam
        self.dm     = np.float(dm)
        self.sigma  = np.float(sigma)
        self.time   = np.float(time)
        self.sample = np.int(sample)
        self.dfact  = np.int(dfact)
        self.line   = line
        self.ra     = 0
        self.dec    = 0
        self.cnum   = -1

        self.ndupes = 0
        self.dupe_dms    = []
        self.dupe_sigmas = []
        self.dupe_times  = []
        self.dupe_dfacts = []

    def add_dupe(self, cand):
        self.ndupes += 1
        self.dupe_dms.append(cand.dm)
        self.dupe_sigmas.append(cand.sigma)
        self.dupe_times.append(cand.time)
        self.dupe_dfacts.append(cand.dfact)


    def __str__(self):
        return self.line

    def __repr__(self):
        out_str = "Pulse(T=%.3f, "    %(self.time) +\
                        "DM=%.2f, "   %(self.dm)   +\
                        "beam=%d, "   %(self.beam) +\
                        "width=%d, "  %(self.dfact) +\
                        "sigma=%.2f) " %(self.sigma)
        return out_str


class Cand(object):
    def __init__(self, pulse):
        self.beam   = pulse.beam
        self.dm     = pulse.dm
        self.sigma  = pulse.sigma
        self.time   = pulse.time
        self.sample = pulse.sample
        self.dfact  = pulse.dfact
        self.line   = pulse.line
        self.ra     = pulse.ra
        self.dec    = pulse.dec
        self.cnum   = pulse.cnum

        self.n_group = pulse.ndupes

        self.ndupes = 0
        self.dupe_dms    = []
        self.dupe_sigmas = []
        self.dupe_times  = []
        self.dupe_dfacts = []
        self.dupe_beams  = []

    def add_dupe(self, cand):
        self.ndupes += 1
        self.dupe_dms.append(cand.dm)
        self.dupe_sigmas.append(cand.sigma)
        self.dupe_times.append(cand.time)
        self.dupe_dfacts.append(cand.dfact)
        self.dupe_beams.append(cand.beam)


    def __str__(self):
        return self.line

    def __repr__(self):
        out_str = "Cand(T=%.3f, "    %(self.time) +\
                        "DM=%.2f, "   %(self.dm)   +\
                        "beam=%d, "   %(self.beam) +\
                        "width=%d, "  %(self.dfact) +\
                        "sigma=%.2f) " %(self.sigma)
        return out_str


def cmp_pulse(p1, p2):
    cmp_val = cmp(p1.sample, p2.sample)
    if cmp_val == 0:
        cmp_val = cmp(p1.dm, p2.dm)
        if cmp_val == 0:
            cmp_val = -cmp(p1.sigma, p2.sigma)
        else: pass
    else: pass
    return cmp_val


def cmp_pulse2(p1, p2):
    cmp_val = cmp(p1.sample, p2.sample)
    if cmp_val == 0:
        cmp_val = -cmp(p1.sigma, p2.sigma)
        if cmp_val == 0:
            cmp_val = cmp(p1.beam, p2.beam)
        else: pass
    else: pass
    return cmp_val


def cmp_snr(p1, p2):
    return -cmp(p1.sigma, p2.sigma)


#############################
## READING / WRITING CANDS ##
#############################

def collect_sp_files(spdir, outdir, outname):
    """
    Read in all *singlepulse files from spdir
    and write them to a single file called 
    outname in outdir
    """
    spglob = "%s/*singlepulse" %(spdir)
    spfiles = glob(spglob)

    outfile = "%s/%s" %(outdir, outname)
    if os.path.isfile(outfile):
        print("%s already exists!" %(outfile))
        return
    else:
        pass

    if len(spfiles) == 0:
        print("No single pulse files found in %s" %(spdir))
        return
    else:
        pass

    with open(outfile, 'w') as fout:
        for ii, spfile in enumerate(spfiles):
            with open(spfile, 'r') as fin:
                for jj, line in enumerate(fin):
                    if "#" in line[0] and ii > 0:
                        continue
                    else:
                        fout.write(line)
    
    return


def collect_sp_files_beams(bnums, topdir, outdir):
    """
    Collect all *singlepulse files for each 
    beam and write them to a single file (per beam)
    in outdir

    topdir is where the beam directories live

    Assume that *singlepulse files are in 

      topdir/beamXXXX/single_pulse/*singlepulse

    Will write beamXXX.cands in outdir
    """
    for bnum in bnums:
        spdir = "%s/beam%04d/single_pulse" %(topdir,  bnum)
        outfile = "beam%04d.cands" %(bnum)
        collect_sp_files(spdir, outdir, outfile)

    return


def cands_from_file(infile, beam):
    """
    Return list of Pulse objects from the file "infile".
    Assumes this file is in the same format as a PRESTO
    *.singlepulse file.
    """
    plist = []
    for line in file(infile):
        if line[0] == "#":
            continue
        else: pass
        plist.append(Pulse(beam, line))
    return plist


def cands_from_many_files(beam_nums, cands_dir):
    """
    Assumes all the files have the name format "beam%03d.cands"
    and are located in cands_dir
    """
    plist = []
    for bnum in beam_nums:
        print(bnum)
        sp_file = "%s/beam%04d.cands" %(cands_dir, bnum)
        # check that file exists
        if not os.path.isfile(sp_file):
            sp_file = "%s/beam%03d.cands" %(cands_dir, bnum)
            if not os.path.isfile(sp_file):
                print("File not found: %s" %sp_file)
                continue
            else: pass
        else: pass
        plist += cands_from_file(sp_file, bnum)

    return plist


def combine_sp_files(sp_files, outfile):
    sp_files.sort()
    fout = open(outfile, 'w')
    for ii, spfile in enumerate(sp_files):
        jj = 0
        for line in file(spfile):
            if (ii==0 and jj==0) or (jj>0):
                fout.write(line)
            else:
                pass
            jj += 1
    fout.close()
    return


def combine_sp_files_multi(beam_list, beams_dir):
    # Make sure beams_dir exists
    if not os.path.exists(beams_dir):
        print("%s does not exist!  Exiting..." %beams_dir)
        return
    else: pass

    for ii, bnum in enumerate(beam_list):
        # Check to see if dat dir exists 
        dat_dir = "%s/beam%04d" %(beams_dir, bnum)
        if not os.path.exists(dat_dir):
            print("dat dir not found: %s" %dat_dir)
            continue
        else: pass

        sp_base = "beam%04d" %bnum
        sp_files = glob("%s/%s*singlepulse" %(dat_dir, sp_base))
        sp_files.sort()

        outfile = "%s/%s.cands" %(dat_dir, sp_base)
        combine_sp_files(sp_files, outfile)
    return


def pulse_groups_to_cands(bnums, pdir):
    """
    Read the (large size) pulse group *npy files 
    and convert them to the (smaller) cand class 
    objects for inter-beam sifting

    assume pulse group npy files have form:

      pdir/beam%04d_dmsift.npy

    """
    pfiles = [ "%s/beam%04d_dmsift.npy" %(pdir, bnum) for bnum in bnums ]
    candlist = []
    for pfile in pfiles:
        if not os.path.isfile(pfile):
            continue
        else: 
            pass
        pdat = np.load(pfile)
        for pulse in pdat:
            candlist.append(Cand(pulse))
        del pdat

    return candlist


######################
## BEAM COORDINATES ##
######################


def beam_coords2offsets(beams, coords, ra0, dec0):
    ras = np.array([ cc[0] for cc in coords ])
    decs = np.array([ cc[1] for cc in coords ])

    dras = ras - ra0
    ddecs = decs - dec0
    dras *= np.cos(dec0 * np.pi / 180.0)

    coords = np.vstack( (dras, ddecs) ).T
    return beams, coords



def get_beam_coords(beamfile):
    beams = []
    coords = []
    for line in file(beamfile):
        if line[0] in ["#", "", "\n"]:
            continue
        else: pass

        bnum = int(line.split()[0])
        ra_deg  = float(line.split()[1])
        dec_deg = float(line.split()[2])

        beams.append(bnum)
        coords.append( (ra_deg, dec_deg) )

    beams = np.array(beams)
    coords = np.array(coords)
    return beams, coords


def get_coords_from_beamfile(beamfile):
    bb = np.load(beamfile)
    beams = np.arange(len(bb), dtype='int')
    ra_strs  = bb[:, 1]
    dec_strs = bb[:, 2]
    dec_strs = np.array([ dd.replace('.', ':', 2) for dd in dec_strs ])
    cc = SkyCoord(ra=ra_strs, dec=dec_strs, unit=(u.hour, u.deg))
    ra_deg = cc.ra.deg
    dec_deg = cc.dec.deg

    coords = zip(ra_deg, dec_deg)
    return beams, coords


#def get_skycoords(beamfile):
#    bb = np.load(beamfile)
#    beams = np.arange(len(bb), dtype='int')
#    ra_strs  = bb[:, 1]
#    dec_strs = bb[:, 2]
#    dec_strs = np.array([ dd.replace('.', ':', 2) for dd in dec_strs ])
#    coords = SkyCoord(ra=ra_strs, dec=dec_strs, unit=(u.hour, u.deg))
#    return coords


def get_skycoords(beams):
    ra_strs  = beams[:, 1]
    dec_strs = beams[:, 2]
    dec_strs = np.array([ dd.replace('.', ':', 2) for dd in dec_strs ])
    coords = SkyCoord(ra=ra_strs, dec=dec_strs, unit=(u.hour, u.deg))
    return coords


def get_beamlist(candfile):
    beams = []
    with open(candfile, 'r') as fin:
        for line in fin:
            if line[0] in ["#", '\n', ' ']:
                continue
            else:
                pass
            beams.append(int(line.split()[1]))
    beams = np.array(beams)
    return beams



#####################
## PARSE CANDLISTS ##
#####################

def attrarr(obj_list, attr):
    if hasattr(obj_list[0], attr):
        out_arr = np.array([ getattr(bb, attr) for bb in obj_list ])
        return out_arr
    else:
        print("List has no attribute \"%s\" " %attr)
        return



#####################
## CAND SELECTION  ##
#####################


def dm_delay(freq_lo, freq_hi, DM):
    return 4.15e3 * (freq_lo**-2.0 - freq_hi**-2.0) * DM


def avg_dm_delay(freqs, DM, ref='last'):
    """
    ref = last -> delay rel to last freq
    ref = inf  -> delay rel to infinity
    """
    if ref=='last':
        tdms = np.array([ dm_delay(ff, freqs[-1], DM) for ff in freqs ])
    else:
        tdms = np.array([ dm_delay(ff, np.inf, DM) for ff in freqs ])

    return np.mean(tdms)


def avg_dm_tsc_delay(freqs, DM, W=0.01, tsc1=1.3, ref='last'):
    """
    ref = last -> delay rel to last freq
    ref = inf  -> delay rel to infinity
    """
    if ref=='last':
        tdms = np.array([ dm_delay(ff, freqs[-1], DM) for ff in freqs ])
    else:
        tdms = np.array([ dm_delay(ff, np.inf, DM) for ff in freqs ])

    tscs = tsc1 * (freqs / 1000.0)**-4.0
    dtdm = 8.3e3 * (np.abs(freqs[1]-freqs[0]) / freqs**3.0) * DM
    wf = W / np.sqrt( W**2.0 + tscs**2.0 + dtdm**2.0)

    return np.sum((tdms - tscs) * wf) / np.sum(wf)


def snr_fac(bw_mhz, Wms, fghz, ddm):
    xparam = 6.91e-3 * np.abs(ddm) * bw_mhz / (Wms * fghz**3.0)
    fac = 0.5 * np.sqrt(np.pi) * erf(xparam) / xparam
    return fac


def filter_dm_cands(cands, W, freqs, twin):
    cands_snr = sorted(cands, cmp=cmp_snr)
    N = len(cands_snr)
    # enumerate snr sorted order
    for ii, cc in enumerate(cands_snr):
        cc.cnum = ii
    # get params in orig order
    toas = attrarr(cands_snr, 'time')
    dms  = attrarr(cands_snr, 'dm')
    # Get first remaining index list
    xx_r = np.arange(len(cands_snr))

    # Begin loop to filter candidates
    out_cands = []
    max_iter = len(cands_snr)
    counter = 0
    print("-- BEFORE SIFTING = %d" %(N))
    sys.stdout.flush()

    ppc_list = range(0, 100, 5)
    while( len(xx_r) and counter < max_iter ):
        #print(xx_r[0])
        counter += 1
        ppc = int((1 - len(xx_r) / float(N)) * 100)
        if ppc in ppc_list:
            ppc_list.remove(ppc)
            print(ppc)

        cand = cands_snr[xx_r[0]]
        dm0 = cand.dm
        t0  = cand.time
        idx = cand.cnum

        #print t0

        tt_x = toas[xx_r]
        dm_x = dms[xx_r]

        yy = np.where( np.abs(tt_x - t0) < twin )[0]
        tt_y = tt_x[yy]
        dm_y = dm_x[yy]

        #tdm = np.array([ avg_dm_delay(freqs, dmii - dm0) for dmii in dm_y ])
        #yy_bool = np.greater(tt_y, t0 - tdm - W) & np.less(tt_y, t0 - tdm + W)
        tdm = np.array([ avg_dm_delay(freqs, dmii - dm0) for dmii in dm_y ])
        tdm2 = np.array([ avg_dm_tsc_delay(freqs, dmii - dm0, W=0.001, tsc1=1.5) for dmii in dm_y ])
        yy_bool = (np.greater(tt_y, t0 - tdm - W)  | \
                   np.greater(tt_y, t0 - tdm2 - W) | \
                   np.greater(tt_y, t0 - W))         \
                   & \
                  (np.less(tt_y, t0 - tdm + W)     | \
                   np.less(tt_y, t0 - tdm2 + W)    | \
                   np.less(tt_y, t0 + W))
        xx_c = xx_r[yy[yy_bool]]

        # Add dupe cands to cand
        for xval in xx_c:
            dcand = cands_snr[xval]
            cand.add_dupe(dcand)

        # Remove cand indices from xx_r
        xx_r = np.setdiff1d(xx_r, xx_c)
        #print(xx_r)

        # Add cand to out list
        out_cands.append( cand )

    print("-- AFTER SIFTING = %d" %(len(out_cands)))
    sys.stdout.flush()

    return out_cands



def filter_dm_cands_dumb(cands, freqs, twin):
    cands_snr = sorted(cands, cmp=cmp_snr)
    N = len(cands_snr)
    # enumerate snr sorted order
    for ii, cc in enumerate(cands_snr):
        cc.cnum = ii
    # get params in orig order
    toas = attrarr(cands_snr, 'time')
    dms  = attrarr(cands_snr, 'dm')
    # Get first remaining index list
    xx_r = np.arange(len(cands_snr))

    # Begin loop to filter candidates
    out_cands = []
    max_iter = len(cands_snr)
    counter = 0
    ppc_list = range(0, 100, 5)
    while( len(xx_r) and counter < max_iter ):
        counter += 1
        ppc = int((1 - len(xx_r) / float(N)) * 100)
        if ppc in ppc_list:
            ppc_list.remove(ppc)
            print(ppc)

        cand = cands_snr[xx_r[0]]
        dm0 = cand.dm
        t0  = cand.time
        idx = cand.cnum

        tt_x = toas[xx_r]
        dm_x = dms[xx_r]

        yy = np.where( np.abs(tt_x - t0) < twin )[0]
        tt_y = tt_x[yy]
        dm_y = dm_x[yy]

        xx_c = xx_r[yy]

        # Add dupe cands to cand
        for xval in xx_c:
            dcand = cands_snr[xval]
            cand.add_dupe(dcand)

        # Remove cand indices from xx_r
        xx_r = np.setdiff1d(xx_r, xx_c)

        # Add cand to out list
        out_cands.append( cand )

    return out_cands


def filter_cands_from_cands_dumb(cands, ref_cands, ddm, freqs):
    N0 = len(ref_cands)
    N = len(cands)
    
    # Remaining index list
    xx_r = np.arange(len(cands))
    
    # target cands / dms
    toas = attrarr(cands, 'time')
    dms  = attrarr(cands, 'dm')
    snrs = attrarr(cands, 'sigma')

    tdm = dm_delay(freqs[0], freqs[-1], ddm)
    print("tdm = %.2f s" %tdm)
    
    # Begin loop to filter candidates
    ppc_list = range(0, 100, 5)
    for ii, cc in enumerate(ref_cands):
        ppc = int((ii / float(N0)) * 100)
        if ppc in ppc_list:
            print(ppc)
            ppc_list.remove(ppc)
        
        dm0  = cc.dm
        t0   = cc.time
        snr0 = cc.sigma

        tt_x  = toas[xx_r]
        dm_x  = dms[xx_r]
        snr_x = snrs[xx_r]

        yy = np.where( (np.abs(tt_x - t0) < tdm) & (np.abs(dm_x - dm0) < ddm) )[0] 
    
        if len(yy) == 0:
            continue
        else: pass

        if np.max(snr_x[yy]) < snr0:
            xx_c = xx_r[yy]

            # Remove cand indices from xx_r
            xx_r = np.setdiff1d(xx_r, xx_c)
        else:
            pass

    out_cands = [cands[xxi] for xxi in xx_r] 
    return out_cands


def filter_cands_from_cands(cands, ref_cands, W, freqs, twin):
    N0 = len(ref_cands)
    N = len(cands)
    
    # Remaining index list
    xx_r = np.arange(len(cands))
    
    # target cands / dms
    toas = attrarr(cands, 'time')
    dms  = attrarr(cands, 'dm')
    
    # Begin loop to filter candidates
    ppc_list = range(0, 100, 5)
    for ii, cc in enumerate(ref_cands):
        ppc = int((ii / float(N0)) * 100)
        if ppc in ppc_list:
            print(ppc)
            ppc_list.remove(ppc)
        
        dm0 = cc.dm
        t0  = cc.time

        tt_x = toas[xx_r]
        dm_x = dms[xx_r]

        yy = np.where( np.abs(tt_x - t0) < twin )[0]
        tt_y = tt_x[yy]
        dm_y = dm_x[yy]

        tdm = np.array([ avg_dm_delay(freqs, dmii - dm0) for dmii in dm_y ])
        tdm2 = np.array([ avg_dm_tsc_delay(freqs, dmii - dm0, W=0.001, tsc1=1.5) for dmii in dm_y ])
        yy_bool = (np.greater(tt_y, t0 - tdm - W)  | \
                   np.greater(tt_y, t0 - tdm2 - W) | \
                   np.greater(tt_y, t0 - W))         \
                   & \
                  (np.less(tt_y, t0 - tdm + W)     | \
                   np.less(tt_y, t0 - tdm2 + W)    | \
                   np.less(tt_y, t0 + W))
        xx_c = xx_r[yy[yy_bool]]

        # Remove cand indices from xx_r
        xx_r = np.setdiff1d(xx_r, xx_c)

    out_cands = [cands[xxi] for xxi in xx_r] 
    return out_cands


def plot_dupes(cands):
    for cc in cands:
        if len(cc.dupe_dms):
            plt.plot(cc.dupe_times, cc.dupe_dms, 'o')
        else:
            pass
    return


def filter_known_psr(cands, dm0, ddm, P, tpk0, window=(0.25, 0.5)):
    tt = attrarr(cands, 'time')
    dm = attrarr(cands, 'dm')
    phase = (tt - tpk0 + P/2.0) % P
    phase_lo = P/2.0 - window[0]
    phase_hi = P/2.0 + window[1]
    condition = (phase > phase_lo) & (phase < phase_hi) & (np.abs(dm - dm0) < ddm) 
    xx_in = np.where( condition )[0]
    xx_out = np.setdiff1d(np.arange(len(tt)), xx_in)
    mag_cands  = [cands[xxi] for xxi in xx_in]
    cands = [cands[xxo] for xxo in xx_out]
    return mag_cands, cands


def check_phases(bnums, P):
    fstr = "beam%04d_sift_cands.npy"
    pks = []
    for bb in bnums:
        infile = fstr % bb
        cands = np.load(infile)
        tt = attrarr(cands, 'time')
        n, b = np.histogram(tt % P, bins=np.linspace(0, P, 100))
        pk_phase = b[np.argmax(n)]
        pks.append(pk_phase)
    return np.array(pks)


def make_plot(cands, max_w=1000, fig_title=None, out_name=None):
    dt = 10 #ms
    tt = attrarr(cands, 'time')
    dm = attrarr(cands, 'dm')
    snr = attrarr(cands, 'sigma')
    ww = attrarr(cands, 'dfact') * dt

    xx = np.where( ww <= max_w )[0]
    
    bpad = 0.1 
    wpad = 0.01
    hpad = 0.05 

    h1 = h4 = 0.5

    w2 = w3 = (1.0 - 2 * bpad - wpad) / 2.0
    h2 = h3 = 1.0 - 2 * bpad - hpad - h1

    w4 = 0.15 
    w1 = 1.0 - 2 * bpad - w4 - wpad
    
    x1 = bpad
    y1 = y4 = bpad

    y2 = y3 = 1.0 - bpad - h2
    x2 = bpad
    x3 = bpad + w2 + wpad 
    x4 = bpad + w1 + wpad 

    if out_name is not None:
        plt.ioff()
    else: pass
    
    fig = plt.figure(figsize=(12, 8))

    ax1 = fig.add_axes([x1, y1, w1, h1])
    ax2 = fig.add_axes([x2, y2, w2, h2])
    ax3 = fig.add_axes([x3, y3, w3, h3])
    ax4 = fig.add_axes([x4, y4, w4, h4])

    print(w1, w4, h2, h3)
    print(x1, x4)

    # T - DM - SNR PLOT
    dmlim = (-100, 5500)
    ax1.scatter(tt[xx], dm[xx], c=ww[xx], s=2 * snr[xx]**2.0, 
                cmap=plt.cm.magma)
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("DM")
    ax1.set_xlim(-100, 7000)
    ax1.set_ylim(dmlim[0], dmlim[1])

    # CANDS SNR HIST
    snrbins = np.arange(int(np.min(snr[xx])), np.max(snr[xx])+1, 1)
    n2, b2 = np.histogram(snr[xx], bins=snrbins)
    ax2.bar(b2[:-1], n2, b2[1]-b2[0], log=True, color='0.4')
    ax2.set_ylim(0.1)
    ax2.set_xlabel("SNR")
    ax2.set_ylabel("Number")

    # DM - SNR PLOT
    cax = ax3.scatter(dm[xx], snr[xx], c=ww[xx], s=30, 
                      cmap=plt.cm.magma)
    ax3.set_xlim(dmlim[0], dmlim[1])
    ax3.yaxis.tick_right()
    cbar = fig.colorbar(cax, fraction=0.10)
    
    # CANDS SNR HIST
    dm_bins = np.arange(0, 5500, 200)
    n4, b4 = np.histogram(dm[xx], bins=dm_bins)
    ax4.barh(b4[:-1], n4, b4[1]-b4[0], color='0.4')
    ax4.set_ylim(0.1)
    ax4.set_ylim(dmlim[0], dmlim[1])
    ax4.set_xlabel("Number")
    ax4.yaxis.tick_right()
    ax4.set_yticklabels([])


    if fig_title is not None:
        fig.suptitle(fig_title)

    if out_name is not None:
        plt.savefig(out_name, dpi=100, bbox_inches='tight')
        plt.close()
        plt.ion()
    else: 
        plt.show()

    return


def make_mult_plots(bdir, bnums, dm0, ddm, P, tpk, max_w=100, 
                    make_plots=True, write_file=False):
    """
    Make plots and write cands

    (max_w only for plot, all cands written to file)
    """
    bstr = "beam%04d_sift_cands.npy"
    for bb in bnums:
        infile = bdir + "/" + bstr %bb
        all_cands = np.load(infile)
        in_cands, cands = filter_known_psr(all_cands, dm0, ddm, P, tpk)
        wmag_plot = "beam%04d_mag.png" %bb
        nomag_plot = "beam%04d_cand.png" %bb
        if make_plots:
            make_plot(in_cands, max_w=max_w, fig_title="beam%04d" %bb, out_name=wmag_plot)
            make_plot(cands, max_w=max_w, fig_title="beam%04d" %bb, out_name=nomag_plot)
        else: pass

        if write_file:
            out_txt = "beam%04d_filter.cands" %bb
            write_cands(out_txt, cands)
        else: pass
    return


def write_cands(outfile, cands):
    """
    Output cand info sorted by sigma
    """
    sigmas = attrarr(cands, 'sigma')
    idx = np.argsort(sigmas)
    idx = idx[::-1]

    fout = open(outfile, 'w')
    hdr = '#{:<10}{:^10}{:^10}{:^15}{:^10}{:^10}{:^10}{:^10}'.format( \
        '  Beam', 'DM', 'Sigma', 'Time', 'Sample', 'Downfact', 'Nhits', 'Group')
    fout.write(hdr + '\n')
    for xx in idx:
        out_str = '{:^10d}{:>10.2f}{:>10.2f}{:>12.3f}{:>10d}{:>10d}{:>10d}{:>10d}'.format( \
            cands[xx].beam, cands[xx].dm, cands[xx].sigma, cands[xx].time, 
            cands[xx].sample, cands[xx].dfact, cands[xx].ndupes, xx)
        fout.write(out_str + '\n')
    fout.close()
    return


def write_cands_mult(bnum):
    for bb in bnum:
        infile  = "beam%04d_filter_cands.npy" %bb
        outfile = "beam%04d_filter.cands" %bb
        cands = np.load(infile)
        write_cands(outfile, cands)
    return


#################
#     MAIN      #
#################

if __name__ == "__main__":
    # START / STOP BEAM NUMBERS 
    start = int(sys.argv[1])
    stop  = int(sys.argv[2])

    print(start, stop)

    collect_sp      = 0
    per_beam_filter = 0
    all_beam_filter = 0
    prelim_cut      = 0
    calc_stats      = 0
    read_stats      = 1
    make_cut        = 1
    make_midx       = 1
    make_dspec      = 1
    make_bplot      = 1

    # BEAM FILE
    #beam_file = '/Users/rsw/work/fastvis/16A-329/tiling/gc_beamlist.npy'
    beam_file = '/hercules/results/rwharton/fastvis_gc/proc/gc_beamlist.npy'
    
    beam_locs = np.load(beam_file)
    beam_locs = beam_locs[ start : stop ]
    beam_coords = get_skycoords(beam_locs)
    bnums = np.arange(start, stop)

    # TOP DIR 
    top_dir = "/hercules/results/rwharton/fastvis_gc/proc/57519"
    #top_dir = "/hercules/results/rwharton/fastvis_gc/proc/sp_test"

    # SEARCH DIR
    search_dir = "%s/search" %(top_dir)
    
    # SINGLEPULSE DIRECTORY
    cands_dir = "%s/sp_cand_files" %(search_dir)
    
    # File paths
    cands_file = "%s/sift_beams_simple.txt" %(cands_dir)
    stats_file = "%s/sifted_stats.txt" %(cands_dir)
    top_file   = "%s/sifted_top.txt" %(cands_dir)
    pulse_groups_npy = "%s/pgroups.npy" %(cands_dir)

    # OBS PARAMETERS
    nchan = 512
    lofreq = 1977.0
    chan_wid = 4.0
    freqs = np.arange(nchan) * chan_wid + lofreq

    # SEARCH/SIFT PARAMS
    P = 3.7661859
    t0 = -0.041
    
    W = 0.15
    twin = 2.0

    # BARY 2 TOPO
    b2t_coeff = np.array([ 9.99941374e-01,  -4.41724462e-04])

    # SEARCH / WIN / SHOW TIMES
    tsearch = 4.0
    tpad    = 0    # Use dfacts instead

    # CUT PARAMS
    snr_cut  = 6.0 
    midx_cut = np.sqrt(nchan) / snr_cut
  
    if collect_sp:
        collect_sp_files_beams(bnums, search_dir, cands_dir)

    if per_beam_filter:
        for bnum in bnums:
            print("\n== BEAM%04d ==" %(bnum))
            sp_file = "%s/beam%04d.cands" %(cands_dir, bnum)
            if not os.path.isfile(sp_file):
                print("FILE NOT FOUND: %s" %(sp_file))
                print("...skipping")
                continue
            else: 
                pass

            cands = cands_from_file(sp_file, bnum)

            # MAKE SMALL FOR TESTING
            #tt = attrarr(cands, 'time')
            #xxt = np.where( tt < 100 )[0]
            #cands = [cands[ii] for ii in xxt ]

            cands_sift_dm = filter_dm_cands(cands, W, freqs, twin)
            npy_path = "%s/beam%04d_dmsift.npy" %(cands_dir, bnum)
            np.save(npy_path, cands_sift_dm)

            del cands
            del cands_sift_dm

    if all_beam_filter:
        candlist = pulse_groups_to_cands(bnums, cands_dir)
        ocands1 = filter_dm_cands_dumb(candlist, freqs, 0.5)
        #ocands2 = filter_dm_cands(candlist, W, freqs, twin)
        write_cands(cands_file, ocands1)
        np.save(pulse_groups_npy, ocands1)
        
    if prelim_cut:
        candlist = dspec.cands_from_candfile(cands_file, b2t_coeff=b2t_coeff)
        dfacts = attrarr(candlist, 'dfact')
        xx = np.where( dfacts < 15 )[0]
        candlist = [ candlist[ii] for ii in xx ]
        cands_file = cands_file.strip('.txt') + "_cut.txt"
        candlist = dspec.write_select_cands(candlist, cands_file)

    if calc_stats:
        candlist = dspec.stats_from_candfile_parts(cands_file, top_dir, 
                                    tsearch, tpad,  b2t_coeff=b2t_coeff)
        stats_file = "%s/sifted_stats.txt" %(cands_dir)
        candlist = dspec.write_select_cands(candlist, stats_file)

    if read_stats:
        candlist = dspec.cands_from_candfile(stats_file, b2t_coeff=b2t_coeff)

    if make_cut:
        snrs  = attrarr(candlist, 'fit_snr')
        midxs = attrarr(candlist, 'fit_midx')
        beams = attrarr(candlist, 'beam')
        dms   = attrarr(candlist, 'dm')
        #xx = np.where( (snrs >= snr_cut) & (midxs < midx_cut) & (beams != 2) & (dms > 300))[0]
        xx = np.where( (snrs >= snr_cut) & (midxs < midx_cut) )[0]
        candlist = [ candlist[ii] for ii in xx ]
        candlist = dspec.write_select_cands(candlist, top_file)

    if make_dspec:
        candlist = dspec.cands_from_candfile(top_file, b2t_coeff=b2t_coeff)
        #dspec.plots_from_candlist_parts(candlist, top_dir, 
        #            tsearch, tpad, tshow, snr_cut=snr_cut, tsamp_avg=None)
        #candlist = candlist[:5]
        sp_plt.summary_pulse_from_candlist2(candlist, top_dir, tsearch, 
                                            use_dfact=False, tsamp_avg=1, nchan_avg=4)

    if make_bplot:
        candlist = dspec.cands_from_candfile(top_file, b2t_coeff=b2t_coeff)
        #candlist = candlist[:5]
        groups = attrarr(candlist, 'group')
        pcands = np.load(pulse_groups_npy)
        pcands = [ pcands[gg] for gg in groups ]  
        #cplt.make_many_snr_plots(beam_coords, bnums, pcands, 3.0, add_label=False)
        cplt.summary_beam_from_candlist(cands_dir, beam_coords, bnums, pcands, 3.0)
