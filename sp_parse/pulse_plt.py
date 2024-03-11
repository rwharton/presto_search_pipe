from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import sys
import pickle
from collections import namedtuple
from astropy.io import fits
import simple_dm as dedisp
import get_dspec as dspec
import glob

########################
###  SETUP PARAMS    ###
########################

try:        
    params = {
              "font.size": 18,
              "font.family": 'serif',
              "font.weight": "normal",
              "xtick.major.size": 7,
              "xtick.minor.size": 4,
              "ytick.major.size": 7,
              "ytick.minor.size": 4,
              "xtick.major.width": 2,
              "xtick.minor.width": 1,
              "ytick.major.width": 2,
              "ytick.minor.width": 1,
              "xtick.major.pad": 8,
              "xtick.minor.pad": 8,
              "ytick.major.pad": 8,
              "ytick.minor.pad": 8,
              "lines.linewidth": 2,
              "lines.markersize": 10,
              "axes.linewidth": 2,
              "legend.loc": "best",
              "text.usetex": False,
              "xtick.labelsize" : 14,
              "ytick.labelsize" : 14
             }
    matplotlib.rcParams.update(params)
except:
    print("Couldn't update rc parameters, ignoring...")
    pass

###########################

def cands_from_candfile(candfile, fixdm=None, b2t_coeff=None):
    """
    Read candidates from candfile

    fixdm forces DM to a given value

    b2t_coeff is polyfit coefficients for function to 
    convert barycentric times (reported in *singlepulse files)
    to topocentric times needed to find pulse in *FITS
    """
    dat = np.loadtxt(candfile)
    beams   = dat[:, 0].astype('int')
    dms_det = dat[:, 1]
    snrs    = dat[:, 2]
    times   = dat[:, 3]
    samps   = dat[:, 4].astype('int')
    dfacts  = dat[:, 5].astype('int')
    nhits   = dat[:, 6].astype('int')

    N = len(dat)

    if dat.shape[-1] == 10:
        tt_topos = dat[:, 7]
        fit_snr  = dat[:, 8]
        fit_midx = dat[:, 9]
    else:
        tt_topos = dspec.bary_to_topo(times, b2t_coeff)
        fit_snr = np.zeros(len(snrs))
        fit_midx = np.zeros(len(snrs))

    Cand = namedtuple('Cand', ['beam', 'dm', 'dm_det', 'snr', 'tt', 
                               'samp', 'dfact', 'nhits', 'tt_topo', 
                               'fit_snr', 'fit_midx'])
    candlist = []

    if fixdm is None:
        dms = dms_det
    else:
        dms = np.ones(len(dms_det)) * fixdm

    for ii in xrange(N):
        cc = Cand(beam=beams[ii], dm=dms[ii], dm_det=dms_det[ii], 
                  snr=snrs[ii], tt=times[ii], samp=samps[ii], 
                  dfact=dfacts[ii], nhits=nhits[ii], tt_topo=tt_topos[ii], 
                  fit_snr=fit_snr[ii], fit_midx=fit_midx[ii])
        candlist.append(cc)
    return candlist


def array_attr(clist, attr):
    if hasattr(clist[0], attr):
        return np.array([ getattr(cc, attr) for cc in clist ])
    else:
        print("Object has no attribute: %s" %attr)
        return

def log_output(text, logfile):
    if logfile is None:
        print(text)
        sys.stdout.flush()
    else:
        fout = open(logfile, 'a+')
        fout.write("%s\n" %text)
        fout.close()
    return



#####################
## DEBUG + TESTING ##
#####################


def dat_from_cand(cand, topdir, tsearch, tpad):
    time  = cand.tt_topo
    dm    = cand.dm
    beam  = cand.beam
    dfact = cand.dfact
    
    fitsfiles = glob.glob("%s/*beam%04d*fits" %(topdir, beam))
    fitsfiles.sort()
    off_starts = dspec.get_starts_for_fits(fitsfiles)
    
    otime, infile, xxf = dspec.get_fits_name_time(time, fitsfiles, off_starts)
    print(otime)

    tmid = otime
    tstart = tmid - tsearch * 0.25
    tt, freqs, dd = dspec.get_data(infile, tstart, tsearch)

    dt0 = tt[1] - tt[0]
    dt = tt[1] - tt[0]
    
    if tpad is not None:
        pk_win = (tmid - tpad, tmid + tpad)
    else:
        pk_win = None

    xx = np.where( np.sum(dd, axis=0) == 0 )[0]
    chan_weights = np.ones(len(freqs))
    chan_weights[xx] = 0

    tpk, ddm, spec, midx = dspec.get_ddm_spec_mod(tt, freqs, dd, dm, chan_weights=chan_weights, 
                                            pk_win=pk_win)
    
    #tpk0, ddm0, spec0, midx0 = dspec.get_ddm_spec_mod(tt, freqs, dd, 0.0, pk_win=pk_win, 
    #                                            chan_weights=chan_weights)
    
    return tpk, tt, ddm, spec, midx, dd


def summary_pulse(cand, topdir, tsearch, tpad=0, tsamp_avg=1, nchan_avg=1,  
                  use_dfact=False, outfile=None):
    ## PARSE CAND PARAMS ##
    time  = cand.tt_topo
    dm    = cand.dm
    beam  = cand.beam
    dfact = cand.dfact


    ## HACKY AND BAD
    if dm == 0:
        dm = 20.0 
    else:
        pass
   
    ## GET FITS ##
    fitsfiles = glob.glob("%s/*beam%04d*fits" %(topdir, beam))
    #fitsfiles = glob.glob("%s/part*/psrfits/*beam%04d*fits" %(topdir, beam))
    fitsfiles.sort()
    off_starts = dspec.get_starts_for_fits(fitsfiles)
    otime, infile, xxf = dspec.get_fits_name_time(time, fitsfiles, off_starts)
    #print(otime)

    tmid = otime
    tstart = tmid - tsearch * 0.25
    tt, freqs, dd = dspec.get_data(infile, tstart, tsearch)

    dt0 = tt[1] - tt[0]
    dt = tt[1] - tt[0]

    df0 = freqs[1] - freqs[0]

    if tpad <= 0:
        tpad = dfact * dt
    else:
        pass

    if use_dfact:
        max_chan_avg = 16
        tsamp_avg = dfact
        df_equiv  = dedisp.equiv_chan(tsamp_avg * dt, np.median(freqs), dm)
        nchan_avg_tmp = int( np.round( df_equiv / df0 ) )
        print(nchan_avg_tmp)
        nchan_avg = max( (nchan_avg_tmp, 1) )
        nchan_avg = min( (nchan_avg,  max_chan_avg) )
        print("tsamp_avg = %d" %tsamp_avg)
        print("nchan_avg = %d" %nchan_avg)

    if nchan_avg > 1 or tsamp_avg > 1:
        #dd = avg_chan(dd, nchan=nchan_avg)
        #freqs, dd = dedisp.dspec_avg_chan_dm(dd, freqs, freqs[-1], dt, 
        #                                     avg_chan=nchan_avg, dm0=DM)
        tt, freqs, dd = dedisp.dspec_avg_dm(dd, tt, freqs, freqs[-1], dt, 
                            avg_tsamp=tsamp_avg, avg_chan=nchan_avg, dm0=dm)
        dt = dt * tsamp_avg

    else: pass
    
    if tpad is not None:
        pk_win = (tmid - tpad, tmid + tpad)
    else:
        pk_win = None

    xx = np.where( np.sum(dd, axis=0) == 0 )[0]
    chan_weights = np.ones(len(freqs))
    chan_weights[xx] = 0
    fac = np.mean( chan_weights.astype('float') )

    ## CALC DDM + SPEC + MIDX ## 
    tpk, ddm, spec, midx = dspec.get_ddm_spec_mod(\
                tt, freqs, dd, dm, chan_weights=chan_weights, pk_win=pk_win)

    
    # OFF PULSE NORM ? #
    #yy = np.where( (np.abs(tt - tpk) > tsearch * 0.1))[0]
    #dd_sig = np.std(ddm[yy])
    #ddm /= dd_sig
    #dd /= dd_sig

    toffs = tt - tpk
    tshift = tt - tmid

    ## GET dDM - dt DATA ##
    delta_DM = 2 * dt / dedisp.dm_delay(freqs[0], freqs[-1], 1.0)
    print("DDM = %.2f" %delta_DM)
    dDMs = np.linspace(-5 * delta_DM, + 5 * delta_DM, 201)
    dms = dDMs + dm
    DM_dt = dedisp.dedisperse_dspec_multi(dd.T, dms, freqs, freqs[-1], dt)
    DM_dt /= fac
    ddt = -0.5 * dedisp.dm_delay(freqs[0], freqs[-1], dDMs)

    ##################
    ## SETUP FIGURE ##
    ##################
    if outfile is not None:
        plt.ioff()
    
    # Make Figure
    fig = plt.figure(figsize=(11, 8.5))

    ax11 = fig.add_subplot(2, 3, 1)
    ax12 = fig.add_subplot(2, 3, 2)
    ax13 = fig.add_subplot(2, 3, 3)
    ax21 = fig.add_subplot(2, 3, 4)
    ax22 = fig.add_subplot(2, 3, 5)
    ax23 = fig.add_subplot(2, 3, 6)

    # Some params 
    dspec_snr = 50 #* np.sqrt(nchan_avg * tsamp_avg)
    min_frac  = 0.25   # dspec_min = dspec_snr * min_frac
    #cmap = plt.cm.magma_r
    cmap = plt.cm.Greys
    im_kwargs = {'interpolation' : 'nearest', 
                 'aspect'        : 'auto', 
                 'origin'        : 'lower'}

    ######################
    ###   dDM - dt     ###
    ######################
    ax_dd = ax21

    dd_cmap = plt.cm.jet
    ext_dd = [toffs[0]-0.5*dt, toffs[-1] + 0.5*dt, dms[0]-0.5*delta_DM, dms[-1] +0.5*delta_DM]
    ax_dd.imshow(DM_dt, extent=ext_dd, cmap=dd_cmap, **im_kwargs) 
    ax_dd.plot(ddt, dms, ls='--', lw=2, color='r')
    ax_dd.axvline(x=0, ls='--', lw=2, c='r')

    dd_xlim = ( 1.2 * np.min(ddt), 1.5 * np.max(ddt) )
    ax_dd.set_xlim(dd_xlim)
    ax_dd.set_ylim(dms[0], dms[-1])
    ax_dd.set_xlabel("Time Offset (s)")
    ax_dd.set_ylabel("DM (pc/cc)")
    ax_dd.tick_params(axis='both', color='w')


    ######################
    ###  DM SNR CURVE  ###
    ######################
    ax_dm = ax11

    xx_pk = np.argmin( np.abs(toffs) ) 
    ax_dm.plot( dms, DM_dt[:, xx_pk], ls='-', c='k', lw=3)
    ax_dm.axvline(x = dm, ls='--', c='r')

    ax_dm.tick_params(labeltop='on', labelbottom='off')
    ax_dm.xaxis.set_label_position('top')

    ax_dm.tick_params(axis='x', labelsize=8)
    
    ax_dm.set_xlabel("DM (pc/cc)")
    ax_dm.set_ylabel("Norm S/N")
   
    
    ###############
    ###  PULSE  ###
    ###############
    tplim = (-0.5, 0.5)
    ax_ts = ax12
    ax_ts.plot(toffs, ddm, c='k')
    ax_ts.axhline(y=0, ls='-', c='k', alpha=0.2)
    ax_ts.axvline(x=0, ls='-', c='r', alpha=0.2)
    
    ax_ts.set_xlim(tplim)
    #ax_ts.set_ylabel("SNR")
    ax_ts.tick_params(labeltop='on', labelbottom='off',
                      labelleft='off')
    ax_ts.xaxis.set_label_position("top")
    ax_ts.set_xlabel("Time Offset (s)")
    pulse_ylim = ax_ts.get_ylim()

    ######################
    ### MOD - INDEX    ###
    ######################
    ax_mi = ax13
    snr_cut = 5.0

    ax_mi.plot(midx, ddm, ls='', marker='o', ms=8)
    xpk = np.where( tt == tpk )[0]
    ax_mi.plot([midx[xpk]], [ddm[xpk]], marker='o', c='r', ms=8)

    if chan_weights is None:
        midx_cut = np.sqrt(len(freqs) / nchan_avg) / snr_cut
    else:
        midx_cut = np.sqrt(np.sum(chan_weights) / nchan_avg) / snr_cut


    info_str = r"$m_{\rm I} = %.2f$" %midx[xpk] + "\n" +\
               r"${\rm SNR} = %.1f$" %ddm[xpk]  + "\n" +\
               r"${\rm DM} = %d$" %dm

    ax_mi.text(0.9, 0.9, info_str, ha='right', va='top', 
               transform=ax_mi.transAxes, fontsize=20)

    ax_mi.axvline(x=midx_cut, c='b', lw=3, alpha=0.2)
    ax_mi.axhline(y=snr_cut, c='b', lw=3, alpha=0.2)

    ax_mi.axvline(x=0, c='k', lw=2, ls='-', alpha=0.1)
    ax_mi.axhline(y=0, c='k', lw=2, ls='-', alpha=0.1)

    ax_mi.set_xlim(-1 , 20 * midx_cut)
    #ax_mi.set_ylim(-1)
    ax_mi.set_ylim(pulse_ylim)

    ax_mi.tick_params(labeltop='on', labelbottom='off', 
                      labelleft='off', labelright='on')
    ax_mi.xaxis.set_label_position("top")
    ax_mi.set_xlabel("Mod Index")
    ax_mi.yaxis.set_label_position("right")
    ax_mi.set_ylabel("SNR")

    ###################
    ### DSPEC SWEEP ###
    ###################
    ax_ds = ax23

    ## Some Useful Limits ##
    df1 = freqs[1] - freqs[0]
    dt1 = tt[1] - tt[0]

    tshow = 2.0
    tsweep = dedisp.dm_delay(freqs[0], freqs[-1], dm)
    tlim_min = max( -0.5 * tsweep, np.min(toffs) )
    tlim_max = min( 1.5 * tsweep, np.max(toffs) )
    tlim = (tlim_min, tlim_max)
    flim = (freqs[0] - 0.5 * df1, freqs[-1] + 0.5 * df1)

    if chan_weights is None:
        dd_sig = np.std(dd)
    else:
        yy = np.where(chan_weights)[0]
        dd_sig = np.std(dd[:, yy])

    print( np.sum(chan_weights) )
    print( len(freqs) )
    print(dd_sig)
    print(dd.shape)
    print( np.min(dd / dd_sig), np.max(dd / dd_sig ))

    im_ext = [tt[0] - tpk - 0.5 * dt1, tt[-1] - tpk + 0.5 * dt1, 
              freqs[0], freqs[-1]]

    ax_ds.imshow(dd.T / dd_sig, interpolation='nearest', origin='lower', 
                 aspect='auto', extent = im_ext, cmap=cmap, 
                 vmin = -1, vmax=3) #vmin= -1 * min_frac * dspec_snr, vmax = dspec_snr)

    # Add Dispersion Sweep
    offset = dt0 * 10
    dm_curve = dedisp.dm_delay(freqs, freqs[-1], dm)
    ax_ds.plot(dm_curve - offset, freqs, lw=3, c='r')
    #ax_ds.plot(dm_curve + offset, freqs, lw=2, c='w')

    # Set limits
    ax_ds.set_xlim(tlim)
    ax_ds.set_ylim(flim)

    ax_ds.tick_params(labelleft='off', labelright='on')
    ax_ds.yaxis.set_label_position('right')

    ax_ds.set_xlabel("Time Offset (s)")
    ax_ds.set_ylabel("Frequency (MHz)")



    ##################
    ### FREQ PHASE ###
    ##################
    ax_fp = ax22
    dd2 = dedisp.dedisperse_dspec(dd.T, dm, freqs, freqs[-1], dt1)

    ax_fp.imshow(dd2 / dd_sig, interpolation='nearest', origin='lower', 
                 aspect='auto', extent = im_ext, cmap=cmap, 
                 vmin=-1, vmax=3) #vmin= -1 * min_frac * dspec_snr, vmax = dspec_snr)

    ax_fp.set_xlim(tplim)
    ax_fp.set_xlabel("Time Offset (s)")
    #ax_fp.set_ylabel("Frequency (GHz)")
    ax_fp.set_yticklabels([])

    #plt.subplots_adjust(left=0.05, right=0.95)
    plt.subplots_adjust(hspace=0.05, wspace=0.05)

    ################
    ## FINISH UP  ##
    ################

    if outfile is not None:
        plt.savefig(outfile, dpi=150, bbox_inches='tight')
        plt.close()
        plt.ion()

    else:
        plt.show()

    return 
    #return tpk, tt, ddm, spec, midx, dd



def summary_pulse2(cand, topdir, tsearch, tpad=0, tsamp_avg=1, nchan_avg=1,  
                  use_dfact=False, outfile=None):
    ## PARSE CAND PARAMS ##
    time  = cand.tt_topo
    dm    = cand.dm
    beam  = cand.beam
    dfact = cand.dfact


    ## HACKY AND BAD
    if dm == 0:
        dm = 20.0 
    else:
        pass

   
    ## GET FITS ##
    #fitsfiles = glob.glob("%s/*beam%04d*fits" %(topdir, beam))
    fitsfiles = glob.glob("%s/part*/psrfits/*beam%04d*fits" %(topdir, beam))
    fitsfiles.sort()
    off_starts = dspec.get_starts_for_fits(fitsfiles)
    otime, infile, xxf = dspec.get_fits_name_time(time, fitsfiles, off_starts)
    #print(otime)

    tmid = otime
    tstart = tmid - tsearch * 0.25
    tt, freqs, dd = dspec.get_data(infile, tstart, tsearch)

    dt0 = tt[1] - tt[0]
    dt = tt[1] - tt[0]

    df0 = freqs[1] - freqs[0]
    
    max_chan_avg = 8
    tsamp_avg = dfact
    df_equiv  = dedisp.equiv_chan(tsamp_avg * dt, np.median(freqs), dm)
    nchan_avg_tmp = int( np.round( df_equiv / df0 ) )
    print(nchan_avg_tmp)
    nchan_avg = max( (nchan_avg_tmp, 1) )
    nchan_avg = min( (nchan_avg,  max_chan_avg) )
    print("tsamp_avg = %d" %tsamp_avg)
    print("nchan_avg = %d" %nchan_avg)
    
    if tpad <= 0:
        tpad = dfact * dt
    else:
        pass
    
    if tpad is not None:
        pk_win = (tmid - tpad, tmid + tpad)
   
    else:
        pk_win = None

    xx = np.where( np.sum(dd, axis=0) == 0 )[0]
    chan_weights = np.ones(len(freqs))
    chan_weights[xx] = 0
    fac = np.mean( chan_weights.astype('float') )

    ## CALC DDM + SPEC + MIDX ## 
    tt1, freqs1, tpk, ddm, spec, midx = dspec.get_ddm_spec_mod2(\
                tt, freqs, dd, dm, avg_tsamp=tsamp_avg, avg_chan=1, 
                chan_weights=chan_weights, pk_win=pk_win)

    toffs = tt - tpk
    tshift = tt - tmid

    ## GET dDM - dt DATA ##
    delta_DM = 2 * dt * tsamp_avg / dedisp.dm_delay(freqs[0], freqs[-1], 1.0)
    print("DDM = %.2f" %delta_DM)
    dDMs = np.linspace(-5 * delta_DM, + 5 * delta_DM, 201)
    dms = dDMs + dm
    DM_dt = dedisp.dedisperse_dspec_multi(dd.T, dms, freqs, freqs[-1], dt)
    DM_dt /= fac
    ddt = -0.5 * dedisp.dm_delay(freqs[0], freqs[-1], dDMs)

    ##################
    ## SETUP FIGURE ##
    ##################
    if outfile is not None:
        plt.ioff()
    
    # Make Figure
    fig = plt.figure(figsize=(11, 8.5))

    ax11 = fig.add_subplot(2, 3, 1)
    ax12 = fig.add_subplot(2, 3, 2)
    ax13 = fig.add_subplot(2, 3, 3)
    ax21 = fig.add_subplot(2, 3, 4)
    ax22 = fig.add_subplot(2, 3, 5)
    ax23 = fig.add_subplot(2, 3, 6)

    # Some params 
    dspec_snr = 50 #* np.sqrt(nchan_avg * tsamp_avg)
    min_frac  = 0.25   # dspec_min = dspec_snr * min_frac
    #cmap = plt.cm.magma_r
    cmap = plt.cm.Greys
    im_kwargs = {'interpolation' : 'nearest', 
                 'aspect'        : 'auto', 
                 'origin'        : 'lower'}

    ######################
    ###   dDM - dt     ###
    ######################
    ax_dd = ax21

    dd_cmap = plt.cm.jet
    #ext_dd = [toffs[0]-0.5*dt, toffs[-1] + 0.5*dt, dms[0]-0.5*delta_DM, dms[-1] +0.5*delta_DM]
    ext_dd = [tshift[0]-0.5*dt, tshift[-1] + 0.5*dt, dms[0]-0.5*delta_DM, dms[-1] +0.5*delta_DM]
    ax_dd.imshow(DM_dt, extent=ext_dd, cmap=dd_cmap, **im_kwargs) 
    ax_dd.plot(ddt, dms, ls='--', lw=2, color='r')
    ax_dd.axvline(x=0, ls='--', lw=2, c='r')

    dd_xlim = ( 1.2 * np.min(ddt), 1.5 * np.max(ddt) )
    ax_dd.set_xlim(dd_xlim)
    ax_dd.set_ylim(dms[0], dms[-1])
    ax_dd.set_xlabel("Time Offset (s)")
    ax_dd.set_ylabel("DM (pc/cc)")
    ax_dd.tick_params(axis='both', color='w')


    ######################
    ###  DM SNR CURVE  ###
    ######################
    ax_dm = ax11

    #xx_pk = np.argmin( np.abs(toffs) ) 
    xx_pk = np.argmin( np.abs(tshift) ) 
    ax_dm.plot( dms, DM_dt[:, xx_pk], ls='-', c='k', lw=3)
    ax_dm.axvline(x = dm, ls='--', c='r')

    ax_dm.tick_params(labeltop='on', labelbottom='off')
    ax_dm.xaxis.set_label_position('top')

    ax_dm.tick_params(axis='x', labelsize=8)
    
    ax_dm.set_xlabel("DM (pc/cc)")
    ax_dm.set_ylabel("Norm S/N")
   
    
    ###############
    ###  PULSE  ###
    ###############
    tp_val = 0.5 * tsamp_avg
    tp_lo  = max( -1.0 * tp_val , tt1[0] - tpk )
    tp_hi  = min( +1.0 * tp_val , tt1[-1] - tpk)
    tplim = (tp_lo, tp_hi)
    ax_ts = ax12
    ax_ts.plot(tt1-tpk, ddm, c='k')
    ax_ts.axhline(y=0, ls='-', c='k', alpha=0.2)
    ax_ts.axvline(x=0, ls='-', c='r', alpha=0.2)
    ax_ts.axvline(x=tmid-tpk, ls='--', c='b', alpha=0.2)

    p_str = r"dfact = %d" %dfact + "\n" +\
            r"t = %.2f" %time 

    ax_ts.text(0.9, 0.9, p_str, ha='right', va='top', 
               transform=ax_ts.transAxes, fontsize=16)
    
    ax_ts.set_xlim(tplim)
    #ax_ts.set_ylabel("SNR")
    ax_ts.tick_params(labeltop='on', labelbottom='off',
                      labelleft='off')
    ax_ts.xaxis.set_label_position("top")
    ax_ts.set_xlabel("Time Offset (s)")
    pulse_ylim = ax_ts.get_ylim()

    ######################
    ### MOD - INDEX    ###
    ######################
    ax_mi = ax13
    snr_cut = 5.0

    ax_mi.plot(midx, ddm, ls='', marker='o', ms=8)
    xpk = np.where( tt1 == tpk )[0]
    ax_mi.plot([midx[xpk]], [ddm[xpk]], marker='o', c='r', ms=8)

    if chan_weights is None:
        midx_cut = np.sqrt(len(freqs) / 1.0) / snr_cut
    else:
        midx_cut = np.sqrt(np.sum(chan_weights) / 1.0) / snr_cut


    info_str = r"$m_{\rm I} = %.2f$" %midx[xpk] + "\n" +\
               r"${\rm SNR} = %.1f$" %ddm[xpk]  + "\n" +\
               r"${\rm DM} = %d$" %dm

    ax_mi.text(0.9, 0.9, info_str, ha='right', va='top', 
               transform=ax_mi.transAxes, fontsize=20)

    ax_mi.axvline(x=midx_cut, c='b', lw=3, alpha=0.2)
    ax_mi.axhline(y=snr_cut, c='b', lw=3, alpha=0.2)

    ax_mi.axvline(x=0, c='k', lw=2, ls='-', alpha=0.1)
    ax_mi.axhline(y=0, c='k', lw=2, ls='-', alpha=0.1)

    ax_mi.set_xlim(-1 , 20 * midx_cut)
    #ax_mi.set_ylim(-1)
    ax_mi.set_ylim(pulse_ylim)

    ax_mi.tick_params(labeltop='on', labelbottom='off', 
                      labelleft='off', labelright='on')
    ax_mi.xaxis.set_label_position("top")
    ax_mi.set_xlabel("Mod Index")
    ax_mi.yaxis.set_label_position("right")
    ax_mi.set_ylabel("SNR")

    ###################
    ### DSPEC SWEEP ###
    ###################
    ax_ds = ax23

    ## Some Useful Limits ##
    df1 = freqs[1] - freqs[0]
    dt1 = tt[1] - tt[0]

    tshow = 2.0
    tsweep = dedisp.dm_delay(freqs[0], freqs[-1], dm)
    tlim_min = max( -0.5 * tsweep, np.min(toffs) )
    tlim_max = min( 1.5 * tsweep, np.max(toffs) )
    tlim = (tlim_min, tlim_max)
    flim = (freqs[0] - 0.5 * df1, freqs[-1] + 0.5 * df1)

    if chan_weights is None:
        dd_sig = np.std(dd)
    else:
        yy = np.where(chan_weights)[0]
        dd_sig = np.std(dd[:, yy])

    print( np.sum(chan_weights) )
    print( len(freqs) )
    print(dd_sig)
    print(dd.shape)
    print( np.min(dd / dd_sig), np.max(dd / dd_sig ))

    im_ext = [tt[0] - tpk - 0.5 * dt1, tt[-1] - tpk + 0.5 * dt1, 
              freqs[0], freqs[-1]]

    ax_ds.imshow(dd.T / dd_sig, interpolation='nearest', origin='lower', 
                 aspect='auto', extent = im_ext, cmap=cmap, 
                 vmin = -1, vmax=3) #vmin= -1 * min_frac * dspec_snr, vmax = dspec_snr)

    # Add Dispersion Sweep
    offset = dt0 * 20
    dm_curve = dedisp.dm_delay(freqs, freqs[-1], dm)
    ax_ds.plot(dm_curve - offset, freqs, lw=2, c='r')
    ax_ds.plot(dm_curve + offset, freqs, lw=2, c='r')

    # Set limits
    ax_ds.set_xlim(tlim)
    ax_ds.set_ylim(flim)

    ax_ds.tick_params(labelleft='off', labelright='on')
    ax_ds.yaxis.set_label_position('right')

    ax_ds.set_xlabel("Time Offset (s)")
    ax_ds.set_ylabel("Frequency (MHz)")



    ##################
    ### FREQ PHASE ###
    ##################
    ax_fp = ax22
    
    tt2, freqs2, dd2 = dedisp.dspec_avg_dm(dd, tt, freqs, freqs[-1], dt, 
                                           avg_tsamp=tsamp_avg, avg_chan=nchan_avg, 
                                           dm0=dm, final_reverse=False)
    #dd2 = dedisp.dedisperse_dspec(dd.T, dm, freqs, freqs[-1], dt1)
    
    dt2 = tt2[1] - tt2[0]
    df2 = freqs2[1] - freqs2[0]
    im_ext2 = [tt2[0] - tpk - 0.5 * dt2, tt2[-1] - tpk + 0.5 * dt2, 
              freqs2[0] - 0.5 * df2, freqs2[-1] + 0.5 * df2]

    dd_sig2 = dd_sig / np.sqrt(tsamp_avg * nchan_avg)
    ax_fp.imshow(dd2.T / dd_sig2 , interpolation='nearest', origin='lower', 
                 aspect='auto', extent = im_ext2, cmap=cmap, 
                 vmin=-1, vmax=3) #vmin= -1 * min_frac * dspec_snr, vmax = dspec_snr)

    ax_fp.set_xlim(tplim)
    ax_fp.set_xlabel("Time Offset (s)")
    #ax_fp.set_ylabel("Frequency (GHz)")
    ax_fp.set_yticklabels([])

    #plt.subplots_adjust(left=0.05, right=0.95)
    plt.subplots_adjust(hspace=0.05, wspace=0.05)

    ################
    ## FINISH UP  ##
    ################

    if outfile is not None:
        plt.savefig(outfile, dpi=150, bbox_inches='tight')
        plt.close()
        plt.ion()

    else:
        plt.show()

    return 
    #return tpk, tt, ddm, spec, midx, dd



def summary_pulse_from_candlist(candlist, topdir, tsearch, tpad=0, 
                                use_dfact=False, tsamp_avg=1, nchan_avg=1):
    for ii, cand in enumerate(candlist):
        outfile = "Cand%04d.png" %(ii)
        summary_pulse(cand, topdir, tsearch, use_dfact=use_dfact, 
                      tsamp_avg=tsamp_avg, nchan_avg=nchan_avg, 
                      outfile=outfile)
    return 


def summary_pulse_from_candlist2(candlist, topdir, tsearch, tpad=0, 
                                use_dfact=False, tsamp_avg=1, nchan_avg=1):
    for ii, cand in enumerate(candlist):
        outfile = "Cand%04d.png" %(ii)
        summary_pulse2(cand, topdir, tsearch, use_dfact=use_dfact, 
                      tsamp_avg=tsamp_avg, nchan_avg=nchan_avg, 
                      outfile=outfile)
    return 



if __name__ == '__main__':
    # MAIN
    #p1 = np.array([ 9.99941374e-01,  -4.41724462e-04])
    #z = np.poly1d(p1)
    #fitsfiles = glob.glob('beam0002/psrfits/*fits')
    #off_starts = get_starts_for_fits(fitsfiles)
    ##off_time, fitsname = get_fits_name_time(z(3091.970), fitsfiles, off_starts)
    pass



