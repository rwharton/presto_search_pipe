from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import pickle
from collections import namedtuple
from astropy.io import fits
import simple_dm as dedisp
import glob


def get_data(infile, tstart, tdur, apply_corr=False, chan_weights=None):
    hdulist = fits.open(infile)
    hdu1, hdu2 = hdulist

    tstop = tstart + tdur 
    
    dt = hdu2.header['tbin']
    freqs = hdu2.data[0]['dat_freq']
    df = freqs[1] - freqs[0]
    nchan = len(freqs)
    tsub = hdu2.data[0]['TSUBINT']
    nsblk = hdu2.header['nsblk']

    row_start = max(int(tstart / tsub), 0)
    row_stop  = max(int(np.ceil(tstop / tsub)), row_start + 1)
    trow_start = row_start * tsub

    dd = hdu2.data[row_start : row_stop]['data'][:, :, 0, :]
    #print(row_stop-row_start)
    #print(dd.shape)
    
    if not apply_corr:
        pass

    else:
        dd_scl = hdu2.data[row_start : row_stop]['dat_scl']
        dd_off = hdu2.data[row_start : row_stop]['dat_offs']
        
        dd_scl = np.reshape(dd_scl, (-1, 1, nchan))
        dd_off = np.reshape(dd_off, (-1, 1, nchan))
        #print(dd_scl.shape)
        #print(dd_off.shape)
        #print(dd.shape)
        dd = dd * dd_scl + dd_off

    #dd = np.reshape(dd, ( (row_stop - row_start) * nsblk, nchan ))
    dd = np.reshape(dd, (-1, nchan ))
    
    idx = max(int( (tstart - trow_start) / dt ), 0)
    ntx = int(tdur / dt)

    dd = dd[idx : idx + ntx, :]

    tt = np.arange(dd.shape[0]) * dt + trow_start + idx * dt

    hdulist.close()

    if chan_weights is not None:
        dd = dd * chan_weights

    # make sure even number of rows
    if dd.shape[0] % 2:
        tt = tt[:-1]
        dd = dd[:-1]
    else: 
        pass
    
    return tt, freqs, dd


def get_obs_info(infile):
    hdulist = fits.open(infile)
    hdu1, hdu2 = hdulist
    
    obs_info = {}
    hdr = hdu1.header
    keys = ['RA', 'DEC', 'STT_IMJD', 'STT_SMJD', 'STT_OFFS']
    for kk in keys:
        obs_info[kk] = hdr.get(kk)

    offs_sub = hdu2.data[0].field('OFFS_SUB')
    tsubint  = hdu2.data[0].field('TSUBINT')
    start_offs = offs_sub - 0.5 * tsubint
    obs_info['SEC_OFFS'] = start_offs

    stt_imjd = obs_info['STT_IMJD']
    stt_smjd = obs_info['STT_SMJD']
    stt_offs = obs_info['STT_OFFS']

    mjd_obs = np.float64(stt_imjd) + (stt_smjd + stt_offs) / (24. * 3600.0)
    mjd_dat = mjd_obs + start_offs / (24.0 * 3600.0)

    obs_info['MJD_OBS'] = mjd_obs
    obs_info['MJD_DAT'] = mjd_dat

    hdulist.close()

    return obs_info


def save_obs_info(fitsfile, basename):
    obs_info = get_obs_info(fitsfile)
    outfile = '%s_obs_info.p' %basename
    pickle.dump(obs_info, open(outfile, 'wb'))
    return


def avg_chan(dd, nchan=4):
    nsteps = int(dd.shape[1] / nchan)
    for nn in xrange(nsteps):
        dd[:,nn*nchan:(nn+1)*nchan] = \
            np.outer(np.mean(dd[:, nn*nchan:(nn+1)*nchan], axis=1), np.ones(nchan))
    return dd


def get_starts_for_fits(fitslist):
    mjd_starts = []

    for fitsfile in fitslist:
        obs_info = get_obs_info(fitsfile)
        mjd_starts.append( obs_info['MJD_DAT'] )
    mjd_starts = np.array( mjd_starts )
    off_starts = (mjd_starts - mjd_starts[0]) * 3600.0 * 24

    return off_starts 


def get_fits_name_time(time, fits_names, fits_times):
    """
    Assumes times and fits names are same, time-sortd order
    """
    xx0 = np.where( time > fits_times )[0][-1]
    off_time = time - fits_times[xx0]
    return off_time, fits_names[xx0], xx0


def get_part_time(time, part_times):
    """
    Assumes times and fits names are same, time-sortd order
    """
    parts = np.arange(1, len(part_times) + 1)
    xx0 = np.where( time > part_times )[0][-1]
    off_time = time - part_times[xx0]
    return off_time, parts[xx0]


def get_dedispersed_timeseries(tt, freqs, dat, dm):
    dt = tt[1] - tt[0]
    df = freqs[1] - freqs[0]
    ddm = dedisp.dedisperse_one(dat.T, dm, dt, df, freqs[0])
    return ddm


def get_meansub_normed(dd, xpk, W, ignore_zero=True):
    N = len(dd)
    xlo = max(0, xpk - W)
    xhi = min(N, xpk + W)
    xon = np.arange(xlo, xhi)
    xoff = np.setdiff1d(np.arange(N), xon)
    ddoff = dd[xoff]

    if ignore_zero:
        xx_nz = np.where( np.abs(ddoff) > 0 )[0]
        if len(xx_nz) == 0:
            avg = 0
            sig = 1
        else:
            avg = np.mean( ddoff[ xx_nz ] )
            sig = np.std(  ddoff[ xx_nz ] )

    else:
        avg  = np.mean(dd[xoff])
        sig  = np.std(dd[xoff])

    return (dd-avg) / sig


def get_ddm_spec_mod(tt, freqs, dat, dm, pk_win=None, chan_weights=None, Wx=10):
    dt = tt[1] - tt[0]
    df = freqs[1] - freqs[0]
    #dd_dm = dedisp.dedisperse_dspec(dat.T, dm, dt, df, freqs[0]).T
    dd_dm = dedisp.dedisperse_dspec(dat.T, dm, freqs, freqs[-1], dt).T
    
    if chan_weights is None:
        ddm = np.mean(dd_dm, axis=1)
        if pk_win is not None:
            xx = np.where( (tt > pk_win[0]) & (tt < pk_win[1]) )[0]
            if len(xx) > 0:
                xpk = xx[np.argmax(ddm[xx])]
            else:
                print("INCORRECT PK_WIN -- USING FULL ARRAY")
                xpk = np.argmax(ddm)
        else:
            xpk = np.argmax(ddm)
        ddm = get_meansub_normed(ddm, xpk, Wx)

        spec = dd_dm[xpk]
        tpk = tt[xpk]
        mod_idx = np.std(dd_dm, axis=1) / np.mean(dd_dm, axis=1)

    else:
        ddm = np.sum(dd_dm * chan_weights, axis=1) / np.sum(chan_weights)
        if pk_win is not None:
            xx = np.where( (tt > pk_win[0]) & (tt < pk_win[1]) )[0]
            if len(xx) > 0:
                xpk = xx[np.argmax(ddm[xx])]
            else:
                print("INCORRECT PK_WIN -- USING FULL ARRAY")
                xpk = np.argmax(ddm)
        else:
            xpk = np.argmax(ddm)
        ddm = get_meansub_normed(ddm, xpk, Wx)
        spec = dd_dm[xpk] * chan_weights
        tpk = tt[xpk]
        mom2 = np.sum(dd_dm**2.0 * chan_weights, axis=1) / np.sum(chan_weights)
        mom1 = np.sum(dd_dm * chan_weights, axis=1) / np.sum(chan_weights)
        mod_idx = np.sqrt((mom2 - mom1**2.0) / mom1**2.0)

    return tpk, ddm, spec, mod_idx


def get_ddm_spec_mod2(tt, freqs, dat, dm, avg_tsamp=1, avg_chan=1, 
                      pk_win=None, chan_weights=None, Wx=10):
    dt0 = tt[1] - tt[0]
    df0 = freqs[1] - freqs[0]
    #dd_dm = dedisp.dedisperse_dspec(dat.T, dm, dt, df, freqs[0]).T
    #dd_dm = dedisp.dedisperse_dspec(dat.T, dm, freqs, freqs[-1], dt).T
    tt1, freqs1, dd_dm = dedisp.dspec_avg_dm(dat, tt, freqs, freqs[-1], dt0, 
                                avg_tsamp=avg_tsamp, avg_chan=avg_chan, dm0=dm, 
                                final_reverse=False)
    
    if chan_weights is None:
        ddm = np.mean(dd_dm, axis=1)
        if pk_win is not None:
            xx = np.where( (tt1 > pk_win[0]) & (tt1 < pk_win[1]) )[0]
            if len(xx) > 0:
                xpk = xx[np.argmax(ddm[xx])]
            else:
                print("INCORRECT PK_WIN -- USING FULL ARRAY")
                xpk = np.argmax(ddm)
        else:
            xpk = np.argmax(ddm)
        ddm = get_meansub_normed(ddm, xpk, Wx)

        spec = dd_dm[xpk]
        tpk = tt1[xpk]
        mod_idx = np.std(dd_dm, axis=1) / np.mean(dd_dm, axis=1)

    else:
        ddm = np.sum(dd_dm * chan_weights, axis=1) / np.sum(chan_weights)
        if pk_win is not None:
            xx = np.where( (tt1 > pk_win[0]) & (tt1 < pk_win[1]) )[0]
            if len(xx) > 0:
                xpk = xx[np.argmax(ddm[xx])]
            else:
                print("INCORRECT PK_WIN -- USING FULL ARRAY")
                xpk = np.argmax(ddm)
        else:
            xpk = np.argmax(ddm)
        ddm = get_meansub_normed(ddm, xpk, Wx)
        spec = dd_dm[xpk] * chan_weights
        tpk = tt1[xpk]
        mom2 = np.sum(dd_dm**2.0 * chan_weights, axis=1) / np.sum(chan_weights)
        mom1 = np.sum(dd_dm * chan_weights, axis=1) / np.sum(chan_weights)
        mod_idx = np.sqrt((mom2 - mom1**2.0) / mom1**2.0)

    return tt1, freqs1, tpk, ddm, spec, mod_idx


def dspec_stats(infile, tmid, tsearch, tpad, tshow, DM, tsamp_avg=1, 
                chan_weights=None):
    # Grab raw data
    tstart = tmid - tsearch * 0.5
    tt, freqs, dd = get_data(infile, tstart, tsearch, chan_weights=chan_weights)
    dt0 = tt[1] - tt[0]
    dt  = tt[1] - tt[0]
    
    # Average in time if necessary
    if tsamp_avg > 1:
        tt, freqs, dd = dedisp.dspec_avg_dm(dd, tt, freqs, freqs[-1], dt, 
                                  avg_tsamp=tsamp_avg, avg_chan=1, dm0=DM)
        dt = tt[1] - tt[0]
    
    else:
        dt = dt0

    # Make sure dd even
    if len(tt) % 2:
        tt = tt[:-1]
        dd = dd[:-1]
    
    if tpad is not None:
        pk_win = (tmid - tpad, tmid + tpad)
    else:
        pk_win = None

    xx = np.where( np.sum(dd, axis=0) == 0 )[0]
    chan_weights = np.ones(len(freqs))
    chan_weights[xx] = 0

    #print(dd.shape)
    tpk, ddm, spec, midx = get_ddm_spec_mod(tt, freqs, dd, DM, 
                           chan_weights=chan_weights, pk_win=pk_win)

    xpk = np.where(tt == tpk)[0][0]
    
    return tpk, ddm[xpk], midx[xpk]


def dspec_stats2(infile, tmid, tsearch, tpad, dm, tsamp_avg=1, 
                chan_weights=None):
    # Grab raw data
    tstart = tmid - tsearch * 0.25
    tt, freqs, dd = get_data(infile, tstart, tsearch)
    dt  = tt[1] - tt[0]

    if tpad <= 0:
        tpad = tsamp_avg * dt
    else:
        pass

    if tpad is not None:
        pk_win = (tmid - tpad, tmid + tpad)
    else:
        pk_win = None

    xx = np.where( np.sum(dd, axis=0) == 0 )[0]
    chan_weights = np.ones(len(freqs))
    chan_weights[xx] = 0

    ## CALC DDM + SPEC + MIDX ## 
    tt1, freqs1, tpk, ddm, spec, midx = get_ddm_spec_mod2(\
                tt, freqs, dd, dm, avg_tsamp=tsamp_avg, avg_chan=1,
                chan_weights=chan_weights, pk_win=pk_win)

    xpk = np.where( tt1 == tpk )[0][0]

    return tpk, ddm[xpk], midx[xpk]


def get_file_name_old(indir, basename, beamnum, nzeros=4):
    return "{}/{}_beam{:0{}}.fits".format(indir, basename, beamnum, nzeros)


def get_file_name(topdir, basename, part, beamnum, nzeros=4):
    path1 = "%s/part%d/psrfits" %(topdir, part) 
    path2 = "{}_part{}_beam{:0{}}.fits".format(basename, part, beamnum, nzeros)
    path = "%s/%s" %(path1, path2)
    return path


def stats_from_candfile_parts(candfile, topdir, tsearch, tpad,
                              fixdm=None, b2t_coeff=None, logfile=None):
    dat = np.loadtxt(candfile)
    beams   = dat[:, 0].astype('int')
    dms_det = dat[:, 1].astype('float')
    snrs    = dat[:, 2].astype('float')
    times   = dat[:, 3]
    samps   = dat[:, 4].astype('int')
    dfacts  = dat[:, 5].astype('int')
    nhits   = dat[:, 6].astype('int')
    groups  = dat[:, 7].astype('int')
    
    tt_topos = bary_to_topo(times, b2t_coeff)
    
    N = len(dat)

    Cand = namedtuple('Cand', ['beam', 'dm', 'dm_det', 'snr', 'tt', 
                               'samp', 'dfact', 'group', 'nhits', 
                               'tt_topo', 'fit_tt', 'fit_snr', 'fit_midx'])
    
    candlist = []

    if fixdm is None:
        dms = dms_det
    else:
        dms = np.ones(len(dms_det)) * fixdm

    for ii in xrange(N):
        log_output("%d / %d" %(ii, N), logfile)
        fitsfiles = glob.glob("%s/part*/psrfits/*beam%04d*fits" %(topdir, beams[ii]))
        #print(fitsfiles)
        fitsfiles.sort()
        off_starts = get_starts_for_fits(fitsfiles)
        
        otime, infile, xxf = get_fits_name_time(tt_topos[ii], fitsfiles, off_starts)
        
        #ftt, fss, fmm = dspec_stats(infile, otime, tsearch, tpad, tshow, dms[ii],
        #                            tsamp_avg=dfacts[ii], chan_weights=chan_weights)
        
        ftt, fss, fmm = dspec_stats2(infile, otime, tsearch, tpad, dms[ii],
                                     tsamp_avg=dfacts[ii])
        
        ftt += off_starts[xxf]

        cc = Cand(tt=times[ii], beam=beams[ii], dm=dms[ii], dm_det=dms_det[ii], 
                  snr=snrs[ii], samp=samps[ii], dfact=dfacts[ii], nhits=nhits[ii], 
                  group=groups[ii], tt_topo=tt_topos[ii], fit_tt = ftt, 
                  fit_snr=fss, fit_midx=fmm)
        candlist.append(cc)
    
    return candlist


def stats_from_candfile(candfile, fitsfile, basename, tsearch, tpad, 
                        fixdm=None, logfile=None):
    dat = np.loadtxt(candfile)
    times = dat[:, 0]
    dms_det   = dat[:, 1]
    snrs  = dat[:, 2]
    N = len(dat)

    Cand = namedtuple('Cand', ['tt', 'dm', 'dm_det', 'snr', 'fit_tt', 
                               'fit_snr', 'fit_midx'])
    candlist = []

    if fixdm is None:
        dms = dms_det
    else:
        dms = np.ones(len(dms_det)) * fixdm

    for ii in xrange(N):
        log_output("%d / %d" %(ii, N), logfile)
        ftt, fss, fmm = dspec_stats2(fitsfile, times[ii], tsearch, tpad, dms[ii])
        cc = Cand(tt=times[ii], dm=dms[ii], dm_det=dms_det[ii], 
                  snr=snrs[ii], fit_tt = ftt, fit_snr=fss, fit_midx=fmm)
        candlist.append(cc)
    return candlist


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
    dms_det = dat[:, 1].astype('float')
    snrs    = dat[:, 2].astype('float')
    times   = dat[:, 3].astype('float')
    samps   = dat[:, 4].astype('int')
    dfacts  = dat[:, 5].astype('int')
    nhits   = dat[:, 6].astype('int')
    groups  = dat[:, 7].astype('int')

    N = len(dat)

    if dat.shape[-1] == 11:
        tt_topos = dat[:, 8].astype('float')
        fit_snr  = dat[:, 9].astype('float')
        fit_midx = dat[:, 10].astype('float')
    else:
        tt_topos = bary_to_topo(times, b2t_coeff)
        fit_snr = np.zeros(len(snrs))
        fit_midx = np.zeros(len(snrs))

    Cand = namedtuple('Cand', ['beam', 'dm', 'dm_det', 'snr', 'tt', 
                               'samp', 'dfact', 'nhits', 'group', 
                               'tt_topo', 'fit_snr', 'fit_midx'])
    candlist = []

    if fixdm is None:
        dms = dms_det
    else:
        dms = np.ones(len(dms_det)) * fixdm

    for ii in xrange(N):
        cc = Cand(beam=beams[ii], dm=dms[ii], dm_det=dms_det[ii], 
                  snr=snrs[ii], tt=times[ii], samp=samps[ii], 
                  dfact=dfacts[ii], nhits=nhits[ii], group=groups[ii], 
                  tt_topo=tt_topos[ii], fit_snr=fit_snr[ii], fit_midx=fit_midx[ii])
        candlist.append(cc)
    return candlist


def array_attr(clist, attr):
    if hasattr(clist[0], attr):
        return np.array([ getattr(cc, attr) for cc in clist ])
    else:
        print("Object has no attribute: %s" %attr)
        return


def write_select_cands(cands, outfile):
    fout = open(outfile, 'w')
    hdr = "#{:^12}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}".format(\
        "Beam", "DM", "SNR", "Time_b", "Sample", "Downfact", "Nhits", "Group", "Time_t", 
        "fit_SNR", "fit_mI")
    fout.write(hdr + "\n")

    fit_snrs = array_attr(cands, 'fit_snr')
    idx = np.argsort(fit_snrs)
    sorted_cands = [ cands[ii] for ii in idx[::-1] ]
    
    for cc in sorted_cands:
        print(cc)

        outstr = "{:^12d}{:>10.2f}{:>10.2f}{:>10.3f}{:>10d}{:^10d}{:^10d}{:^10d}{:>10.3f}{:>10.2f}{:>10.3f}".format(\
            cc.beam, cc.dm, cc.snr, cc.tt, cc.samp, cc.dfact, cc.nhits, cc.group, cc.tt_topo, cc.fit_snr, cc.fit_midx)
        fout.write(outstr + "\n")
    fout.close()
    return sorted_cands


def log_output(text, logfile):
    if logfile is None:
        print(text)
        sys.stdout.flush()
    else:
        fout = open(logfile, 'a+')
        fout.write("%s\n" %text)
        fout.close()
    return


def bary_to_topo(tt_bary, b2t_coeff):
    """
    Convert tt_bary barycentric time to 
    topocentric using coeffecients given 
    by b2t_coeff

    if b2t_coeff is None, just return tt_bary
    """
    if b2t_coeff is None:
        b2t_coeff = np.array([1, 0])
    else:
        pass
    
    z = np.poly1d(b2t_coeff)

    return z(tt_bary)
    



if __name__ == '__main__':
    # MAIN
    #p1 = np.array([ 9.99941374e-01,  -4.41724462e-04])
    #z = np.poly1d(p1)
    #fitsfiles = glob.glob('beam0002/psrfits/*fits')
    #off_starts = get_starts_for_fits(fitsfiles)
    #off_time, fitsname = get_fits_name_time(z(3091.970), fitsfiles, off_starts)
    pass
