import numpy as np
import matplotlib.pyplot as plt
import sys

kdm = 4148.808 # MHz^2 / (pc cm^-3)

def chan_smear(fc, df, DM, kdm=kdm):
    return 2 * kdm * DM * df / fc**3.0

def equiv_chan(tdm, fc, DM, kdm=kdm):
    return tdm * fc**3.0 / (2 * kdm * DM)

def dm_delay(f1, f2, DM, kdm=kdm):
    return kdm * DM * (1.0 / (f1 * f1) - 1.0 / (f2 * f2))


def deltaT_old(ichan, dt, df, f0, kdm=kdm):
    return (kdm / dt) * ((f0 + ichan * df)**-2.0 - f0**-2)


def dmdt_old(DM, ichan, dt, df, f0, kdm=kdm):
    return int(np.round( DM * deltaT_old(ichan, dt, df, f0, kdm=kdm)))


def dmdt_float_old(DM, ichan, dt, df, f0, kdm=kdm):
    return DM * deltaT(ichan, dt, df, f0, kdm=kdm)


def deltaT(freqs, f0, dt, kdm=kdm):
    return (kdm / dt) * (freqs**-2.0 - f0**-2)


def dmdt(DM, freqs, f0, dt, kdm=kdm):
    return DM * deltaT(freqs, f0, dt, kdm=kdm)


def fake_dspec(tt, t0, W, f0, df, nchan, DM):
    pp = np.exp(-(tt - t0)**2.0 / (2 * W**2.0))
    dt = tt[1] - tt[0]
    nt = len(tt)
    dsamp = np.array([ dmdt_old(DM, ichan, dt, df, f0) for ichan in xrange(nchan) ])
    offs = dsamp - np.min(dsamp)
    dspec = np.zeros(shape=(nchan, len(tt)))
    for jj in xrange(nchan):
        dspec[jj, offs[jj]:] = pp[:nt-offs[jj]]
    return dspec


def dedisperse_one(dspec, dm, dt, df, f0, kdm=kdm):
    nchans = dspec.shape[0]
    nt = dspec.shape[1]
    dsamps = np.array([ dmdt_old(dm, ichan, dt, df, f0, kdm=kdm) for ichan in xrange(nchans) ])
    dsamps -= np.min(dsamps)
    tpad = np.max(dsamps)
    outarr = np.zeros( nt + tpad )
    for ii in xrange(nchans):
        osl = slice(tpad - dsamps[ii], nt + tpad - dsamps[ii])
        outarr[osl] += dspec[ii]
    return outarr[tpad:nt + tpad] / float(nchans)


def dedisperse_dspec_old(dspec, dm, dt, df, f0, kdm=kdm):
    nchans = dspec.shape[0]
    nt = dspec.shape[1]
    dsamps = np.array([ dmdt_old(dm, ichan, dt, df, f0, kdm=kdm) for ichan in xrange(nchans) ])
    dsamps -= np.min(dsamps)
    tpad = np.max(dsamps)
    outarr = np.zeros( (nchans, nt + tpad) )
    for ii in xrange(nchans):
        osl = slice(tpad - dsamps[ii], nt + tpad - dsamps[ii])
        outarr[ii, osl] = dspec[ii]
    return outarr[:, tpad:nt + tpad]


def dedisperse_dspec(dspec, dm, freqs, f0, dt, kdm=kdm, reverse=False, use_pad=False):
    nchans = dspec.shape[0]
    nt = dspec.shape[1]
    dsamps = dmdt(dm, freqs, f0, dt, kdm=kdm)
    #dsamps -= np.min(dsamps)
    if reverse:
        sgn = -1.0
    else:
        sgn = +1.0 
   
    if use_pad:
        pad = np.zeros( dspec.shape )
        dds  = np.hstack( (pad, dspec, pad) )
    else:
        dds = dspec

    dout = np.zeros( dspec.shape )
    for ii, dd in enumerate(dds):
        ak = np.fft.rfft(dd)
        bfreq = np.arange(len(ak)) / (1.0 * len(dd))
        shift = np.exp(sgn * 1.0j * 2 * np.pi * bfreq * dsamps[ii])
        dd_shift = np.fft.irfft( ak * shift )
        if use_pad: 
            dout[ii] = dd_shift[nt : 2 * nt]
        else:
            dout[ii] = dd_shift[:]

    return dout


def dedisperse_dspec_multi(dspec, dms, freqs, f0, dt, kdm=kdm):
    Ndms = len(dms)
    Nt   = dspec.shape[1]
    out_arr = np.zeros( (Ndms, Nt) )
    for ii, dm in enumerate(dms):
        dd = dedisperse_dspec(dspec, dm, freqs, freqs[-1], dt)
        out_arr[ii] = np.mean( dd, axis=0 )

    return out_arr


def dspec_avg_chan(dspec, freqs, avg_chan=1):
    Nchan = dspec.shape[1]
    n = int(Nchan / avg_chan)

    freq_out = np.zeros(n)
    dd_out = np.zeros( (dspec.shape[0], n) )

    for ii in xrange(n):
        sl = slice(ii * avg_chan, (ii+1) * avg_chan)
        freq_out[ii] = np.mean(freqs[sl])
        dd_out[:, ii] = np.mean(dspec[:, sl], axis=1)

    # If odd number of channels, make even
    if len(freq_out) % 2:
        freq_out = freq_out[:-1]
        dd_out = dd_out[:, :-1]

    return freq_out, dd_out


def dspec_avg_tsamp(dspec, tt, avg_tsamp=1):
    Nt= dspec.shape[0]
    n = int(Nt / avg_tsamp)

    tt_out = np.zeros(n)
    dd_out = np.zeros( (n, dspec.shape[1]) )

    for ii in xrange(n):
        sl = slice(ii * avg_tsamp, (ii+1) * avg_tsamp)
        tt_out[ii] = np.mean(tt[sl])
        dd_out[ii, :] = np.mean(dspec[sl, :], axis=0)
    
    # If odd number of samples, make even
    if len(tt_out) % 2:
        tt_out = tt_out[:-1]
        dd_out = dd_out[:-1, :]

    return tt_out, dd_out


def dspec_avg_chan_dm(dspec, freqs, f0, dt, avg_chan=1, dm0=0.0):
    dd_dm0 = dedisperse_dspec(dspec.T, dm0, freqs, f0, dt)
    avg_freqs, davg_dm0 = dspec_avg_chan(dd_dm0.T, freqs, avg_chan=avg_chan)
    davg = dedisperse_dspec(davg_dm0.T, dm0, avg_freqs, f0, dt, reverse=True)
    return avg_freqs, davg.T


def dspec_avg_dm_old(dspec, tt, freqs, f0, dt, avg_tsamp=1, avg_chan=1, dm0=0.0):
    dd_dm0 = dedisperse_dspec(dspec.T, dm0, freqs, f0, dt)
    #print(dd_dm0.shape)
    avg_freqs, davg_dm0 = dspec_avg_chan(dd_dm0.T, freqs, avg_chan=avg_chan)
    #print(davg_dm0.shape)
    avg_tt, davg_dm0 = dspec_avg_tsamp(davg_dm0, tt, avg_tsamp=avg_tsamp)
    #print(davg_dm0.shape)
    dt_avg = dt * avg_tsamp
    davg = dedisperse_dspec(davg_dm0.T, dm0, avg_freqs, f0, dt_avg, reverse=True)
    return avg_tt, avg_freqs, davg.T


def dspec_avg_dm(dspec, tt, freqs, f0, dt, avg_tsamp=1, avg_chan=1, dm0=0.0, 
                 final_reverse=True):
    dd_dm0 = dedisperse_dspec(dspec.T, dm0, freqs, f0, dt)
    #print(dd_dm0.shape)
    avg_freqs, davg_dm0 = dspec_avg_chan(dd_dm0.T, freqs, avg_chan=avg_chan)
    #print(davg_dm0.shape)
    avg_tt, davg_dm0 = dspec_avg_tsamp(davg_dm0, tt, avg_tsamp=avg_tsamp)
    #print(davg_dm0.shape)
    dt_avg = dt * avg_tsamp
    if final_reverse:
        davg = dedisperse_dspec(davg_dm0.T, dm0, avg_freqs, f0, dt_avg, reverse=True)
    else:
        davg = davg_dm0.T
    return avg_tt, avg_freqs, davg.T


def dspec_avg_chan_dm_GHz(dspec, freqs, f0, dt, avg_chan=1, dm0=0.0):
    freqs_MHz = freqs * 1e3
    f0_MHz = f0 * 1e3

    dd_dm0 = dedisperse_dspec(dspec.T, dm0, freqs_MHz, f0_MHz, dt)
    avg_freqs_MHz, davg_dm0 = dspec_avg_chan(dd_dm0.T, freqs_MHz, avg_chan=avg_chan)
    davg = dedisperse_dspec(davg_dm0.T, dm0, avg_freqs_MHz, f0_MHz, dt, reverse=True)

    avg_freqs = avg_freqs_MHz / 1e3

    return avg_freqs, davg.T


