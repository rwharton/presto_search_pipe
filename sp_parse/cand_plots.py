import numpy as np
from collections import namedtuple, defaultdict
from glob import glob
from scipy.misc import factorial

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import Normalize
from matplotlib.colors import BoundaryNorm
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import LinearSegmentedColormap

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u

def get_attr(clist, attr):
    if hasattr(clist[0], attr):
        return np.array([ getattr(cc, attr) for cc in clist ])
    else:
        print("Object has no attribute: %s" %attr)
        return


def bin_beams(candlist, beams):
    cand_beams = get_attr(candlist, 'beam')
    d = defaultdict(int)
    for cb in cand_beams:
        d[cb] += 1
    counts = np.array([ d.get(bb, 0) for bb in beams ])
    return counts


def expand_bin_beams(coords, counts, beams, sum_radius=3.0):
    d = defaultdict(int)
    n = []
    for ii, cc in enumerate(coords):
        print(ii)
        dsep = cc.separation(coords).arcsec
        xx = np.where( dsep < sum_radius )[0]
        n.append( len(xx) )
        d[beams[ii]] += np.sum(counts[xx])
    sum_counts = np.array([ d.get(bb, 0) for bb in beams ])    
    return sum_counts, np.array(n)


def get_beam_offsets(coords, beams, coord0=None):
    cc = coords[beams]
    if coord0 is None:
        coord0 = cc[0]
    else:
        pass
    ra_offset = ((cc.ra - coord0.ra) * np.cos(coord0.dec.to('radian'))).arcsec
    dec_offset = (cc.dec - coord0.dec).arcsec
    return ra_offset, dec_offset


def get_beam_vals(beams, bnums, bvals):
    """
    beams -- all beams

    bnums, bvals -- beam vals for beam bnums
    """
    Nb = len(beams)
    out_vals = np.zeros(Nb)

    for ii, bb in enumerate(beams):
        xx = np.where( bnums == bb )[0]
        if len(xx) == 0:
            continue
        else: 
            pass
        out_vals[ii] = np.max( bvals[xx] )

    return out_vals


def make_beam_snr_plot(coords, beams, cand, fwhm, title=None, add_label=True, 
                       add_outline=True, mark_beams=None, outname=None):
    if outname is not None:
        plt.ioff()

    print("beam  = %d" %(cand.beam))
    print("DM    = %.2f" %(cand.dm))
    print("t     = %.3f" %(cand.time))
    print("sigma = %.3f" %(cand.sigma))
    

    dupe_beams = np.array(cand.dupe_beams)
    dupe_sigs = np.array(cand.dupe_sigmas)

    bvals = get_beam_vals(beams, dupe_beams, dupe_sigs) 
    
    ra_offs, dec_offs = get_beam_offsets(coords, beams)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    maxsep = max( np.max(ra_offs), np.max(dec_offs) )

    pad = 0.6 * fwhm
    ax.set_xlim( maxsep + pad, -maxsep - pad )
    ax.set_ylim(-maxsep - pad,  maxsep + pad)

    # define the colormap
    cmap = plt.cm.magma_r
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # create the new map
    cmap = LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)

    # define the bins and normalize
    #bounds = np.arange(-0.5, Ncbins + 0.5, 1)
    #norm = BoundaryNorm(bounds, cmap.N)
    norm = Normalize(vmin=0,vmax=np.max(bvals))

    for ii, bb in enumerate(beams):
        cc = (ra_offs[ii], dec_offs[ii])
        bcolor = cmap(norm(bvals[ii]))
        bcirc = Circle(cc, 0.5 * fwhm, fill=True, color=bcolor, alpha=0.9, lw=0)
        ax.add_artist(bcirc)
        if add_label:
            ax.text(cc[0], cc[1], "%d" %(bb), va='center', ha='center', fontsize=12)

    if add_outline:
        for ii, bb in enumerate(beams):
            cc = (ra_offs[ii], dec_offs[ii])
            bcirc = Circle(cc, 0.5 * fwhm, fill=False, color='0.5', lw=1)
            ax.add_artist(bcirc)

    if mark_beams is not None:
        for ii in mark_beams:
            cc = (ra_offs[ii], dec_offs[ii])
            bcirc = Circle(cc, 0.5 * fwhm, fill=False, color='r', lw=2)
            ax.add_artist(bcirc)        

    cbax = fig.add_axes([0.88, 0.15, 0.05, 0.7])
    #cb1 = ColorbarBase(cbax, cmap=cmap, norm=norm, ticks=np.arange(0, Ncbins, 1))
    cb1 = ColorbarBase(cbax, cmap=cmap, norm=norm)

    ax.set_aspect('equal')
    ax.set_xlabel("RA Offset (arcsec)")
    ax.set_ylabel("Dec Offset (arcsec)")

    if title is not None:
        ax.set_title(title)

    if outname is not None:
        plt.savefig(outname, dpi=150, bbox_inches='tight')
        plt.close()
        plt.ion()

    else:
        plt.show()
    return


def make_many_snr_plots(coords, beams, candlist, fwhm, 
                        title=None, add_label=True, 
                        add_outline=True, outname=None):
    for ii, cc in enumerate(candlist):
        #outname = "cand%04d_beams.png" %(cc.cnum)
        outname = "Cand%04d_beams.png" %(ii)
        bnum = cc.beam
        dm   = cc.dm
        tt   = cc.time
        sig  = cc.sigma

        tstr = "Cand %04d, DM=%.2f, tt=%.2f, sigma=%.1f" %(\
                cc.cnum, dm, tt, sig)

        make_beam_snr_plot(coords, beams, cc, fwhm, title=tstr, 
                add_label=add_label, add_outline=add_outline, 
                mark_beams=[bnum], outname=outname)

    return


def get_skycoords(beams):
    ra_strs  = beams[:, 1]
    dec_strs = beams[:, 2]
    dec_strs = np.array([ dd.replace('.', ':', 2) for dd in dec_strs ])
    coords = SkyCoord(ra=ra_strs, dec=dec_strs, unit=(u.hour, u.deg))
    return coords


#######################
### TESTING / DEBUG ###
#######################


def beam_summary(cand, cands, bcands, beams, coords, fwhm, add_label=False,
                 add_outline=True, mark_beams=None, outfile=None):
    ## HIDE PLT IF NEC ##
    if outfile is not None:
        plt.ioff()

    ## SETUP FIGURE ##
    fig_x = 14.
    fig_y = 14.
    fig = plt.figure(figsize=(fig_x, fig_y))

    # Setting up axes
    bpad = 0.05
    tpad = 0.05 
    lpad = 0.1
    rpad = 0.05
    hspace = 0.075
    vspace = 0.075 
    cpad = 0.075

    #= Widths =#
    ww_t0 = 1.0 - (lpad + rpad)
    ww_ti = (1.0 - (lpad + rpad + 2 * hspace)) / 3.0
    ww_bp = (1.0 - (lpad + rpad + hspace + cpad)) * (0.6)
    ww_dm = (1.0 - (lpad + rpad + hspace + cpad)) * (1-0.6)

    #= Heights =#
    hh_bp = ww_bp / (fig_y / fig_x)
    hh_t0 = (1.0 - (tpad + bpad + 2 * vspace + hh_bp)) / 2.0
    hh_ti = hh_t0
    hh_dm = hh_bp

    #= TOP PLOT =#
    x_t0 = lpad
    y_t0 = 1.0 - (tpad + hh_t0)
    ax_t0 = fig.add_axes([x_t0, y_t0, ww_t0, hh_t0])

    #= Next row =#
    x_t1 = lpad 
    x_t2 = lpad + hspace + ww_ti
    x_t3 = lpad + 2 * hspace + 2 * ww_ti 
    y_ti = 1.0 - (tpad + hh_t0 + vspace + hh_ti)
    ax_t1 = fig.add_axes([x_t1, y_ti, ww_ti, hh_ti])
    ax_t2 = fig.add_axes([x_t2, y_ti, ww_ti, hh_ti])
    ax_t3 = fig.add_axes([x_t3, y_ti, ww_ti, hh_ti])

    #= Beam Plot =#
    x_bp = 1.0 - (rpad + ww_bp + cpad)
    y_bp = bpad 
    ax_bp = fig.add_axes([x_bp, y_bp, ww_bp, hh_bp])
    
    #= Color Bar =#
    x_cb = 1.0 - (rpad + cpad) + 0.01
    y_cb = y_bp
    cbax = fig.add_axes([x_cb, y_cb, cpad * 0.5, hh_bp])

    #= time DM dupes =#
    x_dm = lpad
    y_dm = bpad
    ax_dm = fig.add_axes([x_dm, y_dm, ww_dm, hh_dm])


    #################################
    ###    BEAM PLOT    (ax_bp)   ###
    #################################

    dupe_beams = np.array(cand.dupe_beams)
    dupe_sigs = np.array(cand.dupe_sigmas)

    bvals = get_beam_vals(beams, dupe_beams, dupe_sigs) 
    
    ra_offs, dec_offs = get_beam_offsets(coords, beams)

    maxsep = max( np.max(ra_offs), np.max(dec_offs) )

    pad = 0.6 * fwhm
    ax_bp.set_xlim( maxsep + pad, -maxsep - pad )
    ax_bp.set_ylim(-maxsep - pad,  maxsep + pad)

    # define the colormap
    cmap = plt.cm.magma_r
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # create the new map
    cmap = LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)

    # define the bins and normalize
    #bounds = np.arange(-0.5, Ncbins + 0.5, 1)
    #norm = BoundaryNorm(bounds, cmap.N)
    norm = Normalize(vmin=0,vmax=np.max(bvals))

    for ii, bb in enumerate(beams):
        cc = (ra_offs[ii], dec_offs[ii])
        bcolor = cmap(norm(bvals[ii]))
        bcirc = Circle(cc, 0.5 * fwhm, fill=True, color=bcolor, alpha=0.9, lw=0)
        ax_bp.add_artist(bcirc)
        if add_label:
            ax_bp.text(cc[0], cc[1], "%d" %(bb), va='center', ha='center', 
                       fontsize=12)

    if add_outline:
        for ii, bb in enumerate(beams):
            cc = (ra_offs[ii], dec_offs[ii])
            bcirc = Circle(cc, 0.5 * fwhm, fill=False, color='0.5', lw=1)
            ax_bp.add_artist(bcirc)

    mark_beams = [ cand.beam ]
    if mark_beams is not None:
        for ii in mark_beams:
            cc = (ra_offs[ii], dec_offs[ii])
            bcirc = Circle(cc, 0.5 * fwhm, fill=False, color='r', lw=2)
            ax_bp.add_artist(bcirc)        

    ### MARK SGRA + MAG ? ###
    xx_sgr = np.where( beams == 0 )[0]
    xx_mag = np.where( beams == 2 )[0]

    if len(xx_sgr):
        ra_sgr  = ra_offs[  xx_sgr[0] ]
        dec_sgr = dec_offs[ xx_sgr[0] ]
        ax_bp.plot( ra_sgr, dec_sgr, ls='', marker='+', c='0.5', mew=2)

    if len(xx_mag):
        ra_mag  = ra_offs[  xx_mag[0] ]
        dec_mag = dec_offs[ xx_mag[0] ]
        ax_bp.plot( ra_mag, dec_mag, ls='', marker='x', c='0.5', mew=2)


    #cb1 = ColorbarBase(cbax, cmap=cmap, norm=norm, ticks=np.arange(0, Ncbins, 1))
    cb1 = ColorbarBase(cbax, cmap=cmap, norm=norm)

    ax_bp.set_aspect('equal')
    ax_bp.set_xlabel("RA Offset (arcsec)")
    ax_bp.set_ylabel("Dec Offset (arcsec)")


    ################################ 
    ### TIME - DM - SNR  (ax_t0) ###
    ################################
    
    bb_dms  = get_attr(bcands, 'dm')
    bb_snrs = get_attr(bcands, 'sigma')
    bb_tts  = get_attr(bcands, 'time')

    cc_dms  = get_attr(cands, 'dm')
    cc_snrs = get_attr(cands, 'sigma')
    cc_tts  = get_attr(cands, 'time')
    cc_beams = get_attr(cands, 'beam')

    xx_beam = np.where( cc_beams == cand.beam )[0]
    print(len(xx_beam))

    # check for infs....
    xx_inf = np.isinf( bb_snrs )
    if len(xx_inf):
        bb_snrs[xx_inf] = -1
    else:
        pass

    ax_t0.axvline(x=cand.time, ls='-', c='r', alpha=0.3)
    ax_t0.axhline(y=cand.dm, ls='-', c='r', alpha=0.3)

    ax_t0.scatter(bb_tts, bb_dms, marker='o', s=(bb_snrs/3.0)**2.0, 
                  color='none', edgecolor='k')

    ax_t0.scatter(cc_tts[xx_beam], cc_dms[xx_beam], marker='o', 
                  s=(cc_snrs[xx_beam]/3.0)**2.0, color='r', edgecolor='r')
    
    ax_t0.scatter(cc_tts[xx_beam], cc_dms[xx_beam], marker='o', 
                  s=20**2.0, color='none', edgecolor='r', 
                  linewidths=2)
    
    t0_trange = np.max(bb_tts) - np.min(bb_tts)
    t0_tpad   = 0.05 * t0_trange
    tlims = ( np.min(bb_tts) - t0_tpad, np.max(bb_tts) + t0_tpad )

    t0_dmrange = np.max(bb_dms) - np.min(bb_dms) 
    t0_dmpad   = 0.05 * t0_dmrange
    dmlims = ( np.min(bb_dms) - t0_dmpad, np.max(bb_dms) + t0_dmpad )

    ax_t0.set_xlim(tlims)
    ax_t0.set_ylim(dmlims)

    ax_t0.set_xlabel("Time (s)")
    ax_t0.set_ylabel("DM (pc/cc)")

    tstr = "Cand %04d, DM=%d pc/cc, tt=%.2f s, SNR=%.1f, Beam=%04d" %(\
                    cand.cnum, cand.dm, cand.time, cand.sigma, cand.beam)
    ax_t0.set_title(tstr)


    ################################### 
    ### BEAM SNR HISTOGRAM  (ax_t1) ###
    ###################################
    sig_max = min( np.max(bb_snrs), 500 )

    sig_range = (3, sig_max * 1.1)
    sig_bins  = np.arange(sig_range[0], sig_range[1], 1)

    #print(sig_bins)

    sig_hh, sig_be = np.histogram(bb_snrs, bins=sig_bins)
    ax_t1.plot(sig_be[:-1], sig_hh, ls='-', drawstyle='steps-post', c='k')
    ax_t1.fill_between(sig_be[:-1], y1=sig_hh, y2=0.01, step='post', 
                       color='k', alpha=0.3)
    ax_t1.axvline(x=cand.sigma, ls='-', c='r', alpha=0.3)
    ax_t1.semilogy(nonposy='clip')
    ax_t1.set_ylim(0.3)
    ax_t1.set_xlim(sig_range)

    ax_t1.set_xlabel("SNR")
    ax_t1.set_ylabel("Counts")


    ################################### 
    ### BEAM DM HISTOGRAM  (ax_t2)  ###
    ###################################
    
    dm_min  = 0
    dm_max  = 5000
    dm_step = 20
    dm_bins  = np.arange(dm_min, dm_max + dm_step, dm_step)

    dm_hh, dm_be = np.histogram(bb_dms, bins=dm_bins)
    dm_bb = dm_be[:-1] - cand.dm
    ax_t2.plot(dm_bb, dm_hh, ls='-', drawstyle='steps-post', c='k')
    ax_t2.fill_between(dm_bb, y1=dm_hh, y2=0.01, step='post', 
                       color='k', alpha=0.3)
    ax_t2.axvline(x=0.0, ls='-', c='r', alpha=0.3)
    #ax_t2.semilogy(nonposy='clip')
    ax_t2.set_ylim(0)
    #ax_t2.set_xlim( (dm_min, dm_max))
    ax_t2.set_xlim( -20 * dm_step, +20 * dm_step )
    ax_t2.tick_params(axis='x', labelsize=10)

    ax_t2.set_xlabel("DM - %d" %cand.dm)
    ax_t2.set_ylabel("Counts")


    ################################### 
    ### BEAM DM vs SNR  (ax_t3)     ###
    ###################################
    
    dm_min  = 0
    dm_max  = 5000
   
    ax_t3.scatter(bb_dms, bb_snrs, marker='o', color='none', 
                  edgecolor='k')
    ax_t3.scatter(cand.dm, cand.sigma, marker='o', color='r', 
                  edgecolor='r', s=100)
    ax_t3.axvline(x=cand.dm, ls='-', c='r', alpha=0.3)
    ax_t3.axhline(y=cand.sigma, ls='-', c='r', alpha=0.3)
    #ax_t3.semilogy(nonposy='clip')
    #ax_t3.set_ylim(0.3)
    #ax_t3.set_xlim( (dm_min, dm_max))
    ax_t3.set_xlim( ( cand.dm - 20 * dm_step, cand.dm + 20 * dm_step) )

    ax_t3.set_xlabel("DM")
    ax_t3.set_ylabel("SNR")


    ################################### 
    ### DUPE TIME vs DM (ax_dm)     ###
    ###################################
    
    xx_p = np.where( get_attr(bcands, 'line') == cand.line )[0][0]
    
    bcand  = bcands[xx_p]
    print(bcand)
    bd_tts  = np.array( bcand.dupe_times )
    bd_dms  = np.array( bcand.dupe_dms )
    bd_snrs = np.array( bcand.dupe_sigmas )

    cc_tts  = np.array( cand.dupe_times ) 
    cc_dms  = np.array( cand.dupe_dms )
    cc_snrs = np.array( cand.dupe_sigmas )

    ax_dm.scatter(cc_tts - cand.time, cc_dms, marker='o', s=(cc_snrs/3.0)**2.0, 
                  color='none', edgecolor='k')
    
    ax_dm.scatter(bd_tts - cand.time, bd_dms, marker='o', s=(bb_snrs/3.0)**2.0, 
                  color='none', edgecolor='r')

    #ax_dm.scatter(0.0, cand.dm, marker='+', s=100, color='k', linewidth=2)

    ax_dm.axvline(x=0, lw=2, c='k', alpha=0.3)
    ax_dm.axhline(y=cand.dm, lw=2, c='k', alpha=0.3)

    ax_dm.set_ylim( (dm_min, dm_max) )
    
    ax_dm.set_xlabel("Time Offset (s)")
    ax_dm.set_ylabel("DM (pc/cc)")

    if outfile is not None:
        plt.savefig(outfile, dpi=150, bbox_inches='tight')
        plt.close()
        plt.ion()

    else:
        plt.show()

    return


def summary_beam_from_candlist(cands_dir, coords, beams, candlist, fwhm, 
                               add_label=False, add_outline=True):
    for ii, cc in enumerate(candlist):
        outname = "Cand%04d_beams.png" %(ii)

        bfile = "%s/beam%04d_dmsift.npy" %(cands_dir, cc.beam)
        bcands = np.load(bfile)
        
        beam_summary(cc, candlist, bcands, beams, coords, fwhm, 
                     add_label=add_label, add_outline=add_outline, 
                     outfile=outname)

    return




## MAIN ##

"""
beam_locs = np.load(beamfile)

N = len(beam_locs) #- 24 * 7
#N = 2000
beam_locs = beam_locs[:N]
beam_nums = np.arange(len(beam_locs), dtype='int')
coords    = get_skycoords(beam_locs)
"""
