import numpy as np
from glob import glob
import copy
import itertools
from shutil import copyfile
import sys
import os
import subprocess as sp

#########################
## CANDIDATE CLASS DEF ##
#########################

class Candidate(object):
    def __init__(self, line, beam):
        self.beam     = beam
        cols          = line.split()
        self.filename = cols[0].split(":")[0]
        self.filecand = int( cols[0].split(":")[1] )
        self.DM       = float(cols[1])
        self.SNR      = float(cols[2])
        self.sigma    = float(cols[3])
        self.numharm  = int(cols[4])
        self.ipow     = float(cols[5])
        self.cpow     = float(cols[6])
        self.P        = float(cols[7])
        self.r        = float(cols[8])
        self.z        = float(cols[9])
        self.numhits  = int( cols[10].strip("()") )
        self.line     = line
        self.note     = ''

        self.h_hits  = 0
        self.h_harms = []
        self.h_Ps    = []
        self.h_DMs   = []
        self.h_SNRs  = []
        self.h_beams = []

    def add_harm(self, cand):
        self.h_hits += 1
        self.h_harms.append( self.P / cand.P )
        self.h_Ps.append( cand.P )
        self.h_DMs.append( cand.DM )
        self.h_SNRs.append( cand.SNR )
        self.h_beams.append( cand.beam )

    def __repr__(self):
        repr_string = "Cand(beam=%d, SNR=%.0f, P=%.2f, DM=%1f)" \
            %(self.beam, self.SNR, self.P, self.DM)
        return repr_string


def cmp_snr(a, b):
    return -1 * cmp(a.SNR, b.SNR)


######################
## HARMONIC TESTING ##
######################

def check_harm(P1, P2, harm_pairs, df=1e-3, verbose=0):
    nums = np.array([ hh[0] for hh in harm_pairs ])
    dens = np.array([ hh[1] for hh in harm_pairs ])
    fracs = nums / (1.0 * dens) % 1

    f1 = 1.0 / P1 
    f2 = 1.0 / P2

    hratio = f2 / f1    

    #frac_err = hratio * np.sqrt( (df/f1)**2.0 + (df/f2)**2.0 )
    frac_err = 0.01

    hint  = int(hratio)
    hfrac = hratio - hint

    nmax = np.max(dens)
    if hfrac > (1 - 1.0 / (2 * nmax)):
        hfrac = np.abs(1 - hfrac)
    else:
        pass

    res = np.abs( fracs - hfrac )

    xx = np.argmin(res)
    
    if verbose:
        print (f1, f2)
        print hint 
        print hfrac
        print res[xx]
        print frac_err
    else: pass

    if dens[xx] > 1:
        harm_pair = (hint * dens[xx] + nums[xx], dens[xx])
    else:
        harm_pair = (nums[xx], dens[xx])

    if res[xx] <= frac_err:
        #return harm_pair
        return True

    else:
        #return None
        return False


def check_harm_r(r1, r2, harm_pairs, r_err=1.1, verbose=0, grow=False):
    fac1 = np.array([ hh[1] for hh in harm_pairs ])
    fac2 = np.array([ hh[0] for hh in harm_pairs ])

    drs = np.ones(len(fac1)) * r_err
    if grow:
        drs = drs * np.sqrt(fac1 + fac2)

    xx = np.where( np.abs(fac1 * r1 - fac2 * r2) <= drs )[0]
    
    if verbose:
        print (r1, r2)
        print xx
    else: pass

    if len(xx):
        harmpair = (fac1[xx[0]], fac2[xx[0]])
        return harmpair

    else:
        return None


def check_dupes_r(r1, r2, r_err=1.1):
    dr = np.abs(r1 - r2)
    if dr  <= frac_err:
        return True

    else:
        return False


def gcd(x, y):
    while y != 0:
        (x, y) = (y, x % y)
    return x    


def generate_harmonic_pairs(n):
    pairs = [(1, 1)]
    vals = np.arange(1, n+1)

    pair_iter = itertools.combinations_with_replacement(vals, 2) 
    for pp in pair_iter:
        if pp[0] >= pp[1] or gcd(pp[0], pp[1]) > 1:
            continue
        else:
            pairs.append( pp )

    return pairs



##################
## CAND PARSING ##
##################

def cands_from_file(candfile, beam):
    candlist = []
    with open(candfile, 'r') as fin:
        for line in fin:
            if line[0] in ["#", " ", "\n"]:
                continue
            else:
                pass
            candlist.append( Candidate(line, beam) )
    return candlist


def cands_from_many_files(indir, name_tmp, beams):
    candlist = []
    for bb in beams:
        fname  = name_tmp % bb
        infile = "%s/%s" %(indir, fname)
        cc = cands_from_file(infile, bb)
        candlist += cc
    
    candlist.sort(cmp=cmp_snr)
    return candlist


def parse_cands_from_many_files(indir, name_tmp, beams, n=2):
    candlist = []
    for bb in beams:
        fname  = name_tmp % bb
        infile = "%s/%s" %(indir, fname)
        cc = cands_from_file(infile, bb)
        out_idx = parse_harmonics(cc, n=n)
        candlist += [ cc[ii] for ii in out_idx ]
 
    candlist.sort(cmp=cmp_snr)
    return candlist


def parse_harmonics(candlist, n=2, dP=1.0):
    # Generate pairs
    pairs = generate_harmonic_pairs(n)

    # Index list
    idxs = range(len(candlist))

    # Output index list
    out_idxs = []

    # While still indices to check, check them
    while(len(idxs)):
        ii = idxs.pop(0)
        cand1 = candlist[ii]
        out_idxs.append(ii)
        
        idxs2 = idxs[:]
        for idx in idxs2:
            cand2 = candlist[idx]
            if check_harm(cand1.P, cand2.P, pairs, dP):
                candlist[ii].add_harm(cand2)
                idxs.remove(idx)
            else:
                pass
    return out_idxs


def parse_dupes_r(candlist, r_err=1.1):
    # Index list
    idxs = range(len(candlist))

    # Output index list
    out_idxs = []

    # While still indices to check, check them
    while(len(idxs)):
        ii = idxs.pop(0)
        cand1 = candlist[ii]
        out_idxs.append(ii)
        
        idxs2 = idxs[:]
        for idx in idxs2:
            cand2 = candlist[idx]
            if np.abs(cand1.r - cand2.r) <= r_err:
                candlist[ii].add_harm(cand2)
                idxs.remove(idx)
            else:
                pass
    return out_idxs


def parse_harmonics_r(candlist, n=2, r_err=1.1, grow=False, verbose=False):
    # Generate pairs
    pairs = generate_harmonic_pairs(n)
    #pairs = [pp for pp in pairs if pp[0] == 1]

    # Index list
    idxs = range(len(candlist))

    # Output index list
    out_idxs = []

    # While still indices to check, check them
    while(len(idxs)):
        ii = idxs.pop(0)
        cand1 = candlist[ii]
        out_idxs.append(ii)
        
        idxs2 = idxs[:]
        for idx in idxs2:
            cand2 = candlist[idx]
            hh = check_harm_r(cand1.r, cand2.r, pairs, r_err=r_err,
                              verbose=verbose, grow=grow)
            if hh is not None:
                candlist[ii].add_harm(cand2)
                idxs.remove(idx)
                dr = hh[0] * cand1.r - hh[1] * cand2.r
                cand2.note = "    -> (%d, %d)  dr = %.3f   (%s:%d)\n"\
                             %(hh[0], hh[1], dr, cand1.filename, cand1.filecand)
            else:
                pass
    return out_idxs

####################
## FILE SHUTTLING ##
####################

def filename_from_cand(cand, suffix='png'):
    fname = cand.filename.rstrip('_100')
    fnum  = cand.filecand
    return "%s_Cand_%d.pfd.%s" %(fname, fnum, suffix)


def copy_cand_plots(candlist, topdir, outdir, suffix='png'):
    for cc in candlist:
        fname = filename_from_cand(cc, suffix=suffix)
        infile  = "%s/beam%04d/cands_presto/%s" %(topdir, cc.beam, fname)
        outfile = "%s/%s" %(outdir, fname)
        copyfile(infile, outfile)
    return



################
## WRITE FILE ##
################

def write_filtered_cands(cands, outfile, use_note=False):
    hdr = '#                           ' +\
          'file:candnum                               ' +\
          'DM     SNR    sigma   numharm    ipow     cpow      P(ms) ' +\
          '         r         z     numhits\n'
    fout = open(outfile, 'w')
    fout.write(hdr)
    for cc in cands:
        fout.write(cc.line)
        if use_note:
            fout.write(cc.note)
        else:
            pass
    fout.close()
    return



####################
## FOLD NEW CANDS ##
####################

def try_cmd(cmd, stdout=None, stderr=None):
    """
    Run the command in the string cmd using sp.check_call().  If there
    is a problem running, a CalledProcessError will occur and the
    program will quit.
    """
    print "\n\n %s \n\n" %cmd
    try:
        retval = sp.check_call(cmd, shell=True, stdout=stdout, stderr=stderr)
    except sp.CalledProcessError:
        print("The command:\n %s \ndid not work, quitting..." %cmd)
        sys.exit(0)


def run_prepfold_topcands(filenm, fits_dir, search_dir, zmax, outfile, errfile, 
                          searchtype='none', start=None, stop=None):
    """
    This function will run prepfold on the candidate files produced
    by the presto accelsearch.  This essentially replaces the
    gotocand.py script from the older version of the pipeline
    """
    # Open candfile
    f = open(filenm, 'r')

    if searchtype == 'fine':
        sflag = "-fine"
    elif searchtype =='coarse':
        sflag = "-coarse"
    elif searchtype == 'regular':
        sflag = ""
    else:
        sflag = "-nosearch"

    nsub = 38

    # FITS FILES
    fitsfiles = "%s/*fits" %(fits_dir)

    # CAND DIR
    canddir = "%s/cands_presto" %(search_dir)

    # MASKFILE
    maskfile = "%s/rfi_products/*.mask" %(search_dir)
   
    ii = 0
    for line in f:
        # This just skips down to where the files are
        if line.startswith('#'):
            continue
        else:
            pass
        ii += 1
        if start is not None:
            if ii >= start:
                pass
            else: 
                continue
        if stop is not None:
            if ii < stop:
                pass
            else:
                continue
        namecand = line.split()[0]
        #namesplit = namecand.split("_"+str(zmax)+":")
        namesplit = namecand.rsplit("_", 1)
        if len(namesplit) != 2:
            continue
        else:
            bstr = namesplit[0].split('_')[0]
            zval, cnum = namesplit[1].split(':')
            dm = float( line.split()[1] )
            outname = "cand%03d_"%ii + namecand.split('_ACCEL_')[0] + "_z%s" %zval
            psname  = "%s_ACCEL_Cand%s.pfd.ps" %(outname, cnum)
            if os.path.exists(psname):
                print "File "+psname+" already made, skipping"
            else:
                # Need to copy cand file from canddir to fitsdir
                cand_fname = "%s.cand" %(namecand.split(':')[0])
                orig_candfile = "%s/%s" %(canddir, cand_fname)
                candfile = "%s/%s" %(search_dir, cand_fname)
                if os.path.exists(orig_candfile):
                    copyfile(orig_candfile, candfile)
                    print(orig_candfile)
                    print(candfile)
                    # Now we can run prepfold command
                    cmd = "prepfold -noxwin %s -dm %.2f -nsub %d -accelcand %s -accelfile %s -mask %s -o %s %s"\
                          %(sflag, dm, nsub, cnum, candfile, maskfile, outname, fitsfiles)
                    try_cmd(cmd, stdout=outfile, stderr=errfile)
                    print cmd
                    os.remove(candfile)
                else:
                    print "Could not find %s" %candfile

    # Close file, we're done
    f.close()
    return

############
### MAIN ###
############
#sys.exit()
cands_dir = "/hercules/results/rwharton/phased_gc/day1/search_red_nomag/cands_presto"
best_dir = "%s/top_plots" %(cands_dir)
cand_file = "%s/sgra_day1.ACCEL_300.sifted.cands" %cands_dir
good_file = "%s/good_cands.txt" %(cands_dir)

start = int(sys.argv[-2])
stop  = int(sys.argv[-1])

#candlist = cands_from_file(cand_file, 0)

"""
#out_idx = parse_harmonics(candlist, n=1)
#out_idx = parse_dupes_r(candlist, r_err=1.5)
out_idx = parse_harmonics_r(candlist, n=1, r_err=100, grow=False, verbose=False)
bad_idx = np.setdiff1d(np.arange(len(candlist)), out_idx)

good_cc = [candlist[ii] for ii in out_idx]
bad_cc  = [candlist[ii] for ii in bad_idx]

print("Good Candidates: %d" %len(good_cc))
write_filtered_cands(good_cc, good_file, use_note=False)
"""

fits_dir   = "/hercules/results/rwharton/phased_gc/day1/fix_fits"
search_dir = "/hercules/results/rwharton/phased_gc/day1/search_red_nomag"

run_prepfold_topcands(good_file, fits_dir, search_dir, 300, None, None, searchtype='none', start=start, stop=stop)

