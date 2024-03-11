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

        self.h_hits  = 0
        self.h_harms = []
        self.h_Ps    = []
        self.h_DMs   = []
        self.h_SNRs  = []

    def add_harm(self, cand):
        self.h_hits += 1
        self.h_harms.append( self.P / cand.P )
        self.h_Ps.append( cand.P )
        self.h_DMs.append( cand.DM )
        self.h_SNRs.append( cand.SNR )

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



####################
## FILE SHUTTLING ##
####################

def filename_from_cand(cand):
    fname = cand.filename.rstrip('_100')
    fnum  = cand.filecand
    return "%s_Cand_%d.pfd.png" %(fname, fnum)


def copy_cand_plots(candlist, topdir, outdir):
    for cc in candlist:
        fname = filename_from_cand(cc)
        infile  = "%s/beam%04d/cands_presto/%s" %(topdir, cc.beam, fname)
        outfile = "%s/%s" %(outdir, fname)
        copyfile(infile, outfile)
    return



################
## WRITE FILE ##
################

def write_filtered_cands(cands, outfile):
    hdr = '#                           ' +\
          'file:candnum                               ' +\
          'DM     SNR    sigma   numharm    ipow     cpow      P(ms) ' +\
          '         r         z     numhits\n'
    fout = open(outfile, 'w')
    fout.write(hdr)
    for cc in cands:
        fout.write(cc.line)
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


def run_prepfold_topcands(filenm, bdir, zmax, outfile, errfile):
    """
    This function will run prepfold on the candidate files produced
    by the presto accelsearch.  This essentially replaces the
    gotocand.py script from the older version of the pipeline
    """
    # Open candfile
    f = open(filenm, 'r')
    
    i = 0
    for line in f:
        # This just skips down to where the files are
        if line.startswith('#'):
            i = 1
            continue
        if i==1:
            namecand = line.split()[0]
            namesplit = namecand.split("_"+str(zmax)+":")
            if len(namesplit) != 2:
                continue
            else:
                bstr = namesplit[0].split('_')[0]
                psname = namesplit[0] + "_Cand_" + namesplit[1] + ".pfd.ps"
                if os.path.exists(psname):
                    print "File "+psname+" already made, skipping"
                else:
                    fitsdir = "%s/%s" %(bdir, bstr)
                    fitslist = glob("%s/*fits" %fitsdir)
                    fitslist.sort()
                    fitsfiles = ' '.join(fitslist)
                    canddir = "%s/cands_presto" %fitsdir
                    maskfile = glob("%s/rfi_products/*mask" %fitsdir)[0]
                    candfile = canddir + '/' + namecand.split(':')[0] + '.cand'
                    outname  = namecand.split('_ACCEL_')[0]
                    if ( os.path.exists(candfile) ):
                        cmd = "prepfold -noxwin -accelcand %s -accelfile %s -mask %s -o %s %s"\
                              %(namesplit[1], candfile, maskfile, outname, fitsfiles)
                        try_cmd(cmd, stdout=outfile, stderr=errfile)
                    else:
                        print "Could not find %s" %candfile
                        print "and/or         %s" %datfile
    # Close file, we're done
    f.close()
    return

############
### MAIN ###
############
#sys.exit()
cand_dir = "/hercules/results/rwharton/fastvis_gc/proc/57519/search/cand_files"
candfile_template = "beam%04d.ACCEL_100.sifted.cands"
beamlist = np.arange(300)
candlist = cands_from_many_files(cand_dir, candfile_template, beamlist)

out_idx = parse_harmonics(candlist, n=1)
good_cc = [candlist[ii] for ii in out_idx]

good_cc = parse_cands_from_many_files(cand_dir, candfile_template, beamlist, n=2)

print("Good Candidates: %d" %len(good_cc))

top_dir = "/hercules/results/rwharton/fastvis_gc/proc/57519/search"
out_dir = "/hercules/results/rwharton/fastvis_gc/proc/57519/search/top_plots"

#copy_cand_plots(good_cc, top_dir, out_dir)
