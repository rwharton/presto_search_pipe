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
        if os.path.isfile(infile):
            print("Beam %d" %(bb))
            cc = cands_from_file(infile, bb)
            candlist += cc
        else:
            print("   Beam %d File Missing" %(bb))
    
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
    for ii, cc in enumerate(candlist):
        fname = filename_from_cand(cc, suffix=suffix)
        out_fname = "cand%04d_%s" %(ii, fname)
        infile  = "%s/beam%04d/cands_presto/%s" %(topdir, cc.beam, fname)
        outfile = "%s/%s" %(outdir, out_fname)
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


def parse_name(namecand):
    namesplit = namecand.rsplit("_", 1)
    beamstr = namesplit[0].split('_')[0]
    zval, candnum = namesplit[1].split(':')
    zval = int(zval)
    candnum = int(candnum)
    candfile = "%s.cand" %(namecand.split(':')[0])
    return beamstr, zval, candnum, candfile


def read_top_cands(infile):
    """
    """
    num = 0
    cand_list = []
    with open(infile) as fin:
        for line in fin:
            if line[0] == "#":
                continue
            else: pass

            cols = line.split()
            namecand = cols[0]
            bstr, zval, cnum, cfile = parse_name(namecand)
            dm = float(cols[1])
            cdict = {'bstr' : bstr,
                     'zval' : zval, 
                     'cnum' : cnum, 
                     'cfile': cfile,
                     'dm'   : dm, 
                     'num'  : num}
            cand_list.append(cdict)
            num += 1
    return cand_list


def run_prepfold_topcands(clist, bdir, zmax, outfile, errfile, 
                          searchtype='none', extra_args='', 
                          use_mask=False):
    """
    This function will run prepfold on the candidate files produced
    by the presto accelsearch.  This essentially replaces the
    gotocand.py script from the older version of the pipeline
    """
    if searchtype == 'fine':
        sflag = "-fine"
    elif searchtype =='coarse':
        sflag = "-coarse"
    elif searchtype == 'regular':
        sflag = ""
    else:
        sflag = "-nosearch"

    for cc in clist:
        bstr  = cc['bstr']
        zval  = cc['zval']
        cnum  = cc['cnum']
        cfile = cc['cfile']
        dm    = cc['dm']
        num   = cc['num']

        psname = "cand%04d_%s_DM%.2f_ACCEL_Cand_%d.pfd.ps"\
                 %(num, bstr, dm, cnum) 
        if os.path.exists(psname):
            print "File "+psname+" already made, skipping"
            continue
        else:
            pass

        fitsdir = "%s/%s" %(bdir, bstr)
        fitslist = glob("%s/*fits" %fitsdir)
        fitslist.sort()
        fitsfiles = ' '.join(fitslist)
        fitsfiles = "%s/*fits" %(fitsdir)
        canddir = "%s/cands_presto" %fitsdir
        if use_mask:
            maskfile = "%s/rfi_products/*.mask" %fitsdir
        else:
            maskfile = ''
        
        # Need to copy cand file from canddir to fitsdir
        orig_candfile = "%s/%s" %(canddir, cfile)
        candfile = "%s/%s" %(fitsdir, cfile)
        copyfile(orig_candfile, candfile)
        print(orig_candfile)
        print(candfile)
        # Now we can run prepfold command
        outname = "cand%04d_%s_z%d" %(num, bstr, zval)
        if ( os.path.exists(candfile) ):
            if use_mask:
                cmd = "prepfold -noxwin %s -dm %.2f " %(sflag, dm) +\
                      "-accelcand %d -accelfile %s " %(cnum, candfile) +\
                      "%s " %(extra_args) +\
                      "-noweights -noscales -nooffsets " +\
                      "-mask %s -o %s %s" %(maskfile, outname, fitsfiles)
            else:
                cmd = "prepfold -noxwin %s -dm %.2f " %(sflag, dm) +\
                      "-accelcand %d -accelfile %s " %(cnum, candfile) +\
                      "%s " %(extra_args) +\
                      "-noweights -noscales -nooffsets " +\
                      "-o %s %s" %(outname, fitsfiles)
        
            try_cmd(cmd, stdout=outfile, stderr=errfile)

        else:
            print "Could not find %s" %candfile
            print "and/or         %s" %datfile

        # Now we remove the copied candfile
        os.remove(candfile)

    return



############
### MAIN ###
############
#sys.exit()
group = '57519'
cand_dir = "/hercules/results/rwharton/fastvis_gc/proc/%s/search/cand_files" %(group)
candfile_template = "beam%04d.ACCEL_300.sifted.cands"

beamlist = np.arange(0, 1000)

#"""
candlist = cands_from_many_files(cand_dir, candfile_template, beamlist)
#out_idx = parse_harmonics(candlist, n=1)
#out_idx = parse_dupes_r(candlist, r_err=1.5)
out_idx = parse_harmonics_r(candlist, n=1, r_err=100, grow=False, verbose=False)
bad_idx = np.setdiff1d(np.arange(len(candlist)), out_idx)

good_cc = [candlist[ii] for ii in out_idx]
bad_cc  = [candlist[ii] for ii in bad_idx]

print("Good Candidates: %d" %len(good_cc))
#"""

top_dir = "/hercules/results/rwharton/fastvis_gc/proc/%s/search" %(group)
out_dir = "/hercules/results/rwharton/fastvis_gc/proc/%s/search/top_plots" %(group)

#copy_cand_plots(good_cc, top_dir, out_dir, suffix='ps')

top_name = "%s/good_cands.txt" %(out_dir)
write_filtered_cands(good_cc, top_name)
clist = read_top_cands(top_name)

run_prepfold_topcands(clist, top_dir, 300, None, None, searchtype='nosearch')

