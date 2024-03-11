import os
import sys
import time
import remote_phased_params as params
import sifting
import re
import shlex
import psr_utils
import readhdr
import subprocess as sp
import multiprocessing as mp
from glob import glob
import shutil

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
    return


def accelsift(filenm, zmax, wmax):
    """
    This function is a translation of the PRESTO code ACCEL_sift.py
    so that it can be more easily incorporated into our code.  It 
    sifts through the ACCEL cands, making cuts on various parameters
    and removing duplicates and harmonics.
    """
    # Set a bunch of parameters from our params.py file
    min_num_DMs             = params.min_num_DMs
    low_DM_cutoff           = params.low_DM_cutoff
    sifting.sigma_threshold = params.sigma_threshold
    sifting.c_pow_threshold = params.c_pow_threshold
    sifting.known_birds_p   = params.known_birds_p
    sifting.known_birds_f   = params.known_birds_f
    sifting.r_err           = params.r_err
    sifting.short_period    = params.short_period
    sifting.long_period     = params.long_period
    sifting.harm_pow_cutoff = params.harm_pow_cutoff

    # Try to read the .inf files first, as _if_ they are present, all of    
    # them should be there.  (if no candidates are found by accelsearch     
    # we get no ACCEL files... 
    inffiles = glob("*.inf")

    # Check to see if this is from a short search 
    if len(re.findall("_[0-9][0-9][0-9]M_" , inffiles[0])):
        dmstrs = [x.split("DM")[-1].split("_")[0] for x in candfiles]
    else:
        dmstrs = [x.split("DM")[-1].split(".inf")[0] for x in inffiles]
    dms = map(float, dmstrs)
    dms.sort()
    dmstrs = ["%.2f"%x for x in dms]

    candfiles = glob("*ACCEL_%d_JERK_%d" %(zmax, wmax))

    # Read in candfiles
    cands = sifting.read_candidates(candfiles)
    # Remove candidates that are duplicated in other ACCEL files
    if len(cands):
        cands = sifting.remove_duplicate_candidates(cands)
    # Remove candidates with DM problems
    if len(cands):
        cands = sifting.remove_DM_problems(cands, min_num_DMs, dmstrs, low_DM_cutoff)

    # Remove candidates that are harmonically related to each other
    # Note:  this includes only a small set of harmonics
    if len(cands):
        cands = sifting.remove_harmonics(cands)

    # Write candidates to STDOUT
    if len(cands):
        cands.sort(sifting.cmp_snr)
        sifting.write_candlist(cands, candfilenm=filenm)

    Ncands = len(cands)
    return Ncands


def run_prepfold(filenm, outfile, errfile):
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
            namesplit = namecand.rsplit(":")
            if len(namesplit) != 2:
                continue
            else:
                bname = namesplit[0]
                cnum  = namesplit[1]
                psname = bname + "_Cand_" + cnum + ".pfd.ps"
                if os.path.exists(psname):
                    print "File "+psname+" already made, skipping"
                else:
                    candfile = namecand.split(':')[0] + '.cand'
                    datfile  = namecand.split('_ACCEL_')[0] + '.dat'
                    outname  = namecand.split('_ACCEL_')[0]
                    if ( os.path.exists(candfile) and os.path.exists(datfile) ):
                        cmd = "prepfold -noxwin -accelcand %s -accelfile %s -o %s %s"\
                                %(cnum, candfile, outname, datfile)
                        try_cmd(cmd, stdout=outfile, stderr=errfile)
                    else:
                        print "Could not find %s" %candfile
                        print "and/or         %s" %datfile
    f.close()
    return





basename = "sgra_day1"
zmax     = 300
wmax     = 900

accel_out = file('accelsearch.out', 'w')
accel_err = file('accelsearch.err', 'w')

candsfilenm = "%s.ACCEL_%d_JERK_%d.sifted.cands" %(basename, zmax, wmax)

Ncands = accelsift(candsfilenm, zmax, wmax)

if Ncands:
    run_prepfold(candsfilenm, accel_out, accel_err)

