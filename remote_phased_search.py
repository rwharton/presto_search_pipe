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

class Timer:
    def __init__(self):
        self.rfifind = 0.0
        self.prepsubband = 0.0
        self.realfft = 0.0
        self.zap = 0.0
        self.accelsearch = 0.0
        self.presto_sp = 0.0
        self.total = 0.0

    def print_summary(self):
        print "****************************************************"
        print "                  TIME SUMMARY                      "
        print "****************************************************"
        print "\n"
        print "Program:                         Running Time (min): "
        print "--------                         -----------------  "
        print "rfifind                              %.2f" %(self.rfifind/60.)
        print "prepsubband                          %.2f" %(self.prepsubband/60.)
        print "realfft                              %.2f" %(self.realfft/60.)
        print "zap                                  %.2f" %(self.zap/60.)
        print "accelsearch                          %.2f" %(self.accelsearch/60.)
        print "single_pulse_search                  %.2f" %(self.presto_sp/60.)
        print "\n"
        print "Total Runtime = %.2f min" %(self.total/60.)

    def write_summary(self, outfile):
        fout = open(outfile, 'w')
        fout.write( "****************************************************\n")
        fout.write( "                  TIME SUMMARY                      \n)")
        fout.write( "****************************************************\n)")
        fout.write( "\n"                                                     )
        fout.write( "Program:                         Running Time (min): \n")
        fout.write( "--------                         -----------------  \n")
        fout.write( "rfifind                              %.2f\n" %(self.rfifind/60.))
        fout.write( "prepsubband                          %.2f\n" %(self.prepsubband/60.))
        fout.write( "realfft                              %.2f\n" %(self.realfft/60.))
        fout.write( "zap                                  %.2f\n" %(self.zap/60.))
        fout.write( "accelsearch                          %.2f\n" %(self.accelsearch/60.))
        fout.write( "single_pulse_search                  %.2f\n" %(self.presto_sp/60.))
        fout.write( "\n"                                                   )
        fout.write( "Total Runtime = %.2f min\n" %(self.total/60.))
        fout.close()


def check_dependencies(work_dir, fits_dir, fitsbase):
    """
    Since we allow for some processing steps to be skipped, we need
    to check to make sure that all depedencies are still in place.

    This function checks for these dependencies and exits with a 
    descriptive error if they do not exist
    """
    # Print to screen what processing steps have been selected
    print "The following processing steps have been selected:\n"
    if params.do_rfifind:
        print "   - PRESTO rfifind (RFI mitigation tools)"
    if params.do_prepsub:
        print "   - PRESTO prepsubband (dedispersion)"
    if params.do_candsearch:
        print "   - PRESTO acceleration search and candidate sifting"
    if params.do_presto_sp:
        print "   - PRESTO singlepulse search (singlepulse.py)"
    # Print to screen what processing steps are being skipped
    print "\nThe following processing steps are being skipped:\n"
    if params.do_rfifind == 0:
        print "   - PRESTO rfifind (RFI mitigation tools)"
    if params.do_prepsub == 0:
        print "   - PRESTO prepsubband (dedispersion)"
    if params.do_candsearch == 0:
        print "   - PRESTO acceleration search and candidate sifting"
    if params.do_presto_sp == 0:
        print "   - PRESTO singlepulse search (singlepulse.py)"
    print "\nChecking dependencies...\n"
    # There must be at least one .fits file in the fits directory
    # if we want to do rfifind or prepsubband
    fl = glob(fits_dir + '/%s*.fits' %fitsbase)
    if len(fl):
        print "  Found %d file(s) in %s:\n" %(len(fl), fits_dir)
        for i in fl:
            print "    %s\n" %(i.split('/')[-1])
    else:
        print "  No %s*.fits files found in %s !\n    Exiting...\n" %(fitsbase, fits_dir)
        if (params.do_rfifind or params.do_prepsub):
            print "  ... Need FITS files for selected mode(s) \n"
            print "  ... Exiting .... \n"
            sys.exit(0)
        else:
            print "  ... but we don\'t need them! \n"
            print "  ... full steam ahead! \n"
    
    # If skipping the RFIFIND step in processing but want to do
    # processing steps further down the line, then there must be a
    # rfi_products folder in the results directory with a .mask file
    # in it
    if params.do_rfifind == 0 and params.use_mask and \
            (params.do_prepsub or params.do_candsearch or params.do_presto_sp):
        mlist = glob(work_dir + '/rfi_products/*.mask')
        if len(mlist):
            print "  Using RFI .mask:\n    %s\n" %(mlist[0])
        else:
            print "  No RFI .mask found in %s/rfi_products!\n    Exiting...\n"\
                %(work_dir)
            sys.exit(0)
   
    # If skipping the PREPSUBBAND step in processing but want to
    # do processing steps further down the line, then there must be
    # de-dispersed time series files in the results directory of
    # the form basename*DM*.dat and basename*DM*.inf
    if params.do_prepsub == 0 and (params.do_candsearch or 
                                  params.do_presto_sp):
        dats = glob(work_dir + '/*DM*dat')
        infs = glob(work_dir + '/*DM*inf')
        if not (len(dats) and len(infs)):
            print "  No .dat and/or .inf files in %s!\n    Exiting...\n" %(work_dir)
            sys.exit(0)
    # If we haven't exited by now, then things should be good
    print "\nLooks good...\n\n"
    # Pause for a few seconds so you can actually read the output
    time.sleep(5)

def format_name(name_dir):
    """
    Remove trailing '/' on path names (for consistency)
    """
    if(name_dir.endswith('/')):
        name_dir = name_dir.rstrip('/')
    return(name_dir)

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

def run_rfifind(fitslist, fitsname, work_dir):
    print("Running rfifind on the psrfits files")
    t_rfi_start = time.time()

    fitsfiles = ' '.join(fitslist) 

    # Get flag values from params file
    rfi_time = params.rfi_time
    tsig     = params.time_sig
    fsig     = params.freq_sig
    chanfrac = params.chan_frac
    intfrac  = params.int_frac
    other_flags = params.rfi_otherflags

    cmd = 'rfifind -psrfits -o %s -time %d -timesig %f -freqsig %f '\
          '-chanfrac %f -intfrac %f %s %s' %(fitsname, rfi_time, tsig, fsig,
                                             chanfrac, intfrac, other_flags,
                                             fitsfiles)
    try_cmd(cmd)
    # If no directory exists for the rfi products, make one
    # and move them over
    rfi_dir = work_dir+'/rfi_products'
    if not os.path.exists(rfi_dir):
        os.makedirs(rfi_dir)

    cmd = 'mv ./*rfifind* '+rfi_dir
    try_cmd(cmd)

    t_rfi_end = time.time()
    rfi_time = (t_rfi_end-t_rfi_start)
    print("RFI Flagging took %f hours" %(rfi_time/3600.))
    return rfi_time

def run_prepsubband(basename, maskname, fitslist, dmlow=params.dmlow, \
                        ddm=params.ddm, ndm=params.dmspercall, \
                        downsample=params.downsample, nsub=params.nsub):
    t_prep_start = time.time()
    fitsfiles = ' '.join(fitslist) 
    print("Dedispersing with 1st batch of DMs")
    orig_N = readhdr.get_samples(fitslist, params.dat_type)
    numout = psr_utils.choose_N(orig_N)
    print(orig_N, numout)

    other_flags = params.prep_otherflags
    if params.use_mask:
        cmd = 'prepsubband -o %s -psrfits -nsub %d -numout %d -lodm %.6f -dmstep %.6f '\
            '-numdms %d -downsamp %d %s -mask %s %s' %(basename, nsub, numout/downsample, 
                                                                 dmlow, ddm, ndm, downsample,
                                                                 other_flags, maskname, fitsfiles)
    else:
        cmd = 'prepsubband -o %s -psrfits -nsub %d -numout %d -lodm %.6f -dmstep %.6f '\
            '-numdms %d -downsamp %d  %s %s' %(basename, nsub, numout/downsample,
                                                        dmlow, ddm, ndm, downsample,
                                                        other_flags, fitsfiles)

    try_cmd(cmd)

    t_prep_end = time.time()
    dt = t_prep_end - t_prep_start
    print "De-dispersion took %.2f hours.\n" %(dt/3600.)
    return dt

def multi_call_prepsubband(basename, maskname, fitslist, dmlow=params.dmlow, \
                               ddm=params.ddm, downsample=params.downsample, \
                               dmcalls=params.dmcalls, nsub=params.nsub,  \
                               dsubDM=params.dsubDM, \
                               dmspercall=params.dmspercall):
    t_prep_start = time.time()
    fitsfiles = ' '.join(fitslist)

    orig_N = readhdr.get_samples(fitslist, params.dat_type)
    numout = psr_utils.choose_N(orig_N)
    other_flags = params.prep_otherflags

    # Downsample organization as in PRESTO dedisp.py (why?)
    sub_downsample = downsample / 2
    dat_downsample = 2
    if downsample < 2: sub_downsample = dat_downsample = 1

    print("Dedispersing using %d calls on %d subbands\n" %(dmcalls, nsub))
    for ii in xrange(dmcalls):
        subDM = dmlow + (ii+0.5)*dsubDM
        # Make subband
        if params.use_mask:
            cmd_sub = "prepsubband -o %s -sub -subdm %.2f -nsub %d -downsamp %d %s -mask %s %s" \
                %(basename, subDM, nsub, sub_downsample, other_flags, maskname, fitsfiles)
        else:
            cmd_sub = "prepsubband -o %s -sub -subdm %.2f -nsub %d -downsamp %d %s %s" \
                %(basename, subDM, nsub, sub_downsample, other_flags, fitsfiles)
        try_cmd(cmd_sub)
        
        # Get dedispersed time series
        sub_dmlow = dmlow + ii*dsubDM
        subfiles =  basename+"_DM%.2f.sub[0-9]*" %subDM
        if params.use_mask:
            cmd_dat = "prepsubband -o %s -numout %d -lodm %.2f -dmstep %.2d "\
                "-numdms %d -downsamp %d %s -mask %s %s" \
                %(basename, numout/downsample, sub_dmlow, ddm, dmspercall, dat_downsample, other_flags, maskname, subfiles)
        else:
            cmd_dat = "prepsubband -o %s -numout %d -lodm %.2f -dmstep %.2d "\
                "-numdms %d -downsamp %d %s %s" \
                %(basename, numout/downsample, sub_dmlow, ddm, dmspercall, dat_downsample, other_flags, subfiles)
        try_cmd(cmd_dat)
    
    t_prep_end = time.time()
    dt = t_prep_end - t_prep_start
    print "De-dispersion took %.2f hours.\n" %(dt/3600.)
    return dt


def accelsift_old(filenm):
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
    candfiles = glob("*ACCEL_" + str(params.zmax))
    # Check to see if this is from a short search                                         
    if len(re.findall("_[0-9][0-9][0-9]M_" , inffiles[0])):
        dmstrs = [x.split("DM")[-1].split("_")[0] for x in candfiles]
    else:
        dmstrs = [x.split("DM")[-1].split(".inf")[0] for x in inffiles]
    dms = map(float, dmstrs)
    dms.sort()
    dmstrs = ["%.2f"%x for x in dms]

    # Read in all the candidates
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


def accelsift(filenm):
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
    
    lo_candfiles = glob("*ACCEL_0")
    hi_candfiles = glob("*ACCEL_" + str(params.zmax))

    # Read in the lo-z candidates
    lo_cands = sifting.read_candidates(lo_candfiles)
    # Remove candidates that are duplicated in other ACCEL files
    if len(lo_cands):
        lo_cands = sifting.remove_duplicate_candidates(lo_cands)
    # Remove candidates with DM problems
    if len(lo_cands):
        lo_cands = sifting.remove_DM_problems(lo_cands, min_num_DMs, dmstrs, low_DM_cutoff)
    
    # Read in the hi-z candidates
    hi_cands = sifting.read_candidates(hi_candfiles)
    # Remove candidates that are duplicated in other ACCEL files
    if len(hi_cands):
        hi_cands = sifting.remove_duplicate_candidates(hi_cands)
    # Remove candidates with DM problems
    if len(hi_cands):
        hi_cands = sifting.remove_DM_problems(hi_cands, min_num_DMs, dmstrs, low_DM_cutoff)

    all_cands = lo_cands + hi_cands
    # Remove candidates that are harmonically related to each other
    # Note:  this includes only a small set of harmonics
    if len(all_cands):
        all_cands = sifting.remove_harmonics(all_cands)
   
    # Write candidates to STDOUT
    if len(all_cands):
        all_cands.sort(sifting.cmp_snr)
        sifting.write_candlist(all_cands, candfilenm=filenm)

    Ncands = len(all_cands)
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
            #namesplit = namecand.split("_"+str(params.zmax)+":")
            namesplit = namecand.rsplit("_", 1)
            if len(namesplit) != 2:
                continue
            else:
                bname = namesplit[0]
                cnum  = namesplit[1].split(":")[1]
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
    # Close file, we're done
    f.close()


def run_one_fft(datfile, fft_out, fft_err):
    cmd = "realfft %s" %datfile
    fout = file(fft_out, 'a+')
    ferr = file(fft_err, 'a+')
    try_cmd(cmd, stdout=fout, stderr=ferr)


def run_realfft(workdir, basename):
    """
    Take FFT of all the dedispersed .dat files
    """
    t_fft_start = time.time()
    # Grab the .dat files
    datfiles = glob('*.dat')
    datfiles.sort()
    
    print("Will take FFTs of %d *.dat files\n" %len(datfiles))

    # Take FFTs using multiple cores
    ncores = params.accel_cores
    pool_accel = mp.Pool(processes=ncores)
    fft_out = 'fft.out'
    fft_err = 'fft.err'
    for datfile in datfiles:
        current_nm = datfile.split('.dat')[0] + '.fft'
        if os.path.exists(current_nm):
            print "File " + current_nm + "already made, skipping"
        else:
            arg_tup = (datfile, fft_out, fft_err)
            pool_accel.apply_async(func=run_one_fft, args=arg_tup)
    pool_accel.close()
    pool_accel.join()

    t_fft_end = time.time()
    dt = t_fft_end - t_fft_start
    return dt


def run_zap(workdir, basename):
    """
    Zap all *fft files with zaplist

    Needs zapfile and baryv (set in params)
    """
    t_zap_start = time.time()
    # Set params
    baryv = params.baryv
    zap_path = params.zapfile

    print("HERE!")
    sys.stdout.flush()

    # Grab the .fft files
    fftfiles = glob('*.fft')
    fftfiles.sort()
    print("Will zap the %d *.fft files\n" %len(fftfiles))
    sys.stdout.flush()

    # Check if zapfile exists
    if os.path.exists(zap_path):
        print("Using zapfile:  %s" %(zap_path))
        # Copy it over
        zapfile = zap_path.split('/')[-1]
        shutil.copyfile(zap_path, zapfile)

        for fftfile in fftfiles:
            cmd = "zapbirds -baryv %.12f -zap -zapfile %s %s" \
                    %(baryv, zapfile, fftfile)
            try_cmd(cmd)

    else:
        print("ZAP FILE NOT FOUND:  %s" %(zap_path))
        sys.stdout.flush()

    t_zap_end = time.time()
    dt = t_zap_end - t_zap_start
    return dt


def run_one_accel(datfile, zmax, nharm, flo, fhi, sigma_min, zap_str, suffix):
    fname = datfile.split('/')[-1].rstrip(suffix)
    accel_out = file(fname+'_accel.out', 'w')
    accel_err = file(fname+'_accel.err', 'w')
    cmd = 'accelsearch -zmax %d -numharm %d -flo %.6f -fhi %.6f -sigma %.2f %s %s' \
          %(zmax, nharm, flo, fhi, sigma_min, zap_str, datfile)
    try_cmd(cmd, stdout=accel_out, stderr=accel_err)
    return


def run_accelsearch(work_dir, basename):
    """
    This function will run PRESTO's accelsearch on all the *.dat files
    in the current directory. The candidates ``sifted'' to remove
    candidates with obvious problems, harmonics, and duplicates. The
    final sifted list is written to a *.sifted.cands file.

    This function will also run prepfold on all the ACCEL cands

    The maximum accel bin, zmax, can be set in the params.py file.
    """
    t_accel_start = time.time()

    # Set data file suffix (.dat or .fft)
    if params.use_fft:
        suffix = '.fft'
    else:
        suffix = '.dat'
    datfiles = glob('*%s' %suffix)
    datfiles.sort()
    
    # Create files for the output and errors from accelsearch
    accel_out = file('accelsearch.out', 'w')
    accel_err = file('accelsearch.err', 'w')

    # Read in needed parameters from params file
    zmax  = params.zmax
    nharm = params.numharm
    flo   = params.freq_lo
    fhi   = params.freq_hi
    zap_str = params.zap_str
    ncores = params.accel_cores
    sigma_min = params.sigma_threshold
    
    # Run lo_accel accelsearch on each .dat file if we have not already done so
    print "Running lo-accel (z=0) accelsearch on %d %s files...\n" %(len(datfiles), suffix)
    pool_accel = mp.Pool(processes=ncores)
    for datfile in datfiles:
        current_nm = datfile.split(suffix)[0]
        current_nm = current_nm + "_ACCEL_" + str(0) + ".cand"
        if os.path.exists(current_nm):
            print "File " + current_nm + "already made, skipping"
        else:
            arg_tup = (datfile, 0, nharm, flo, fhi, sigma_min, zap_str, suffix)
            pool_accel.apply_async(func=run_one_accel, args=arg_tup)
    pool_accel.close()
    pool_accel.join()
    
    # Run hi_accel accelsearch on each .dat file if we have not already done so
    print "Running hi-accel (z=%d) accelsearch on %d %s files...\n" %(zmax, len(datfiles), suffix)
    pool_accel = mp.Pool(processes=ncores)
    for datfile in datfiles:
        current_nm = datfile.split(suffix)[0]
        current_nm = current_nm + "_ACCEL_" + str(params.zmax) + ".cand"
        if os.path.exists(current_nm):
            print "File " + current_nm + "already made, skipping"
        else:
            arg_tup = (datfile, zmax, nharm, flo, fhi, sigma_min, zap_str, suffix)
            pool_accel.apply_async(func=run_one_accel, args=arg_tup)
    pool_accel.close()
    pool_accel.join()
    
    # Make candsfile (if it doesn't exist), run the accel sifting function
    # on the cands, which sifts and outputs to file
    candsfilenm = basename + '.ACCEL_' + str(zmax) + '.sifted.cands'
    print "Sifting through the candidates...\n"
    Ncands = accelsift(candsfilenm)

    # Run prepfold on the ACCEL files
    if Ncands:
        print "Running prepfold on the ACCEL cands...\n"
        run_prepfold(candsfilenm, accel_out, accel_err)

    # Make an output directory for the PRESTO cands
    outcand_dir = work_dir + "/cands_presto/"
    if not os.path.exists(outcand_dir):
        os.makedirs(outcand_dir)
    
    # Move the candidates over
    print "Moving candidates over to %s" %outcand_dir
    cmd = "mv %s/*ACCEL* %s" %(work_dir, outcand_dir)
    try_cmd(cmd)
    
    t_accel_end = time.time()
    dt = t_accel_end - t_accel_start
    return dt


def run_singlepulse_search(work_dir):
    sp_exe = params.singlepulse
    w_max = params.max_width 
    dfac  = params.dtrend 
    flags = params.sp_otherflags
    
    print("Looking for single pulses...\n")
    t_sp_start = time.time()
    #cmd = sp_exe+' -m '+str(w_max)+' '+work_dir+'/*.dat'
    cmd = "%s -m %.6f -d %d %s %s/*.dat" %(sp_exe, w_max, dfac, flags, work_dir)
    try_cmd(cmd)

    sp_dir = work_dir+'/single_pulse/'
    if not os.path.exists(sp_dir):
        os.makedirs(sp_dir)
    cmd = 'mv '+work_dir+'/*singlepulse* '+work_dir+'/single_pulse/'
    try_cmd(cmd)
    t_sp_end = time.time()
    dt = t_sp_end - t_sp_start
    return dt


def organize_results(work_dir):
    """
    Put *out + *err files into a directory

    Put *dat + *fft + *inf files into a directory
    """
    out_dir = "%s/output_files" %(work_dir)

    # Collect all errs and outs
    ostr_list = ["out", "err"]
    for ostr in ostr_list:
        ofiles = glob("%s/*%s" %(work_dir, ostr))
        if len(ofiles):
            # If files exist and out_dir doesnt, make it
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            for ofile in ofiles:
                shutil.move(ofile, out_dir)
        else: pass

    # Collect all infs, dats, and ffts
    dmstr_list = ["inf", "dat", "fft"]
    for dmstr in dmstr_list:
        dmfiles = glob("%s/*%s" %(work_dir, dmstr))
        dm_dir = "%s/dm_%s" %(work_dir, dmstr)
        if len(dmfiles):
            # If files exist and dm_dir doesnt, make it
            if not os.path.exists(dm_dir):
                os.makedirs(dm_dir)
            for dmfile in dmfiles:
                shutil.move(dmfile, dm_dir)
        else: pass

    return
    

def search_beam(fitsname, fits_dir, work_dir):
    tt = Timer()
    t_start = time.time()
    
    print("Results Directory: %s\n" %work_dir)
    print("FITS Directory: %s\n" %fits_dir)
    print("File Prefix: %s\n" %fitsname)
    
    # Check to see if results directory exists. If not, create it.
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    
    # Check dependencies for planned processing steps.
    check_dependencies(work_dir, fits_dir, fitsname)
    
    # Copy over param file if so desired.
    if params.do_param_cp:
        cp_cmd = 'cp params.py %s/params.txt' %work_dir
        try_cmd(cp_cmd)
    
    # If we haven't done so already, go to results directory
    os.chdir(work_dir)
    
    # Need the following if we are doing rfifind or prepsubband
    if params.do_rfifind or params.do_prepsub:
        fitslist = glob('%s/%s*.fits' %(fits_dir, fitsname))
        fitslist.sort()
        #fitsfiles = ' '.join(fitsfiles_arr)
    
    # Run rfifind on the .fits files and put rfifind products in a 
    # folder called rfi_products
    if params.do_rfifind:
        tt.rfifind = run_rfifind(fitslist, fitsname, work_dir)
    rfi_dir = work_dir+'/rfi_products'
    if os.path.exists(rfi_dir):
        maskname = glob(rfi_dir+'/*.mask')[0]
    
    # Run prepsubband on the .fits files
    if params.do_prepsub:
        if params.dmcalls > 1:
            tt.prepsubband = multi_call_prepsubband(fitsname, maskname, fitslist)
        else:
            tt.prepsubband = run_prepsubband(fitsname, maskname, fitslist)
    
    # Run realfft on the dedispersed time series, if selected
    if params.do_fft:
        tt.realfft = run_realfft(work_dir, fitsname)

    # Zap FFT files
    if params.do_zap:
        tt.zap = run_zap(work_dir, fitsname)
    
    # Search dedispersed time series for pulsar candidates
    # This includes an acceleration search
    if params.do_candsearch:
        tt.accelsearch = run_accelsearch(work_dir, fitsname)
    
    # Do a single pulse search and move the results into
    # a singlepulse directory
    if params.do_presto_sp:
        tt.presto_sp = run_singlepulse_search(work_dir)

    # Organize stray files into folders
    organize_results(work_dir)

    # Finish up time profiling and print summary to screen
    t_finish = time.time()
    tt.total = t_finish - t_start
    tt.print_summary()
    tt.write_summary("%s.log" %fitsname)

    return


def setup(basename, local_fits, local_results, host_fits, host_results):
    """
    Copy over files from LOCAL to HOST as necessary
    """
    resume = params.resume
    copy_fits = params.copy_fits

    sub_dirs = ["rfi_products", "dm_dat", "dm_inf", "dm_fft"]

    # Make HOST fits dir
    if not os.path.exists(host_fits):
        os.makedirs(host_fits)

    # Make HOST results dir 
    if not os.path.exists(host_results):
        os.makedirs(host_results)

    # Copy FITS if desired
    if copy_fits:
        print("COPYING OVER FITS FILES")
        fitslist = glob("%s/%s*fits" %(local_fits, basename))
        if len(fitslist):
            for fitsfile in fitslist:
                fname = fitsfile.split('/')[-1]
                dest_file = "%s/%s" %(host_fits, fname)
                shutil.copyfile(fitsfile, dest_file)
        else:
            print("No FITS files found in %s" %(local_fits))
    
    # Copy last results if desired
    if resume:
        print("RESUMING PREVIOUS SEARCH")
        for sub_dir in sub_dirs:
            if os.path.exists("%s/%s" %(local_results, sub_dir)):
                print("COPYING %s" %sub_dir)
                local_dir = "%s/%s" %(local_results, sub_dir)
                host_dir  = "%s/%s" %(host_results, sub_dir)
                shutil.copytree(local_dir, host_dir)
            else: pass

        dat_files = glob("%s/dm_dat/*" %(host_results))
        if len(dat_files):
            for dat_file in dat_files:
                shutil.move(dat_file, host_results)
            shutil.rmtree("%s/dm_dat" %(host_results))
        
        inf_files = glob("%s/dm_inf/*" %(host_results))
        if len(inf_files):
            for inf_file in inf_files:
                shutil.move(inf_file, host_results)
            shutil.rmtree("%s/dm_inf" %(host_results))

        fft_files = glob("%s/dm_fft/*" %(host_results))
        if len(fft_files):
            for fft_file in fft_files:
                shutil.move(fft_file, host_results)
            shutil.rmtree("%s/dm_fft" %(host_results))

    return


def get_results(local_results, host_results):
    """
    Copy over files from HOST to LOCAL
    """
    # Check for and copy directories
    sub_dirs = ["cands_presto", "dm_dat", "dm_inf", "dm_fft", 
                "output_files", "rfi_products", "single_pulse"]
    for sub_dir in sub_dirs:
        local_path = "%s/%s" %(local_results, sub_dir)
        host_path  = "%s/%s" %(host_results, sub_dir)
        if not os.path.exists(local_path):
            if os.path.exists(host_path):
                shutil.copytree(host_path, local_path)
            else: pass
        else: pass

    # Check for and copy log file
    log_files = glob("%s/*log" %(host_results))
    if len(log_files):
        for log_file in log_files:
            fname = log_file.split('/')[-1]
            shutil.copyfile(log_file, "%s/%s" %(local_results, fname))
    else: pass
    
    return


def cleanup(host_fits, host_results):
    """
    Remove remaining files from HOST
    """
    if os.path.exists(host_fits):
        print("REMOVING HOST FITS FILES in %s" %host_fits)
        shutil.rmtree(host_fits)
    else: pass
    
    if os.path.exists(host_results):
        print("REMOVING HOST RESULTS in %s" %host_results)
        shutil.rmtree(host_results)
    else: pass

    return 


def print_dirs(h_top, h_fits, h_results, l_fits, l_results):
    print("========== HOST ===========")
    print(" TOP     = %s" %h_top)
    print(" FITS    = %s" %h_fits)
    print(" RESULTS = %s" %h_results)
    print("===========================")
    print("\n\n")
    print("========== LOCAL ===========")
    print(" FITS    = %s" %l_fits)
    print(" RESULTS = %s" %l_results)
    print("===========================")
    print("\n\n")
    return 



####################
##     MAIN       ##
####################


if __name__ == "__main__":
    # Relevant directories on compute node
    top_HOST    = sys.argv[1]
    fits_HOST   = "%s/psrfits" %(top_HOST)
    results_HOST = "%s/search" %(top_HOST)

    # Relevant directories on local space
    fits_LOCAL    = params.fits_dir
    results_LOCAL = params.results_dir
   
    # basename for fits files
    basename = params.basename

    # print path info
    print_dirs(top_HOST, fits_HOST, results_HOST, 
               fits_LOCAL, results_LOCAL)

    # TRY SETTING UP AND SEARCHING
    try:
        # SET-UP on HOST 
        setup(basename, fits_LOCAL, results_LOCAL, 
              fits_HOST, results_HOST)

        print(glob("%s/*" %fits_HOST))
        print(glob("%s/*" %results_HOST))

        # Actually do the search
        search_beam(basename, fits_HOST, results_HOST)

        print(glob("%s/*" %results_HOST))
    
        # Copy back results
        get_results(results_LOCAL, results_HOST)

    except:
        print("Something failed!!!!")

    finally:
        # Delete everything from compute node
        cleanup(fits_HOST, results_HOST)

