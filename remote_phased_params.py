# Important paths
fits_dir = "/hercules/results/rwharton/phased_gc/day2/fix_fits"
#results_dir = "/hercules/results/rwharton/phased_gc/day2/search_fix"
results_dir = "/hercules/results/rwharton/phased_gc/day2/search_red_nomag"
basename = "sgra_day2"

# Resume previous processing?
resume    = 1
copy_fits = 0

# Processing steps to do
do_rfifind    = 0      # Run PRESTO rfifind and generate a mask
do_prepsub    = 0      # Run PRESTO prepsubband dedispersion
do_fft        = 0      # FFT *dat files before accelsearch
do_zap        = 0      # Zap the *fft files
do_candsearch = 1      # Run PRESTO accelsearch on the data
do_presto_sp  = 0      # Run PRESTO singlepulse.py
do_param_cp   = 0      # copy parameter file to output directory 

# data format type
# Use PRESTO naming conventions (ie, SIGPROC filterbank => 'filterbank'
# and PSRFITS => 'psrfits')
dat_type = 'psrfits'

# Locations of all the required scripts
singlepulse = 'single_pulse_search.py'

# rfifind params
use_mask  = 1     # use the rfi mask?
rfi_time  = 5.0   # sec, time over which stats are calc'd in rfifind (d: 30)
time_sig  = 15    # The +/-sigma cutoff to reject time-domain chunks (d: 10)
freq_sig  = 6     # The +/-sigma cutoff to reject freq-domain chunks (d: 4)
chan_frac = 0.3   # The fraction of bad channels that will mask a full interval (d: 0.7)
int_frac  = 0.3   # The fraction of bad intervals that will mask a full channel (d: 0.3)
rfi_otherflags = '-ncpus 8 -zapchan 512:704 -noweights ' #

# Dedispersion parameters
# NOTE: Will do multi-pass sub-band dedispersion
#       if and only if dmcalls > 1
dmlow = 20 
ddm   = 20
dmcalls    = 1
dmspercall = 250    # 250
nsub       = 608    # 608
dsubDM     = 0.0
downsample = 1
prep_otherflags = '-ncpus 8 -noweights ' # other flags for prepsubband

# Zapping params
baryv = 0.0 # 9.841459e-05  # Needs to change for each obs
zapfile = "" #"/hercules/results/rwharton/phased_gc/day1/zaplist/phased_gc_mag.zaplist"

# accelsearch parameters
use_fft = 1        # Search FFTs (1) instead of dats (0)
zmax = 300         # Max acceleration in fourier bins
numharm = 16       # Number of harmonics to sum
freq_lo = 0.1      # Hz, lowest freq to consider real
freq_hi = 2500.0   # Hz, highest freq to consider real
zap_str = '' # Put zap stuff here, if desired
accel_cores = 8   # number of processing cores

# Single pulse
max_width = 1.0   # Max pulse width (seconds)
dtrend    = 32    # Detrend factor 1-32 powers of two 
sp_otherflags = "-f -p " # Other flags


# Sifting parameters (copied mostly from PRESTO's ACCEL_sift.py):
#--------------------------------------------------------------
min_num_DMs = 1      # Min number of DM bins a candidate must show up in to be "good"
low_DM_cutoff = 100.0  # Lowest DM to be considered a "real" pulsar

# The following show up in the sifting.py PRESTO module:

# Ignore candidates with a sigma (from incoherent power summation) less than this
sigma_threshold = 2.0   # was 3.0    
# Ignore candidates with a coherent power less than this
c_pow_threshold = 100.0
# If the birds file works well, the following shouldn't
# be needed at all...  If they are, add tuples with the bad
# values and their errors.
#        (ms, err)
known_birds_p = []
#        (Hz, err)
known_birds_f = []
# The following are all defined in the sifting module.
# But if we want to override them, uncomment and do it here.
# You shouldn't need to adjust them for most searches, though. 

# How close a candidate has to be to another candidate to
# consider it the same candidate (in Fourier bins)
r_err = 1.1
# Shortest period candidates to consider (s)
short_period = 1./freq_hi
# Longest period candidates to consider (s)
long_period = 1./freq_lo
# Ignore any candidates where at least one harmonic does exceed this power
harm_pow_cutoff = 2.0
#--------------------------------------------------------------
