# presto_search_pipe

A very simple PRESTO search pipeline that takes in 
filterbank files 

## Usage

Running the pipeline is pretty straightforward.  You give 
an input directory containing the fil or FITS files, the 
outdir where you want the search results to go, and the 
basename for the input data files:

     usage: presto_search.py [-h] indir outdir basename
     
     PRESTO Search Pipeline
     
     positional arguments:
       indir       Directory containing *.fil or *fits files
       outdir      Output directory for results
       basename    Basename for input data files
     
     optional arguments:
       -h, --help  show this help message and exit

For example, if I have the data file `gc-zap-xlcp_final.fil` in my 
current working directory and want to have the output go to a directory 
called `results`, I would just do:

    python -u /path/to/src/presto_search.py . results gc-zap-xlcp_final 

The results directory will be created if it doesn't already exist.

## Parameter File

The parameter file `search_params.py` is how we set the parameters 
of the search.  We can decide which steps to do:

    # Processing steps to do
    do_rfifind    = 1      # Run PRESTO rfifind and generate a mask
    do_prepsub    = 1      # Run PRESTO prepsubband dedispersion
    do_fft        = 1      # FFT *dat files before accelsearch
    do_zap        = 0      # Zap the *fft files
    do_candsearch = 1      # Run PRESTO accelsearch on the data
    do_presto_sp  = 1      # Run PRESTO singlepulse.py
    do_param_cp   = 0      # copy parameter file to output directory 

where there will be some check to make sure the requested steps 
can be done.



