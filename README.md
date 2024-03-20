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




