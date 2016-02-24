README.txt for pal_finder_v0.02.04

Copyright (c) 2012 Alex Poole
Univeristy of Colorado
All rights reserved

-----------
DESCRIPTION
-----------

This file serves as documentation for pal_finder_v0.02.04.pl The purpose
of this script is to find microsatellite repeat elements directly from raw 454
or Illumina paired-end sequencing reads, and design PCR primers to amplify
these repeat loci in an automated fashion. Exact matches to repeats or 2-, 3-,
4-, 5-, and/or 6-mers are located and primer3 is then used to generate primer
pairs to amplify regions containing microsatellite loci.  The minimun number
of repeat units for each n-mer size is specified by the user.

-------------
PREREQUISITES
-------------

The script has been tested on Mac OS X and linux systems. It will not work on
Windows (but may if running cygwin). Perl version 5.10 or above is recommended
although it may work on previous versions.

A working primer3 executable (version 2.0.0) is also required. This must be
downloaded from the primer3 project site and compiled on the system on which
you will be running pal_finder_v0.02.02.pl. Newer versions of primer3
will not work with this script so it is imperative to get version 2.0.0. The
correct version can be downloaded from here:

http://primer3.sourceforge.net/releases.php

Once primer3 has been compiled move the executables (primer3_core, oligotm,
ntdpal, long_seq_tm_test) to a directory into the system $PATH or modify the
$PATH variable to include the directory with the primer3 executables.

-------
RUNNING
-------

In addition to the perl script there are two more files that need to be in
place before running the program. The most important of which is the
configuration file. This file contains all the settings and input parameters
for the script. The file "config.txt" in the top level directory is an exaple
configuration file. The config.txt file also serves as the primary
documentation for all the input parameters and settings and should be read
thoroughly before running pal_finder_v0.02.03.pl.

The other file is a list of sequences that primer3 will use to prevent primers
being designed to repeat, low complexity, or any other sequences specified.
The file "simple.ref" is such a file and contains low complexity sequence. 
This file can be modified to include species-specific repeat sequences which
can be obtained from the repBase website:

http://www.girinst.org/repbase/update/index.html

The name and location of the file can also be altered by changing the value of
the PRIMER_MISPRIMING_LIBRARY variable in the configuration file. This file is
not necessary if you are not designing primers and are only interested in the
microsatellite content of your reads.

To run the script:
Once all the setting have be adjusted in the configuration file open a
terminal window and at the command prompt type:

perl pal_finder_v0.02.04.pl <config.txt>

where <config.txt> is the name of your configuration file.

----
TEST
----

To test if the script is running correctly on your system, the first thing you
should do is move into the pal_finder_v0.02.02 directory and type:

/scratch/lfs/chensc/Trans_ea_ena/20150409/transcript_analysis/scripts

If the script finishes and outputs "done!!" and there are two files in the
test/output directory named "test_microsat_summary.txt" and
"test/output/test_PAL_summary.txt".

This will test the scripts' functionality for finding designing primers for
di-nucleotide repeat microsatellites on Illumina Paired-end reads, which can
be found in the test/data directory. There is also a 454 data set in the same
directory to test on 454 reads.

