BEAR: Better Emulation for Artificial Reads
------------------------------------------
Forked from https://github.com/sej917/BEAR.

BEAR is intended to be an easy-to-use collection of scripts for generating simulated WGS metagenomic reads with read lengths, quality scores, error profiles, and species abundances derived from real user-supplied WGS data.

This fork of BEAR implements normalization of reference sequences by length and supports scaffolds as a reference. It also whipes out some minor bugs of the original script.

BEAR is implemented as a collection of Perl and Python scripts, and is known to work with Perl v5.14.2 and Python v2.7.3. BEAR has the following dependencies:

    -Perl Getopt::Long, Bio::SeqIO (part of BioPerl), List::Util modules
    -Python Bio and Numpy packages
    -Python sys, csv, StringIO, random, decimal, argparse modules
    -DRISEE, which can be downloaded at https://github.com/MG-RAST/DRISEE


