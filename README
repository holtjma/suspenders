Introduction
============

As of version 0.2.2, Suspenders merges multiple alignments of the same reads under different pretenses.  The 
classic case is to two different in silico genomes which typically represent the mother's and father's genomes.
However, this can be expanded to multiple alignments that take place with different alignment settings such as
variations in gap length, mismatches allowed, etc.

For the classic case, two files are taken as input:
	(1) a BAM file that contains the first alignment to the mother's in silico genome
	(2) a BAM file that contains the second alignment to the father's in silico genome
These two alignments must be pre-processed by Lapels before running Suspenders so they can be compared in the 
same reference coordinates system.  Additionally, both input files need to be sorted by read name (as opposed 
to coordinates).

For more information on Suspenders, refer to "Read Annotation Pipeline for High-Throughput Sequencing Data" by
James Holt and Shunping Huang et al.

For detailed usage, please type after installation:

	pysuspenders -h

System Requirements
===================

Suspenders and its modules have been tested under Python 2.7.

Several python modules are required to run the code.

*[pysam] - Tested with pysam 0.7.4

As a wrapper of Samtools, the pysam module facilitates the manipulation of SAM/BAM files in Python. Its latest 
package can be downloaded from:

	http://code.google.com/p/pysam/


*[argparse] - Tested with argparse 1.2.1

The argparse module is used to parse the command line arguments of the module. It has been maintained in Python 
Standard Library since Python 2.7.  Its latest package can be downloaded from:

	http://code.google.com/p/argparse/

*[matplotlib]

The matplotlib module helps generate some statistical output charts for the pileup data.  It will also allow future
modification to the code for generating more informative statistical images.  Its latest package can be downloaded
from:

	http://matplotlib.org

* others
Reads that have multiple alignments are required to have the 'HI' tag to specify the hit index.  Recent 
aligners (eg. bowtie >= 0.12.8 and tophat >= 1.4.0) will create this tag.  It is recommended to use a recent 
version for read alignment.

Lapels is a pre-processing requirement after the two alignments.  This will put both input files into the same 
coordinate system for comparison while merging.  Its latest package can be downloaded from:

	http://code.google.com/p/lapels/


Installation
============

It is recommended to use easy-install (http://packages.python.org/distribute/easy_install.html) for the 
installation.

	easy_install suspenders

Alternatively, users can download the tarball of source from

	http://code.google.com/p/suspenders/

and then type:

	easy_install suspenders-<version>.tar.gz

By default, the package will be installed under the directory of Python dist-packages, and the executable of 
pysuspenders can be found under '/usr/local/bin/'.

If you don't have permission to install it in the system-owned directory, you can install it in locally following 
the next steps:

(1) Create a local package directory for python:

	mkdir -p <local_dir>

(2) Add the absolute path of <local_dir> to the environment variable PYTHONPATH:

	export PYTHONPATH=$PYTHONPATH:<local_dir>

(3) Use easy_install to install the package in that directory:

	easy_install -d <local_dir> suspenders-<version>.tar.gz

For example, if you want to install the package under the home directory in
a Linux system, you can type:

	mkdir -p /home/$USER/.local/lib/python/dist-packages/
	export PYTHONPATH=$PYTHONPATH:/home/$USER/.local/lib/python/dist-packages/
	easy_install -d /home/$USER/.local/lib/python/dist-packages/ suspenders-<version>.tar.gz

After installation, pysuspenders will be located in '/home/$USER/.local/lib/python/dist-packages/'.


Merge Types
===========
The different types are based on a set of filters to pull out reads that have already been successfully merged.  
Each filter catches a set of reads and marks them with a 'ct' tag denoting which filter caught that particular 
read.  See the 'ct' tag in examples for specific on how the filter works.

Union: Keep all the read alignments from both files, but if the alignments are identical (position, cigar string, 
edit distance), only store one copy of the read.  For example, if read A aligns to positions 1 and 2 in the mother 
and positions 2 and 3 in the father, the result will be a single read at each of the positions 1, 2, and 3.  
Filter order: Unique->Kept-All

Quality: Keep the single best alignment based on the quality score from 
'http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#local-alignment-score-example'.  If two or more alignments 
have the same score, they will be passed on to the next filter.  
Filter order: Unique->Quality->[Random, Kept-All]

Pileup: Keep the single best alignment based on the pileup heights.  The pileup heights are calculated using only 
alignments that have already been filtered by Unique or Quality.  Note: storing pileup data requires more memory 
than the other two merge types.  
Filter order: Unique->Quality->Pileup->[Random, Kept-All]

Examples
========

Examples can be downloaded from

	http://code.google.com/p/suspenders/

To run on the example files, navigate to the examples folder and run:

	pysuspenders -t ./merged.bam ./mother.bam ./father.bam

The following are snippets of the content in the input and output BAM file in the examples (gathered using 
'samtools view merged.bam').

mother.bam:
...
UNC9-SN296_0254:7:1101:1482:32883#CGATGT 153 chr2 22809107 50 25M3406N75M * 0 0 <SEQ> <MAPQ> 
NH:i:1 MD:Z:30C69 OM:i:1 XO:i:0 XM:i:1 i0:i:0 s0:i:0 OC:Z:25M3405N75M XG:i:0 AS:i:-1 XS:A:- d0:i:0
...
UNC9-SN296_0254:7:1101:1483:94304#CGATGT 137 chr1 81336983 50 100M * 0 0 <SEQ> <MAPQ> 
NH:i:1 MD:Z:69A13C16 OM:i:2 XN:i:0 XO:i:0 XM:i:2 i0:i:0 s0:i:0 OC:Z:100M XG:i:0 AS:i:-7 YT:Z:UU d0:i:0
...

father.bam:
...
UNC9-SN296_0254:7:1101:1482:32883#CGATGT 153 chr2 22809107 50 25M3406N75M * 0 0 <SEQ> <MAPQ> 
NH:i:1 MD:Z:30C69 OM:i:1 XO:i:0 XM:i:1 i0:i:0 s0:i:0 OC:Z:25M3407N75M XG:i:0 AS:i:-1 XS:A:- d0:i:0
...
UNC9-SN296_0254:7:1101:1483:94304#CGATGT 137 chr1 81336983 50 100M * 0 0 <SEQ> <MAPQ> 
NH:i:1 MD:Z:69A30 OM:i:1 XN:i:0 XO:i:0 XM:i:1 i0:i:0 s0:i:1 OC:Z:100M XG:i:0 AS:i:-1 YT:Z:UU d0:i:0
...

merged.bam:
...
UNC9-SN296_0254:7:1101:1482:32883#CGATGT 153 chr2 22809107 50 25M3406N75M * 0 0 <SEQ> <MAPQ> 
NH:i:1 MD:Z:30C69 OM:i:1 XO:i:0 XM:i:1 i0:i:0 s0:i:0 OC:Z:25M3405N75M XG:i:0 AS:i:-1 XS:A:- d0:i:0 
YA:A:3 ms:i:0 mi:i:0 md:i:0 ps:i:0 pi:i:0 pd:i:0 pc:Z:25M3407N75M pm:i:1 po:A:3 ct:A:U
...
UNC9-SN296_0254:7:1101:1483:94304#CGATGT 137 chr1 81336983 50 100M * 0 0 <SEQ> <MAPQ> 
NH:i:1 MD:Z:69A30 OM:i:1 XN:i:0 XO:i:0 XM:i:1 i0:i:0 s0:i:1 OC:Z:100M XG:i:0 AS:i:-1 YT:Z:UU d0:i:0 
YA:A:3 ms:i:0 mi:i:0 md:i:0 ps:i:1 pi:i:0 pd:i:0 mc:Z:100M mm:i:2 po:A:2 ct:A:Q
...

In the output, reads are merged if they are determined to be 'equal' based on a variety of filters.

Additionally, new tags have been added for each read to identify how the reads were merged:

Major tags:
po : an integer flag set representing the Parent of Origin of this particular read.  If multiple inputs
	 have identical positions, cigar strings, and number of mismatches, then we consider it to be an 
	 identical alignment.  In this case, the bit for each input will be set.  For example, given three
	 inputs and an identical read in the second and third input, our flag would be 0b110 which would be
	 stored as integer value 6.
ct : the Choice Type for this read, eg the filter that determined which read to save
	'U' for unique filter: only one possible alignment available, so it was kept
	'Q' for quality score filter: chose the possible alignment with the highest score (calculated using the Bowtie schema)
	'P' for pileup height filter: chose the possible alignment with the highest pileup height
	'R' for random filter: chose a random possible alignment from a set of 'equal' alignments based on the previous filters
	'K' for kept-all filter: keep all possible alignments from a set of 'equal' alignments based on the previous filters