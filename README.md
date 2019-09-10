# IDU
Illumina Demux Utilities

Scripts for automating the demux process for Illumina sequencing platform:

- auto_bcl2fq_script.py   
  Given a samplesheet, the script automatically detects different types of barcodes (libraries) in it, and generates sub-samplesheet files as well as a shell script for demux.
```
	usage: auto_bcl2fq_script.py [-h] [-v] [-r [RUNFOLDER]] -i [samplesheet_file]
        	                     [-s [lanes_to_skip [lanes_to_skip ...]]]
	                             [-b [bcl2fastq]] [-o [script_file]] [-f] [-n]

	Given a samplesheet, automatically generate shell script for demux.

	optional arguments:
	  -h, --help            show this help message and exit
	  -v, --version         show program's version number and exit
	  -r [RUNFOLDER], --runfolder [RUNFOLDER]
	                        sequencing run folder
	  -i [samplesheet_file], --samplesheet [samplesheet_file]
	                        samplesheet file
	  -s [lanes_to_skip [lanes_to_skip ...]], --skiplanes [lanes_to_skip [lanes_to_skip ...]]
	                        skip lanes
	  -b [bcl2fastq], --bcl2fastq [bcl2fastq]
	                        location of bcl2fastq program
	  -o [script_file], --script [script_file]
	                        output script name
	  -f, --force           whether or not to overwrite output script if existing
	  -n, --no-lane-splitting
	                        whether or not to add --no-lane-splitting option to
	                        shell script
	
	Examples:
	auto_bcl2fq_script.py -i /gc4/HiSeq4000/flowcellA/190828_ST-K00128_0425_AHF2NMBBXY/Data/Intensities/BaseCalls/Hiseq4000FCA_08282019.csv  -s 2  -r /gc4/HiSeq4000/flowcellA/190828_ST-K00128_0425_AHF2NMBBXY
```

