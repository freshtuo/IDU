#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
# auto_bcl2fq_script.py
# automatically generate script for running demux given a samplesheet file
# split samples based on index

Created on Tue Sep  3 12:55:03 2019

@author: taz2008
@version: 1.01
"""

import sys
import argparse
import os
from re import search
from string import atoi

# functions
def diff(idx1, idx2):
        # calculate Hamming distance between two indexes
        return len([i for i in xrange(len(idx1)) if idx1[i] != idx2[i]])

def minDiff(tidxList):
        # calculate the minimum distance between each pair of indexes with in a index group
        # only one index? then return the length of that index
        if len(tidxList) == 1:
                return len(tidxList[0])
        # more than one index, check distance between each pair of indexes
        tdis = []
        for ti in xrange(len(tidxList)):
                for tj in xrange(ti+1, len(tidxList)):
                        tdis.append(diff(tidxList[ti], tidxList[tj]))
        return min(tdis)

def splitSampleSheet(tfile):
        # split samplesheet into two parts: header info and sample info
        with open(tfile, 'rU') as fin:
                tdata = fin.readlines()
                tk = 0
                while tk != len(tdata):
                        if search("^\[Data\]", tdata[tk]):
                                break
                        tk += 1
                if tk == len(tdata):# end of file
                        print "Failed to locate sample info section within samplesheet file: %s"%tfile
                        sys.exit(99)
                return (tdata[:tk+1], tdata[tk+1:])# (header, sample)

def scanSampleInfo(tdata, tskip):
        # scan sample info in the samplesheet
        # skip unwanted lanes
        # split columns by comma
        tinfo = [tx.strip().split(",") for tx in tdata]
        # header dictionary
        thdic = {}
        for tk,tx in enumerate(tinfo[0]):
                thdic[tx] = tk
        # remove samples from unwanted lanes
        if "Lane" in thdic and tskip:
                return (thdic, [tx for tx in tinfo[1:] if tx[thdic["Lane"]] not in tskip])
        else:
                return (thdic, tinfo[1:])# header dictionary, split sample info

def splitByLaneIndex(thdic, tinfo):
        # split sample info entries by lane & index
        # get group type per sample: Lane & index size & index2 size
        tlabels = []
        for tk,tx in enumerate(tinfo):
                if "Lane" in thdic:# lane column exists
                        if "index2" in thdic:# dual index column exists
                                tlabels.append(":".join([tx[thdic["Lane"]], "%d"%(len(tx[thdic["index"]])), "%d"%(len(tx[thdic["index2"]]))]))
                        else:# only a single index available
                                tlabels.append(":".join([tx[thdic["Lane"]], "%d"%(len(tx[thdic["index"]])), "0"]))
                else:# no lane column found
                        if "index2" in thdic:# dual index column exists
                                tlabels.append(":".join(["1", "%d"%(len(tx[thdic["index"]])), "%d"%(len(tx[thdic["index2"]]))]))
                        else:# only a single index available
                                tlabels.append(":".join(["1", "%d"%(len(tx[thdic["index"]])), "0"]))
        # summarize possible groups
        tgrpdic = {}
        for tk,tx in enumerate(tlabels):
                if tx not in tgrpdic:
                        tgrpdic[tx] = [tk]
                else:
                        tgrpdic[tx].append(tk)
        return tgrpdic

def allowMismatch(thdic, tinfo, tgrpdic):
        # check per group: allow mismatch? True of False
        tmisdic = {}
        for tgrp in tgrpdic:# each group
                if "index2" in thdic:# dual index column exists
                        if tgrp.split(":")[1] != "0" and tgrp.split(":")[2] != "0":# dual index!!!
                                tidxes = [tinfo[tk][thdic["index"]]+tinfo[tk][thdic["index2"]] for tk in tgrpdic[tgrp]]
                                if minDiff(tidxes) < 3:# no mismatch allowed
                                        tmisdic[tgrp] = False
                                else:
                                        tmisdic[tgrp] = True
                        elif tgrp.split(":")[1] != "0":# single index!!!
                                tidxes = [tinfo[tk][thdic["index"]] for tk in tgrpdic[tgrp]]
                                if minDiff(tidxes) < 3:# no mismatch allowed
                                        tmisdic[tgrp] = False
                                else:
                                        tmisdic[tgrp] = True
                        else:# lane without any index
                                tmisdic[tgrp] = True
                else:# single index
                        if tgrp.split(":")[1] != "0":# single index confirmed !!!
                                tidxes = [tinfo[tk][thdic["index"]] for tk in tgrpdic[tgrp]]
                                if minDiff(tidxes) < 3:# no mismatch allowed
                                        tmisdic[tgrp] = False
                                else:
                                        tmisdic[tgrp] = True
                        else:# lane without any index
                                tmisdic[tgrp] = True
        # any groups from the same lane, do NOT allow mismatch then!!!
        tlanedic = {}
        for tgrp in tgrpdic:# each group
                tlane = tgrp.split(":")[0]
                if tlane not in tlanedic:
                        tlanedic[tlane] = [tgrp]
                else:
                        tlanedic[tlane].append(tgrp)
        for tlane in tlanedic:
                if len(tlanedic[tlane]) > 1:
                        for tgrp in tlanedic[tlane]:
                                tmisdic[tgrp] = False
        return tmisdic

def mergeGroups(tgrpdic, tmisdic):
        # merge groups by Indexes & Mismatch
        tmergedic = {}
        for tgrp in tgrpdic:# each group
                # ignore lanes, only takes indexes
                # mismatch: True --> 1; False --> 0
                tid = ":".join(tgrp.split(":")[1:] + ["%d"%(tmisdic[tgrp])])# index 1 size + index 2 size + allow mismatch
                if tid not in tmergedic:
                        tmergedic[tid] = [tgrp]
                else:
                        tmergedic[tid].append(tgrp)
        return tmergedic

def getVal(patname,targetStrList):
        # search for <patname>val</patname> in xml file
        for targetstr in targetStrList:
                tpat = search("<%s>(.*)</%s>"%(patname,patname),targetstr)
                if tpat:
                        return (True,atoi(tpat.groups()[0]))
        return (False,None)

def fetchSeqInfo(trunparfile, tr1id, ti1id, ti2id, tr2id):
        # extract sequencing length info
        with open(trunparfile, 'r') as fseq:
                rp = fseq.readlines()
                # get read 1 length
                mark,r1 = getVal(tr1id, rp)
                if not mark:
                        print "Failed to locate <Read1>?</Read1>"
                        sys.exit(98)
                # get index 1 length
                mark,i1 = getVal(ti1id, rp)
                if not mark:
                        print "Failed to locate <IndexRead1>?</IndexRead1>"
                        sys.exit(97)
                # get index 2 length
                mark,i2 = getVal(ti2id, rp)
                if not mark:
                        print "Failed to locate <IndexRead2>?</IndexRead2>"
                        sys.exit(96)
                # get read 2 length
                mark,r2 = getVal(tr2id, rp)
                if not mark:
                        print "Failed to locate <Read2>?</Read2>"
                        sys.exit(95)
                return (r1,i1,i2,r2)

def guessBaseMask(tmergedic, r1, i1, i2, r2):
        # guess base mask for each merged group
        tmaskdic = {}
        for tid in tmergedic:
                gi1 = atoi(tid.split(":")[0])
                gi2 = atoi(tid.split(":")[1])
                # index size longer than sequencing read length?
                if gi1 > i1 or gi2 > i2:
                        print "Error: index in the samplesheet is longer than that was actually sequenced!"
                        print "       please manually check!!!"
                        sys.exit(94)
                # guess base mask
                tmask = []
                # add read 1
                if r1 > 0:
                        tmask.append("Y*")
                # add index 1
                if i1 > 0:
                        if gi1 == 0:# no index 1
                                tmask.append("N*")
                        else:
                                tmask.append("I*" + "n" * (i1 - gi1))
                # add index 2
                if i2 > 0:
                        if gi2 == 0:# no index 2
                                tmask.append("N*")
                        else:
                                tmask.append("I*" + "n" * (i2 - gi2))
                # add read 2
                if r2 > 0:
                        tmask.append("Y*")
                tmaskdic[tid] = ",".join(tmask)
        return tmaskdic

def makeScript(theaders, tsamples, tinfo, tgrpdic, tmergedic, tmaskdic, toutfile, tssdir, tworkdir, tprogram, tnosplit):
        # generate sub-samplesheet files and shell script
        with open(toutfile, 'w') as fout:
                for tj,tid in enumerate(tmergedic):# each merged group
                        # sub-samplesheet file
                        with open(tssdir + "samplesheet_grp%d.csv"%(tj+1), 'w') as fss:
                                fss.write("".join(theaders) + tsamples[0])# header part
                                for tgrp in tmergedic[tid]:# each group
                                        for tk in tgrpdic[tgrp]:# each sample
                                                fss.write(",".join(tinfo[tk])+"\n")
                        # shell script
                        fout.write("nohup %s  --runfolder-dir %s/  --output-dir %s/Unaligned_%d/"%(tprogram, tworkdir, tworkdir, tj+1))
                        fout.write("  --sample-sheet %ssamplesheet_grp%d.csv"%(tssdir,tj+1))
                        fout.write("  --use-bases-mask %s"%(tmaskdic[tid]))
                        if tid.split(":")[2] == "0":# no allow mismatch
                                fout.write("  --barcode-mismatches 0")
                        if tnosplit:
                                fout.write("  --no-lane-splitting")
                        fout.write(" >run_bcl2fq.%d.log\n\n"%(tj+1))

ver = "1.01"

# main
if __name__ != "__main__":
        #print "Do nothing, import functions only."
        sys.exit(0)

# command line arguments
parser = argparse.ArgumentParser(description="Given a samplesheet, automatically generate shell script for demux.")
parser.add_argument("-v", "--version", action="version", version="Ver %s"%ver)
parser.add_argument("-r", "--runfolder", nargs="?", required=False, default=os.getcwd(), help="sequencing run folder", metavar="run_folder", dest="runfolder")
parser.add_argument("-i", "--samplesheet", nargs="?", required=True, help="samplesheet file", metavar="samplesheet_file", dest="infile")
parser.add_argument("-s", "--skiplanes", nargs="*", required=False, default=[], type=int, choices=[1,2,3,4,5,6,7,8], help="skip lanes", metavar="lanes_to_skip", dest="skip")
parser.add_argument("-b", "--bcl2fastq", nargs="?", required=False, default="/usr/local/bin/bcl2fastq", help="location of bcl2fastq program", metavar="bcl2fastq", dest="bcl2fastq")
parser.add_argument("-o", "--script", nargs="?", required=False, default="run_bcl2fq.sh", help="output script name", metavar="script_file", dest="outfilename")
parser.add_argument("-f", "--force", action="store_true", required=False, default=False, help="whether or not to overwrite output script if existing", dest="force")
parser.add_argument("-n", "--no-lane-splitting", action="store_true", required=False, default=False, help="whether or not to add --no-lane-splitting option to shell script", dest="nosplit")
parser.add_argument("-p", "--platform", nargs="?", required=False, default="hiseq", choices=["hiseq","nextseq","novaseq"], help="sequencing platform", dest="platform")

# parse arguments
args = parser.parse_args()
#print args
workdir = args.runfolder
infile = args.infile
program = args.bcl2fastq
skip = ["%d"%x for x in args.skip]
outfile = workdir + "/" + args.outfilename
overwrite = args.force
nosplit = args.nosplit
samplesheetdir = workdir + "/Data/Intensities/BaseCalls/"
runparfile = workdir + "/runParameters.xml"

# labels for acquiring read length for read1, read2, index1, index2 (default sequencing platform: HiSeq)
r1id = "Read1"
r2id = "Read2"
i1id = "IndexRead1"
i2id = "IndexRead2"

# modify for NextSeq/NovaSeq runs
if args.platform == "nextseq":
        # runParameters.xml file
        runparfile = workdir + "/RunParameters.xml"
        # labels for index reads
        i1id = "Index1Read"
        i2id = "Index2Read"
elif args.platform == "novaseq":
        # runParameters.xml file
        runparfile = workdir + "/RunParameters.xml"
        # labels for reads
        r1id = "Read1NumberOfCycles"
        r2id = "Read2NumberOfCycles"
        i1id = "IndexRead1NumberOfCycles"
        i2id = "IndexRead2NumberOfCycles"

# samplesheet file exists?
if not os.path.exists(infile):
        print "Cannot find samplesheet file: %s"%infile
        sys.exit(1)

# runParameters.xml exists?
if not os.path.exists(runparfile):
        print "Cannot find runParameters.xml file: %s"%runparfile
        sys.exit(2)

# output file already exists?
if os.path.exists(outfile) and not overwrite:
        print "Shell script %s already exists, quit without doing anything."%outfile
        print "Remove the target file or use -f/--force option to force overwriting the file."
        sys.exit(3)

# get sequencing length info
r1,i1,i2,r2 = fetchSeqInfo(runparfile, r1id, i1id, i2id, r2id)
#print r1,i1,i2,r2

# split header and sample info
headers, samples = splitSampleSheet(infile)

# scan samplesheet file
hdic, sinfo = scanSampleInfo(samples, skip)
#print len(sinfo)

# group samples by Lane & Indexes
grpdic = splitByLaneIndex(hdic, sinfo)
#print len(grpdic), grpdic.keys()

# check if allows mismatch for each group
misdic = allowMismatch(hdic, sinfo, grpdic)
#print misdic

# merge groups by Indexes & Mismatch
mergedic = mergeGroups(grpdic, misdic)
#print mergedic

# guess base mask
maskdic = guessBaseMask(mergedic, r1, i1, i2, r2)
#print maskdic

# prepare sub-samplesheet files and shell script for demux
makeScript(headers, samples, sinfo, grpdic, mergedic, maskdic, outfile, samplesheetdir, workdir, program, nosplit)

print "All complete!"
