#!/usr/bin/env python

import os, sys, re
import argparse

Description="Program to get common or alternate SNPs between two VCFs"
usage="""
python {script} --r1 sample1_R1.fastq sample2_R1.fastq sample3_R1.fastq --r2 sample1_R2.fastq sample2_R2.fastq sample3_R2.fastq -d projectdir --Ssample Sbulksample --Rsample Resistantsample --Rreference resistantparent.fasta --Sreference Susceptibleparent.ffasta --s1 sample1_R1.fastq sample2_R1.fastq sample3_R1.fast --s2 sample1_R2.fastq sample2_R2.fastq sample3_R2.fastq --Sreference susceptibleparent.fasta --b1 sample1_R1.fastq sample2_R1.fastq sample3_R1.fastq --b2 sample1_R2.fastq sample2_R2.fastq sample3_R2.fastq
""".format(script=sys.argv[0])

parser=argparse.ArgumentParser(description=Description, version="1.0", epilog=usage,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-d", "--projectdir", action="store", dest="projectdir", help="Project Directory name. This will be created if does not exist.")
parser.add_argument("--r1", nargs="+", action="store", dest="Rread1", help="Forward Reads for Resistant sample. Multiple forward reads can be provided(comma or space separated)")
parser.add_argument("--r2", nargs="+", action="store", dest="Rread2", help="Reverse Reads for Resistant sample. Multiple reverse reads can be provided(comma or space separated). Reverse reads files must be in the same order as forward reads.")
parser.add_argument("--s1", nargs="+", action="store", dest="Sread1", help="Forward Reads for Susceptible sample. Multiple forward reads can be provided(comma or space separated)")
parser.add_argument("--s2", nargs="+", action="store", dest="Sread2", help="Reverse Reads for Susceptible sample. Multiple reverse reads can be provided(comma or space separated). Reverse reads files must be in the same order as forward reads.")
parser.add_argument("--b1", nargs="+", action="store", dest="sbulkread1", help="Susceptible bulk forward reads")
parser.add_argument("--b2", nargs="+", action="store", dest="sbulkread2", help="Susceptible bulk reverse reads")

parser.add_argument("--sampleid", action="store", dest="sampleid", help="A unique name for the sample")
parser.add_argument("--Rsample", action="store", dest="resistantsample", help="Resistant sample name. Use this name if resistant reads are provided in --r1/--r2")
parser.add_argument("--Ssample", action="store", dest="Susceptiblename", help="Susceiptible samplename. Use this name if susceiptible reads are provided in --s1/--s2")
parser.add_argument("--sbulksample", action="store", dest="sbulksample", help="Susceptible Bulk population sample name. Use this name if susceptible reads are povided in --b1/--b2")
parser.add_argument("--Rreference", action="store", dest="Rreference", help="Resistant parent reference sequence (fasta format)")
parser.add_argument("--Sreference", action="store", dest="Sreference", help="Susceptible parent reference sequence (fasta format)")

options=parser.parse_args()

#print options


####### Test whether the tools are installed in the system or not. If not, exit the pipeline and print the message ##########

tools_not_found = []
from subprocess import PIPE, Popen
def get_system_cmdline(cmd):
    toolname=cmd[1]
    p=Popen(cmd, stdout=PIPE)
    out,err = p.communicate()
    if "no " + toolname + " in" in out or "not found" in out or out == "":
        tools_not_found.append(toolname)



get_system_cmdline(["which", "rake"])
get_system_cmdline(["which", "fastqc"])
get_system_cmdline(["which", "trimmomatic"])
get_system_cmdline(["which", "bowtie2-build"])
get_system_cmdline(["which", "bowtie2"])
get_system_cmdline(["which", "samtools"])


if len(tools_not_found) > 0:
    print ",".join(tools_not_found) + " are not found in your system. Please install these tools before running the pipeline."
    if "rake" in tools_not_found:
        print "To get rake, install the ruby programming language. rake comes with Ruby."
    print "\nThe presence or abasence of tools are checked using the command : which. e.g. which samtools\nPlease test the tools yourself using \"which\" command in the command line."
    exit(0)
else:
    print "OK. It looks good. All prerequisite tools are installed in this system. Let's start the pipeline."

##############################################################################################################################


########## Create projectdir #####

if not os.path.exists(options.projectdir):
    os.makedirs(options.projectdir)

os.chdir(options.projectdir)
#### run rake command from here. One rake command for each pair of reads ####################
if options.Rread1 and options.Rread2:
    for R1,R2 in zip(options.Rread1, options.Rread2):
        print "rake -f Rakefile projectdir=" + options.projectdir + " R1=" + R1 + " R2=" + R2 + " sampleid=ResistantSample Rreference=" + options.Rreference

if options.Sread1 and options.Sread2:
    for R1,R2 in zip(options.Sread1, options.Sread2):
        print "rake -f Rakefile projectdir=" + options.projectdir + " R1=" + R1 + " R2=" + R2 + " sampleid=SusceptibleSample Sreference=" + options.Sreference

if options.sbulkread1 and options.sbulkread2:
    for R1, R2 in zip(options.sbulkread1, options.sbulkread2):
        print "rake -f Rakefile projectdir=" + options.projectdir + " R1=" + R1 + " R2=" + R2 + " sampleid=SbulkSusSample Sreference=" + options.Sreference


print "rake -f Rakefile Rreference=" + options.Rreference + " Sreference=" + options.Sreference + " sbulksample=" + options.sbulksample + "filtersnps:run get_snps_subseq:run"

##############################################################################################

exit(0)
