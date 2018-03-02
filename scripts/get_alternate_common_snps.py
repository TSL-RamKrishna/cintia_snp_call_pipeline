#!/usr/bin/env python

import os, sys, re
import argparse

from natsort import natsorted

Description="Program to get common or alternate SNPs between two VCFs"
usage="""
python {script} --vcf1 VCF1.vcf -2 VCF2.vcf --alternate #finds positions where alternate snps are different
python {script} --vcf1 VCF1.vcf -2 VCF2.vcf --common    #finds positions where alternate snps are common
python {script} --vcf1 VCF1.vcf -2 VCF2.vcf --12        #finds positions in VCF1 that are not present in VCF2
python {script} --vcf1 VCF1.vcf -2 VCF2.vcf --21        #finds positions in VCF2 that are not present in VCF1
python {script} --vcf1 VCF1.vcf -2 VCF2.vcf --common12
python {script} --vcf1 VCF1.vcf -2 VCF2.vcf --common21
python {script} --vcf1 VCF1.vcf -2 VCF2.vcf --alternate12
python {script} --vcf1 VCF1.vcf -2 VCF2.vcf --alternate21
""".format(script=sys.argv[0])

parser=argparse.ArgumentParser(description=Description, epilog=usage)
parser.add_argument("-1", "--vcf1", action="store", dest="vcf1", help="First VCF input file")
parser.add_argument("-2", "--vcf2", action="store", dest="vcf2", help="Second VCF input file")
parser.add_argument("--vcfout", action="store_true", dest="vcfout", help="Output in VCF format")
parser.add_argument("--common", action="store_true", dest="common", help="Gets positions with common SNPs")
parser.add_argument("--alternate", action="store_true", dest="alternate", help="Gets positions with alternate SNPs")
parser.add_argument("--12", action="store_true", dest="onetwo", help="Gets positions with SNPS in first VCF but not in second")
parser.add_argument("--21", action="store_true", dest="twoone", help="Gets positions with SNPS in second VCF but not in first")
parser.add_argument("--common12", action="store_true", dest="common12", help="Gets positions with SNPs that are common in both VCFs and those that are unique to VCF1")
parser.add_argument("--common21", action="store_true", dest="common21", help="Gets positons with SNPs that are common in both VCFs and those that are unique ot VCF2")
parser.add_argument("--alternate12", action="store_true", dest="alternate12", help="Gets positions with SNPs that are alternate in both VCFs and those that are unique to VCF1")
parser.add_argument("--alternate21", action="store_true", dest="alternate21", help="Gets positions with SNPs that are alternate in both VCFs and those that are unique to VCF2")
parser.add_argument("--refaltsame", action="store_true", dest="refaltsame", help="Output if reference and alternate bases are same")
parser.add_argument("--refaltdiff", action="store_true", dest="refaltdiff", help="Output if reference and alternate bases are different")
parser.add_argument("--homo", action="store_true", dest="homozygous", help="Consider only the homozygous SNPs")
parser.add_argument("--hetero", action="store_true", dest="heterozygous", help="Consider only the heterozygous SNPs")
parser.add_argument("--out", action="store", dest="output", help="Output file name")

options=parser.parse_args()

def sort_2_chromosomes(vcf1chr, vcf2chr):
    return natsorted([vcf1chr, vcf2chr])




class compare_2_vcfs():


    def __init__(self, vcf1, vcf2, output):

        self.vcf1=vcf1
        self.vcf2=vcf2
        self.readlinevcf1=""
        self.readlinevcf2=""
        self.outputfh=""
        if options.output:
            print "ouptut file specified"
            self.output = options.output
            self.outputfh = open(self.output, "w")
        else:
            self.output = None

    def define_mode(self, mode):
        self.mode=mode

    def open_input_files(self):
        self.vcf1fh=open(self.vcf1)
        self.vcf2fh=open(self.vcf2)

    def close_input_files(self):
        self.vcf1fh.close()
        self.vcf2fh.close()

    def write_to_output(self, vcfline):
        if self.outputfh:
            if options.vcfout:
                self.outputfh.write(vcfline + "\n")
            else:
                self.outputfh.write(altsnp + "\n")
        else:
            if options.vcfout:
                print vcfline
            else:
                print " ".join(vcfline.split()[:5])

    def close_output_file(self):
        if self.output:
            self.outputfh.close()

    def get_ref_alt_same_base(self, vcfline):

        array=vcfline.split()
        if (len(array[3])==1 and array[4] == ".") or array[3] == array[4]:
            return True
        else:
            return False

    def get_ref_alt_diff_base(self, vcfline):

        array=vcfline.split()
        if array[3] != array[4]:
            return True
        else:
            return False

    def check_ref_alt_bases_and_write(self, vcfline):
        if options.refaltsame:
            if self.get_ref_alt_same_base(vcfline) == True:
                self.write_to_output(vcfline)
        elif options.refaltdiff:
            if self.get_ref_alt_diff_base(vcfline) == True:
                self.write_to_output(vcfline)
        else:
            self.write_to_output(vcfline)

    def read_vcf1_line(self):
        return self.vcf1fh.readline()

    def read_vcf2_line(self):
        return self.vcf2fh.readline()

    def add_header(self, header):
        if self.output:
            self.outputfh.write(header + "\n")
        else:
            print header

    def get_line_data(self, line):
        linearray=line.rstrip().split()
        if len(linearray) >= 5:
            chromosome=linearray[0]
            position=linearray[1]
            ref=linearray[3]
            alt=linearray[4]
            return chromosome, position, ref, alt
        else:
            return None, None, None, None

    def get_higher_dp_vcfline(self,vcf1line, vcf2line):
        dp1array=re.findall('DP\d*=\d{1,10}', vcf1line)
        dp1=sorted(map(lambda x: int(x.split("=")[1]), dp1array), reverse=True)[0]
        dp2array=re.findall('DP\d*=\d{1,10}', vcf2line)
        dp2=sorted(map(lambda x: int(x.split("=")[1]), dp2array), reverse=True)[0]

        if dp1 > dp2:
            return vcf1line
        else:
            return vcf2line

    def get_vcf1_linearray(self):

        while True:
            vcf1line=self.read_vcf1_line().rstrip()
            self.readlinevcf1=vcf1line

            if vcf1line.startswith("#"):
                if options.onetwo or options.common12 or options.alternate12:
                    if options.vcfout:
                        self.outputfh.write(self.readlinevcf1 + "\n")
                    else:
                        continue
            elif not vcf1line:
                return None, None, None, None
            elif vcf1line == '':
                continue
            elif options.heterozygous == True:
                if vcf1line.split()[7].split(";")[2] == "HET=1":
                    return self.get_line_data(vcf1line)
                else:
                    continue
            elif options.homozygous == True:
                if vcf1line.split()[7].split(";")[3] == "HOM=1":
                    return self.get_line_data(vcf1line)
                else:
                    continue

            else:
                return self.get_line_data(vcf1line)

    def get_vcf2_linearray(self):
        while True:
            vcf2line=self.read_vcf2_line().rstrip()
            self.readlinevcf2=vcf2line
            if vcf2line.startswith("#"):
                if options.twoone or options.common21 or options.alternate21:
                    if options.vcfout:
                        self.outputfh.write(self.readlinevcf2 + "\n")
                    else:
                        continue
            elif not vcf2line:
                return None, None, None, None
            elif vcf2line == '':
                continue
            elif options.heterozygous == True:
                if vcf2line.split()[7].split(";")[2] == "HET=1":
                    return self.get_line_data(vcf2line)
                else:
                    continue
            elif options.homozygous == True:
                if vcf2line.split()[7].split(";")[3] == "HOM=1":
                    return self.get_line_data(vcf2line)
                else:
                    continue

            else:
                return self.get_line_data(vcf2line)

    def get_snps_positions(self):
        self.open_input_files()

        vcf1chr, vcf1pos, vcf1ref, vcf1alt = self.get_vcf1_linearray()
        vcf2chr, vcf2pos, vcf2ref, vcf2alt = self.get_vcf2_linearray()

        while True:
            if vcf1alt == ".":
                vcf1alt = vcf1ref

            if vcf2alt == ".":
                vcf2alt = vcf2ref

            #print vcf1chr, vcf1pos, vcf1alt, vcf2chr, vcf2pos, vcf2alt
            if vcf1chr == None and vcf2chr == None:
                break
            elif vcf1chr and vcf2chr and vcf1chr == vcf2chr:
                if int(vcf1pos) > int(vcf2pos):
                    if self.mode == "twoone" or self.mode == "common21" or self.mode == "alternate21":
                        self.check_ref_alt_bases_and_write(self.readlinevcf2)

                    vcf2chr, vcf2pos, vcf2ref, vcf2alt = self.get_vcf2_linearray()

                elif int(vcf1pos) == int(vcf2pos):
                    #print vcf1chr, vcf1pos, vcf1alt, vcf2chr, vcf2pos, vcf2alt
                    vcfline_with_greater_dp = self.get_higher_dp_vcfline(self.readlinevcf1, self.readlinevcf2)
                    if vcf1alt != vcf2alt and (self.mode == "alternate" or self.mode == "alternate12" or self.mode == "alternate21") :
                        if self.mode == "alternate":
                            self.write_to_output(self.readlinevcf1)
                        elif self.mode == "alternate12":
                            self.write_to_output(self.readlinevcf1)
                        elif self.mode == "alternate21":
                            self.write_to_output(self.readlinevcf2)
                            #self.write_to_output(vcf2chr + " " + str(vcf2pos) + " " + vcf2ref + " " + vcf2alt, self.readlinevcf2)

                    elif vcf1alt == vcf2alt and (self.mode == "common" or self.mode == "common12" or self.mode == "common21") :
                        self.check_ref_alt_bases_and_write(vcfline_with_greater_dp)
                        #self.write_to_output(vcf2chr + " " + str(vcf2pos) + " " + vcf2ref + " " + vcf2alt, self.readlinevcf2)

                    else:
                        pass

                    vcf1chr, vcf1pos, vcf1ref, vcf1alt = self.get_vcf1_linearray()
                    vcf2chr, vcf2pos, vcf2ref, vcf2alt = self.get_vcf2_linearray()

                elif int(vcf1pos) < int(vcf2pos):

                    if self.mode == "onetwo" or self.mode == "common12" or self.mode == "alternate12":
                        #self.write_to_output(vcf1chr + " " + str(vcf1pos) + " " + vcf1ref + " " + vcf1alt, self.readlinevcf1)
                        self.check_ref_alt_bases_and_write(self.readlinevcf1)
                    vcf1chr, vcf1pos, vcf1ref, vcf1alt = self.get_vcf1_linearray()


            else:
                chr_sorted=sort_2_chromosomes(vcf1chr, vcf2chr)
                #print chr_sorted, vcf1chr, vcf2chr, vcf1pos, vcf1alt, vcf2pos, vcf2alt

                if vcf1chr == None and vcf2chr != None:
                    if self.mode == "twoone" or self.mode=="common21" or self.mode == "alternate21":
                        #self.write_to_output(vcf2chr + " " + str(vcf2pos) + " " + vcf2ref + " " + vcf2alt, self.readlinevcf1)
                        self.check_ref_alt_bases_and_write(self.readlinevcf2)
                    vcf2chr, vcf2pos, vcf2ref, vcf2alt = self.get_vcf2_linearray()
                elif vcf1chr != None and vcf2chr == None:
                    if self.mode == "onetwo" or self.mode == "common12" or self.mode == "alternate12":
                        #self.write_to_output(vcf1chr + " " + str(vcf1pos) + " " + vcf1ref + " " + vcf1alt, self.readlinevcf1)
                        self.check_ref_alt_bases_and_write(self.readlinevcf1)
                    vcf1chr, vcf1pos, vcf1ref, vcf1alt = self.get_vcf1_linearray()
                elif vcf1chr != None and vcf2chr != None:
                    if vcf1chr == chr_sorted[0]:
                        if self.mode == "onetwo" or self.mode == "common12" or self.mode == "alternate12":
                            #self.write_to_output(vcf1chr + " " + str(vcf1pos) + " " + vcf1ref + " " + vcf1alt, self.readlinevcf1)
                            self.check_ref_alt_bases_and_write(self.readlinevcf1)
                        vcf1chr, vcf1pos, vcf1ref, vcf1alt = self.get_vcf1_linearray()
                    elif vcf2chr == chr_sorted[0]:
                        if self.mode == "twoone" or self.mode == "common21" or self.mode == "alternate21":
                            #self.write_to_output(vcf2chr + " " + str(vcf2pos) + " " + vcf2ref + " " + vcf2alt, self.readlinevcf1)
                            self.check_ref_alt_bases_and_write(self.readlinevcf2)
                        vcf2chr, vcf2pos, vcf2ref, vcf2alt = self.get_vcf2_linearray()

                else:
                    pass

        self.close_input_files()

def func_main():
    if not options.vcf1:
        print "You must supply first VCF file. Use option -1 or --vcf1 to provide first vcf file"
        exit(1)
    elif not options.vcf2:
        print "You must supply second VCF file. Use option -2 or --vcf2 to provide second vcf file"
        exit(1)
    elif not options.common and not options.alternate and not options.onetwo and not options.twoone and not options.common12 and not options.common21 and not options.alternate12 and not options.alternate21:
        print "You must specify at least one mode. The valid modes are common, alternate, onetwo, twoone, common12, common21, alternate12, alternate21"
        print "Please check the help page using --help option"
        exit(1)
    else:

        result = compare_2_vcfs(options.vcf1, options.vcf2, options.output)

        if options.common == True:
            result.define_mode("common")
            #result.add_header("Chromosome Position Ref Alt")
        elif options.alternate == True:
            result.define_mode("alternate")
            #result.add_header("Chromosome Position Ref Alt-1 Alt-2")
        elif options.onetwo == True:
            result.define_mode("onetwo")
            #result.add_header("Chromosome Position Ref Alt")
        elif options.twoone == True:
            result.define_mode("twoone")
            #result.add_header("Chromosome Position Ref Alt")
        elif options.common12:
            result.define_mode("common12")
            #result.add_header("Chromosome Position Ref Alt")
        elif options.common21:
            result.define_mode("common21")
            #result.add_header("Chromosome Position Ref Alt")
        elif options.alternate12:
            result.define_mode("alternate12")
            #if not options.vcfout:
            #    result.add_header("Chromosome Position Ref Alt-1 Alt-2")
        elif options.alternate21:
            result.define_mode("alternate21")
            #if not options.vcfout:
            #    result.add_header( "Chromosome Position Ref Alt-1 Alt-2")
        else:
            pass

        result.get_snps_positions()
        result.close_input_files()
        result.close_output_file()

if __name__ == '__main__':
    func_main()
    exit(0)
