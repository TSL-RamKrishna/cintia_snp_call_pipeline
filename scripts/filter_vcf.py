#!/usr/bin/env python
import sys, re
#############################
### program to filter VCF snps with following filter criteria
### 1) snp position should not have another snp 200 bp forward or backward
### 2) min depth = 100
#############################
### usage : python scriptname refseq.lengths vcf.filename outputfilename

min_distance_between_snps = 200
min_dp_threshold=100

def filtersnps_by_position(positions, reflength):

    # min array size should be 3.
    selected_positions = []
    positions=sorted(map(lambda x:int(x), positions)) # incase there are string numbers
    if len(positions) ==1:
        if positions[0] >= min_distance_between_snps and reflength - min_distance_between_snps >= positions[0]:
            selected_positions.append(positions[0])
    elif len(positions) == 2:
        if positions[1] - positions[0] > 200:
            if positions[0] >= min_distance_between_snps:
                selected_positions.append(positions[0])
            if reflength - min_distance_between_snps >= positions[1]:
                selected_positions.append(positions[1])
    elif len(positions) >=3:

        #if positions[0] > min_distance_between_snps and positions[1] - positions[0] >= min_distance_between_snps:
        if positions[1] - positions[0] >= min_distance_between_snps:
            selected_positions.append(positions[0])

        for pos in range(1, len(positions)-1):
            if positions[pos] - positions[pos-1] > min_distance_between_snps and positions[pos+1] - positions[pos] > min_distance_between_snps:
                selected_positions.append(positions[pos])

        #if reflength - min_distance_between_snps > positions[-1] and positions[-1] - positions[-2] > min_distance_between_snps:
        if positions[-1] - positions[-2] > min_distance_between_snps:
            selected_positions.append(positions[-1])

    return selected_positions

def get_highest_dp(vcfline):
    dparray=re.findall('DP\d*=\d{1,10}', vcfline)
    highest_dp=sorted(map(lambda x: int(x.split("=")[1]), dparray), reverse=True)[0]
    return highest_dp


def main(ref_lengths, vcf_file, filtered_vcf_file):

    readlengthfh=open(ref_lengths)      # File with reference sequences and their lengths
    readlengths=dict()
    for line in readlengthfh:
        line=line.rstrip()
        if line=="":
            continue
        else:
            linearray=line.split()
            readlengths[linearray[0]] = int(linearray[1])

    readlengthfh.close()

    vcfhandle=open(vcf_file)     #VCF file after comparing Sparent and Sbulk with Rparent

    filtersnps=dict()
    for line in vcfhandle:
        line=line.rstrip()
        if line=="":
            continue
        else:
            reference = line.split()[0]
            position = int(line.split()[1])
            if reference in filtersnps.keys():
                filtersnps[reference]["positions"].append(position)
            else:
                filtersnps[reference]={"positions":[position]}

    vcfhandle.close()

    # lets select positions with snps that has no snps around 200 positions

    for refseqkey in filtersnps.keys():
        #print filtersnps[refseqkey]["positions"]
        #print readlengths[refseqkey]
        selected_positions = filtersnps_by_position(filtersnps[refseqkey]["positions"], readlengths[refseqkey])
        filtersnps[refseqkey]["selected_positions"]=selected_positions

    # now we have selected the positions, lets get the snps from vcf file if read depth is above minimum threshold
    vcfhandle=open(vcf_file)
    output = open(filtered_vcf_file, "w")
    for vcfline in vcfhandle:
        vcfline=vcfline.rstrip()
        if line=="":
            continue
        else:
            reference = vcfline.split()[0]
            position = int(vcfline.split()[1])
            if get_highest_dp(vcfline) >=min_dp_threshold:
                if position in filtersnps[reference]["selected_positions"]:
                    output.write(vcfline + "\n")
            else:
                continue

    vcfhandle.close()
    output.close()

if __name__ ==  "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
exit(0)
