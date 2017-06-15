ENV["R1"] ? @read1 = ENV["R1"] : nil
@read1 ? @R1_basename  = "#{@read1}".split(".")[0].split("/")[-1] : nil

ENV["R2"] ? @read2 = ENV["R2"] : nil
@read2 ? @R2_basename = "#{@read2}".split(".")[0].split("/")[-1] : nil

ENV["samplename"] ? @sample=ENV["samplename"] : nil
ENV["sampleid"] ? @sampleid=ENV["sampleid"] : nil
ENV["reference"] ? @reference=ENV["reference"] : nil


directory  "results"


namespace :fastqc  do
  desc "Do fastqc quality check of the input data reads"
  directory  "results/#{@sample}"
  desc "do fastqc for R1"

  file "results/#{@sample}/#{@R1_basename}_fastqc.html" => [ "results/#{@sample}", "#{@read1}"] do
    sh "source fastqc-0.11.5; fastqc -outdir results/#{@sample} -extract #{@read1}"
  end
  file "results/#{@sample}/#{@R1_basename}_fastqc.zip" => [ "results/#{@sample}", "#{@read1}"] do
  end
  file "results/#{@sample}/#{@R1_basename}_fastqc" => [ "results/#{@sample}", "#{@read1}"] do
  end

  task :R1 => ["results/#{@sample}/#{@R1_basename}_fastqc.html", "results/#{@sample}/#{@R1_basename}_fastqc.zip", "results/#{@sample}/#{@R1_basename}_fastqc"] do
    puts "R1 FASTQC completed"
  end

  desc "do fastqc for R2"
  file "results/#{@sample}/#{@R2_basename}_fastqc.html" => [ "results/#{@sample}", "#{@read2}"] do
    sh "source fastqc-0.11.5;fastqc -outdir results/#{@sample} -extract #{@read2}"
  end
  file "results/#{@sample}/#{@R2_basename}_fastqc.zip" => [ "results/#{@sample}", "#{@read2}"] do
  end
  file "results/#{@sample}/#{@R2_basename}_fastqc" => [ "results/#{@sample}", "#{@read2}"] do
  end

  task :R2 => ["results/#{@sample}/#{@R2_basename}_fastqc.html", "results/#{@sample}/#{@R2_basename}_fastqc.zip", "results/#{@sample}/#{@R2_basename}_fastqc"] do
    puts "R2 FASTQC completed"
  end

  task :run => ["R1", "R2"] do
    puts "FASTQC completed"
  end


end


namespace :trimmomatic do
  desc "Runs Trimmomatic quality trimming tool"

  file "results/#{@sample}/#{@R1_basename}_paired.fastq" => ["#{@read1}", "#{@read2}"] do
    sh "source trimmomatic-0.36; source jre-1.7.0.11; trimmomatic PE -threads 2 -phred33 -trimlog results/#{@sample}/trimmomatic.log -quiet -validatePairs  #{@read1} #{@read2} results/#{@sample}/#{@R1_basename}_paired.fastq results/#{@sample}/#{@R1_basename}_unpaired.fastq results/#{@sample}/#{@R2_basename}_paired.fastq results/#{@sample}/#{@R2_basename}_unpaired.fastq ILLUMINACLIP:TruSeq2-PE_CKedited.fa:2:30:10 LEADING:15 SLIDINGWINDOW:4:20 TRAILING:15 MINLEN:65"
  end

  file "results/#{@sample}/#{@R2_basename}_paired.fastq" => ["#{@read1}", "#{@read2}"] do
  end

  task :run =>  ["fastqc:run", "results/#{@sample}/#{@R1_basename}_paired.fastq", "results/#{@sample}/#{@R2_basename}_paired.fastq"] do
    puts "Trimmomatic completed"
  end

end


namespace :BWAindex do
  desc "Creates BWA reference index"

  file "#{@reference}.amb" => ["#{@reference}"] do
    sh "source bwa-0.7.15; bwa index #{@reference}"
  end

  file "#{@reference}.ann" => ["#{@reference}"] do
  end

  file "#{@reference}.bwt" => ["#{@reference}"] do
  end

  file "#{@reference}.pac" => ["#{@reference}"] do
  end

  file "#{@reference}.sa" => ["#{@reference}"] do
  end

  task :run => ["#{@reference}.amb", "#{@reference}.ann", "#{@reference}.bwt", "#{@reference}.pac", "#{@reference}.sa"] do
    puts "Reference indexing completed"
  end
end

#namespace :BWA do
#  desc "Runs BWA alignment of reads to the reference"
#  file "results/#{@sample}/#{@sampleid}_aligned.sam" => ["results/#{@sample}/#{@R1_basename}_paired.fastq", "results/#{@sample}/#{@R2_basename}_paired.fastq"] do
#    sh "source bwa-0.7.15; bwa mem -t 3 #{@reference} results/#{@sample}/#{@R1_basename}_paired_trimmed.fastq results/#{@sample}/#{@R2_basename}_paired_trimmed.fastq 1>results/#{@sample}/#{@sampleid}_aligned.sam 2>results/#{@sample}/#{@sampleid}_bwa_align.log;"
#  end

#  task :run => [ "trimmomatic:run", "BWAindex:run", "results/#{@sample}/#{@sampleid}_aligned.sam"] do
#    puts "bwa alignment completed"
#  end
#end

namespace :BowtieIndex do
 file "#{@reference}.1.bt2" => ["#{@reference}"] do
   sh "source bowtie2-2.1.0;  bowtie2-build -f   #{@reference}  #{@reference} "
 end
 file "#{@reference}.2.bt2" => ["#{@reference}"] do
 end
 file "#{@reference}.3.bt2" => ["#{@reference}"] do
 end
 file "#{@reference}.4.bt2" => ["#{@reference}"] do
 end
 file "#{@reference}.rev.1.bt2" => ["#{@reference}"] do
 end
 file "#{@reference}.rev.2.bt2" => ["#{@reference}"] do
 end

 task :run => ["#{@reference}.1.bt2", "#{@reference}.2.bt2", "#{@reference}.3.bt2", "#{@reference}.4.bt2", "#{@reference}.rev.1.bt2", "#{@reference}.rev.2.bt2" ] do
   puts "Bowtie reference indexing completed"
 end
end

namespace :Bowtie do
 file "results/#{@sample}/#{@sampleid}_aligned.sam" => ["#{@reference}", "results/#{@sample}/#{@R1_basename}_paired.fastq", "results/#{@sample}/#{@R2_basename}_paired.fastq" ] do
 sh "source bowtie2-2.1.0; bowtie2 -q --phred33 -k 1 --reorder --no-mixed --no-discordant --very-sensitive-local --no-unal --rg-id #{@sampleid} --rg \"platform:Illumina\" --rg \"sequencer:EI\" --un results/#{@sample}/#{@R1_basename}_unaligned_unpaired.fastq --un-conc results/#{@sample}/#{@R1_basename}_unconc.fastq  -x #{@reference} -1  results/#{@sample}/#{@R1_basename}_paired.fastq -2 results/#{@sample}/#{@R2_basename}_paired.fastq  -S results/#{@sample}/#{@sampleid}_aligned.sam 2> results/#{@sample}/#{@sampleid}_aligned.log; "
 end

 task :run => ["BowtieIndex:run", "results/#{@sample}/#{@sampleid}_aligned.sam"] do
   puts "Bowtie mapping completed"
 end
end

namespace :samtools do
  desc "create dict of reference"

  file "#{@reference}.dict" => ["#{@reference}"] do
    system "source samtools-1.3.1; samtools dict #{@reference} > #{@reference}.dict"
  end
  desc "Index the reference sequence"
  file "#{@reference}.fai" => ["#{@reference}"] do
    system "source samtools-1.3.1; samtools faidx #{@reference} > #{@reference}.fai"
  end
  desc "convert sam to bam"
  file "results/#{@sample}/#{@sampleid}_aligned.bam" => ["#{@reference}.dict"] do
    system "source samtools-1.3.1; samtools view -bS -t #{@reference}.dict -o results/#{@sample}/#{@sampleid}_aligned.bam results/#{@sample}/#{@sampleid}_aligned.sam "
  end

  desc "sort BAM file"
  file "results/#{@sample}/#{@sampleid}_alignedSorted.bam" => ["results/#{@sample}/#{@sampleid}_aligned.bam"] do
    system "source samtools-1.3.1; samtools sort results/#{@sample}/#{@sampleid}_aligned.bam > results/#{@sample}/#{@sampleid}_alignedSorted.bam"
  end

  desc "Index BAM file"
  file "results/#{@sample}/#{@sampleid}_alignedSorted.bam.bai" => ["results/#{@sample}/#{@sampleid}_alignedSorted.bam"] do
    system "source samtools-1.3.1; samtools index -b results/#{@sample}/#{@sampleid}_alignedSorted.bam > results/#{@sample}/#{@sampleid}_alignedSorted.bam.bai"
  end

  desc "Convert  BAM TO BED"
  file "results/#{@sample}/#{@sampleid}_alignedSorted.bed" => ["results/#{@sample}/#{@sampleid}_alignedSorted.bam"] do
    sh "source bedtools-2.20.1; bedtools bamtobed -i results/#{@sample}/#{@sampleid}_alignedSorted.bam > results/#{@sample}/#{@sampleid}_alignedSorted.bed"
  end

  #desc "Get bedCov from bed file"
  #file "results/#{@sample}/alignedSorted.bedCov" => ["results/#{@sample}/alignedSorted.bed", "results/#{@sample}/alignedSorted.bam"] do
  #  sh "source samtools-1.3.1; samtools bedcov -Q 20 results/#{@sample}/alignedSorted.bed results/#{@sample}/alignedSorted.bam > results/#{@sample}/alignedSorted.bedCov"
  #end

  file "results/#{@sample}/#{@sampleid}_alignedSorted.genomecov.bed" => ["results/#{@sample}/#{@sampleid}_alignedSorted.bam"] do
    sh "source bedtools-2.20.1; bedtools genomecov -d -ibam results/#{@sample}/#{@sampleid}_alignedSorted.bam | awk '{if($3>=1){print}}' > results/#{@sample}/#{@sampleid}_alignedSorted.genomecov.bed"
  end


  desc "create mpileup from the bam file"
  file "results/#{@sample}/#{@sampleid}_mpileup.vcf" => ["results/#{@sample}/#{@sampleid}_alignedSorted.bam", "#{@reference}"] do
    #sh "samtools mpileup -q 20 -Q 15 -d #{@depth} -Bf #{@reference} results/#{@sample}/alignSorted.bam > #{@outdir}/alignedSorted.pileup 2>> #{@outdir}/rake_log.txt"
    sh "source samtools-1.3.1; samtools mpileup -f #{@reference} results/#{@sample}/#{@sampleid}_alignedSorted.bam -C 20 -d 250 -q 20 -Q 13 -O -s --output results/#{@sample}/#{@sampleid}_mpileup.vcf"
    #samtools mpileup -B -f #{@reference}   sample1/sample1.realigned.bam sample2/sample2.realigned.bam sample3/sample3.realigned.bam | varscan mpileup2cns --min-coverage 2 --min-reads2 3 --min-avg_qual 20 --min-var-freq 0.8 --p-value 0.005 --variants --output-vcf > variants.vcf
 end


  task :sambam => ["results/#{@sample}/#{@sampleid}_aligned.bam", "results/#{@sample}/#{@sampleid}_alignedSorted.bam", "results/#{@sample}/#{@sampleid}_alignedSorted.bam.bai"] do
    puts "samtools conversion SAM -> BAM, Sort BAM and Index BAM completed"
  end

  task :bedCoverage => ["results/#{@sample}/#{@sampleid}_alignedSorted.bed", "results/#{@sample}/#{@sampleid}_alignedSorted.genomecov.bed"] do
    puts "BedCoverage completed"
  end

  task :mpileup => ["results/#{@sample}/#{@sampleid}_alignedSorted.vcf"] do
    puts "mpileup completed"
  end

  task :merge  => do

  task :run => [:sambam, :bedCoverage, :mpileup, :merge ] do
    puts "Generating sambam and bedcoverage completed"
  end

end


namespace :VCF do
  file "results/#{@sample}/#{@sampleid}_snps_and_indels.vcf" => "results/#{@sample}/#{@sampleid}_mpileup.vcf" do
    sh "varscan mpileup2cns results/#{@sample}/#{@sampleid}_mpileup.vcf --min-coverage 2 --min-reads2 3 --min-avg_qual 20 --min-var-freq 0.8 --p-value 0.005 --variants --output-vcf > #{@sampleid}_snps_and_indels.vcf"
  end
  task :run => "results/#{@sample}/#{@sampleid}_snps_and_indels.vcf" do
    puts "Consensus and Variants (SNPs and indels) call completed "
  end
end

task :run_pipeline => [ "fastqc:run", "trimmomatic:run", "Bowtie:run", "samtools:run", "VCF:run"] do
  puts "Pipeline completed"
end

task :default => [ "fastqc:run", "trimmomatic:run", "Bowtie:run"] do
  puts "Pipeline completed"
end
