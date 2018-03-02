ENV["R1"] ? @read1 = ENV["R1"] : nil
#@read1 ? @R1_dir=ENV["R1"].pathmap("%d") : nil
#@read1 ? @R1_basename=ENV["R1"].pathmap("%n") : nil
@read1 ? @R1_basename  = "#{@read1}".split(".")[0].split("/")[-1] : nil
ENV["R2"] ? @read2 = ENV["R2"] : nil
#@read2 ? @R2_dir=ENV["R2"].pathmap("%d") :nil
#@read2 ? @R2_basename=ENV["R2"].pathmap("%n") : nil
@read2 ? @R2_basename = "#{@read2}".split(".")[0].split("/")[-1] : nil
#ENV["samplename"] ? @sample=ENV["samplename"] : nil
ENV["sampleid"] ? @sampleid=ENV["sampleid"] : nil
ENV["sbulksample"] ? @sbulksample=ENV["sbulksample"] : nil
ENV["reference"] ? @reference=ENV["reference"] : nil

ENV["Rreference"] ? @Rreference=ENV["Rreference"] : nil
ENV["Sreference"] ? @Sreference=ENV["Sreference"] : nil

ENV["projectdir"] ? @projectdir=ENV["projectdir"] : nil


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
  task :runR1 => ["R1"] #these tasks are in case if you are running single end reads
  task :runR2 => ["R2"]

  multitask :run => [:runR1, :runR2] do
    puts "FASTQC completed"
  end


end

namespace :trimmomatic do
  desc "Runs Trimmomatic quality trimming tool"

  file "results/#{@sample}/#{@R1_basename}_trimmed_paired.fastq" => ["#{@read1}", "#{@read2}"] do
    sh "source trimmomatic-0.36; source jre-1.7.0.11; trimmomatic PE -threads 2 -phred33 -trimlog results/#{@sample}/trimmomatic.log -quiet -validatePairs  #{@read1} #{@read2} results/#{@sample}/#{@R1_basename}_trimmed_paired.fastq results/#{@sample}/#{@R1_basename}_unpaired.fastq results/#{@sample}/#{@R2_basename}_trimmed_paired.fastq results/#{@sample}/#{@R2_basename}_unpaired.fastq ILLUMINACLIP:TruSeq2-PE_CKedited.fa:2:30:10 LEADING:15 SLIDINGWINDOW:4:20 TRAILING:15 MINLEN:65"
  end

  file "results/#{@sample}/#{@R2_basename}_trimmed_paired.fastq" => ["#{@read1}", "#{@read2}"] do
  end

  multitask :run =>  ["fastqc:run", "results/#{@sample}/#{@R1_basename}_trimmed_paired.fastq", "results/#{@sample}/#{@R2_basename}_trimmed_paired.fastq"] do
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



 file "#{@Rreference}.1.bt2" => ["#{@Rreference}"] do
   sh "source bowtie2-2.1.0;  bowtie2-build -f   #{@Rreference}  #{@Rreference} "
 end
 file "#{@Rreference}.2.bt2" => ["#{@Rreference}"] do
 end
 file "#{@Rreference}.3.bt2" => ["#{@Rreference}"] do
 end
 file "#{@Rreference}.4.bt2" => ["#{@Rreference}"] do
 end
 file "#{@Rreference}.rev.1.bt2" => ["#{@Rreference}"] do
 end
 file "#{@Rreference}.rev.2.bt2" => ["#{@Rreference}"] do
 end


 task :resistance_run => ["#{@Rreference}.1.bt2", "#{@Rreference}.2.bt2", "#{@Rreference}.3.bt2", "#{@Rreference}.4.bt2", "#{@Rreference}.rev.1.bt2", "#{@Rreference}.rev.2.bt2" ] do
   puts "Bowtie Resistance reference indexing completed"
 end

 file "#{@Sreference}.1.bt2" => ["#{@Sreference}"] do
   sh "source bowtie2-2.1.0;  bowtie2-build -f   #{@Sreference}  #{@Sreference} "
 end
 file "#{@Sreference}.2.bt2" => ["#{@Sreference}"] do
 end
 file "#{@Sreference}.3.bt2" => ["#{@Sreference}"] do
 end
 file "#{@Sreference}.4.bt2" => ["#{@Sreference}"] do
 end
 file "#{@Sreference}.rev.1.bt2" => ["#{@Sreference}"] do
 end
 file "#{@Sreference}.rev.2.bt2" => ["#{@Sreference}"] do
 end

 task :susceptible_run => ["#{@Sreference}.1.bt2", "#{@Sreference}.2.bt2", "#{@Sreference}.3.bt2", "#{@Sreference}.4.bt2", "#{@Sreference}.rev.1.bt2", "#{@Sreference}.rev.2.bt2" ] do
   puts "Bowtie Susceptible reference indexing completed"
 end

end

namespace :Bowtie do
  directory "results/#{@sample}/mappedToRparent"
  directory "results/#{@sample}/mappedToSusparent"

  file "#{@Rreference}.dict" => ["#{@Rreference}"] do
    sh "source samtools-1.3.1; samtools dict -o #{@Rreference}.dict #{@Rreference}"
  end
  file "#{@Sreference}.dict" => ["#{@Sreference}"] do
    sh "source samtools-1.3.1; samtools dict -o #{@Sreference}.dict #{@Sreference}"
  end

 file "results/#{@sample}/mappedToRparent/#{@sampleid}_paired_aligned.sam" => ["BowtieIndex:resistance_run", "results/#{@sample}/mappedToRparent", "#{@Rreference}", "results/#{@sample}/#{@R1_basename}_trimmed_paired.fastq",  "results/#{@sample}/#{@R2_basename}_trimmed_paired.fastq"] do
  sh "source bowtie2-2.1.0; bowtie2 -q --phred33 -k 1 --reorder --very-sensitive-local --no-unal --rg-id #{@sampleid} --rg \"platform:Illumina\" --no-unal -x #{@Rreference} -1 results/#{@sample}/#{@R1_basename}_trimmed_paired.fastq -2 results/#{@sample}/#{@R2_basename}_trimmed_paired.fastq -S results/#{@sample}/mappedToRparent/#{@sampleid}_paired_aligned.sam 2> results/#{@sample}/mappedToRparent/#{@sampleid}_aligned.log; "
 end

 file "results/#{@sample}/mappedToRparent/#{@sampleid}_paired_aligned.bam" => [ "#{@Rreference}.dict", "results/#{@sample}/mappedToRparent/#{@sampleid}_paired_aligned.sam" ] do
   sh "source samtools-1.3.1; samtools view -bS -t #{@Rreference}.dict -o results/#{@sample}/mappedToRparent/#{@sampleid}_paired_aligned.bam results/#{@sample}/mappedToRparent/#{@sampleid}_paired_aligned.sam"
 end

 file "results/#{@sample}/mappedToSusparent/#{@sampleid}_paired_aligned.sam" => ["BowtieIndex:susceptible_run", "results/#{@sample}/mappedToSusparent", "#{@Sreference}", "results/#{@sample}/#{@R1_basename}_trimmed_paired.fastq",  "results/#{@sample}/#{@R2_basename}_trimmed_paired.fastq"] do
  sh "source bowtie2-2.1.0; bowtie2 -q --phred33 -k 1 --reorder --very-sensitive-local --no-unal --rg-id #{@sampleid} --rg \"platform:Illumina\" --no-unal -x #{@Sreference} -1 results/#{@sample}/#{@R1_basename}_trimmed_paired.fastq -2 results/#{@sample}/#{@R2_basename}_trimmed_paired.fastq -S results/#{@sample}/mappedToSusparent/#{@sampleid}_paired_aligned.sam 2> results/#{@sample}/mappedToSusparent/#{@sampleid}_aligned.log; "
 end

 file "results/#{@sample}/mappedToSusparent/#{@sampleid}_paired_aligned.bam" => ["#{@Sreference}.dict", "results/#{@sample}/mappedToSusparent/#{@sampleid}_paired_aligned.sam" ] do
   sh "source samtools-1.3.1; samtools view -bS -t #{@Sreference}.dict -o results/#{@sample}/mappedToSusparent/#{@sampleid}_paired_aligned.bam results/#{@sample}/mappedToSusparent/#{@sampleid}_paired_aligned.sam"
 end


 multitask :run => ["results/#{@sample}/mappedToSusparent/#{@sampleid}_paired_aligned.bam", "results/#{@sample}/mappedToRparent/#{@sampleid}_paired_aligned.bam"] do
   puts "Bowtie mapping completed. SAM file converted to BAM. Original SAM file removed."
 end

 #file "results/#{@sample}/#{@sampleid}_aligned.sam" => ["#{@reference}", "results/#{@sample}/#{@R1_basename}_trimmed_paired.fastq", "results/#{@sample}/#{@R2_basename}_trimmed_paired.fastq" ] do
 #sh "source bowtie2-2.1.0; source samtools-1.3.1; bowtie2 -q --phred33 -k 1 --reorder --no-mixed --no-discordant --very-sensitive-local --no-unal --rg-id #{@sampleid} --rg \"platform:Illumina\" --rg \"sequencer:EI\" --un results/#{@sample}/#{@R1_basename}_unaligned_unpaired.fastq --un-conc results/#{@sample}/#{@R1_basename}_unconc.fastq  -x #{@reference} -1  results/#{@sample}/#{@R1_basename}_trimmed_paired.fastq -2 results/#{@sample}/#{@R2_basename}_trimmed_paired.fastq  -S #results/#{@sample}/#{@sampleid}_aligned.sam 2> results/#{@sample}/#{@sampleid}_aligned.log "
 #end

end

namespace :samtools do

  directory "results/#{@sample}/mappedToRparent/mergedbams"
  directory "results/#{@sample}/mappedToSusparent/mergedbams"
  #----------------------------------------------------------------------------------------------------------------------
  rsamfiles=FileList["results/#{@sample}/mappedToRparent/*_paired_aligned.sam"]
  rbamfiles=rsamfiles.pathmap("%X.bam")
  rbamfiles_sorted=rsamfiles.pathmap("%XSorted.bam")
  rindexedfiles=rsamfiles.pathmap("%XSorted.bam.bai")
  Rmerged=rsamfiles.pathmap("%d/merged.bam")
  Rmerged_sorted=rsamfiles.pathmap("%d/merged_sorted.bam")
  Rmerged_sorted_bai=rsamfiles.pathmap("%d/merged_sorted.bam.bai")

  rsamfiles.zip(rbamfiles, rbamfiles_sorted, rindexedfiles).each do |rsam, rbam, rsorted, rbai|
    file rbam => [rsam] do
      sh "source samtools-1.3.1; samtools view -b -h -o #{rbam} -q 30 --threads 4 #{rsam}"
    end
    file rsorted => [rbam] do
      sh "source samtools-1.3.1;  samtools sort --threads 4 -O bam -o #{rsorted} #{rbam} "
    end

    file rbai => [rsorted] do
      sh "source samtools-1.3.1; samtools index #{rsorted} #{rbai}"
    end
  end


  #-----------------------------------------------------------------------------------------------------------------------
  s_samfiles=FileList["results/#{@sample}/mappedToSusparent/*_paired_aligned.sam"]
  s_bamfiles=s_samfiles.pathmap("%X.bam")
  s_bamfiles_sorted=s_samfiles.pathmap("%XSorted.bam")
  s_indexedfiles=s_samfiles.pathmap("%XSorted.bam.bai")
  Smerged=s_samfiles.pathmap("%d/merged.bam")
  Smerged_sorted=s_samfiles.pathmap("%d/merged_sorted.bam")
  Smerged_sorted_bai=s_samfiles.pathmap("%d/merged_sorted.bam.bai")

  s_samfiles.zip(s_bamfiles, s_bamfiles_sorted, s_indexedfiles).each do |s_sam, s_bam, s_sorted, s_bai|
    file s_bam => [s_sam] do
      sh "source samtools-1.3.1; samtools view -b -h -o #{s_bam} -q 30 --threads 4 #{s_sam}"
    end
    file s_sorted => [s_bam] do
      sh "source samtools-1.3.1;  samtools sort --threads 4 -O bam -o #{s_sorted} #{s_bam} "
    end

    file s_bai => [s_sorted] do
      sh "source samtools-1.3.1; samtools index #{s_sorted} #{s_bai}"
    end
  end

  #-------------------------------------------------------------------------------------------------------------------

  file "results/#{@sample}/mappedToRparent/mergedbams/merged.bam" =>  rindexedfiles  do
    sh "source samtools-1.3.1; samtools merge -f -r -O BAM --reference #{@Rreference} results/#{@sample}/mappedToRparent/mergedbams/merged.bam #{rbamfiles_sorted}"
  end
  file "results/#{@sample}/mappedToSusparent/mergedbams/merged.bam" => s_indexedfiles  do
    sh "source samtools-1.3.1; samtools merge -f -r -O BAM --reference #{@Sreference} results/#{@sample}/mappedToSusparent/mergedbams/merged.bam #{s_bamfiles_sorted}"
  end
  #-------------------------------------------------------------------------------------------------------------------

  #file "results/#{@sample}/mergedbams/merged.bam" => [:convert, "results/#{@sample}/mergedbams"] do
  #  sh "source samtools-1.3.1; samtools merge -n -r -f --threads 4 -O BAM --reference #{@reference} results/#{@sample}/mergedbams/merged.bam #{bamfiles_sorted}"
  #end

	file  "results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.bam" =>  "results/#{@sample}/mappedToRparent/mergedbams/merged.bam" do
		sh "source samtools-1.3.1; samtools sort --reference #{@Rreference} --threads 4 -O bam -o results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.bam results/#{@sample}/mappedToRparent/mergedbams/merged.bam"
	end
  file  "results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.bam" =>  "results/#{@sample}/mappedToSusparent/mergedbams/merged.bam" do
		sh "source samtools-1.3.1; samtools sort --reference #{@Sreference} --threads 4 -O bam -o results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.bam results/#{@sample}/mappedToSusparent/mergedbams/merged.bam"
	end

	file "results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.bam.bai" => ["results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.bam"] do
		sh "source samtools-1.3.1; samtools index results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.bam results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.bam.bai"
	end
  file "results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.bam.bai" => ["results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.bam"] do
		sh "source samtools-1.3.1; samtools index results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.bam results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.bam.bai"
	end

  multitask :Rmerge => [ "results/#{@sample}/mappedToRparent/mergedbams", "results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.bam.bai"]

  multitask :Smerge => [ "results/#{@sample}/mappedToSusparent/mergedbams", "results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.bam.bai"]

	multitask merge: [:Rmerge, :Smerge] do
	  puts "Merging done"
	end

  desc "create mpileup from the bam file"
   file "results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.bcf" => [:Rmerge] do
     sh "source samtools-1.3.1; samtools mpileup -d 250 -m 1 -E --BCF -f #{@Rreference} --output results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.bcf results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.bam"
  end
  file "results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.bcf" => [:Smerge] do
    sh "source samtools-1.3.1; samtools mpileup  -d 250 -m 1 -E --BCF -f #{@Sreference} --output results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.bcf results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.bam"
 end

  file "results/#{@sample}/mappedToRparent/mergedbams/mergedSorted_indexed.bcf" => ["results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.bcf"] do
    sh "source bcftools-1.3.1; bcftools index --force results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.bcf > results/#{@sample}/mappedToRparent/mergedbams/mergedSorted_indexed.bcf"
  end
  file "results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted_indexed.bcf" => ["results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.bcf"] do
    sh "source bcftools-1.3.1; bcftools index --force results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.bcf > results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted_indexed.bcf"
  end
  desc "mapped to Rparent"
  file "results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.vcf" => ["results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.bcf"] do
    sh "source bcftools-1.3.1; bcftools call --multiallelic-caller -O v -o results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.vcf results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.bcf"
  end
  desc "mapped to Sparent"
  file "results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.vcf" => ["results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.bcf"] do
    sh "source bcftools-1.3.1; bcftools call --keep-alts --prior-freqs --multiallelic-caller -O v -o results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.vcf results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.bcf"
  end

  #get VCF mpileup directly
  desc "generate mpileup for Rparent"
  file "results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.mpileup" => [:Rmerge] do
    sh "source samtools-1.3.1; samtools mpileup --fasta-ref #{@Rreference} --VCF -u --output-MQ --output-BP --output-tags DP,AD,ADF,ADR,SP --reference #{@Rreference} --output results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.mpileup results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.bam"
  end
  desc "generate mpileup for Sparent"
  file "results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.mpileup" => [:Smerge] do
    sh "source samtools-1.3.1; samtools mpileup --fasta-ref #{@Sreference} --VCF -u --output-MQ --output-BP --output-tags DP,AD,ADF,ADR,SP --reference #{@Sreference} --output results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.mpileup results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.bam"
  end

  task :get_snps => ["results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.vcf", "results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.vcf"]
  multitask :get_mpileups => ["results/#{@sample}/mappedToRparent/mergedbams/mergedSorted.mpileup",  "results/#{@sample}/mappedToSusparent/mergedbams/mergedSorted.mpileup"]

end

namespace :filtersnps do

  file "results/common_Sparent108_#{@sbulksample}_mappedToRparent.vcf" => ["results/Sparent108/mappedToRparent/mergedbams/mergedSorted.vcf", "results/#{@sbulksample}/mappedToRparent/mergedbams/mergedSorted.vcf"] do
    sh "source python-2.7.11; python scripts/get_alternate_common_snps.py -1 results/Sparent108/mappedToRparent/mergedbams/mergedSorted.vcf -2 results/#{@sbulksample}/mappedToRparent/mergedbams/mergedSorted.vcf --common --vcfout --refaltdiff --out result/common_Sparent108_Sbulk293_mappedToRparent.vcf"
  end
  file "results/alternate_mappedToRparent.vcf" => ["results/Rparent614/mappedToRparent/mergedbams/mergedSorted.vcf", "results/common_Sparent108_#{@sbulksample}_mappedToRparent.vcf"] do
    sh "source python-2.7.11; python scripts/get_alternate_common_snps.py -1 results/Rparent614/mappedToRparent/mergedbams/mergedSorted.vcf -2 results/common_Sparent108_#{@sbulksample}_mappedToRparent.vcf --alternate --vcfout --out results/alternate_mappedToRparent.vcf"
  end
  file "results/alternate_mappedToRparent_filtered.vcf" => ["results/alternate_mappedToRparent.vcf"] do
    sh "python scripts/filter_vcf.py #{@Rreference}.fai results/alternate_mappedToRparent.vcf results/alternate_mappedToRparent_filtered.vcf"
  end
  file "results/common_Sparent108_#{@sbulksample}_mappedToSusparent.vcf" => ["results/Sparent108/mappedToSusparent/mergedbams/mergedSorted.vcf", "results/#{@sbulksample}/mappedToSusparent/mergedbams/mergedSorted.vcf"] do
    sh "source python-2.7.11; python scripts/get_alternate_common_snps.py -1 results/Sparent108/mappedToSusparent/mergedbams/mergedSorted.vcf -2 results/#{@sbulksample}/mappedToSusparent/mergedbams/mergedSorted.vcf --common --vcfout --refaltsame --out results/common_Sparent108_Sbulk293_mappedToSusparent.vcf"
  end
  file "results/alternate_mappedToSusparent.vcf" => ["results/Rparent614/mappedToSusparent/mergedbams/mergedSorted.vcf", "results/common_Sparent108_#{@sbulksample}_mappedToSusparent.vcf"] do
    sh "source python-2.7.11; python scripts/get_alternate_common_snps.py -1 results/Rparent614/mappedToSusparent/mergedbams/mergedSorted.vcf -2  results/common_Sparent108_#{@sbulksample}_mappedToSusparent.vcf --alternate --vcfout --out results/alternate_mappedToSusparent.vcf"
  end
  file "results/alternate_mappedToSusparent_filtered.vcf" => ["results/alternate_mappedToSusparent.vcf"] do
    sh "python scripts/filter_vcf.py #{@Sreference}.fai results/alternate_mappedToSusparent.vcf results/alternate_mappedToSusparent_filtered.vcf"
  end

  task :mappedToRparent => ["results/alternate_mappedToRparent_filtered.vcf"]
  task :mappedToSusparent => ["results/alternate_mappedToSusparent_filtered.vcf"]
  task :run => [:mappedToRparent, :mappedToSusparent]
end

namespace :get_snps_subseq do
  file "results/alternate_mappedToRparent_filtered_snp_subseq.fasta" => ["filtersnps:mappedToRparent"] do
    sh "python scripts/get_snp_subsequence.py results/alternate_mappedToRparent_filtered.vcf #{@Rreference} results/alternate_mappedToRparent_filtered_snp_subseq.fasta"
  end
  file "results/alternate_mappedToSusparent_filtered_snp_subseq.fasta" => ["filtersnps:mappedToSusparent"] do
    sh "python scripts/get_snp_subsequence.py results/alternate_mappedToSusparent_filtered.vcf #{@Sreference} results/alternate_mappedToSusparent_filtered_snp_subseq.fasta"
  end

  task :get_Rparent_snp_subseq => ["results/alternate_mappedToRparent_filtered_snp_subseq.fasta"]
  task :get_Susparent_snp_subseq => ["results/alternate_mappedToSusparent_filtered_snp_subseq.fasta"]

  multitask :run => [:get_Rparent_snp_subseq, :get_Susparent_snp_subseq]
end

task :run_pipeline => [ "fastqc:run", "trimmomatic:run", "Bowtie:run", "samtools:run", "VCF:run"] do
  puts "Pipeline completed"
end

task :default => [ "samtools:get_snps"] do
  puts "Pipeline completed"
end
