## This WDL pipeline implements SV calling of mobile elements with Scramble by Rebecca I. Torene.
## https://github.com/GeneDx/scramble
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in mapped BAM format

version 1.0

import "Structs.wdl"

workflow Scramble {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    String sample_name
    File reference_fasta
    Boolean detect_deletions
    String scramble_docker
    RuntimeAttr? runtime_attr_scramble
  }
    
  parameter_meta {
    bam_or_cram_file: "A .bam or .cram file to search for SVs. crams are preferable because they localise faster and use less disk."
    bam_or_cram_index: "Index for bam_or_cram_file."
    sample_name: "A sample name. Outputs will be sample_name+'.scramble.insertions.vcf.gz' and sample_name+'.scramble.deletions.vcf.gz'."
    reference_fasta: "A FASTA file with the reference used to align bam or cram file."
    detect_deletions: "Run deletion detection as well as mobile element insertion."
  }
  
  meta {
    author: "Ted Sharpe, et al"
    email: "tsharpe@broadinstitute.org"
  }

  call RunScramble {
    input:
      bam_or_cram_file = bam_or_cram_file,
      bam_or_cram_index = bam_or_cram_index,
      sample_name = sample_name,
      reference_fasta = reference_fasta,
      detect_deletions = detect_deletions,
      scramble_docker = scramble_docker,
      runtime_attr_override = runtime_attr_scramble
  }

  output {
    File insertions_vcf = RunScramble.insertions_vcf
    File insertions_index = RunScramble.insertions_index
    File? deletions_vcf = RunScramble.deletions_vcf
    File? deletions_index = RunScramble.deletions_index
  }
}

task RunScramble {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    String sample_name
    File reference_fasta
    Boolean detect_deletions
    String scramble_docker
    RuntimeAttr? runtime_attr_override
    File? NOT_A_FILE
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 16,
    disk_gb: 250,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File insertions_vcf = "~{sample_name}.scramble.insertions.vcf.gz"
    File insertions_index = "~{sample_name}.scramble.insertions.vcf.gz.tbi"
    File? deletions_vcf = if detect_deletions then "~{sample_name}.scramble.deletions.vcf.gz" else NOT_A_FILE
    File? deletions_index = if detect_deletions then "~{sample_name}.scramble.deletions.vcf.gz.tbi" else NOT_A_FILE
  }
  command <<<
    set -euo pipefail

    xDir=$PWD
    scrambleDir="/app"
    meiRef=$scrambleDir/cluster_analysis/resources/MEI_consensus_seqs.fa

    # create a blast db from the reference
    cat ~{reference_fasta} | makeblastdb -in - -parse_seqids -title ref -dbtype nucl -out ref

    # Scramble first step: identify clusters of split reads
    $scrambleDir/cluster_identifier/src/build/cluster_identifier ~{bam_or_cram_file} > clusters.txt

    # split the file of clusters to keep memory bounded
    # The awk script removes lines where field 4 (the left consensus) contains nothing but 'n's
    #   because the deletion hunter in Scramble barfs on these.
    awk 'length(gensub(/n/,"","g",$4))' clusters.txt | split -a3 -l2500 - xyzzy

    # Scramble 2nd step: produce a *_MEIs.txt file for each split
    for fil in xyzzy???
      do Rscript --vanilla $scrambleDir/cluster_analysis/bin/SCRAMble.R --out-name $xDir/$fil \
           --cluster-file $xDir/$fil --install-dir $scrambleDir/cluster_analysis/bin \
           --mei-refs $meiRef \
           --ref $xDir/ref --no-vcf --eval-meis ~{true='--eval-dels' false='' detect_deletions}
    done

    # create a header for the output vcf(s)
    echo \
    '##fileformat=VCFv4.3
    ##reference=~{reference_fasta}
    ##source=scramble' > hdr

    grep '^>' $meiRef | awk \
    '{mei=toupper(substr($0,2)); if (mei=="L1") mei="LINE1"
      print "##ALT=<ID=INS:ME:" mei ",Description=\"" mei " element insertion\">"}' >> hdr

    echo \
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
    ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
    ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
    ##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description="Source algorithms">
    ##INFO=<ID=STRANDS,Number=1,Type=String,Description="Breakpoint strandedness [++,+-,-+,--]">
    ##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate">
    ##INFO=<ID=MEI_START,Number=1,Type=Integer,Description="Start of alignment to canonical MEI sequence">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> hdr

    blastdbcmd -db ref -entry all -outfmt '##contig=<ID=%a,length=%l>' >> hdr
    echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	~{sample_name}" >> hdr

    cp hdr insertsVCF.tmp

    # use awk to write the first part of an awk script that initializes an array mapping MEI name
    #   onto its consensus sequence length
    awk \
    'BEGIN { FS=OFS="\t"; print "BEGIN{" }
     /^>/  {if ( seq != "" ) print "seqLen[\"" seq "\"]=" len; seq = substr($0,2); len = 0}
     !/^>/ {len += length($0)}
     END   {if ( seq != "" ) print "seqLen[\"" seq "\"]=" len; print "}"}' $meiRef > awkScript.awk

    # write the rest of the awk script that transforms the contents of the *_MEIs.txt files into a VCF
    echo \
    'BEGIN{ FS=OFS="\t" }
    { if(FNR<2)next
      split($1,loc,":")
      start=loc[2]+1
      end=start+1
      len=seqLen[$2]-$10
      mei=toupper($2); if (mei=="L1") mei="LINE1"
      print loc[1],start,".","N","<INS:ME:" mei ">",int($6),"PASS",\
            "END=" end ";SVTYPE=INS;SVLEN=" len ";MEI_START=" $10 ";STRANDS=+-;CHR2=" loc[1] ";ALGORITHMS=scramble",\
            "GT","0/1" }' >> awkScript.awk

    # transform the MEI descriptions into VCF lines
    awk -f awkScript.awk xyzzy???_MEIs.txt >> insertsVCF.tmp

    # sort and index the output VCF for insertions
    bcftools sort -Oz <insertsVCF.tmp >"~{sample_name}.scramble.insertions.vcf.gz"
    bcftools index -ft "~{sample_name}.scramble.insertions.vcf.gz"

    if [ ~{detect_deletions} == "true" ]
    then
      cp hdr deletesVCF.tmp

      # transform the *_PredictedDeletions.txt files into an output VCF for deletions
      awk \
      'BEGIN{ FS=OFS="\t" }
      { if(FNR<2)next
        Q= $11=="NA" ? ($15=="NA"?".":$15) : ($15=="NA"?$11:($11+$15)/2)
        print $1,$2+1,".","N","<DEL>",Q=="."?".":int(Q),"PASS",\
              "END=" $3+1 ";SVTYPE=DEL;SVLEN=" $5 ";STRANDS=+-;CHR2=" $1 ";ALGORITHMS=scramble",\
              "GT","0/1" }' xyzzy???_PredictedDeletions.txt >> deletesVCF.tmp

      # sort and index the output VCF for deletions
      bcftools sort -Oz <deletesVCF.tmp > "~{sample_name}.scramble.deletions.vcf.gz"
      bcftools index -ft "~{sample_name}.scramble.deletions.vcf.gz"
    fi
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: scramble_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
