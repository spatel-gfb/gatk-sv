version 1.0

import "Structs.wdl"
import "TasksClusterBatch.wdl" as tasks_cluster
import "TasksMakeCohortVcf.wdl" as tasks_cohort
import "Utils.wdl" as util

workflow ClusterPESR {
  input {
    Array[File] vcfs

    File ploidy_table
    String batch
    String caller

    File? svtk_to_gatk_script
    File? gatk_to_svtk_script
    Float? java_mem_fraction

    Int min_size
    File exclude_intervals
    File contig_list

    Float pesr_interval_overlap
    Int pesr_breakend_window
    String? clustering_algorithm

    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_multi_svtk_to_gatk_vcf
    RuntimeAttr? runtime_attr_exclude_intervals_pesr
    RuntimeAttr? runtime_attr_length_filter
    RuntimeAttr? runtime_attr_svcluster
    RuntimeAttr? runtime_override_concat_vcfs_pesr
    RuntimeAttr? runtime_attr_multi_gatk_to_svtk_vcf
    RuntimeAttr? runtime_attr_index_vcfs
  }

  call PreparePESRVcfs {
    input:
      vcfs=vcfs,
      exclude_intervals=exclude_intervals,
      exclude_intervals_index=exclude_intervals + ".tbi",
      ploidy_table=ploidy_table,
      min_size=min_size,
      output_suffix="formatted",
      script=svtk_to_gatk_script,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_multi_svtk_to_gatk_vcf
  }

  scatter (maybe_vcf in PreparePESRVcfs.out) {
    if (basename(maybe_vcf, ".vcf.gz") + ".vcf.gz" == basename(maybe_vcf)) {
      File indexed_vcfs_ = maybe_vcf
    }
  }
  Array[File] indexed_vcfs = select_all(indexed_vcfs_)

  Array[String] contigs = transpose(read_tsv(contig_list))[0]
  scatter (contig in contigs) {
    call tasks_cluster.SVCluster {
      input:
        vcfs=indexed_vcfs,
        ploidy_table=ploidy_table,
        output_prefix="~{batch}.~{caller}.~{contig}.clustered",
        contig=contig,
        fast_mode=true,
        algorithm=clustering_algorithm,
        pesr_sample_overlap=0,
        pesr_interval_overlap=pesr_interval_overlap,
        pesr_breakend_window=pesr_breakend_window,
        reference_fasta=reference_fasta,
        reference_fasta_fai=reference_fasta_fai,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        variant_prefix="~{batch}_~{caller}_~{contig}_",
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_svcluster
    }
  }

  call tasks_cohort.ConcatVcfs {
    input:
      vcfs=SVCluster.out,
      vcfs_idx=SVCluster.out_index,
      naive=true,
      outfile_prefix="~{batch}.~{caller}.clustered",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_vcfs_pesr
  }

  call tasks_cluster.GatkToSvtkVcf {
    input:
      vcf=ConcatVcfs.concat_vcf,
      output_prefix="~{batch}.~{caller}.clustered.svtk_formatted",
      script=gatk_to_svtk_script,
      source=caller,
      contig_list=contig_list,
      remove_formats="CN",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_multi_gatk_to_svtk_vcf
  }

  output {
    File clustered_vcf = GatkToSvtkVcf.out
    File clustered_vcf_index = GatkToSvtkVcf.out_index
  }
}


task PreparePESRVcfs {
  input {
    Array[File] vcfs
    File exclude_intervals
    File exclude_intervals_index
    File ploidy_table
    Int min_size
    File? script
    String? remove_infos
    String? remove_formats
    String output_suffix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(10 + size(vcfs, "GB") * 2),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    Array[File] out = glob("out/*")
  }
  command <<<
    set -euo pipefail
    mkdir out/
    i=0
    while read VCF; do
      NAME=$(basename $VCF .vcf.gz)
      SAMPLE_NUM=`printf %05d $i`
      # Convert format
      python ~{default="/opt/sv-pipeline/scripts/format_svtk_vcf_for_gatk.py" script} \
        --vcf $VCF \
        --out tmp.vcf.gz \
        --ploidy-table ~{ploidy_table} \
        ~{"--remove-infos " + remove_infos} \
        ~{"--remove-formats " + remove_formats}

      # Interval and size filtering
      bcftools query -f '%CHROM\t%POS\t%POS\t%ID\t%SVTYPE\n%CHROM\t%END\t%END\t%ID\t%SVTYPE\n%CHR2\t%END2\t%END2\t%ID\t%SVTYPE\n' tmp.vcf.gz \
        | awk '$1!="."' \
        > ends.bed
      bedtools intersect -wa -a ends.bed -b ~{exclude_intervals} | cut -f4 | sort | uniq | cut -f2 \
        > excluded_vids.list
      bcftools view -i '%ID!=@excluded_vids.list' -e 'INFO/SVLEN<~{min_size}' tmp.vcf.gz -Oz -o out/$SAMPLE_NUM.$NAME.~{output_suffix}.vcf.gz
      tabix out/$SAMPLE_NUM.$NAME.~{output_suffix}.vcf.gz
      i=$((i+1))
    done < ~{write_lines(vcfs)}
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}