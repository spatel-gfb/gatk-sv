##########################################################################################

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

## Copyright Broad Institute, 2020
## 
## This WDL pipeline implements Duphold 
##
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

version 1.0

import "Structs.wdl"
import "TasksBenchmark.wdl" as mini_tasks


workflow SplitPerSampleGTGQPerVcf{
    input{
        File cleanVcf
        File cleanVcfIdx

        Array[File] SampleLists
        String prefix

        String rdpesr_benchmark_docker

        RuntimeAttr? runtime_split_per_sample_gtgq
    }

    scatter (SampleList in SampleLists){
        call split_per_sample_gtgq{
            input:
                vcf = cleanVcf,
                vcf_idx = cleanVcfIdx,
                sample_list = SampleList,
                chr = prefix,
                rdpesr_benchmark_docker = rdpesr_benchmark_docker,
                runtime_attr_override = runtime_split_per_sample_gtgq
        }
    }

    output{
        Array[Array[File]] gtgq = split_per_sample_gtgq.gtgq_file
    }
}


task split_per_sample_gtgq {
  input {
    File vcf
    File vcf_idx
    File sample_list
    String chr
    String rdpesr_benchmark_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 1.5, 
    disk_gb: ceil(2.0 +  size(vcf, "GB")),
    boot_disk_gb: 30,
    preemptible_tries: 1,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    mkdir per_sample_GTGQ/
    python  /src/split_gtgq_per_sample.py ~{sample_list} ~{vcf} ~{chr}

  >>>

  output {
    Array[File] gtgq_file =  glob("per_sample_GTGQ/*~{chr}.gtgq.gz")
  }
  
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: rdpesr_benchmark_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
