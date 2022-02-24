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

workflow AnnotateInheritancePerSamplePerVcf{
    input{
        File cleanVcf
        String contig

        String probands_id
        String father_id
        String mother_id

        String sv_pipeline_docker
        String sv_base_mini_docker
        String rdpesr_benchmark_docker

        RuntimeAttr? runtime_attr_split_vcf
        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_bed_comp
    }

    call mini_tasks.split_per_sample_vcf as split_proband {
        input:
            vcf = cleanVcf,
            sample = probands_id,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_split_vcf
    }

    call mini_tasks.split_per_sample_vcf as split_father {
        input:
            vcf = cleanVcf,
            sample = father_id,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_split_vcf
    }

    call mini_tasks.split_per_sample_vcf as split_mother {
        input:
            vcf = cleanVcf,
            sample = mother_id,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_split_vcf
    }

    call mini_tasks.vcf2bed as vcf2bed_proband {
        input:
            vcf = split_proband.vcf_file,
            prefix = probands_id,
            sample = probands_id,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_vcf2bed
    }

    call mini_tasks.vcf2bed as vcf2bed_father {
        input:
            vcf = split_father.vcf_file,
            prefix = father_id,
            sample = father_id,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_vcf2bed
    }

    call mini_tasks.vcf2bed as vcf2bed_mother {
        input:
            vcf = split_mother.vcf_file,
            prefix = mother_id,
            sample = mother_id,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_vcf2bed
    }

    call mini_tasks.InheritanceComparison as Inheritance_proband_vs_father {
        input:
            query_bed = vcf2bed_father.bed,
            ref_bed = vcf2bed_proband.bed,
            prefix = "${probands_id}_vs_father",
            sv_pipeline_docker=sv_pipeline_docker
    }

    call mini_tasks.InheritanceComparison as Inheritance_proband_vs_mother {
        input:
            query_bed = vcf2bed_mother.bed,
            ref_bed = vcf2bed_proband.bed,
            prefix = "${probands_id}_vs_mother",
            sv_pipeline_docker=sv_pipeline_docker
    }

    call IntegrateInheritance {
        input:
             proband_bed = vcf2bed_proband.bed,
             proband_vs_father = Inheritance_proband_vs_father.comparison,
             proband_vs_mother = Inheritance_proband_vs_mother.comparison,
             prefix = "~{probands_id}.~{contig}",
             rdpesr_benchmark_docker = rdpesr_benchmark_docker
    }

    output{
        File inheritance = IntegrateInheritance.inheri
    }
}

task IntegrateInheritance{
    input{
        File proband_bed
        File proband_vs_father
        File proband_vs_mother
        String prefix
        String rdpesr_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 7.5, 
        disk_gb: 15,
        boot_disk_gb: 15,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File inheri = "~{prefix}.inheri.gz"
    }

    command <<<
        Rscript /src/integrate_inheritance.R \
            --bed  ~{proband_bed} \
            --pb_vs_fa ~{proband_vs_father} \
            --pb_vs_mo ~{proband_vs_mother} \
            --output ~{prefix}.inheri 

        bgzip ~{prefix}.inheri 
    >>>

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
