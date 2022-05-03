nextflow.enable.dsl = 2

// Task description
// Develop a nextflow pipeline that includes the following analysis steps:
// o Trimming the data using fastp (with default parameters)
// o Generation of a quality report on the trimmed data using fastqc
// o Resistance prediction using srst2
//
// o The pipeline should meet the following requirements:
// o All external tools (fastp, fastqc, srst2) should be run in
// singularity images
// o The pipeline should get a whole folder with fastq files as an input
// o Analysis should be performed on all fastq files (parallelization
// where possible)
// o Output results should be summarized in a single file
// o The FASTA-file containing the resistance prediction reference
// database should be supplied as a parameter on command line
// o To test and apply the pipeline please use CARD_v3.0.8_SRST2 database
// from SRST2 repository
// (https://github.com/katholt/srst2/blob/master/data/CARD_v3.0.8_SRST2.
// fasta)

params.outdir = "LEK_results"
params.indir = "/home/cq/CQ-Bildung-ABI2022-1/Next_Generation_Sequencing/LEK/rawdata"
params.gendb = false

log.info """\
 SAHMORIE CAMERON - ANTIBIOTIC RESISTANCE PIPELINE
 ===================================
 outdir       : ${params.outdir}
 indir        : ${params.indir}
 Gene Database  : ${params.gendb}

 est. 2022
 """

 process fastp {
   container "https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0"
   publishDir "${params.outdir}/fastp", mode: "copy", overwrite: true
   input:
   path infastq
   output:
   path "${infastq.getSimpleName()}", emit: fastq_out
   path "${infastq.getSimpleName()}.html"
   path "${infastq.getSimpleName()}.json", emit: fastpreport
   script:
     """
     fastp -i ${infastq} -o ${infastq.getSimpleName()}.fastq -h ${infastq.getSimpleName()}.html -j ${infastq.getSimpleName()}.json
     """
 }

 process fastqc {
   publishDir "${params.outdir}/qc", mode: "copy", overwrite: true
   container "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"
   input:
     path fastq
    output:
     path "${fastq}"
    script:
     """
     fastqc ${fastq}
     """
 }

 process srst2 {
   publishDir "${params.outdir}/srst", mode: "copy", overwrite: true
   container "https://depot.galaxyproject.org/singularity/srst2%3A0.2.0--py27_2"
   input:
     path fastq
    output:
     path "*"
   script:
     """
      srst2 --input_se INPUT --output . --log --gene_db ${params.gendb}
     """
 }

 workflow {
   fastqs = channel.fromPath("/home/cq/Antibiotic-Resistance-Pipeline-Nextflow-/rawdata/*.fasta").collect()
   fastqs.view()
   fastpout = fastp(fastqs)
   qcied = fastqc(all_fastqs.flatten(), params.patient1)
   srst2_input_channel = fastpout.fastpreport.concat(qcied.zipped)
   srstresults = srst2(srst2_input_channel.collect())
  }



//nextflow /home/cq/CQ-Bildung-ABI2022-1/Next_Generation_Sequencing/LEK.nf --patient1 /home/cq/CQ-Bildung-ABI2022-1/Next_Generation_Sequencing/LEK/rawdata/patient1.fastq --outdir LEK_Results -profile singularity

//nextflow /home/cq/CQ-Bildung-ABI2022-1/Next_Generation_Sequencing/LEK.nf --patient1 /home/cq/CQ-Bildung-ABI2022-1/Next_Generation_Sequencing/LEK/rawdata/patient1.fastq --outdir LEK_Results --reference /home/cq/CQ-Bildung-ABI2022-1/Next_Generation_Sequencing/LEK/rawdata/card.fasta  -profile singularity
