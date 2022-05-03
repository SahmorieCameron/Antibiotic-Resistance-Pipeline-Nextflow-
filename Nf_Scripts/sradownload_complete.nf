nextflow.enable.dsl = 2

//  container "https://depot.galaxyproject.org/singularity/"

// everyting
// the main difference between storeDir and publishDir:
// store dir SKIPS the process if the file exists at destination
// publishDir RUNS the process, and overwrites the file to the destination

// container: if mentioned
// AND -profile singularity is given on a command line
// AND we have nextflow.config with singularity instruction
// if singularity image is not in the singularity cache directory
// it is downloaded to the singularity cacheDir (mentioned in nextflow.config)
// if singularity image is in the cacheDir, it will be started when starting a process

process prefetch {
  storeDir "${params.outdir}" //storing moves the output of the process to store directory from the work folder
  container "https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5262h314213e_0"
  input:
    val sra_id
  output:
    path "${sra_id}/${sra_id}.sra" // this folder structure is the special one produced by prefeth!
  script:
    """
    prefetch ${sra_id} # produces the SRR*/SRR*.sra
    """
}

process fastqdump {
  container "https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5262h314213e_0"
  publishDir "${params.outdir}", mode: "copy", overwrite: true // copies output to the mentioned dir.
  input:
    path infile_sra // this file path will be linked into a working dir of this process!
  output:
    path "${infile_sra.getSimpleName()}*.fastq" // SRR*.fastq
  script:
    """
    fasterq-dump "${infile_sra}" # translates to fastq-dump SRR*.sra
    """
}


// ${infastq[0].getSimpleName()} vs ${infastq.getSimpleName()} depending on the input!
process fastp {
  container "https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0"
  publishDir "${params.outdir}/fastp", mode: "copy", overwrite: true
  input:
    path infastq
    val accession
  output:
    path "${accession}*fastp.fastq", emit: fastq_cut
    path "${accession}*fastp.html"
    path "${accession}*fastp.json", emit: fastpreport
  script:
  if(infastq instanceof List){
    """
    fastp -i ${infastq[0]} -o ${infastq[0].getSimpleName()}_fastp.fastq -I ${infastq[1]} -O ${infastq[1].getSimpleName()}_fastp.fastq -h ${accession}_fastp.html -j ${accession}_fastp.json
    """
  } else{
    """
    fastp -i ${infastq} -o ${infastq.getSimpleName()}_fastp.fastq -h ${infastq.getSimpleName()}_fastp.html -j ${accession}_fastp.json
    """
  }

}

process fastqc {
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  container "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"
  input:
    path fastq
    val accession
  output:
    path "${accession}*_fastqc.html"
    path "${accession}*_fastqc.zip", emit: zipped
  script:
    """
    fastqc ${fastq}
    """
}

process multiqc {
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  container "https://depot.galaxyproject.org/singularity/multiqc:1.9--py_1"
  input:
    path report_files
    val accession
  output:
    path "${accession}_multiqc_report.html"
  script:
    """
    multiqc . -n ${accession}_multiqc_report.html
    """
}

workflow {
  // starting process prefetch with a value channel params.accession from command line (--accession)
  // throw the prefetch output path into a channel called sraresult
  // pipe the sraresult channel into the process fastqdump
  sraresult = prefetch(params.accession)
  fastqs = fastqdump(sraresult)
  fastpout = fastp(fastqs, params.accession)
  // this channel has both trimmed and original fastqs!
  all_fastqs = fastqs.concat(fastpout.fastq_cut)
  qcied = fastqc(all_fastqs.flatten(), params.accession)
  multiqc_input_channel = fastpout.fastpreport.concat(qcied.zipped)
  multiqc(multiqc_input_channel.collect(), params.accession)
}
