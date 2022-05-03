nextflow.enable.dsl = 2

params.hashlen = 51
params.outdir = "results"

params.with_quastref = false

process velvet{
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  container "https://depot.galaxyproject.org/singularity/velvet:1.2.10--h7132678_5"
  input:
    path fastq
    val hashlen
  output:
    path "velvetdir/LastGraph"
    path "velvetdir/contigs.fa", emit: velvetcontigs
  script:
    if(fastq instanceof List){
      """
      velveth velvetdir ${hashlen} -shortPaired -fastq -separate ${fastq}
      velvetg velvetdir
      """
    } else {
      """
      velveth velvetdir ${hashlen} -fastq -short ${fastq}
      velvetg velvetdir
      """
    }
}

process quast{
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  container "https://depot.galaxyproject.org/singularity/quast:5.0.2--py37pl526hb5aa323_2"
  input:
    path infile
  output:
    path "quast_results"
  script:
    """
    quast ${infile}
    """
}

process quast_ref{
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  container "https://depot.galaxyproject.org/singularity/quast:5.0.2--py37pl526hb5aa323_2"
  input:
    path infile
    path ref
  output:
    path "quast_results"
  script:
    """
    quast ${infile} -r ${ref}
    """
}

//params.indir = null
workflow {
  fastqchannel = channel.fromPath("${params.indir}/*.fastq").collect()
  fastqchannel.view()
  vout = velvet(fastqchannel, params.hashlen)
  if(params.with_quastref){
    reference = channel.fromPath(params.with_quastref)
    quast_ref(vout.velvetcontigs, reference)
  } else {
    quast(vout.velvetcontigs)
  }

}
