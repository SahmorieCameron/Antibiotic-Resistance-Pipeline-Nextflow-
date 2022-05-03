nextflow.enable.dsl = 2

process bowtie2_index{
  publishDir "${params.outdir}", mode: 'copy', overwrite: true
  container "https://depot.galaxyproject.org/singularity/bowtie2%3A2.4.5--py39hd2f7db1_2"
  input:
    path reference
  output:
    path "index*"
  script:
    """
    bowtie2-build ${reference} index
    """

}

process bowtie2_map{
  publishDir "${params.outdir}", mode: 'copy', overwrite: true
  container "https://depot.galaxyproject.org/singularity/bowtie2%3A2.4.5--py39hd2f7db1_2"
  input:
    path fastq
    path index
  output:
    path "mapped.sam"
  script:
    """
    bowtie2 -x index -U ${fastq} -S mapped.sam
    """
}

workflow{
  reffile = channel.fromPath(params.ref)
  fastqs = channel.fromPath("${params.indir}/*.fastq").collect()
  indexchannel = bowtie2_index(reffile)
  //indexchannel.view()
  bowtie2_map(fastqs, indexchannel)
}
