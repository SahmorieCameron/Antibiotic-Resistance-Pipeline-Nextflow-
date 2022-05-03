nextflow.enable.dsl = 2

params.kmerlen = 51
params.covcutoff = 0
params.mincontiglength = 300
params.quastref = null

process velvet {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  container "https://depot.galaxyproject.org/singularity/velvet:1.2.10--hed695b0_3"

  input:
    path fastqfiles
    val kmerlen
    val covcutoff
    val mincontiglength

  output:
    path "velvet_assembly", emit: assemblydir
    path "velvet_assembly/velvetcontigs.fa", emit: contigs

  script:
    if(fastqfiles instanceof List) {
      """
      velveth velvet_assembly/ ${kmerlen} -fastq -shortPaired -separate ${fastqfiles}
      velvetg velvet_assembly/ -cov_cutoff ${covcutoff} -min_contig_lgth ${mincontiglength}
      mv velvet_assembly/contigs.fa velvet_assembly/velvetcontigs.fa
      """
    } else {
      """
      velveth velvet_assembly/ ${kmerlen} -fastq -short ${fastqfiles}
      velvetg velvet_assembly/ -cov_cutoff ${covcutoff} -min_contig_lgth ${mincontiglength}
      mv velvet_assembly/contigs.fa velvet_assembly/velvetcontigs.fa
      """
    }

}


process spades {
  storeDir "${params.outdir}"

  container "https://depot.galaxyproject.org/singularity/spades:3.15.3--h95f258a_0"

  input:
    path fastqfiles

  output:
    path "spades_assembly", emit: assemblydir
    path "spades_assembly/spadescontigs.fasta", emit: contigs

  script:
    if(fastqfiles instanceof List) {
      """
      spades.py --only-assembler -o spades_assembly --careful -1 ${fastqfiles[0]} -2 ${fastqfiles[1]}
      mv spades_assembly/contigs.fasta spades_assembly/spadescontigs.fasta
      """
    } else {
      """
      spades.py --only-assembler -o spades_assembly --careful -s ${fastqfiles}
      mv spades_assembly/contigs.fasta spades_assembly/spadescontigs.fasta
      """
    }
}

process quast {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  container "https://depot.galaxyproject.org/singularity/quast:5.0.2--py37pl526hb5aa323_2"

  input:
    path contigs
  output:
    path "quast"

  script:
    """
    quast.py ${contigs} -o quast
    """
}

process quast_ref {
  storeDir "${params.outdir}"

  container "https://depot.galaxyproject.org/singularity/quast:5.0.2--py37pl526hb5aa323_2"

  input:
    path contigs
    path quastref

  output:
    path "quast"

  script:
    """
    quast.py ${contigs} -r ${quastref} -o quast
    """
}

workflow assembly {
  take:
    kmerlen
    covcutoff
    mincontiglength
    quastref
    fastqfiles
  main:
    velvet_out = velvet(fastqfiles, kmerlen, covcutoff, mincontiglength)
    spades_out = spades(fastqfiles)
    all_contigs = velvet_out.contigs.concat(spades_out.contigs)
    //all_contigs.view()
    reference = quastref
    if(quastref) {
      reference = Channel.fromPath(quastref)
      quast_out = quast_ref(all_contigs.collect(), reference)
    } else {
      quast_out = quast(all_contigs.collect())
    }
  emit:
    quastfolder = quast_out
}


workflow {


  fastqfiles = channel.fromPath("${params.indir}/*.fastq").collect()
  outchannel = assembly(params.kmerlen, params.covcutoff, params.mincontiglength, params.quastref, fastqfiles)


}
