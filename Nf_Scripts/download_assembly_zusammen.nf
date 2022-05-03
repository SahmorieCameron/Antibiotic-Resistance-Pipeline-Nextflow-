nextflow.enable.dsl = 2

include { sradownload } from "./sradownload_module_zusammen"
include { assembly } from "./assembly_module_zusammen"

params.outdir = null
params.accession = null

params.kmerlen = 51
params.covcutoff = 0
params.mincontiglength = 300
params.quastref = null

workflow {
  sraresults = sradownload(params.accession)
  fastqfiles = sraresults.fastq_trimmed
  fastqfiles.view()
  assemblyres = assembly(params.kmerlen, params.covcutoff, params.mincontiglength, params.quastref, fastqfiles)
  //multiqc(params.outdir, assemblyres.quastout)
}
