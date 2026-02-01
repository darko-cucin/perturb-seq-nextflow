import groovy.json.JsonSlurper

include { paramsSummaryMap                                            } from 'plugin/nf-schema'
include { samplesheetToList                                           } from 'plugin/nf-schema'
include { STAR_GENOMEGENERATE as    STAR_GENOMEGENERATE_GEX           } from '../modules/star/genome_generate'
include { STAR_GENOMEGENERATE as    STAR_GENOMEGENERATE_CRISPR        } from '../modules/star/genome_generate'
include { GUNZIP as                 GUNZIP_GTF_CRISPR                 } from '../modules/gunzip'
include { GUNZIP as                 GUNZIP_GTF_GEX                    } from '../modules/gunzip'
include { GUNZIP as                 GUNZIP_FASTA_CRISPR               } from '../modules/gunzip'
include { GUNZIP as                 GUNZIP_FASTA_GEX                  } from '../modules/gunzip'
include { FASTQC                                                      } from '../modules/fastqc'
include { STARSOLO as               STARSOLO_CRISPR                   }   from '../modules/star/solo_align'      
include { STARSOLO as               STARSOLO_GEX                      }   from '../modules/star/solo_align'

workflow PERTURBSEQ{

  main:
  reads_ch =
    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { meta, r1, r2 ->
            tuple(
                tuple(meta.id, meta.lyb_type),   // grouping key
                meta,
                [r1, r2]
            )
        }
        .groupTuple()
        .map { key, metas, fastqs ->

            assert metas.unique().size() == 1 :
                "Conflicting meta for ${key}"

            def meta = metas[0]

            def all_fastqs = fastqs.flatten()
            tuple(meta, all_fastqs)
        }
    
    // separate CRISPR from GEX reads
    separated = reads_ch
    .branch {
        CRISPR:  it[0].lyb_type == 'CRISPR' 
        GEX:     it[0].lyb_type == 'GEX' 
    }
  
  // load protocol config with potential whitelist
  protocol_config = protocols(params.tenx_version)
  
  if (params.barcode_whitelist) {
        ch_barcode_whitelist = file(params.barcode_whitelist, checkIfExists: true)
    } else if (protocol_config.containsKey("whitelist")) {
        ch_barcode_whitelist = file("$projectDir/${protocol_config['whitelist']}", checkIfExists: true)
    } else {
        ch_barcode_whitelist = []
    }

 FASTQC (
    reads_ch
  )
 
  if (params.index_reference_crispr) {

    crispr_index = Channel.value(
        file(params.index_reference_crispr, checkIfExists: true)
    )

  // code block that will be executed if crispr_index is not provided. In this case pipeline need to generate index
  } else {

    genomeFastaWarn('crispr')

    ch_gtf_crispr =
        params.gtf_reference_crispr.endsWith('.gz')
            ? GUNZIP_GTF_CRISPR([[:], file(params.gtf_reference_crispr, checkIfExists: true)]).gunzip.map { it[1] }
            : Channel.value(file(params.gtf_reference_crispr, checkIfExists: true))

    ch_fasta_crispr =
        params.fasta_reference_crispr.endsWith('.gz')
            ? GUNZIP_FASTA_CRISPR([[:], file(params.fasta_reference_crispr, checkIfExists: true)]).gunzip.map { it[1] }
            : Channel.value(file(params.fasta_reference_crispr, checkIfExists: true))

   STAR_GENOMEGENERATE_CRISPR(
        ch_fasta_crispr,
        ch_gtf_crispr
    )

    crispr_index = STAR_GENOMEGENERATE_CRISPR.out.index

  }
 
  if (params.index_reference_gex) {

    gex_index = Channel.value(
        file(params.index_reference_gex, checkIfExists: true)
    )

  // code block that will be executed if crispr_gex is not provided. In this case pipeline need to generate index
  } else {

    genomeFastaWarn('gex')

    ch_gtf_gex =
        params.gtf_reference_gex.endsWith('.gz')
            ? GUNZIP_GTF_GEX([[:], file(params.gtf_reference_gex, checkIfExists: true)]).gunzip.map { it[1] }
            : Channel.value(file(params.gtf_reference_gex, checkIfExists: true))

    ch_fasta_gex =
        params.fasta_reference_gex.endsWith('.gz')
            ? GUNZIP_FASTA_GEX([[:], file(params.fasta_reference_gex, checkIfExists: true)]).gunzip.map { it[1] }
            : Channel.value(file(params.fasta_reference_gex, checkIfExists: true))
    
    
    STAR_GENOMEGENERATE_GEX (
         ch_fasta_gex,
         ch_gtf_gex
          )
  
    gex_index = STAR_GENOMEGENERATE_GEX.out.index
    
}

// Running starsolo for different types of libraries 
  STARSOLO_CRISPR (
    separated.CRISPR,
    crispr_index,
    ch_barcode_whitelist,
  )

  STARSOLO_GEX (
    separated.GEX,
    gex_index,
    ch_barcode_whitelist,
)
    }  

def genomeFastaWarn(lib_type) {
    log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
        "  '--fasta_$lib_type' parameter has been provided.\n" +
        "  Make sure transcript names in this file match those in the GFF/GTF file.\n\n" +
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
}

    // Retrieve the aligner-specific protocol based on the specified protocol.
    // Returns a map ["protocol": protocol, "extra_args": <extra args>, "whitelist": <path to whitelist>]
    // extra_args and whitelist are optional.
    def protocols (protocol) {
        def jsonSlurper = new JsonSlurper()
        def json = new File("${workflow.projectDir}/assets/protocols.json").text
        def protocols = jsonSlurper.parseText(json)

        if(protocol) {
            return protocols['starsolo'][protocol]
        } else {
            log.warn("Protocol '${protocol}' not recognized by the pipeline. Passing on the protocol to the aligner unmodified.")
            return ["protocol": protocol]
        }
    }
