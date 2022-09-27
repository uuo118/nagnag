//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA            } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_GTF              } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_GFF              } from '../../modules/nf-core/modules/gunzip/main'

include { UNTAR as UNTAR_STAR_INDEX         } from '../../modules/nf-core/modules/untar/main'

include { STAR_GENOMEGENERATE               } from '../../modules/nf-core/modules/star/genomegenerate/main'

workflow PREPARE_GENOME {

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (params.fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], params.fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = file(params.fasta)
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ( [ [:], params.gtf ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = file(params.gtf)
        }
    }
    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = Channel.empty()
    if (params.star_index) {
        if (params.star_index.endsWith('.tar.gz')) {
            ch_star_index = UNTAR_STAR_INDEX ( [ [:], params.star_index ] ).untar.map { it[1] }
            ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
        } else {
            ch_star_index = file(params.star_index)
        }
    } else {
            ch_star_index = STAR_GENOMEGENERATE ( ch_fasta, ch_gtf ).index
            ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    emit:
        fasta            = ch_fasta            //    path: genome.fasta
        gtf              = ch_gtf              //    path: genome.gtf
        star_index       = ch_star_index       //    path: star/index/

        versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]

}
