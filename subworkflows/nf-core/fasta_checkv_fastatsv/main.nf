// Import modules
include { CHECKV_DOWNLOADDATABASE   } from '../../../modules/nf-core/checkv/downloaddatabase/main'
include { CHECKV_ENDTOEND           } from '../../../modules/nf-core/checkv/endtoend/main'

workflow FASTA_CHECKV_FASTATSV {
    take:
    virus_fasta_gz  // channel: [ [ meta ], fasta.gz ]
    checkv_db       // path: checkv_db

    main:
    ch_versions = Channel.empty()

    // if checkv DB exists, skip CHECKV_DOWNLOADDATABASE
    if (checkv_db){
        ch_checkv_db = checkv_db
    } else {
        //
        // MODULE: download standard CheckV database
        //
        CHECKV_DOWNLOADDATABASE(
        )
        ch_checkv_db    = CHECKV_DOWNLOADDATABASE.out.checkv_db
        ch_versions     = ch_versions.mix(CHECKV_DOWNLOADDATABASE.out.versions)
    }

    //
    // MODULE: Run CheckV end to end
    //
    CHECKV_ENDTOEND(
        virus_fasta_gz,
        ch_checkv_db.first()
    )
    ch_versions = ch_versions.mix(CHECKV_ENDTOEND.out.versions)

    emit:
    viruses_fasta_gz    = CHECKV_ENDTOEND.out.viruses           // channel: [ [ meta ], viruses.fna.gz ]
    proviruses_fasta_gz = CHECKV_ENDTOEND.out.proviruses        // channel: [ [ meta ], proviruses.fna.gz ]
    quality_summary_tsv = CHECKV_ENDTOEND.out.quality_summary   // channel: [ [ meta ], quality_summary.tsv ]
    contamination_tsv   = CHECKV_ENDTOEND.out.contamination     // channel: [ [ meta ], contamination.tsv ]
    completeness_tsv    = CHECKV_ENDTOEND.out.completeness      // channel: [ [ meta ], completeness.tsv ]
    versions            = ch_versions                           // path: [ versions.yml ]
}
