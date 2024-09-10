// Import modules
include { PHABLES_INSTALL       } from '../../../modules/nf-core/phables/install/main'
include { PHABLES_RUN           } from '../../../modules/nf-core/phables/run/main'

workflow FASTQGFA_PHABLES_FASTA {
    take:
    fastq_gz            // [ [ meta.id ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ]  , read files (MANDATORY)
    gfa_gz              // [ [ meta.id ], gfa.gz ]                                  , assembly graph (MANDATORY)
    phables_config      // [ phables.config ]                                       , phables config (OPTIONAL)
    phables_db          // [ phables_db ]                                           , Phables database (OPTIONAL)

    main:
    ch_versions = Channel.empty()

    // if genomad_db exists, skip PHABLES_INSTALL
    if (phables_db) {
        ch_phables_db = phables_db
    } else {
        //
        // MODULE: download phables database
        //
        PHABLES_INSTALL(
        )
        ch_phables_db   = PHABLES_INSTALL.out.phables_db
        ch_versions     = ch_versions.mix(PHABLES_INSTALL.out.versions)
    }

    // combine fastq and gfa for phables
    ch_phables_input = fastq_gz
        .join(gfa_gz)
        .multiMap { meta, fastq, gfa ->
            fastq : [ meta, fastq ]
            gfa   : [ meta, gfa ]
        }

    //
    // MODULE: Run phables to extend assemblies
    //
    PHABLES_RUN(
        ch_phables_input.fastq,
        ch_phables_input.gfa,
        phables_config,
        ch_phables_db.first()
    )
    ch_phables_fasta_gz = PHABLES_RUN.out.fasta
    ch_versions = ch_versions.mix(PHABLES_RUN.out.versions)

    emit:
    fasta_gz    = ch_phables_fasta_gz   // [ [ meta ], extended_contigs.fna.gz ]    , FASTA file containing extended contigs
    versions    = ch_versions.unique()  // [ versions.yml ]
}
