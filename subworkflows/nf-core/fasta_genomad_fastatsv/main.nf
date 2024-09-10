// Import modules
include { GENOMAD_DOWNLOAD  } from '../../../modules/nf-core/genomad/download/main'
include { GENOMAD_ENDTOEND  } from '../../../modules/nf-core/genomad/endtoend/main'

workflow FASTA_GENOMAD_FASTATSV {
    take:
    fasta_gz    // [ [ meta ], fasta.gz ]   , assemblies/genomes (mandatory)
    genomad_db  // [ genomad_db ]           , genomad database directory (optional)

    main:
    ch_versions = Channel.empty()

    // if genomad_db exists, skip GENOMAD_DOWNLOAD
    if (genomad_db){
        ch_genomad_db = genomad_db
    } else {
        //
        // MODULE: download geNomad database
        //
        GENOMAD_DOWNLOAD(
        )
        ch_genomad_db   = GENOMAD_DOWNLOAD.out.genomad_db
        ch_versions     = ch_versions.mix(GENOMAD_DOWNLOAD.out.versions)
    }

    //
    // MODULE: Run geNomad end-to-end
    //
    GENOMAD_ENDTOEND(
        fasta_gz,
        ch_genomad_db
    )
    ch_versions = ch_versions.mix(GENOMAD_ENDTOEND.out.versions)

    emit:
    plasmid_fasta_gz    = GENOMAD_ENDTOEND.out.plasmid_fasta    // [ [ meta ], plasmid.fna.gz ]
    plasmid_faa_gz      = GENOMAD_ENDTOEND.out.plasmid_proteins // [ [ meta ], plasmid.faa.gz ]
    plasmid_tsv         = GENOMAD_ENDTOEND.out.plasmid_summary  // [ [ meta ], plasmid_summary.tsv ]
    virus_fasta_gz      = GENOMAD_ENDTOEND.out.virus_fasta      // [ [ meta ], virus.fna.gz ]
    virus_faa_gz        = GENOMAD_ENDTOEND.out.virus_proteins   // [ [ meta ], virus.faa.gz ]
    virus_tsv           = GENOMAD_ENDTOEND.out.virus_summary    // [ [ meta ], virus_summary.tsv ]
    protein_faa_gz      = GENOMAD_ENDTOEND.out.proteins         // [ [ meta ], proteins.faa.gz ]
    genes_tsv           = GENOMAD_ENDTOEND.out.genes            // [ [ meta ], genes.tsv ]
    features_tsv        = GENOMAD_ENDTOEND.out.features         // [ [ meta ], features.tsv ]
    scores_tsv          = GENOMAD_ENDTOEND.out.scores           // [ [ meta ], aggregated_classification.tsv ]
    taxonomy_tsv        = GENOMAD_ENDTOEND.out.taxonomy         // [ [ meta ], taxonomy.tsv ]
    versions            = ch_versions                           // channel: [ versions.yml ]
}
