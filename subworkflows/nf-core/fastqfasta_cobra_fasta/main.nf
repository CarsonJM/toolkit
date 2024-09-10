// Import modules
include { COVERM_CONTIG } from '../../../modules/nf-core/coverm/contig/main'
include { COBRAMETA     } from '../../../modules/nf-core/cobrameta/main'

workflow FASTQFASTA_COBRA_FASTA {

    take:
    fastq_gz            // channel: [ [ meta.id, meta.single_end ], [ path(fastq_1), path(fastq_2) ] ]
    fasta_gz            // channel: [ [ meta.id, meta.single_end, meta.assembler, meta.mink, meta.maxk ], path(fasta)]
    query_contigs_tsv   // channel: [ [ meta.id, meta.single_end, meta.assembler, meta.mink, meta.maxk ], path(queries) ]

    main:

    ch_versions = Channel.empty()

    // join fastq and fasta by meta.id
    ch_coverm_input = fastq_gz.map { meta, fastq -> [ meta.id, fastq ] }
        .join(fasta_gz.map { meta, fasta -> [ meta.id, meta, fasta ] })
        .multiMap { meta_id, fastq, meta_fasta, fasta ->
            fastq: [ meta_fasta, fastq ]
            fasta: [ meta_fasta, fasta ]
        }

    //
    // MODULE: Align reads to their corresponding assembly
    //
    COVERM_CONTIG(
        ch_coverm_input.fastq,
        ch_coverm_input.fasta
    )
    ch_versions = ch_versions.mix(COVERM_CONTIG.out.versions.first())

    // prepare input for cobra
    ch_cobra_input = fasta_gz.map { meta, fasta -> [ meta.id, meta, fasta ] }
        .join(query_contigs_tsv.map { meta, query -> [ meta.id, query ] })
        .join(COVERM_CONTIG.out.tsv.map { meta, coverage -> [ meta.id, coverage ] })
        .join(COVERM_CONTIG.out.bam.map { meta, bam -> [ meta.id, bam ] })
        .multiMap { meta_id, meta_fasta, fasta, query, coverage, bam ->
            fasta:      [ meta_fasta, fasta ]
            coverage:   [ meta_fasta, coverage ]
            query:      [ meta_fasta, query ]
            bam:        [ meta_fasta, bam ]
        }

    //
    // MODULE: Extend query contigs
    //
    COBRAMETA(
        ch_cobra_input.fasta,
        ch_cobra_input.coverage,
        ch_cobra_input.query,
        ch_cobra_input.bam,
        ch_cobra_input.fasta.map { meta, fasta -> meta.assembler },
        ch_cobra_input.fasta.map { meta, fasta -> meta.mink },
        ch_cobra_input.fasta.map { meta, fasta -> meta.maxk },
    )


    emit:
    fasta_gz    = COBRAMETA.out.fasta   // channel: [ val(meta), [ fasta.gz ] ]
    versions    = ch_versions           // channel: [ versions.yml ]
}
