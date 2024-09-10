// Import modules
include { CAT_FASTQ                             } from '../../../modules/nf-core/cat/fastq/main'
include { MEGAHIT_MEGAHIT as MEGAHIT_SINGLE     } from '../../../modules/nf-core/megahit/megahit/main'
include { MEGAHIT_MEGAHIT as MEGAHIT_COASSEMBLY } from '../../../modules/nf-core/megahit/megahit/main'
include { PLASS_PENGUIN as PENGUIN_SINGLE       } from '../../../modules/nf-core/plass/penguin/main'
include { PLASS_PENGUIN as PENGUIN_COASSEMBLY   } from '../../../modules/nf-core/plass/penguin/main'
include { SPADES as SPADES_SINGLE               } from '../../../modules/nf-core/spades/main'
include { SPADES as SPADES_COASSEMBLY           } from '../../../modules/nf-core/spades/main'

// Import subworkflows
include { getWorkDirs; rmEmptyFastAs    } from '../utils_nfmicrobe_functions/main'

// Run workflow
workflow FASTQ_READASSEMBLY_FASTA {

    take:
    fastq_gz                // channel: [ [ meta.id, meta.single_end, meta.group ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ] (MANDATORY)
    run_megahit_single      // boolean: false
    run_megahit_coassembly  // boolean: false
    run_spades_single       // boolean: false
    run_spades_coassembly   // boolean: false
    use_spades_scaffolds    // boolean: false
    run_penguin_single      // boolean: false
    run_penguin_coassembly  // boolean: false

    main:

    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()

    // create empty channels for combining assemblies
    ch_assemblies_prefilt_fasta_gz  = Channel.empty()
    ch_assembly_graph_gz            = Channel.empty()
    ch_spades_logs                  = Channel.empty()

    if (run_megahit_single) {
        //
        // MODULE: Assemble reads individually with MEGAHIT
        //
        MEGAHIT_SINGLE(
            fastq_gz.map { meta, fasta -> [ meta + [ assembler: 'megahit_single' ], fasta ] }
        ).contigs
        ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(MEGAHIT_SINGLE.out.contigs)
            .join(MEGAHIT_SINGLE.out.min_kmer)
            .join(MEGAHIT_SINGLE.out.max_kmer)
            .map { meta, fasta, min_kmer, max_kmer ->
                [ meta + [ mink: min_kmer, maxk: max_kmer ], fasta ]
            }
        ch_assembly_graph_gz            = ch_assembly_graph_gz.mix(MEGAHIT_SINGLE.out.gfa)
        ch_versions                     = ch_versions.mix(MEGAHIT_SINGLE.out.versions)
    }

    if (run_spades_single) {
        // prepare reads for metaspades input
        ch_spades_single_input = fastq_gz
            .map { meta, fastq ->
                [ meta + [ assembler: 'spades_single' ], fastq, [], [] ]
            }
            .branch { meta, fastq, extra1, extra2 ->
                single_end: meta.single_end
                paired_end: true
            }

        //
        // MODULE: Assemble reads individually with metaSPAdes
        //
        SPADES_SINGLE(
            ch_spades_single_input.paired_end,
            [],
            []
        )
        if (use_spades_scaffolds) {
            ch_spades_single_fasta_gz   = SPADES_SINGLE.out.scaffolds
                .join(SPADES_SINGLE.out.min_kmer)
                .join(SPADES_SINGLE.out.max_kmer)
                .map { meta, fasta, min_kmer, max_kmer ->
                    [ meta + [ mink: min_kmer, maxk: max_kmer ], fasta ]
                }
        } else {
            ch_spades_single_fasta_gz   = SPADES_SINGLE.out.contigs
                .join(SPADES_SINGLE.out.min_kmer)
                .join(SPADES_SINGLE.out.max_kmer)
                .map { meta, fasta, min_kmer, max_kmer ->
                    [ meta + [ mink: min_kmer, maxk: max_kmer ], fasta ]
                }
        }
        ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(ch_spades_single_fasta_gz)
        ch_assembly_graph_gz            = ch_assembly_graph_gz.mix(SPADES_SINGLE.out.gfa)
        ch_spades_logs                  = ch_spades_logs.mix(SPADES_SINGLE.out.log)
        ch_versions                     = ch_versions.mix(SPADES_SINGLE.out.versions)
    }

    if (run_penguin_single) {
        //
        // MODULE: Assemble reads individually with PenguiN
        //
        PENGUIN_SINGLE(
            fastq_gz.map { meta, fasta -> [ meta + [ assembler: 'penguin_single' ], fasta ] }
        )
        ch_penguin_single_fasta_gz      = PENGUIN_SINGLE.out.contigs

        ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(ch_penguin_single_fasta_gz)
        ch_versions                     = ch_versions.mix(PENGUIN_SINGLE.out.versions)
    }

    if (run_megahit_coassembly || run_spades_coassembly || run_penguin_coassembly) {
        // group and set group as new id
        ch_cat_coassembly_fastq_gz = fastq_gz
            .map { meta, reads -> [ meta.group, meta, reads ] }
            .groupTuple( by: 0, sort:'deep' )
            .map { group, meta, reads ->
                def meta_new                = [:]
                meta_new.id                 = "group-$group"
                meta_new.group              = group
                meta_new.single_end         = meta.single_end[0]
                if ( meta_new.single_end ) {
                    return [ meta_new, reads.collect { fastqs -> fastqs } ]
                } else {
                    return [ meta_new, reads.flatten() ]
                }
            }
            .branch {
                meta, reads ->
                    coassembly: meta.single_end && reads.size() >= 2 || !meta.single_end && reads.size() >= 4
                    skip_coassembly: true   // Can skip coassembly if there is not multiple samples
            }

        //
        // MODULE: Combine reads within groups for coassembly
        //
        CAT_FASTQ(
            ch_cat_coassembly_fastq_gz.coassembly
        )
        ch_coassembly_fastq_gz  = CAT_FASTQ.out.reads
        ch_versions             = ch_versions.mix(CAT_FASTQ.out.versions)
    }

    if (run_megahit_coassembly) {
        //
        // MODULE: Co-assemble reads with MEGAHIT
        //
        MEGAHIT_COASSEMBLY(
            ch_coassembly_fastq_gz.map { meta, fastq -> [ meta + [ assembler: 'megahit_coassembly' ], fastq ] }
        )
        ch_versions                     = ch_versions.mix(MEGAHIT_COASSEMBLY.out.versions)
        ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(MEGAHIT_COASSEMBLY.out.contigs)
            .join(MEGAHIT_COASSEMBLY.out.min_kmer)
            .join(MEGAHIT_COASSEMBLY.out.max_kmer)
            .map { meta, fasta, min_kmer, max_kmer ->
                [ meta + [ mink: min_kmer, maxk: max_kmer ], fasta ]
            }
        ch_assembly_graph_gz            = ch_assembly_graph_gz.mix(MEGAHIT_COASSEMBLY.out.gfa)

    }

    if (run_spades_coassembly) {
        // prepare reads for metaspades input
        ch_metaspades_coassembly_input = ch_coassembly_fastq_gz
            .map { meta, fastq ->
                [ meta + [ assembler: 'spades_coassembly' ], fastq, [], [] ]
            }
            .branch { meta, fastq, extra1, extra2 ->
                single_end: meta.single_end
                paired_end: true
            }

        //
        // MODULE: Co-assemble reads with metaSPAdes
        //
        SPADES_COASSEMBLY(
            ch_metaspades_coassembly_input.paired_end,
            [],
            []
        )
        if (use_spades_scaffolds) {
            ch_spades_co_fasta_gz   = SPADES_COASSEMBLY.out.scaffolds
                .join(SPADES_COASSEMBLY.out.min_kmer)
                .join(SPADES_COASSEMBLY.out.max_kmer)
                .map { meta, fasta, min_kmer, max_kmer ->
                    [ meta + [ mink: min_kmer, maxk: max_kmer ], fasta ]
                }
        } else {
            ch_spades_co_fasta_gz   = SPADES_COASSEMBLY.out.contigs
                .join(SPADES_COASSEMBLY.out.min_kmer)
                .join(SPADES_COASSEMBLY.out.max_kmer)
                .map { meta, fasta, min_kmer, max_kmer ->
                    [ meta + [ mink: min_kmer, maxk: max_kmer ], fasta ]
                }
        }
        ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(ch_spades_co_fasta_gz)
        ch_assembly_graph_gz            = ch_assembly_graph_gz.mix(SPADES_COASSEMBLY.out.gfa)
        ch_spades_logs                  = ch_spades_logs.mix(SPADES_COASSEMBLY.out.log)
        ch_versions                     = ch_versions.mix(SPADES_COASSEMBLY.out.versions)
    }

    if (run_penguin_coassembly) {
        //
        // MODULE: Co-assemble reads with PenguiN
        //
        PENGUIN_COASSEMBLY(
            ch_cat_coassembly_fastq_gz.coassembly.map { meta, fastq -> [ meta + [ assembler: 'penguin_coassembly' ], fastq ] }
        )
        ch_penguin_co_fasta_gz          = PENGUIN_COASSEMBLY.out.contigs
        ch_versions                     = ch_versions.mix(PENGUIN_COASSEMBLY.out.versions)
        ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(ch_penguin_co_fasta_gz)
    }

    // REMOVE EMPTY FASTA FILES FROM CHANNEL
    ch_assemblies_fasta_gz = rmEmptyFastAs(ch_assemblies_prefilt_fasta_gz, false)

    emit:
    assemblies_fasta_gz = ch_assemblies_fasta_gz    // channel: [ [ meta.id, meta.single_end, meta.group, meta.assembler, meta.mink, meta.maxk ], assembly.fasta.gz ]
    assembly_graph_gz   = ch_assembly_graph_gz      // channel: [ [ meta.id, meta.single_end, meta.group, meta.assembler ], assembly_graph.gfa.gz ]
    spades_logs         = ch_spades_logs            // channel: [ [ meta.id, meta.single_end, meta.group, meta.assembler ], spades.log ]
    multiqc_files       = ch_multiqc_files          // channel: /path.to/multiqc_files
    versions            = ch_versions               // channel: [ path(versions.yml) ]
}
