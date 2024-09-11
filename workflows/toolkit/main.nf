/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Import plugins
include { paramsSummaryMap              } from 'plugin/nf-validation'

// Import modules
include { CAT_FASTQ                     } from '../../modules/nf-core/cat/fastq'
include { CUSTOM_CLEANWORKDIRS          } from '../../modules/nf-core/custom/cleanworkdirs'
include { FASTQC as FASTQC_RAW          } from '../../modules/nf-core/fastqc'
include { FASTQC as FASTQC_PREPROCESSED } from '../../modules/nf-core/fastqc'
include { FASTP                         } from '../../modules/nf-core/fastp'
include { MULTIQC                       } from '../../modules/nf-core/multiqc'

// Import functions from subworkflows
include { getGenomeAttribute            } from '../../subworkflows/local/utils_nfmicrobe_toolkit_pipeline'
include { getWorkDirs; rmEmptyFastQs    } from '../../subworkflows/nf-core/utils_nfmicrobe_functions'
include { methodsDescriptionText        } from '../../subworkflows/local/utils_nfmicrobe_toolkit_pipeline'

// Import subworkflows
include { FASTQ_BOWTIE2_FASTQ           } from '../../subworkflows/nf-core/fastq_bowtie2_fastq'
include { FASTQ_READASSEMBLY_FASTA      } from '../../subworkflows/nf-core/fastq_readassembly_fasta'
include { FASTQFASTA_COBRA_FASTA        } from '../../subworkflows/nf-core/fastqfasta_cobra_fasta'
include { FASTQGFA_PHABLES_FASTA        } from '../../subworkflows/nf-core/fastqgfa_phables_fasta'
include { paramsSummaryMultiqc          } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML        } from '../../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TOOLKIT {

    take:
    input_fastqs    // channel: fastq files read in from --input (samplesheet) and --fastqs
    input_fastas    // channel: fasta files read in from --input (samplesheet) and --fastas
    parameters      // pipeline parameters

    main:

    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_workdirs_to_clean    = Channel.empty()

    /*
    -------------------------------------------------
        RAW READ QUALITY ASSESSMENT
    -------------------------------------------------
    */
    //
    // MODULE: Run FastQC on raw reads
    //
    FASTQC_RAW (
        input_fastqs
    )
    ch_multiqc_files    = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{ zip -> zip[1] })
    ch_versions         = ch_versions.mix(FASTQC_RAW.out.versions.first())


    /*
    -------------------------------------------------
        READ PREPROCESSING
    -------------------------------------------------
    */
    ch_preprocessed_prefilt_fastq_gz    = Channel.empty()

    if (parameters.run_fastp) {
        //
        // MODULE: Run fastp on raw reads
        //
        FASTP(
            input_fastqs,
            [],
            false,
            false,
            false
        )
        ch_preprocessed_prefilt_fastq_gz    = ch_preprocessed_prefilt_fastq_gz.mix(FASTP.out.reads)
        ch_multiqc_files                    = ch_multiqc_files.mix(FASTP.out.json.collect{ json -> json[1] })
        ch_versions                         = ch_versions.mix(FASTP.out.versions)
    }

    // REMOVE EMPTY FASTQ FILES
    ch_preprocessed_fastq_gz = rmEmptyFastQs(ch_preprocessed_prefilt_fastq_gz, false)

    // IDENTIFY WORKDIRS TO CLEAN
    ch_pre_preprocessing_workdirs = getWorkDirs(
        input_fastqs,
        ch_preprocessed_fastq_gz,
        false
    )
    ch_workdirs_to_clean = ch_workdirs_to_clean.mix(
        ch_pre_preprocessing_workdirs.map { meta, dir -> [ meta, dir, 'RAW' ] }
    )

    // Use input fasta if no preprocessing is selected
    if (!parameters.run_fastp) {
        ch_preprocessed_fastq_gz = input_fastqs
    }

    /*
    -------------------------------------------------
        RUN MERGING
    -------------------------------------------------
    */
    ch_runmerged_fastq_gz   = Channel.empty()

    if (parameters.perform_run_merging) {
        // prepare reads for concatenating within runs
        ch_reads_forcat = ch_preprocessed_fastq_gz
            .map { meta, reads -> [ meta - meta.subMap('run'), reads ] }
            .groupTuple(sort: 'deep')
            .branch {
                meta, reads ->
                    cat:      reads.size() >= 2 // SE: [ [ meta ], [ S1_R1, S2_R1 ] ]; PE: [ [ meta ], [ [ S1_R1, S1_R2 ], [ S2_R1, S2_R2 ] ] ]
                    skip_cat: true              // Can skip merging if only single lanes
            }

        //
        // MODULE: Concatenate reads across runs, within samples
        //
        CAT_FASTQ(
            ch_reads_forcat.cat.map { meta, reads -> [ meta, reads.flatten() ] }
        )
        ch_runmerged_fastq_gz   = ch_runmerged_fastq_gz.mix(CAT_FASTQ.out.reads)
        ch_versions             = ch_versions.mix(CAT_FASTQ.out.versions)

        // Ensure we don't have nests of nests so that structure is in form expected for assembly
        ch_reads_forcat_skipped = ch_reads_forcat.skip_cat
            .map { meta, reads ->
                def new_reads = meta.single_end ? reads[0] : reads.flatten()
                [ meta, new_reads ]
            }

        // Combine single run and multi-run-merged data
        ch_runmerged_fastq_gz = ch_runmerged_fastq_gz.mix(ch_reads_forcat_skipped)
    }

    // IDENTIFY WORKDIRS TO CLEAN
    ch_pre_runmerge_workdirs = getWorkDirs(
        ch_preprocessed_fastq_gz.map { meta, reads -> [ meta - meta.subMap('run'), reads[0] ] },
        ch_runmerged_fastq_gz,
        false
    )
    ch_workdirs_to_clean = ch_workdirs_to_clean.mix(
        ch_pre_runmerge_workdirs.map { meta, dir -> [ meta, dir, 'PREPROCESSING' ] }
    )

    // Use preprocessing fasta if no run merging is selected
    if (!parameters.perform_run_merging) {
        ch_runmerged_fastq_gz = ch_preprocessed_fastq_gz
    }

    /*
    -------------------------------------------------
        HOST READ REMOVAL
    -------------------------------------------------
    */
    // Load igenomes fasta and bowtie2 index
    igenomes_fasta   = getGenomeAttribute('fasta')
    igenomes_index   = getGenomeAttribute('bowtie2')

    if (parameters.run_bowtie2_host_removal) {
        // create host fasta and index channels based on input parameters
        if (parameters.genome) {
            ch_host_fasta_gz    = file(igenomes_fasta, checkIfExists: true)
            ch_bowtie2_index    = file(igenomes_index, checkIfExists: true)
        } else if (parameters.host_bowtie2_index) {
            if (!parameters.host_fasta) {
                error "[nf-microbe/toolkit]: --run_bowtie2_host_removal = true but no iGenomes (--genome) or custom fasta (--host_fasta) provided"
            }
            ch_host_fasta_gz    = file(parameters.host_fasta, checkIfExists: true)
            ch_bowtie2_index    = file(parameters.host_bowtie2_index, checkIfExists: true)
        } else {
            if (!parameters.host_fasta) {
                error "[nf-microbe/toolkit]: --run_bowtie2_host_removal = true but no iGenomes (--genome) or custom fasta (--host_fasta) provided"
            }
            ch_host_fasta_gz    = file(parameters.host_fasta, checkIfExists: true)
            ch_bowtie2_index    = null
        }
    }

    ch_hostremoved_prefilt_fastq_gz = Channel.empty()

    if (parameters.run_bowtie2_host_removal) {
        //
        // SUBWORKFLOW: Remove host reads using Bowtie2
        //
        FASTQ_BOWTIE2_FASTQ(
            ch_runmerged_fastq_gz,
            ch_host_fasta_gz,
            ch_bowtie2_index
        )
        ch_hostremoved_prefilt_fastq_gz = ch_hostremoved_prefilt_fastq_gz.mix(FASTQ_BOWTIE2_FASTQ.out.fastq_gz)
        ch_multiqc_files                = ch_multiqc_files.mix(FASTQ_BOWTIE2_FASTQ.out.bt2_log.collect{ bt2_log -> bt2_log[1] })
        ch_versions                     = ch_versions.mix(FASTQ_BOWTIE2_FASTQ.out.versions)
    }

    // REMOVE EMPTY FASTQ FILES FROM CHANNEL
    ch_hostremoved_fastq_gz = rmEmptyFastQs(ch_hostremoved_prefilt_fastq_gz, false)

    // IDENTIFY WORKDIRS TO CLEAN
    ch_pre_hostremove_workdirs = getWorkDirs(
        ch_runmerged_fastq_gz.map { meta, reads -> [ meta, reads[0] ] },
        ch_hostremoved_prefilt_fastq_gz,
        false
    )
    ch_workdirs_to_clean = ch_workdirs_to_clean.mix(
        ch_pre_hostremove_workdirs.map { meta, dir -> [ meta, dir, 'RUNMERGE' ] }
    )

    // Use runmerged fasta if no host removal is selected
    if (!parameters.run_bowtie2_host_removal) {
        ch_hostremoved_fastq_gz = ch_runmerged_fastq_gz
    }

    /*
    -------------------------------------------------
        PREPROCESSED READ QUALITY ASSESSMENT
    -------------------------------------------------
    */
    if (parameters.run_fastp ||
        parameters.perform_run_merging ||
        parameters.run_bowtie2_host_removal) {
        //
        // MODULE: Run FastQC on preprocessed reads
        //
        FASTQC_PREPROCESSED (
            ch_hostremoved_fastq_gz
        )
        ch_multiqc_files    = ch_multiqc_files.mix(FASTQC_PREPROCESSED.out.zip.collect{ zip -> zip[1]})
        ch_versions         = ch_versions.mix(FASTQC_PREPROCESSED.out.versions.first())
    }

    /*
    -------------------------------------------------
        READ-BASED TAXONOMY
    -------------------------------------------------
    */

    /*
    -------------------------------------------------
        READ-BASED FUNCTION
    -------------------------------------------------
    */

    /*
    -------------------------------------------------
        READ ASSEMBLY
    -------------------------------------------------
    */
    //
    // SUBWORKFLOW: Read assembly
    //
    FASTQ_READASSEMBLY_FASTA(
        ch_preprocessed_fastq_gz,           // channel: [ [ meta.id, meta.single_end, meta.run, meta.group ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ]
        parameters.run_megahit_single,      // boolean: false
        parameters.run_megahit_coassembly,  // boolean: false
        parameters.run_spades_single,       // boolean: false
        parameters.run_spades_coassembly,   // boolean: false
        parameters.use_spades_scaffolds,    // boolean: false
        parameters.run_penguin_single,      // boolean: false
        parameters.run_penguin_coassembly,  // boolean: false
    )
    ch_assemblies_fasta_gz  = FASTQ_READASSEMBLY_FASTA.out.assemblies_fasta_gz
    ch_assembly_gfa_gz      = FASTQ_READASSEMBLY_FASTA.out.assembly_graph_gz
    ch_multiqc_files        = ch_multiqc_files.mix(FASTQ_READASSEMBLY_FASTA.out.multiqc_files)
    ch_versions             = FASTQ_READASSEMBLY_FASTA.out.versions

    /*
    -------------------------------------------------
        ASSEMBLY EXTENSION
    -------------------------------------------------
    */
    //
    // SUBWORKFLOW: Extend assemblies
    //
    FASTQFASTAGFA_ASSEMBLYEXTENSION_FASTA(
        ch_preprocessed_fastq_gz,   // channel: [ [ meta.id, meta.single_end, meta.run, meta.group ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ]
        parameters.run_cobra,       // boolean: false
        ch_assemblies_fasta_gz,     // channel: [ [ meta.id, meta.single_end, meta.assembler, meta.mink, meta.maxk ], [ fasta.gz ] ]
        ch_query_contigs_tsv,       // channel: [ [ meta.id, meta.single_end, meta.assembler, meta.mink, meta.maxk ], [ queies.tsv ] ]
        parameters.run_phables,     // boolean: false
        ch_assemblies_gfa_gz,       // channel: [ [ meta.id, meta.single_end, meta.assembler, meta.mink, meta.maxk ], [ gfa.gz ] ]
        phables_config,             // channel: path(phables_config.yml)
        phables_db                  // channel: path(phables_db)
    )
    ch_extended_fasta_gz    = FASTQFASTAGFA_ASSEMBLYEXTENSION_FASTA.out.fasta_gz
    ch_versions             = FASTQFASTAGFA_ASSEMBLYEXTENSION_FASTA.out.versions

    // combine assemblies, extended assemblies, and input fastas
    ch_combined_fasta_gz = input_fastas
        .mix(ch_assemblies_fasta_gz)
        .mix(ch_extended_fasta_gz)

    /*
    -------------------------------------------------
        ASSEMBLY QC/FILTERING
    -------------------------------------------------
    */

    /*
    -------------------------------------------------
        MGE IDENTIFICATION
    -------------------------------------------------
    */
    //
    // SUBWORKFLOW: Identify viruses/plasmids with geNomad
    //
    if (parameters.run_genomad || parameters.run_pyhmmer) {
        FASTA_MGEIDENTIFICATION_FASTATSV(
            ch_assembly_qc_fasta_gz,    // channel: [ [ meta.id, meta.single_end, meta.run, meta.group ], fasta.gz ]
            ch_assembly_qc_faa_gz,      // channel: [ [ meta.id, meta.single_end, meta.run, meta.group ], faa.gz ]
            parameter.run_genomad,      // boolean: false
            parameters.genomad_db,      // val: path(genomad_db)
            parameters.run_pyhmmer      // boolean: false
            parameters.use_mge_fastas   // boolean: false
        )
        ch_mge_fasta_gz         = FASTA_MGEIDENTIFICATION_FASTATSV.out.virus_fasta_gz
            .mix(FASTA_MGEIDENTIFICATION_FASTATSV.out.plasmid_fasta_gz)
            .mix(FASTA_MGEIDENTIFICATION_FASTATSV.out.pyhmmer_virus_fasta_gz)
            .mix(FASTA_MGEIDENTIFICATION_FASTATSV.out.pyhmmer_plasmid_fasta_gz)
        ch_mge_faa_gz           = FASTA_MGEIDENTIFICATION_FASTATSV.out.virus_faa_gz.mix(FASTA_MGEIDENTIFICATION_FASTATSV.out.plasmid_faa_gz)
        ch_genomad_scores_tsv   = FASTA_MGEIDENTIFICATION_FASTATSV.out.scores_tsv
        ch_genomad_taxonomy_tsv = FASTA_MGEIDENTIFICATION_FASTATSV.out.taxonomy_tsv
        ch_genomad_genes_tsv    = FASTA_MGEIDENTIFICATION_FASTATSV.out.genes_tsv
        ch_genomad_features_tsv = FASTA_MGEIDENTIFICATION_FASTATSV.out.features_tsv
        ch_pyhmmer_tsv          = FASTA_MGEIDENTIFICATION_FASTATSV.out.pyhmmer_tsv
    } else {
        ch_mge_fasta_gz         = ch_assembly_qc_fasta_gz
        ch_mge_faa_gz           = ch_assembly_qc_faa_gz
        ch_genomad_scores_tsv   = []
        ch_genomad_taxonomy_tsv = []
        ch_genomad_genes_tsv    = []
        ch_genomad_features_tsv = []
        ch_pyhmmer_tsv          = []
    }

    /*
    -------------------------------------------------
        VIRUS QC
    -------------------------------------------------
    */
    //
    // SUBWORKFLOW: Identify viruses/plasmids with geNomad
    //
    if (parameters.run_genomad || parameters.run_pyhmmer) {
        FASTA_MGEIDENTIFICATION_FASTATSV(
            ch_assembly_qc_fasta_gz,    // channel: [ [ meta.id, meta.single_end, meta.run, meta.group ], fasta.gz ]
            ch_assembly_qc_faa_gz,      // channel: [ [ meta.id, meta.single_end, meta.run, meta.group ], faa.gz ]
            parameter.run_genomad,      // boolean: false
            parameters.genomad_db,      // val: path(genomad_db)
            parameters.run_pyhmmer      // boolean: false
            parameters.use_mge_fastas   // boolean: false
        )
        ch_mge_fasta_gz         = FASTA_MGEIDENTIFICATION_FASTATSV.out.virus_fasta_gz
            .mix(FASTA_MGEIDENTIFICATION_FASTATSV.out.plasmid_fasta_gz)
            .mix(FASTA_MGEIDENTIFICATION_FASTATSV.out.pyhmmer_virus_fasta_gz)
            .mix(FASTA_MGEIDENTIFICATION_FASTATSV.out.pyhmmer_plasmid_fasta_gz)
        ch_mge_faa_gz           = FASTA_MGEIDENTIFICATION_FASTATSV.out.virus_faa_gz.mix(FASTA_MGEIDENTIFICATION_FASTATSV.out.plasmid_faa_gz)
        ch_genomad_scores_tsv   = FASTA_MGEIDENTIFICATION_FASTATSV.out.scores_tsv
        ch_genomad_taxonomy_tsv = FASTA_MGEIDENTIFICATION_FASTATSV.out.taxonomy_tsv
        ch_genomad_genes_tsv    = FASTA_MGEIDENTIFICATION_FASTATSV.out.genes_tsv
        ch_genomad_features_tsv = FASTA_MGEIDENTIFICATION_FASTATSV.out.features_tsv
        ch_pyhmmer_tsv          = FASTA_MGEIDENTIFICATION_FASTATSV.out.pyhmmer_tsv
    } else {
        ch_mge_fasta_gz         = ch_assembly_qc_fasta_gz
        ch_mge_faa_gz           = ch_assembly_qc_faa_gz
        ch_genomad_scores_tsv   = []
        ch_genomad_taxonomy_tsv = []
        ch_genomad_genes_tsv    = []
        ch_genomad_features_tsv = []
        ch_pyhmmer_tsv          = []
    }


    /*
    -------------------------------------------------
        VIRUS TAXONOMY
    -------------------------------------------------
    */

    /*
    -------------------------------------------------
        VIRUS HOST
    -------------------------------------------------
    */

    /*
    -------------------------------------------------
        VIRUS LIFESTYLE
    -------------------------------------------------
    */

    /*
    -------------------------------------------------
        VIRUS FUNCTION
    -------------------------------------------------
    */

    //
    // MODULE: Clean intermediate files
    //
    if (parameters.remove_intermediate_files) {
        //
        // MODULE: Clean up intermediate working directories
        //
        ch_workdirs_to_clean_unique = ch_workdirs_to_clean.unique()
        CUSTOM_CLEANWORKDIRS(ch_workdirs_to_clean_unique)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${parameters.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = parameters.multiqc_config ?
        Channel.fromPath(parameters.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = parameters.multiqc_logo ?
        Channel.fromPath(parameters.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = parameters.multiqc_methods_description ?
        file(parameters.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    // create a channel to contain all output files for testing
    ch_outputs = ch_preprocessed_fastq_gz
        .mix(ch_runmerged_fastq_gz)
        .mix(ch_hostremoved_fastq_gz)
        .mix(ch_assemblies_fasta_gz)

    emit:
    output_files    = ch_outputs
    multiqc_report  = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions        = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
