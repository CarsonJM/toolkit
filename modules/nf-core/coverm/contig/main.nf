process COVERM_CONTIG {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coverm:0.7.0--h07ea13f_1':
        'biocontainers/coverm:0.7.0--h07ea13f_1' }"

    input:
    tuple val(meta) , path(fastq)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("${prefix}_alignment_results.tsv")    , emit: tsv
    tuple val(meta), path("${prefix}.bam")                      , emit: bam
    tuple val(meta), path("${prefix}.log")                      , emit: log
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def reads = !meta.single_end ? "--coupled ${fastq}" : "--single ${fastq}"
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    coverm \\
        contig \\
        ${reads} \\
        --reference ${fasta} \\
        --output-file ${prefix}_alignment_results.tsv \\
        --bam-file-cache-directory ${prefix}_bam_files \\
        --threads ${task.cpus} \\
        ${args} \\
        2> ${prefix}.log

    mv ${prefix}_bam_files/*.bam ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(echo \$(coverm --version 2>&1) | sed 's/^.*coverm //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def reads = !meta.single_end ? "--coupled ${fastq}" : "--single ${fastq}"
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_alignment_results.tsv
    touch ${prefix}.bam
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(echo \$(coverm --version 2>&1) | sed 's/^.*coverm //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
