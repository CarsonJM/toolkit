process COBRAMETA {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cobra-meta:1.2.3--pyhdfd78af_0':
        'biocontainers/cobra-meta:1.2.3--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(coverage)
    tuple val(meta3), path(query)
    tuple val(meta4), path(bam)
    val assembler
    val mink
    val maxk

    output:
    tuple val(meta), path("${prefix}_COBRA_extended.fasta.gz")      , optional: true    , emit: fasta
    tuple val(meta), path("${prefix}_COBRA_joining_summary.txt")    , optional: true    , emit: joining_summary
    tuple val(meta), path("${prefix}.cobra.log")                                        , emit: log
    path "versions.yml"                                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    # reformat coverage/virus files to remove header
    tail ${query} -n +2 | awk '{print \$1}' > ${prefix}_queries.txt
    tail ${coverage} -n +2 > ${prefix}_cobra_coverage.txt

    cobra-meta \\
        --fasta ${fasta_name} \\
        --coverage ${prefix}_cobra_coverage.txt \\
        --query ${prefix}_queries.txt \\
        --mapping ${bam} \\
        --assembler ${assembler} \\
        --mink ${mink} \\
        --maxk ${maxk} \\
        --threads ${task.cpus} \\
        --output ${prefix} \\
        ${args}

    mv ${prefix}/log ${prefix}.cobra.log

    if [ -f ${prefix}/COBRA_category_ii-a_extended_circular_unique.fasta ]; then
        cat ${prefix}/COBRA_category_ii-a_extended_circular_unique.fasta ${prefix}/COBRA_category_ii-b_extended_partial_unique.fasta > ${prefix}_COBRA_extended.fasta
        gzip ${prefix}_COBRA_extended.fasta
        mv ${prefix}/COBRA_joining_summary.txt ${prefix}_COBRA_joining_summary.txt
    else
        echo "" | gzip > ${prefix}_COBRA_extended.fasta.gz
        touch ${prefix}_COBRA_joining_summary.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cobra: \$(echo \$(cobra-meta --version 2>&1) | sed 's/^.*cobra v//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_COBRA_extended.fasta.gz
    touch ${prefix}_COBRA_joining_summary.txt
    touch ${prefix}.cobra.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cobra: \$(echo \$(cobra-meta --version 2>&1) | sed 's/^.*cobra v//' ))
    END_VERSIONS
    """
}
