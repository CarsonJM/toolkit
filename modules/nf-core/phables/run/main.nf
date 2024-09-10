process PHABLES_RUN {
    tag "${meta.id}"
    label 'process_high'
    beforeScript "mkdir -p /tmp/phables_${meta.id}"
    containerOptions workflow.containerEngine == 'singularity' ?
        "--overlay /tmp/phables_${meta.id}":
        "--privileged"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://carsonjm/phables:1.4.0':
        'docker.io/carsonjm/phables:1.4.0' }"

    input:
    tuple val(meta) , path(fastq, stageAs: "reads/*")
    tuple val(meta2), path(graph)
    path phables_config
    path phables_db

    output:
    tuple val(meta), path("${prefix}.phables.fasta.gz")         , emit: fasta
    tuple val(meta), path("${prefix}.phables_complex_info.tsv") , emit: complex_info
    tuple val(meta), path("${prefix}.phables_genome_info.tsv")  , emit: genome_info
    tuple val(meta), path("${prefix}.phables.log")              , emit: log
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def is_compressed = graph.getExtension() == "gz" ? true : false
    def graph_name = is_compressed ? graph.getBaseName() : graph
    def phables_config_cp = phables_config ? "cp ${phables_config} ${meta.id}_config.yml" : ''
    def phables_config_cmd = phables_config ? "--configfile ${meta.id}_config.yml" : ''
    def conda_options = workflow.containerEngine == 'docker' ? "--no-use-conda" : ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # decompress graph if compressed
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${graph} > ${graph_name}
    fi

    mkdir phables_out

    # create a copy of config
    ${phables_config_cp}

    phables \\
        run \\
        ${phables_config_cmd} \\
        --output phables_out \\
        --input ${graph_name} \\
        --reads reads \\
        --threads ${task.cpus} \\
        --prefix ${prefix} \\
        ${conda_options} \\
        ${args} || if [ -f phables_out/phables/resolved_component_info.txt ]; then true; else false; fi

    # check if phables outputs were created
    if [ ! -f phables_out/phables/resolved_paths.fasta ]; then
        touch phables_out/phables/resolved_paths.fasta
    fi

    # gzip resolved genomes
    gzip -c phables_out/phables/resolved_paths.fasta > ${prefix}.phables.fasta.gz
    mv phables_out/phables/resolved_component_info.txt ${prefix}.phables_complex_info.tsv
    mv phables_out/phables/resolved_genome_info.txt ${prefix}.phables_genome_info.tsv
    mv phables_out/phables.log ${prefix}.phables.log

    # remove extra files
    # rm -rf phables_out/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phables: \$(echo \$(phables -v 2>&1) | sed -n 's/phables, version //' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def is_compressed = graph.getExtension() == "gz" ? true : false
    def graph_name = is_compressed ? graph.getBaseName() : graph
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.phables.fasta.gz
    touch ${prefix}.phables_complex_info.tsv
    touch ${prefix}.phables_genome_info.tsv
    touch ${prefix}.phables.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phables: \$(echo \$(phables -v 2>&1) | sed -n 's/phables, version //' ))
    END_VERSIONS
    """
}
