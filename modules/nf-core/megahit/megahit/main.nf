process MEGAHIT_MEGAHIT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/megahit_pigz:657d77006ae5f222' :
        'community.wave.seqera.io/library/megahit_pigz:87a590163e594224' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.contigs.fa.gz")                            , emit: contigs
    tuple val(meta), path("intermediate_contigs/k*.contigs.fa.gz")      , emit: k_contigs
    tuple val(meta), path("intermediate_contigs/k*.addi.fa.gz")         , emit: addi_contigs
    tuple val(meta), path("intermediate_contigs/k*.local.fa.gz")        , emit: local_contigs
    tuple val(meta), path("intermediate_contigs/k*.final.contigs.fa.gz"), emit: kfinal_contigs
    tuple val(meta), path("${prefix}.megahit.log")                      , emit: log
    tuple val(meta), path("${prefix}.gfa.gz")                           , emit: gfa
    tuple val(meta), env(min_kmer)                                      , emit: min_kmer
    tuple val(meta), env(max_kmer)                                      , emit: max_kmer
    path "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def reads_command = meta.single_end ? "-r ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    megahit \\
        ${reads_command} \\
        ${args} \\
        -t ${task.cpus} \\
        --out-prefix ${prefix}

    # identify kmer size used to create assemblies
    kmer_size=\$(grep "^>" megahit_out/*.fa | sed 's/>k//; s/_.*//' | head -n 1)

    # convert contigs to FastG format
    megahit_toolkit \\
        contig2fastg \\
        \$kmer_size \\
        megahit_out/*.fa > ${prefix}.graph.fastg

    # convert FastG to GFA format
    fastg2gfa \\
        ${prefix}.graph.fastg > ${prefix}.gfa

    # gzip megahit contigs
    pigz \\
        --no-name \\
        -p ${task.cpus} \\
        ${args2} \\
        megahit_out/*.fa \\
        megahit_out/intermediate_contigs/*.fa

    # move output files to top level directory
    mv megahit_out/* .
    mv ./*log ${prefix}.megahit.log

    # gzip GFA files
    gzip ${prefix}.gfa ${prefix}.graph.fastg

    # identify megahit min/max kmer size
    kmer_string=\$(grep "k list: " ${prefix}.megahit.log | sed 's/.*k list: //; s/ .*//')
    kmer_array=(\${kmer_string//,/ })
    min_kmer=\$(IFS=\$'\\n'; echo "\${kmer_array[*]}" | sort -nr | tail -n 1)
    max_kmer=\$(IFS=\$'\\n'; echo "\${kmer_array[*]}" | sort -nr | head -n 1)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def reads_command = meta.single_end ? "-r ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    mkdir -p intermediate_contigs
    echo "" | gzip > ${prefix}.contigs.fa.gz
    echo "" | gzip > intermediate_contigs/k21.contigs.fa.gz
    echo "" | gzip > intermediate_contigs/k21.addi.fa.gz
    echo "" | gzip > intermediate_contigs/k21.local.fa.gz
    echo "" | gzip > intermediate_contigs/k21.final.contigs.fa.gz
    touch ${prefix}.megahit.log
    echo "" | gzip > ${prefix}.graph.fastg.gz
    echo "" | gzip > ${prefix}.gfa.gz
    min_kmer=21
    max_kmer=141

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
    END_VERSIONS
    """
}
