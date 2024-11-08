process ALIGN {
    tag { "${id}" }

    input:
    tuple val(lane), val(id), path(tar), path(index)

    output:
    tuple val(lane), val(id), path("${id}.sam"), emit: alignedFiles
    path("${id}.log"), emit: alignResults
    path("${id}.sam.stats"), emit: alignStats
    path("${id}.sam.flagstat"), emit: alignFlagstats

    script:
    """
    tar -x --use-compress-program=pigz -f ${tar}

    bowtie2 \
        --threads \$((${task.cpus})) \
        -x ${index}/index \
        --ignore-quals \
        -L 21 \
        -N 0 \
        --very-sensitive \
        --score-min L,-36,0 \
        --np 25 \
        --rdg 25,25 \
        --rfg 25,25 \
        --no-overlap \
        --no-contain \
        --fr \
        -1 ${id}_R1.fastq.gz \
        -2 ${id}_R2.fastq.gz 2> ${id}.log > ${id}.sam

    samtools stats ${id}.sam > ${id}.sam.stats
    samtools flagstats ${id}.sam > ${id}.sam.flagstat

    """
}