process ALIGN {
    tag { "${id}" }

    input:
    tuple val(lane), val(id), path(fastq_files), path(index)

    output:
    tuple val(lane), val(id), path("${id}.sam"), emit: alignedFiles
    path("${id}.log"), emit: alignResults
    path("${id}.sam.stats"), emit: alignStats
    path("${id}.sam.flagstat"), emit: alignFlagstats

    script:
    R1 = fastq_files[0]
    R2 = fastq_files[1]
    """
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
        -1 ${R1} \
        -2 ${R2} 2> ${id}.log > ${id}.sam

    samtools stats ${id}.sam > ${id}.sam.stats
    samtools flagstats ${id}.sam > ${id}.sam.flagstat

    """
}