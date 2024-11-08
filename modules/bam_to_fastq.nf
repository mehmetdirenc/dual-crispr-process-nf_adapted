process BAM_TO_FASTQ {
    tag { lane }

    input:
    tuple val(lane), path(bam)

    output:
    tuple val(lane), path("${lane}.tar"), emit: fastq

    script:
   """
    samtools fastq -@ ${task.cpus} ${bam} -1 ${lane}_R1.fastq.gz -2 ${lane}_R2.fastq.gz
    tar -c --use-compress-program=pigz -f ${lane}.tar ${lane}_R1.fastq.gz ${lane}_R2.fastq.gz
    rm ${lane}_R1.fastq.gz ${lane}_R2.fastq.gz
    """
}