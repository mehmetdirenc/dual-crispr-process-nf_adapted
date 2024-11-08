process FASTQC {
    tag { id }

    input:
    tuple val(id), path(file)

    output:
    path "*_fastqc.{zip,html}", emit: fastqcResults

    script:
    """
    tar -x --use-compress-program=pigz -f ${file}
    fastqc -t ${task.cpus} -q *.fastq.gz
    """
}