process TRIM_BARCODE_AND_SPACER {
    tag { id }

    input:
    tuple val(lane), val(id), path(tar)

    output:
    tuple val(lane), val(id), path("${id}.tar"), emit: trimmed

    script:
    """
    barcode=\$(printf "%${params.barcode_length}s" | tr ' ' "N")
    barcode_spacer_R1="\${barcode}${params.spacer_seq_R1}"    
    barcode_spacer_R2="\${barcode}${params.spacer_seq_R2}"
    length_barcode_spacer_R1=\${#barcode_spacer_R1}
    length_barcode_spacer_R2=\${#barcode_spacer_R2}
    
    tar -x --use-compress-program=pigz -f ${tar}

    mv ${tar} input.tar
    mv ${id}_R1.fastq.gz input_R1.fastq.gz
    mv ${id}_R2.fastq.gz input_R2.fastq.gz

    cutadapt \
        -j ${task.cpus} \
        -u \${length_barcode_spacer_R1} \
        -U \${length_barcode_spacer_R2} \
        -l ${params.guide_length} \
        --minimum-length ${params.guide_length} \
        -o ${id}_R1.fastq.gz \
        -p ${id}_R2.fastq.gz \
        input_R1.fastq.gz input_R2.fastq.gz

    tar -c --use-compress-program=pigz -f ${id}.tar ${id}_R1.fastq.gz ${id}_R2.fastq.gz
    """
}