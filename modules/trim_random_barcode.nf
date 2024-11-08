process TRIM_RANDOM_BARCODE {
    tag { lane }

    input:
    tuple val(lane), path(tar)

    output:
    tuple val(lane), path("${lane}.tar"), emit: trimmed

    script:
    """
    barcode=\$(printf "%${params.barcode_length}s" | tr ' ' "N")
    barcode_spacer_R1="\${barcode}${params.spacer_seq_R1}"    
    barcode_spacer_R2="\${barcode}${params.spacer_seq_R2}"
    length_barcode_spacer_R1=\${#barcode_spacer_R1}
    length_barcode_spacer_R2=\${#barcode_spacer_R2}

    tar -x --use-compress-program=pigz -f ${tar}
    
    mv ${lane}_R1.fastq.gz input_R1.fastq.gz
    mv ${lane}_R2.fastq.gz input_R2.fastq.gz
    
    cutadapt \
        -O \${length_barcode_spacer_R1} \
        -g \${barcode_spacer_R1} \
        --action=retain \
        -j ${task.cpus} \
        -o ${lane}_R1.fastq.gz \
        input_R1.fastq.gz 
   
    cutadapt \
        -O \${length_barcode_spacer_R2} \
        -g \${barcode_spacer_R2} \
        --action=retain \
        -j ${task.cpus} \
        -o ${lane}_R2.fastq.gz \
        input_R2.fastq.gz 

    tar -c --use-compress-program=pigz -f ${lane}.tar ${lane}_R1.fastq.gz ${lane}_R2.fastq.gz
    """
}