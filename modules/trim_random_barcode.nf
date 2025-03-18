process TRIM_RANDOM_BARCODE {
    tag { lane }
    cpus 2
    input:
    tuple val(lane), path(fastq_files)

    output:
    tuple val(lane), path("output/${lane}*.fastq.gz"), emit: trimmed

    script:
    R1 = fastq_files[0]
    R2 = fastq_files[1]
    """
    echo added_random_barcode
    barcode=\$(printf "%${params.barcode_length}s" | tr ' ' "N")
    barcode_spacer_R1="\${barcode}${params.spacer_seq_R1}"    
    barcode_spacer_R2="\${barcode}${params.spacer_seq_R2}"
    length_barcode_spacer_R1=\${#barcode_spacer_R1}
    length_barcode_spacer_R2=\${#barcode_spacer_R2}
    trimmed_length_R1=\$((length_barcode_spacer_R1 + 8))
    trimmed_length_R2=\$((length_barcode_spacer_R2 + 8))

    mkdir -p output
    
    cutadapt \
        -O \${trimmed_length_R1} \
        -g \${barcode_spacer_R1} \
        --action=retain \
        -j 2 \
        -o output/${lane}_R1.fastq.gz \
        ${lane}_R1.fastq.gz 
   
    cutadapt \
        -O \${trimmed_length_R2} \
        -g \${barcode_spacer_R2} \
        --action=retain \
        -j 2 \
        -o output/${lane}_R2.fastq.gz \
        ${lane}_R2.fastq.gz

    """
}