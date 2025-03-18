process TRIM_BARCODE_AND_SPACER {
    tag { id }

    input:
    tuple val(lane), val(id), path(fastq_files)

    output:
    tuple val(lane), val(id), path("output/${lane}*.fastq.gz"), emit: trimmed

    script:
    R1 = fastq_files[0]
    R2 = fastq_files[1]
    """
    barcode=\$(printf "%${params.barcode_length}s" | tr ' ' "N")
    barcode_spacer_R1="\${barcode}${params.spacer_seq_R1}"    
    barcode_spacer_R2="\${barcode}${params.spacer_seq_R2}"
    length_barcode_spacer_R1=\${#barcode_spacer_R1}
    length_barcode_spacer_R2=\${#barcode_spacer_R2}
    trimmed_length_R1=\$((length_barcode_spacer_R1 + 8))
    trimmed_length_R2=\$((length_barcode_spacer_R2 + 8))
    mkdir -p output

    cutadapt \
        -j 16 \
        -u \${trimmed_length_R1} \
        -U \${trimmed_length_R2} \
        -l ${params.guide_length} \
        --minimum-length ${params.guide_length} \
        -o output/${id}_R1.fastq.gz \
        -p output/${id}_R2.fastq.gz \
        ${R1} ${R2}
    """
}