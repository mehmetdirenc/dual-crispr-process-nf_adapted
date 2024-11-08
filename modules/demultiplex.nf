process DEMULTIPLEX {
    tag { lane }

    input:
    tuple val(lane), path(tar), path(barcode_R1), path(barcode_R2)

    output:
    tuple val(lane), path("*.tar"), emit: demuxed

    script:
    """
    tar -x --use-compress-program=pigz -f ${tar}
    rm ${tar}

    if [[ ${params.barcode_demux_location} == 'forward' ]]
    then
        cutadapt \
            -j ${task.cpus} \
            -e ${params.barcode_demux_mismatches} \
            --no-indels \
            -g file:${barcode_R1} \
            --action=none \
            -o "${lane}#{name}_R1.fastq.gz" \
            -p "${lane}#{name}_R2.fastq.gz" \
            ${lane}_R1.fastq.gz ${lane}_R2.fastq.gz
    fi

    if [[ ${params.barcode_demux_location} == 'reverse' ]]
    then
        cutadapt \
            -j ${task.cpus} \
            -e ${params.barcode_demux_mismatches} \
            --no-indels -g file:${barcode_R2} \
            --action=none \
            -o "${lane}#{name}_R2.fastq.gz" \
            -p "${lane}#{name}_R1.fastq.gz" \
            ${lane}_R2.fastq.gz ${lane}_R1.fastq.gz
    fi

    if [[ ${params.barcode_demux_location} == 'both' ]]
    then
        cutadapt \
            -j ${task.cpus} \
            -e ${params.barcode_demux_mismatches} \
            --no-indels \
            --pair-adapters \
            -g file:${barcode_R1} \
            -G file:${barcode_R2} \
            --action=none \
            -o "${lane}#{name}_R2.fastq.gz" \
            -p "${lane}#{name}_R1.fastq.gz" \
            ${lane}_R2.fastq.gz ${lane}_R1.fastq.gz
    fi

    for file in *_R1.fastq.gz;
    do
      filename="\${file%_R1.fastq.gz}"
      tar -c --use-compress-program=pigz -f \${filename}.tar \${filename}_R1.fastq.gz \${filename}_R2.fastq.gz
    done
    """
}