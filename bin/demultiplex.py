from Bio import SeqIO
import gzip
import re



# Read forward primers from a file
def read_primers(primer_seq_file, sample_primer_file, new_barcodes_file):
    primer_dict = {}
    sample_to_barcode_dict = {}
    with open(primer_seq_file, "r") as f:
        # Skip the header line
        next(f)
        # Read each line and extract the forward primer sequence
        for line in f:
            split_line = line.strip().split("\t")
            primer_dict[split_line[0].strip()] = split_line[1].strip()
    header = "lane	sample_name	barcode_R1	barcode_R2\n"
    with open(new_barcodes_file, "w") as b:
        b.write(header)
        with open(sample_primer_file, "r") as f:
            # Skip the header line
            next(f)
            # Read each line and extract the forward primer sequence
            for line in f:
                split_line = line.strip().split("\t")
                sample = split_line[0]
                primer_name = split_line[1]
                print([primer_dict.keys()][0], primer_name)
                sequence = primer_dict[primer_name]
                match = re.search(r'[a-z]+([A-Z]{4})', sequence)
                barcode = match.group(1) if match else None
                # barcode = split_line[1]
                b.write("%s\t%s\t%s\n" % ("lane", sample, barcode))
    return primers

# Extract sample-specific barcodes from forward primers
def extract_barcodes(primers):
    sample_barcodes = {}
    for i, primer in enumerate(primers):
        # Find the position of the constant region
        constant_region = "ACACGACGCTCTTCCGATCT"
        start = primer.find(constant_region) + len(constant_region)
        # Extract the 4 nucleotides after the staggering nucleotides
        sample_barcode = primer[start+5:start+9]  # Skip the staggering nucleotide(s)
        # Map the barcode to a sample name
        sample_name = f"Sample{i+1}"  # Sample1, Sample2, etc.
        sample_barcodes[sample_barcode] = sample_name
    return sample_barcodes





if __name__ == '__main__':
    # Path to the primers file
    primers_file = "/home/direnc/inputs/kirsten/crispr/primers.tsv"
    sample_to_primers_file = "/home/direnc/inputs/kirsten/crispr/sample_to_primer.tsv"
    new_barcodes_file = "/home/direnc/PycharmProjects/dual-crispr-process-nf/test_data/barcodes_kirsten.txt"
    # Read primers and create barcode mapping
    primers = read_primers(primers_file, sample_to_primers_file, new_barcodes_file)
    sample_barcodes = extract_barcodes(primers)

    # Print the barcode mapping
    print("Sample Barcodes:")
    for barcode, sample in sample_barcodes.items():
        print(f"{sample}: {barcode}")

    # Open output files
    output_files = {sample: open(f"{sample}_R1.fastq", "w") for sample in sample_barcodes.values()}

    # Process forward reads (R1)
    with gzip.open("Undetermined_S0_R1_001.fastq.gz", "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            read_seq = str(record.seq)
            # Check if the constant region is present
            if "ACACGACGCTCTTCCGATCT" in read_seq:
                start = read_seq.find("ACACGACGCTCTTCCGATCT") + len("ACACGACGCTCTTCCGATCT")
                # Extract the 4 nucleotides after the staggering nucleotides
                sample_barcode = read_seq[start+5:start+9]  # Skip the staggering nucleotide(s)
                # Assign to sample
                if sample_barcode in sample_barcodes:
                    sample_name = sample_barcodes[sample_barcode]
                    SeqIO.write(record, output_files[sample_name], "fastq")

    # Close output files
    for file in output_files.values():
        file.close()