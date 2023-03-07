# dual-crispr-process-nf
Process dual-sgRNA CRISPR functional genetic screening data

## Installation

### Nextflow
Install `nextflow` following the instructions at https://www.nextflow.io/docs/latest/getstarted.html

### Singularity
Install `singularity` following the instructions at
https://singularity.lbl.gov/install-linux

### dual-crispr-process-nf pipeline
The most convenient way is to install `dual-crispr-process-nf` is to use `nextflow`'s built-in `pull` command
```bash
nextflow pull zuberlab/dual-crispr-process-nf
```

## Test (not yet implemented)
Before you start, make sure `nextflow` and `singularity` are properly installed on your system.

Clone git repository from Github and run the pipeline using the provided test data.
```bash
git clone https://github.com/ZuberLab/dual-crispr-process-nf.git
cd dual-crispr-process-nf
./test
```

## Documentation
```bash
nextflow run zuberlab/dual-crispr-process-nf --help
```

## Credits
Nextflow:  Paolo Di Tommaso - https://github.com/nextflow-io/nextflow

Singularity: Singularityware - https://singularity.lbl.gov
