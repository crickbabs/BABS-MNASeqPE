# PIPELINE & DOCS UNDER CONSTRUCTION. ETA OCTOBER 2018.

# ![BABS-MNASeqPE](https://raw.githubusercontent.com/crickbabs/BABS-MNASeqPE/master/docs/images/BABS-MNASeqPE_logo.png)

## Introduction

A [Nextflow](https://www.nextflow.io/) pipeline for processing paired-end Illumina MNASeq sequencing data.

The pipeline was written by [The Bioinformatics & Biostatistics Group](https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/) at [The Francis Crick Institute](https://www.crick.ac.uk/), London.

## Pipeline summary

1. Raw read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [`Fastq Screen`](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/))
2. Adapter trimming ([`cutadapt`](http://cutadapt.readthedocs.io/en/stable/installation.html))
3. Alignment ([`BWA`](https://sourceforge.net/projects/bio-bwa/files/))
4. Mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
5. Filtering to remove:
    * reads that are marked as duplicates ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    * reads that arent marked as primary alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    * reads that are unmapped ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    * reads that map to multiple locations ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    * reads containing > 3 mismatches in either read of the pair ([`BAMTools`](https://github.com/pezmaster31/bamtools))
    * reads that have a user-defined insert size ([`BAMTools`](https://github.com/pezmaster31/bamtools))
    * reads that are soft-clipped ([`BAMTools`](https://github.com/pezmaster31/bamtools))
    * reads that map to different chromosomes ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html))
    * reads that arent in FR orientation ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html))
    * reads where only one read of the pair fails the above criteria ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html))
6. Merge alignments at replicate-level ([`picard`](https://broadinstitute.github.io/picard/))
    * Re-mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
    * Remove duplicate reads (optional; [`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    * Create normalised bigWig files scaled to 1 million mapped read pairs ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`wigToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
7. Create IGV session file containing bigWig tracks for data visualisation ([`IGV`](https://software.broadinstitute.org/software/igv/)).
8. Collect and present QC at the raw read and alignment-level ([`MultiQC`](http://multiqc.info/))

## Documentation

The documentation for the pipeline can be found in the `docs/` directory:

1. [Installation](docs/install.md)
2. [Pipeline configuration](docs/config.md)
3. [Reference genome](docs/genome.md)
4. [Design file](docs/design.md)
5. [Running the pipeline](docs/usage.md)
6. [Output and interpretation of results](docs/output.md)
7. [Troubleshooting](docs/troubleshooting.md)

## Pipeline DAG

# ![BABS-MNASeqPE directed acyclic graph](https://raw.githubusercontent.com/crickbabs/BABS-MNASeqPE/master/docs/images/BABS-MNASeqPE_dag.png)

## Credits

The pipeline was written by the [The Bioinformatics & Biostatistics Group](https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/) at [The Francis Crick Institute](https://www.crick.ac.uk/), London.

The pipeline was developed by [Harshil Patel](mailto:harshil.patel@crick.ac.uk).

The [NGI-RNAseq](https://github.com/SciLifeLab/NGI-RNAseq) pipeline developed by Phil Ewels was used a template for this pipeline. Many thanks to Phil and the team at SciLifeLab. The help, tips and tricks provided by Paolo Di Tommaso were also invaluable. Thank you!

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
