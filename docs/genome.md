
## Reference genome

The genome parameters required to run the pipeline are listed below:

| Parameter          | Description                                                                                                                                                                                                       |
|--------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `fasta`            | Path to multi-fasta file containing reference genome assembly.                                                                                                                                                    |
| `gtf`              | Path to [GTF](https://www.ensembl.org/info/website/upload/gff.html) file containing gene annotation which is typically available for download with the reference assembly.                                        |
| `bwa_index`        | Path to [BWA index](http://bio-bwa.sourceforge.net/bwa.shtml) for reference genome assembly. See **Indexing genome** section below.                                                                               |

The parameters can either be specified at the command-line when running the pipeline  

```bash
nextflow run main.nf --design design.csv --fasta <FASTA_FILE> --gtf <GTF_FILE> --bwa_index <BWA_INDEX> -profile babs_modules
```

or you can edit the [genomes.config](https://github.com/crickbabs/BABS-MNASeqPE/blob/master/conf/genomes.config) file to define and store these parameters for multiple genome assemblies. Using this method you will only need to provide the specified shorthand name for the reference genome when running the pipeline.

```bash
nextflow run main.nf --design design.csv --genome hg19 -profile babs_modules
```

### Indexing genome

The fasta file will need to be indexed with SAMtools and BWA before running the pipeline.

```bash
samtools faidx <FASTA_FILE>
bwa index <FASTA_FILE>
```

See [BWA documentation](http://bio-bwa.sourceforge.net/bwa.shtml) for more information.
