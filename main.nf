#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
==============================================================================
    B A B S - M N A S E Q  P A I R E D - E N D   B E S T    P R A C T I C E
==============================================================================
 MNASeq Analysis Pipeline For Paired-End Illumina Samples.
 Started 29th August 2018.
 #### Homepage / Documentation
 https://github.com/crickbabs/BABS-MNASeqPE
 #### Authors
 Harshil Patel <harshil.patel@crick.ac.uk>
-------------------------------------------------------------------------------
*/

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                             PARAMETERS                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

def helpMessage() {
    log.info"""
    ======================================================
    BABS MNASeq Paired-End Pipeline v${params.version}
    ======================================================

    Usage:
    The typical command for running the pipeline is as follows:

    nextflow run main.nf --design design.csv --genome sacCer3 -profile conda

    Mandatory arguments:
      --design                      Comma separted file containing information about the samples in the experiment (see README.md)
      --genome                      Genome shortname (e.g. sacCer3)
      -profile                      Hardware config to use i.e. conda, babs_modules, standard

    References:                     If not specified in the configuration file or you wish to overwrite any of the reference parameters
      --fasta                       Path to Fasta reference file containing all chromosomes/contigs
      --bwa_index                   Path to BWA index

    Trimming options:
      --adapter                     3' adapter sequence trimmed by cutadapt (default: AGATCGGAAGAGC)

    Filtering options:
      --insert_min                  Minimum insert size for filtering of mono-nucleosome paired-end reads (default: 100)
      --insert_max                  Maximum insert size for filtering of mono-nucleosome paired-end reads (default: 200)
      --rm_duplicates               Remove duplicate reads from BAM alignment file

    Output options:
      --outdir                      The output directory where the results will be saved (default: './results')

    Other options:
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}

// SHOW HELP MESSAGE
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                SET UP CONFIGURATION VARIABLES                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////
/* --         DEFAULT PARAMETER VALUES         -- */
////////////////////////////////////////////////////

params.design = false
params.genome = false
params.profile = false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa_index ?: false : false
params.adapter = false
params.rm_duplicates = false
params.insert_min = 100
params.insert_max = 200
params.outdir = './results'
params.outdir_abspath = new File(params.outdir).getCanonicalPath().toString()
params.name = false
params.project = false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.fastqscreen_config = "$baseDir/conf/fastq_screen.conf.txt"
params.bamtools_filter_pe_config = "$baseDir/conf/bamtools_filter_pe_template.json"

// PRESET ADAPTER TRIMMING OPTION
adapter_seq = 'AGATCGGAAGAGC'
if (params.adapter){
    adapter_seq = params.adapter
}

// SPECIFY THE RUN NAME
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
    custom_runName = workflow.runName
}

output_docs = file("$baseDir/docs/output.md")

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

multiqc_config = file(params.multiqc_config)
if( !multiqc_config.exists() ) exit 1, "MultiQC config file not found: ${params.multiqc_config}"
fastqscreen_config = file(params.fastqscreen_config)

if( file(params.bamtools_filter_pe_config).exists() ){
    Channel
      .fromPath(params.bamtools_filter_pe_config)
      .set { bamtools_filter_pe_config_ch }
} else {
    exit 1, "BamTools filter config file not found: ${params.bamtools_filter_pe_config}"
}

////////////////////////////////////////////////////
/* --          CHECK INPUT FILES               -- */
////////////////////////////////////////////////////

if( params.fasta && file(params.fasta).exists() ){
    Channel
      .fromPath(params.fasta)
      .into { fasta_index_ch;
              fasta_markdup_collectmetrics_ch;
              fasta_igv_ch }
} else {
    exit 1, "Genome fasta file not found: ${params.fasta}"
}

if( params.bwa_index && file("${params.bwa_index}.amb").exists() ) {
    Channel
        .fromPath("${params.bwa_index}*.{amb,ann,bwt,pac,sa}")
        .into { bwa_index_bwa_aln_r1_ch;
                bwa_index_bwa_aln_r2_ch;
                bwa_index_bwa_sampe_ch }
}  else {
    exit 1, "BWA indices not found: ${params.bwa_index}"
}

if( params.design && file(params.design).exists() ){
    Channel
      .fromPath(params.design)
      .splitCsv(header:true, sep:',')
      .map { row -> [ [row.sample,"R"+row.replicate,"L"+row.run].join("_"),
                       row.sample,
                       row.replicate,
                       row.run,
                       file(row.fastq_1),
                       file(row.fastq_2) ] }
      .into { design_raw_fastqc_ch;
              design_raw_fastqscreen_ch;
              design_cutadapt_ch }
} else {
    exit 1, "Design file not found: ${params.design}"
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       HEADER LOG INFO                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

log.info "======================================================"
log.info "BABS MNASeq Paired-End Pipeline v${params.version}"
log.info "======================================================"
def summary = [:]
summary['Run name']                   = custom_runName ?: workflow.runName
summary['Design file']                = params.design
summary['Genome version']             = params.genome
summary['Genome fasta file']          = params.fasta
summary['BWA index']                  = params.bwa_index
summary['Adapter sequence']           = adapter_seq
summary['Minimum insert size']        = params.insert_min
summary['Maximum insert size']        = params.insert_max
summary['Remove duplicates']          = params.rm_duplicates
summary['Max memory']                 = params.max_memory
summary['Max CPUs']                   = params.max_cpus
summary['Max time']                   = params.max_time
summary['Output directory']           = params.outdir
summary['Output directory path']      = params.outdir_abspath
summary['Working directory']          = workflow.workDir
summary['Current home']               = "$HOME"
summary['Current user']               = "$USER"
summary['Current path']               = "$PWD"
summary['Script directory']           = workflow.projectDir
summary['Config profile']             = workflow.profile
if(workflow.revision) summary['Pipeline release'] = workflow.revision
if(params.project) summary['BABS project'] = params.project
log.info summary.collect { k,v -> "${k.padRight(26)}: $v" }.join("\n")
log.info "======================================================"

// CHECK THAT NEXTFLOW VERSION IS UP TO DATE ENOUGH
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "============================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                     PREPARE ANNOTATION FILES                        -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// CREATE FAIDX FOR REFERENCE GENOME
// CREATE CHROMOSOME SIZES FILE FOR BEDTOOLS
// CREATE BAMTOOLS CONFIG FILE WITH USER-DEFINED INSERT SIZES FOR FILTERING
process prep_genome {

    tag "$fasta"

    publishDir "${params.outdir}/genome", mode: 'copy'

    input:
    file fasta from fasta_index_ch
    file bamtools from bamtools_filter_pe_config_ch

    output:
    file "*.fai" into prep_genome_index_ch
    file "*.sizes" into prep_genome_sizes_replicate_bedgraph_ch,
                        prep_genome_sizes_replicate_bigwig_ch
    file "*.json" into prep_genome_bamtools_filter_pe_config_ch

    script:
          """
          samtools faidx ${fasta}
          cut -f 1,2 ${fasta}.fai > ${fasta}.sizes
          awk 'BEGIN{OFS="\t"}{print \$1, '0' , \$2}' ${fasta}.sizes > ${fasta}.bed

          sed 's/MIN_INSERT_SIZE/${params.insert_min}/g' <${bamtools} >bamtools_filter_pe.json
          sed -i -e 's/MAX_INSERT_SIZE/${params.insert_max}/g' bamtools_filter_pe.json
          """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        FASTQ QC                                     -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// RUN FASTQC ON RAW FASTQ FILES
process raw_fastqc {

    tag "$sampleid"

    label 'mediumcpu'

    publishDir "${params.outdir}/qc/fastqc/raw", mode: 'copy',
               saveAs: {filename -> filename.endsWith(".zip") ? "zip/$filename" : "$filename"}

    input:
    set val(sampleid), val(sample), val(replicate), val(run), file(fastq_1), file(fastq_2) from design_raw_fastqc_ch

    output:
    set val(sampleid), file("*.{zip,html}") into raw_fastqc_ch

    script:
        """
        ln -s ${fastq_1} ${sampleid}_1.fastq.gz
        ln -s ${fastq_2} ${sampleid}_2.fastq.gz
        fastqc --outdir ./ --threads ${task.cpus} ${sampleid}_1.fastq.gz
        fastqc --outdir ./ --threads ${task.cpus} ${sampleid}_2.fastq.gz
        """
}

// RUN FASTQSCREEN ON RAW FASTQ FILES
process raw_fastqscreen {

    tag "$sampleid"

    label 'mediumcpu'

    publishDir "${params.outdir}/qc/fastq_screen", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".txt")) "txt/$filename"
                            else if (filename.endsWith(".png")) "png/$filename"
                            else "$filename"
                        }

    input:
    set val(sampleid), val(sample), val(replicate), val(run), file(fastq_1), file(fastq_2) from design_raw_fastqscreen_ch

    output:
    set val(sampleid), file("*.{png,txt,html}") into raw_fastqscreen_ch

    when:
    params.fastqscreen_config

    script:
        """
        ln -s ${fastq_1} ${sampleid}_1.fastq.gz
        ln -s ${fastq_2} ${sampleid}_2.fastq.gz
        fastq_screen --outdir ./ \\
                     --conf ${fastqscreen_config} \\
                     --threads ${task.cpus} \\
                     ${sampleid}_1.fastq.gz
        fastq_screen --outdir ./ \\
                     --conf ${fastqscreen_config} \\
                     --threads ${task.cpus} \\
                     ${sampleid}_2.fastq.gz
        """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        ADAPTER TRIMMING                             -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// RUN CUTADAPT ON RAW FASTQ FILES
process cutadapt {

    tag "$sampleid"

    publishDir "${params.outdir}/qc/cutadapt", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".log")) "$filename"
                            else null
                        }

    input:
    set val(sampleid), val(sample), val(replicate), val(run), file(fastq_1), file(fastq_2) from design_cutadapt_ch

    output:
    set val(sampleid), file('*.trim.fastq.gz') into cutadapt_fastqc_ch,
                                                    cutadapt_bwa_aln_r1_ch,
                                                    cutadapt_bwa_aln_r2_ch,
                                                    cutadapt_bwa_sai_to_sam_ch
    set val(sampleid), file('*.log') into cutadapt_log_ch

    script:
        """
        ln -s ${fastq_1} ${sampleid}_1.fastq.gz
        ln -s ${fastq_2} ${sampleid}_2.fastq.gz
        cutadapt -a ${adapter_seq} \\
                 -A ${adapter_seq} \\
                 --minimum-length=25 \\
                 --quality-cutoff=20 \\
                 -o ${sampleid}_1.trim.fastq.gz \\
                 -p ${sampleid}_2.trim.fastq.gz \\
                 ${sampleid}_1.fastq.gz \\
                 ${sampleid}_2.fastq.gz > ${sampleid}.cutadapt.log
        """
}

// RUN FASTQC ON CUTADAPT TRIMMED FASTQ FILES
process trim_fastqc {

    tag "$sampleid"

    label 'mediumcpu'

    publishDir "${params.outdir}/qc/fastqc/trim", mode: 'copy',
               saveAs: {filename -> filename.endsWith(".zip") ? "zip/$filename" : "$filename"}

    input:
    set val(sampleid), file(fastqs) from cutadapt_fastqc_ch

    output:
    set val(sampleid), file("*.{zip,html}") into trim_fastqc_ch

    script:
        """
        fastqc --outdir ./ --threads ${task.cpus} ${fastqs[0]}
        fastqc --outdir ./ --threads ${task.cpus} ${fastqs[1]}
        """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        ALIGN                                        -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// RUN BWA ALN ON CUTADAPT TRIMMED FASTQ FILE FOR READ 1
process bwa_aln_r1 {

    tag "$sampleid"

    label 'highcpu'

    publishDir "${params.outdir}/align", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".sysout")) "sysout/$filename"
                            else null
                        }

    input:
    set val(sampleid), file(fastqs) from cutadapt_bwa_aln_r1_ch
    file index from bwa_index_bwa_aln_r1_ch.collect()

    output:
    set val(sampleid), file("*.sai") into bwa_aln_r1_ch
    set val(sampleid), file("*.sysout") into bwa_aln_r1_sysout_ch

    script:
        """
        bwa aln -t ${task.cpus} ${params.bwa_index} ${fastqs[0]} > ${sampleid}_1.sai 2> ${sampleid}_bwa_aln_r1.sysout
        """
}

// RUN BWA ALN ON CUTADAPT TRIMMED FASTQ FILE FOR READ 2
process bwa_aln_r2 {

    tag "$sampleid"

    label 'highcpu'

    publishDir "${params.outdir}/align", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".sysout")) "sysout/$filename"
                            else null
                        }

    input:
    set val(sampleid), file(fastqs) from cutadapt_bwa_aln_r2_ch
    file index from bwa_index_bwa_aln_r2_ch.collect()

    output:
    set val(sampleid), file("*.sai") into bwa_aln_r2_ch
    set val(sampleid), file("*.sysout") into bwa_aln_r2_sysout_ch

    script:
        """
        bwa aln -t ${task.cpus} ${params.bwa_index} ${fastqs[1]} > ${sampleid}_2.sai 2> ${sampleid}_bwa_aln_r2.sysout
        """
}

// RUN BWA SAMPE FOR SAI TO SAM CONVERSION
process bwa_sampe {

    tag "$sampleid"

    label "sampe"

    publishDir "${params.outdir}/align", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".sysout")) "sysout/$filename"
                            else null
                        }

    input:
    set val(sampleid), file(fastqs), file(sai1), file(sai2) from cutadapt_bwa_sai_to_sam_ch.join(bwa_aln_r1_ch, by: [0]).join(bwa_aln_r2_ch, by: [0])
    file index from bwa_index_bwa_sampe_ch.collect()

    output:
    set val(sampleid), file("*.sam") into bwa_sampe_ch
    set val(sampleid), file("*.sysout") into bwa_sampe_sysout_ch

    script:
        rg="\'@RG\\tID:${sampleid}\\tSM:${sampleid.toString().subSequence(0, sampleid.length() - 3)}\\tPL:illumina\\tLB:1\\tPU:1\'"
        """
        bwa sampe -r $rg ${params.bwa_index} ${sai1} ${sai2} ${fastqs} > ${sampleid}.sam 2> ${sampleid}_bwa_sampe.sysout
        """
}

// CONVERT SAM TO COORDINATE SORTED BAM
process bwa_bam {

    tag "$sampleid"

    label 'lowcpu'

    input:
    set val(sampleid), file(sam) from bwa_sampe_ch

    output:
    set val(sampleid), file("*.sorted.{bam,bam.bai}") into bwa_bam_ch
    set val(sampleid), file("*.flagstat") into bwa_bam_flagstat_ch

    script:
        out_prefix="${sampleid}"
        """
        samtools view -b -h -O BAM -@ ${task.cpus} -o ${out_prefix}.bam ${sam}
        samtools sort -@ ${task.cpus} -o ${out_prefix}.sorted.bam -T ${out_prefix} ${out_prefix}.bam
        samtools index ${out_prefix}.sorted.bam
        samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
        """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                    BAM POST-PROCESSING                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// RUN PICARD MARK DUPLICATES ON COORDINATE SORTED BAM
process markdup {

    tag "$sampleid"

    label 'lowcpu'

    publishDir "${params.outdir}/align", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".flagstat")) "flagstat/$filename"
                            else if (filename.endsWith(".idxstats")) "idxstats/$filename"
                            else if (filename.endsWith(".sysout")) "sysout/$filename"
                            else if (filename.endsWith(".metrics.txt")) "picard_metrics/$filename"
                            else null
                        }

    input:
    set val(sampleid), file(bam) from bwa_bam_ch

    output:
        set val(sampleid), file("*.{bam,bam.bai}") into markdup_filter_bam_ch, markdup_collectmetrics_in_ch
        set val(sampleid), file("*.flagstat") into markdup_flagstat_ch
        set val(sampleid), file("*.idxstats") into markdup_idxstats_ch
        set val(sampleid), file("*.metrics.txt") into markdup_metrics_ch
        set val(sampleid), file("*.sysout") into markdup_sysout_ch

    script:
        out_prefix="${sampleid}.mkD"
        """
        picard -Xmx${task.memory.toString().split(" ")[0]}g \\
               MarkDuplicates \\
               VALIDATION_STRINGENCY=LENIENT \\
               REMOVE_DUPLICATES=false \\
               ASSUME_SORTED=true \\
               TMP_DIR=tmp \\
               INPUT=${bam[0]} \\
               OUTPUT=${out_prefix}.sorted.bam \\
               METRICS_FILE=${out_prefix}.MarkDuplicates.metrics.txt \\
               >> ${out_prefix}.MarkDuplicates.sysout 2>&1
        samtools index ${out_prefix}.sorted.bam
        samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
        samtools idxstats ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.idxstats
        """
}

// RUN PICARD COLLECTMULTIPLE METRICS ON COORDINATE SORTED BAM
process markdup_collectmetrics {

    tag "$sampleid"

    label 'lowcpu'

    publishDir "${params.outdir}/align", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith("_metrics")) "picard_metrics/$filename"
                            else if (filename.endsWith(".pdf")) "picard_metrics/pdf/$filename"
                            else if (filename.endsWith(".sysout")) "sysout/$filename"
                        }

    input:
    set val(sampleid), file(bam) from markdup_collectmetrics_in_ch
    file fasta from fasta_markdup_collectmetrics_ch.collect()

    output:
    set val(sampleid), file("*{_metrics,.pdf}") into markdup_collectmetrics_ch
    set val(sampleid), file("*.sysout") into markdup_collectmetrics_sysout_ch

    script:
        out_prefix="${sampleid}.mkD"
        """
        picard -Xmx${task.memory.toString().split(" ")[0]}g \\
               CollectMultipleMetrics \\
               VALIDATION_STRINGENCY=LENIENT \\
               TMP_DIR=tmp \\
               INPUT=${bam[0]} \\
               OUTPUT=${out_prefix}.CollectMultipleMetrics \\
               REFERENCE_SEQUENCE=${fasta} \\
               >> ${out_prefix}.CollectMultipleMetrics.sysout 2>&1
        """
}

// FILTER BAM FILE TO KEEP (1) UNIQUELY MAPPED, (2) PRIMARY ALIGNMENT, (3) IN FR ORIENTATION ON SAME CHROMOSOME
//                         (4) NON-SOFTCLIPPED (bamtools), (5) INSERT SIZE > 100BP & < 200BP (bamtools), (6) MISMATCH <= 3 (bamtools)
process filter_bam {

    tag "$sampleid"

    label 'lowcpu'

    input:
    set val(sampleid), file(bam) from markdup_filter_bam_ch
    file bamtools_config from prep_genome_bamtools_filter_pe_config_ch.collect()

    output:
    set val(sampleid), file("*.flT.bam") into filter_bam_ch
    set val(sampleid), file("*.sorted.{bam,bam.bai}") into filter_bam_sort_ch
    set val(sampleid), file("*.flagstat") into filter_bam_flagstat_ch

    script:
        // 0x0001 = read paired
        // 0x0004 = read unmapped
        // 0x0008 = mate unmapped
        // 0x0100 = not primary alignment
        // 0x0400 = read is PCR or optical duplicate
        out_prefix="${sampleid}.flT"
        if (params.rm_duplicates) {
            """
            samtools view \\
                     -f 0x001 \\
                     -F 0x004 \\
                     -F 0x0008 \\
                     -F 0x0100 \\
                     -F 0x0400 \\
                     -q 1 \\
                     -b ${bam[0]} \\
                     | bamtools filter \\
                                -out ${out_prefix}.sorted.bam \\
                                -script ${bamtools_config}
            samtools index ${out_prefix}.sorted.bam
            samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
            samtools sort -n -@ ${task.cpus} -o ${out_prefix}.bam -T ${out_prefix} ${out_prefix}.sorted.bam
            """
        } else {
            """
            samtools view \\
                     -f 0x001 \\
                     -F 0x004 \\
                     -F 0x0008 \\
                     -F 0x0100 \\
                     -q 1 \\
                     -b ${bam[0]} \\
                     | bamtools filter \\
                                -out ${out_prefix}.sorted.bam \\
                                -script ${bamtools_config}
            samtools index ${out_prefix}.sorted.bam
            samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
            samtools sort -n -@ ${task.cpus} -o ${out_prefix}.bam -T ${out_prefix} ${out_prefix}.sorted.bam
            """
        }
}

// FILTER BAM FILE TO REMOVE READS WHERE ONE MATE HAS BEEN PHYSICALLY FILTERED OUT OF BAM FILE FROM PREVIOUS STEP.
process rm_orphan {

    tag "$sampleid"

    input:
    set val(sampleid), file(bam) from filter_bam_ch

    output:
    set val(sampleid), file("*.bam") into rm_orphan_ch

    script:
        out_prefix="${sampleid}.clN"
        """
        python $baseDir/bin/bampe_rm_orphan.py ${bam} ${out_prefix}.bam --only_prop_pair
        """
}

// SORT ORPHAN BAM FILE BY COORDINATE
process rm_orphan_sort_bam {

    tag "$sampleid"

    label 'lowcpu'

    publishDir "${params.outdir}/align", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".flagstat")) "flagstat/$filename"
                            else if (filename.endsWith(".idxstats")) "idxstats/$filename"
                            else if (filename.endsWith(".bam")) "$filename"
                            else if (filename.endsWith(".bai")) "$filename"
                        }

    input:
    set val(sampleid), file(orphan_bam) from rm_orphan_ch

    output:
    set val(sampleid), file("*.{bam,bam.bai,flagstat}") into rm_orphan_sort_bam_replicate_ch,
                                                             rm_orphan_sort_bam_replicate_rmdup_ch
    set val(sampleid), file("*.idxstats") into rm_orphan_sort_bam_idxstats_ch

    script:
        out_prefix="${sampleid}.clN"
        """
        samtools sort -@ ${task.cpus} -o ${out_prefix}.sorted.bam -T ${out_prefix} ${out_prefix}.bam
        samtools index ${out_prefix}.sorted.bam
        samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
        samtools idxstats ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.idxstats
        """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                    MERGE REPLICATE BAM                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// CREATE CHANNEL TO MERGE AT REPLICATE LEVEL
rm_orphan_sort_bam_replicate_ch.map { it -> [ it[0].toString().subSequence(0, it[0].length() - 3), it[1] ] }
                               .groupTuple(by: [0])
                               .map { it ->  [ it[0], it[1].flatten() ] }
                               .into { rm_orphan_sort_bam_replicate_merge_ch;
                                       rm_orphan_sort_bam_replicate_markdup_ch;
                                       rm_orphan_sort_bam_replicate_rmdup_ch }

// MERGE FILTERED BAM FILES AT REPLICATE LEVEL USING PICARD MERGESAMFILES
process merge_replicate {

    tag "$sampleid"

    label 'lowcpu'

    publishDir "${params.outdir}/align/replicateLevel", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".sysout")) "sysout/$filename"
                            else null
                        }

    input:
    set val(sampleid), file(bams) from rm_orphan_sort_bam_replicate_merge_ch

    output:
    set val(sampleid), file("*.{bam,bam.bai}") into merge_replicate_ch
    set val(sampleid), file("*.flagstat") into merge_replicate_flagstat_ch
    set val(sampleid), file("*.sysout") into merge_replicate_sysout_ch

    script:
        out_prefix="${sampleid}.RpL"
        bam_files = bams.findAll { it.toString().endsWith('.bam') }.sort()
        flagstat_files = bams.findAll { it.toString().endsWith('.flagstat') }.sort()
        if (bam_files.size() > 1) {
            """
            picard -Xmx${task.memory.toString().split(" ")[0]}g \\
                   MergeSamFiles \\
                   VALIDATION_STRINGENCY=LENIENT \\
                   SORT_ORDER=coordinate \\
                   TMP_DIR=tmp \\
                   ${'INPUT='+bam_files.join(' INPUT=')} \\
                   OUTPUT=${out_prefix}.sorted.bam \\
                   >> ${out_prefix}.MergeSamFiles.sysout 2>&1
            samtools index ${out_prefix}.sorted.bam
            samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
            """
        } else {
            """
            touch ${out_prefix}.sorted.bam
            touch ${out_prefix}.sorted.bam.bai
            touch ${out_prefix}.MergeSamFiles.sysout
            cp ${flagstat_files[0]} ${out_prefix}.sorted.bam.flagstat
            """
        }
}

// RUN PICARD MARK DUPLICATES
process merge_replicate_markdup {

    tag "$sampleid"

    label 'lowcpu'

    publishDir "${params.outdir}/align/replicateLevel", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".sysout")) "sysout/$filename"
                            else if (filename.endsWith(".metrics.txt")) "picard_metrics/$filename"
                            else null
                        }

    input:
    set val(sampleid), file(orphan_bams), file(merged_bam) from rm_orphan_sort_bam_replicate_markdup_ch.join(merge_replicate_ch, by: [0])

    output:
    set val(sampleid), file("*.{bam,bam.bai}") into merge_replicate_markdup_ch
    set val(sampleid), file("*.flagstat") into merge_replicate_markdup_flagstat_ch
    set val(sampleid), file("*.metrics.txt") into merge_replicate_markdup_metrics_ch
    set val(sampleid), file("*.sysout") into merge_replicate_markdup_sysout_ch

    script:
        out_prefix="${sampleid}.RpL.mkD"
        bam_files = orphan_bams.findAll { it.toString().endsWith('.bam') }.sort()
        flagstat_files = orphan_bams.findAll { it.toString().endsWith('.flagstat') }.sort()
        if (bam_files.size() > 1) {
            """
            picard -Xmx${task.memory.toString().split(" ")[0]}g \\
                   MarkDuplicates \\
                   VALIDATION_STRINGENCY=LENIENT \\
                   REMOVE_DUPLICATES=false \\
                   ASSUME_SORTED=true \\
                   TMP_DIR=tmp \\
                   INPUT=${merged_bam[0]} \\
                   OUTPUT=${out_prefix}.sorted.bam \\
                   METRICS_FILE=${out_prefix}.MarkDuplicates.metrics.txt \\
                   >> ${out_prefix}.MarkDuplicates.sysout 2>&1
            samtools index ${out_prefix}.sorted.bam
            samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
            """
        } else {
            """
            touch ${out_prefix}.sorted.bam
            touch ${out_prefix}.sorted.bam.bai
            touch ${out_prefix}.MarkDuplicates.sysout
            touch ${out_prefix}.MarkDuplicates.metrics.txt
            cp ${flagstat_files[0]} ${out_prefix}.sorted.bam.flagstat
            """
        }
}

// REMOVE DUPLICATES USING SAMTOOLS
process merge_replicate_rmdup {

    tag "$sampleid"

    publishDir "${params.outdir}/align/replicateLevel", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".flagstat")) "flagstat/$filename"
                            else if (filename.endsWith(".bam")) "$filename"
                            else if (filename.endsWith(".bai")) "$filename"
                        }

    input:
    set val(sampleid), file(orphan_bams), file(markdup_bam) from rm_orphan_sort_bam_replicate_rmdup_ch.join(merge_replicate_markdup_ch, by: [0])

    output:
    set val(sampleid), file("*.{bam,bam.bai}") into merge_replicate_rmdup_name_bam_ch,
                                                    merge_replicate_rmdup_bedgraph_ch
    set val(sampleid), file("*.flagstat") into merge_replicate_rmdup_flagstat_ch,
                                               merge_replicate_rmdup_flagstat_bedgraph_ch

    script:
        out_prefix="${sampleid}.RpL.rmD"
        bam_files = orphan_bams.findAll { it.toString().endsWith('.bam') }.sort()
        if (bam_files.size() > 1) {
            if (params.rm_duplicates) {
                """
                samtools view -bF 0x400 ${markdup_bam[0]} > ${out_prefix}.sorted.bam
                samtools index ${out_prefix}.sorted.bam
                samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
                """
            } else {
                """
                cp ${markdup_bam[0]} ${out_prefix}.sorted.bam
                samtools index ${out_prefix}.sorted.bam
                samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
                """
            }
        } else {
            """
            cp ${bam_files[0]} ${out_prefix}.sorted.bam
            samtools index ${out_prefix}.sorted.bam
            samtools flagstat ${out_prefix}.sorted.bam > ${out_prefix}.sorted.bam.flagstat
            """
        }
}

// SORT BAM FILE BY NAME
process merge_replicate_name_bam {

    tag "$sampleid"

    label 'lowcpu'

    input:
    set val(sampleid), file(bam) from merge_replicate_rmdup_name_bam_ch

    output:
    set val(sampleid), file("*.bam") into merge_replicate_name_bam_ch

    script:
        out_prefix="${sampleid}.RpL.rmD"
        """
        samtools sort -n -@ ${task.cpus} -o ${out_prefix}.bam -T ${out_prefix} ${out_prefix}.sorted.bam
        """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                         COVERAGE TRACKS                             -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// CREATE NORMALISED BEDGRAPH FILE USING BEDTOOLS GENOMECOVERAGEBED
// CALCULATE SCALE-FACTOR FROM FLAGSTAT FILE
process merge_replicate_bedgraph {

    tag "$sampleid"

    publishDir "${params.outdir}/align/replicateLevel/bigwig", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".txt")) "scale_factor/$filename"
                            else null
                        }

    input:
    set val(sampleid), file(bam), file(flagstat) from merge_replicate_rmdup_bedgraph_ch.join(merge_replicate_rmdup_flagstat_bedgraph_ch, by: [0])
    file chrom_sizes from prep_genome_sizes_replicate_bedgraph_ch.collect()

    output:
    set val(sampleid), file("*.bg") into merge_replicate_bedgraph_ch
    set val(sampleid), file("*.txt") into merge_replicate_bedgraph_scale_ch

    script:
        out_prefix="${sampleid}.RpL.rmD"
        """
        SCALE_FACTOR=\$(grep 'read1' ${flagstat} | awk '{print 1000000/\$1}')
        echo \$SCALE_FACTOR > ${out_prefix}.scale_factor.txt
        genomeCoverageBed -ibam ${bam[0]} -bg -trackline -scale \$SCALE_FACTOR -pc -g ${chrom_sizes} >  ${out_prefix}.bg
        """
}

// CONVERT BEDGRAPH TO BIGWIG USING KENTOOLS
process merge_replicate_bigwig {

    tag "$sampleid"

    label 'bigwig'

    publishDir "${params.outdir}/align/replicateLevel/bigwig", mode: 'copy'

    input:
    set val(sampleid), file(bedgraph) from merge_replicate_bedgraph_ch
    file chrom_sizes from prep_genome_sizes_replicate_bigwig_ch.collect()

    output:
    set val(sampleid), file("*.bigWig") into merge_replicate_bigwig_ch

    script:
        out_prefix="${sampleid}.RpL.rmD"
        """
        wigToBigWig -clip ${bedgraph}  ${chrom_sizes} ${out_prefix}.bigWig
        """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                         DANPOS2                                     -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

process merge_replicate_bam_to_bed {

    tag "$sampleid"

    input:
    set val(sampleid), file(bam) from merge_replicate_name_bam_ch

    output:
    set val(sampleid), file("*.bed") into merge_replicate_bam_to_bed_ch

    script:
        out_prefix="${sampleid}.RpL.rmD"
        """
        bamToBed -i ${bam[0]} > ${out_prefix}.bed
        """
}

process merge_replicate_danpos {

    tag "$sampleid"

    publishDir "${params.outdir}/align/replicateLevel/danpos", mode: 'copy',
                saveAs: {filename ->
                            if (filename.endsWith(".sysout")) "sysout/$filename"
                            else if (filename.endsWith(".wig.gz")) "$filename"
                            else if (filename.endsWith(".xls")) "$filename"
                            else null
                        }

    input:
    set val(sampleid), file(bed) from merge_replicate_bam_to_bed_ch

    output:
    set val(sampleid), file("*.xls") into merge_replicate_danpos_xls_ch
    set val(sampleid), file("*.wig.gz") into merge_replicate_danpos_wig_ch
    set val(sampleid), file("*.sysout") into merge_replicate_danpos_sysout_ch

    script:
        out_prefix="${sampleid}.RpL.rmD"
        """
        danpos.py dpos \\
                  ./${bed} \\
                  --paired 1 \\
                  --span 1 \\
                  --smooth_width 20 \\
                  --width 40 \\
                  --count 1000000 \\
                  --out ./result/ \\
                  > ${out_prefix}.danpos.sysout 2>&1
        mv ./result/*/*.wig .
        mv ./result/*/*.xls .
        gzip *.wig
        """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                             IGV                                     -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

merge_replicate_bigwig_ch.map { it -> [ "igv", it[1] ] }
                         .groupTuple(by: [0])
                         .map { it ->  it[1].flatten().sort() }
                         .set { igv_input_ch }

// CUSTOM SCRIPT TO REORDER TRACKS IN IGV SESSION FILE MORE SENSIBLY.
// PROBABLY POSSIBLE WITH NEXTFLOW BUT MUCH EASIER IN PYTHON!
process igv_session {

    tag "igv_session"

    publishDir "${params.outdir}/igv", mode: 'copy'

    input:
    file igvs from igv_input_ch
    file fasta from fasta_igv_ch.collect()

    output:
    file "*.{xml,txt}" into igv_session_ch

    script:
        """
        [ ! -f ${params.outdir_abspath}/genome/${fasta.getName()} ] && ln -s ${params.fasta} ${params.outdir_abspath}/genome/${fasta.getName()}
        python $baseDir/bin/igv_get_files.py ${params.outdir_abspath} igv_files.txt
        python $baseDir/bin/igv_files_to_session.py igv_session.xml igv_files.txt ${params.outdir_abspath}/genome/${fasta.getName()}
        """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                          MULTIQC                                    -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// RUN THIS AFTER IGV BECAUSE ITS AT THE END OF THE PIPELINE
process multiqc {

    tag "multiqc"

    publishDir "${params.outdir}/qc/multiqc", mode: 'copy'

    input:
    file igvs from igv_session_ch

    output:
    file "*" into multiqc_ch

    script:
        """
        multiqc ${params.outdir_abspath} \\
                -f \\
                --config ${params.multiqc_config} \\
                --filename BABS-MNASeqPE_multiqc_report.html \\
                -m fastqc \\
                -m fastq_screen \\
                -m cutadapt \\
                -m samtools \\
                -m picard
        """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                          QC LOG FILE                                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// RUN THIS AFTER MULTIQC BECAUSE ITS AT THE END OF THE PIPELINE
process qc_to_tsv {

    tag "qc_to_tsv"

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    file multiqcs from multiqc_ch

    output:
    file "*.tsv" into pipeline_qc_to_tsv_ch

    script:
        """
        python $baseDir/bin/pipeline_qc_to_tsv.py ${params.outdir_abspath} BABS-MNASeqPE_pipeline_qc.tsv
        """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        END OF PIPELINE                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
