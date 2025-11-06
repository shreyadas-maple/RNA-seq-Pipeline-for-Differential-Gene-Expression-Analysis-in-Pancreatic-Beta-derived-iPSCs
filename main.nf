#!/usr/bin/env nextflow

include {FASTQC} from './modules/fastqc'
include {GENERATEGENEIDS} from './modules/generateGeneIDs'
include {STAR} from './modules/star'
include {STAR_ALIGNER} from './modules/star_aligner'
include {MULTIQC} from './modules/multiqc'
include {VERSE} from './modules/verse'
include {CONCAT_COUNTS} from './modules/concat_counts'

workflow {

    // Make the align_ch to pair up the replicates -> total 6 rows
    Channel.fromFilePairs(params.reads)
    | set{align_ch}

    // Make the fastqc_channel to be used for quality evaluation of the short reads
    // -> total 12 rows for each replicate
    Channel.fromFilePairs(params.reads)
    | transpose()
    | set {fastqc_channel}

    // Run FASTQC on the paired-end short reads
    FASTQC(fastqc_channel)

    // Make the text file with the delimeted geneIDs, gene_names
    GENERATEGENEIDS(params.gtf)

    // Make the STAR index file using the reference genome and the gtf file
    STAR(params.genome, params.gtf)

    // Make the channel for STAR_ALIGN process
    Channel.fromFilePairs(params.reads)
    | map {sample_id, files -> tuple (sample_id, files[0], files[1])}
    | set {star_align}

    // Call the STAR aligner on the reads and the star index files
    STAR_ALIGNER(star_align, params.index)

    // Make a channel for MultiQC 

    // First we will grab just the log files from STAR_ALIGNER and name the channel ch_1
    STAR_ALIGNER.out.log
                    .map {sample_id, logs -> logs}
                    .set {ch_1}

    // Second we get the zip files from the output of FASTQC, flatten it, then mix it with
    // ch_1 (STAR_ALIGNER log output) and collect it into a list
    // NOTE: this will result in all the files being out of order, but this doesn't matter
    // with MULTIQC, because it is smart ;)
    FASTQC.out.zip
              .map {sample_id, files -> files}
              .flatten()
              .mix(ch_1)
              .collect()
              .set{multiqc_channel}
    
    // Send the channel with the fastqc.html and Star align log files to multiqc
    MULTIQC(multiqc_channel)

    VERSE(STAR_ALIGNER.out.bam, params.gtf)

    // Grab all the .exon.txt files from the results directory
    // and name the channel as counts_ch
    Channel.fromPath("${params.outdir}/*.exon.txt")
           .collect()
           .set{counts_ch}

    // Pass the counts_ch into the CONCAT_COUNTS process to
    // generate a .csv file with the gene counts for each sample
    CONCAT_COUNTS(counts_ch)

}
