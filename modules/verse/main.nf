process VERSE{
    label 'process_medium'

    // The container for the process VERSE
    container 'ghcr.io/bf528/verse:latest'

    // Push the results to the results directory
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample), path(bam)
    path(gtf)

    output:
    path("*.exon.txt")

    script:
    """
    verse -S -a $gtf -o ${sample} $bam
    """
}