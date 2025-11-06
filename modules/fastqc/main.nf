
process FASTQC {
    // Define the number of cpus to be used for this process
    label 'process_low'

    // Define the Docker image to be used for this process
    container 'ghcr.io/bf528/fastqc:latest'

    // Put the output files in the results directory
    publishDir params.outdir, mode: 'copy'

    // These are the inputs from the fastqc_channel
    input:
    tuple val(name), path(fastq)

    // We want the .zip and .html files to be in the results directory
    output:
    tuple val(name), path('*.zip'), emit: zip
    tuple val(name), path('*.html'), emit: html

    // We write the script to call fastqc on the path to the fastq file and with specified cpus
    script:
    """
    fastqc $fastq -t $task.cpus
    """

}