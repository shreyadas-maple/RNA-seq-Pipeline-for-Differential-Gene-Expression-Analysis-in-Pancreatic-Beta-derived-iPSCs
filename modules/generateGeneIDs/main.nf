process GENERATEGENEIDS {
    // Define the number of cpus to be used for this process
    label 'process_low'

    // Define the Docker image to be used for this process
    container 'ghcr.io/bf528/biopython:latest'

    // Put the output files in the results directory
    publishDir params.outdir, mode: 'copy'

    // The input is the gtf file
    input:
    path(gtf)

    // This is the .csv file that is produced by the script
    output:
    path("genes.csv")

    script:
    """
    generateGeneIDs.py -i $gtf -o genes.csv
    """

}