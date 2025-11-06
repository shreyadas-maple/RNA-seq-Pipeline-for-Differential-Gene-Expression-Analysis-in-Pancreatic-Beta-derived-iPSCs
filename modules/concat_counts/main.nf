process CONCAT_COUNTS{

    container 'ghcr.io/bf528/pandas:latest'

    publishDir params.outdir, mode: 'copy'

    input:
    path(exons)

    output:
    path("counts.csv")

    script:
    """
    concat_counts.py -i $exons -o counts.csv
    """
}