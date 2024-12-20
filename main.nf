#!/usr/bin/env nextflow

println( "Running PhyloTwin diversity estimation pipeline" )

// Count number of occurrences per species
process subset_data {

    input:
      path occurrences
      path specieskeys
      path polygon

    output:
      path "aggregated_counts.parquet", emit: occ_counts_long
      path "dataset_keys.tsv",          emit: dataset_keys

    script:
    def spkeysArg = specieskeys ? "--specieskeys ${specieskeys}" : ''
    def polyArg = polygon ? "--polygon ${polygon}" : ''
    """
    echo -e "Subsetting data\n"

    echo "Input directory: " ${occurrences}

    subset_data.R \
      --inpdir      ${occurrences} \
      --output      ./results \
      --tree        ${params.tree} \
      --phylum      ${params.phylum} \
      --class       ${params.class} \
      --order       ${params.order} \
      --family      ${params.family} \
      --genus       ${params.genus} \
      ${spkeysArg} \
      --resolution  ${params.resolution} \
      --country     ${params.country} \
      --latmin      ${params.latmin} \
      --latmax      ${params.latmax} \
      --lonmin      ${params.lonmin} \
      --lonmax      ${params.lonmax} \
      ${polyArg} \
      --data        ${params.data} \
      --threads     ${task.cpus}

    echo "..Done"
    """
}

