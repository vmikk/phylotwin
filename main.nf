#!/usr/bin/env nextflow

println( "Running PhyloTwin diversity estimation pipeline" )

// Count number of occurrences per species
process subset_data {

    cpus 2

    input:
      path occurrences
      path specieskeys

    output:
      path "aggregated_counts.parquet", emit: occ_counts_long
      path "dataset_keys.tsv",          emit: dataset_keys

    script:
    """
    echo -e "Subsetting data\n"

    echo "Input directory: " ${occurrences}

    subset_data.R \
      --input       ${occurrences} \
      --output      ./results \
      --tree        ${params.tree} \
      --phylum      ${params.phylum} \
      --class       ${params.class} \
      --order       ${params.order} \
      --family      ${params.family} \
      --genus       ${params.genus} \
      --specieskeys ${specieskeys} \
      --resolution  ${params.resolution} \
      --country     ${params.country} \
      --latmin      ${params.latmin} \
      --latmax      ${params.latmax} \
      --lonmin      ${params.lonmin} \
      --lonmax      ${params.lonmax} \
      --polygon     ${params.polygon} \
      --data        ${params.data} \
      --threads     ${task.cpus}

    echo "..Done"
    """
}

