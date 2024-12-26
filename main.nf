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

workflow {

  // Channels
  ch_occ  = Channel.fromPath(params.occ, type: 'dir', checkIfExists: true)
  ch_tree = Channel.fromPath(params.tree)
    
  // Optional input files
  ch_poly   = params.polygon ? Channel.fromPath(params.polygon) : Channel.empty()
  ch_spkeys = params.specieskeys ? Channel.fromPath(params.specieskeys) : Channel.empty()

  ch_occ.view()
  ch_tree.view()
  ch_poly.view()
  ch_spkeys.view()

  // Pipeline data
  ch_data = Channel.fromPath(params.data, type: 'dir', checkIfExists: true)
  // ch_taxatables   = params.data + "/TaxonomyTables"
  // ch_phylotrees   = params.data + "/Phylotrees"
  // ch_countries_h3 = params.data + "/Countries_H3"

  // Subset data
  subset_data(ch_occ, ch_spkeys, ch_poly)


// On completion
workflow.onComplete {
    println "Pipeline completed at : $workflow.complete"
    println "Duration              : ${workflow.duration}"
    println "Execution status      : ${workflow.success ? 'All done!' : 'Failed' }"
}


