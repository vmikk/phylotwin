#!/usr/bin/env nextflow

println( "Running PhyloTwin diversity estimation pipeline" )


// Prepare data subset for diversity estimation (based on spatial and taxonomy filters)
process subset_data {

    input:
      path occurrences
      path specieskeys
      path polygon

    output:
      path "aggregated_counts.parquet", emit: occ_counts_long
      path "dataset_keys.tsv",          emit: dataset_keys

    script:
    def spkeysArg = specieskeys.name != 'no_specieskeys' ? "--specieskeys ${specieskeys}" : ''
    def polyArg = polygon.name != 'no_polygon' ? "--polygon ${polygon}" : ''
    """
    echo -e "Subsetting data\n"

    subset_data.R \
      --inpdir      ${occurrences} \
      --output      ./results \
      --tree        ${params.tree} \
      --phylum      ${params.phylum} \
      --class       ${params.classs} \
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



// Estimate diversity (mostly fast indices)
process estimate_diversity {

    input:
      path table
      path tree

    output:
      path "diversity_estimates.txt",  emit: txt
      path "diversity_estimates.qs",   emit: qs
      path "diversity_estimates.gpkg", emit: geopackage

    script:
    """
    echo -e "Estimating diversity\n"

    echo "Species abundances: " ${table}
    echo "Phylogenetic tree: "  ${tree}

    estimate_diversity.R \
      --input       ${table} \
      --tree        ${tree} \
      --output      diversity_estimates \
      --div         ${params.div} \
      --threads     ${task.cpus}

    echo "..Done"
    """
}


// Plot diversity indices (interactive map - Leaflet-based choropleth)
process viz_leaflet {

    input:
      path divestimates  // table with diversity estimates
      val index          // index name

    output:
      path "Choropleth_*.html", emit: choropleth
      path "Leaflet_*.qs",      emit: leaflet

    script:
    """
    Visualization_Leaflet.R \
      --divestimates    ${divestimates} \
      --variables       ${index} \
      --shortid         TRUE \
      --antimeridianfix TRUE

    """
}



workflow {

  // Channels
  ch_occ  = Channel.fromPath(params.occ, type: 'dir', checkIfExists: true)
  ch_tree = Channel.fromPath(params.tree)
    
  // Optional input files
  ch_poly   = params.polygon ? Channel.fromPath(params.polygon) : Channel.fromPath("no_polygon", checkIfExists: false)
  ch_spkeys = params.specieskeys ? Channel.fromPath(params.specieskeys) : Channel.fromPath("no_specieskeys", checkIfExists: false)

  // Pipeline data
  ch_data = Channel.fromPath(params.data, type: 'dir', checkIfExists: true)
  // ch_taxatables   = params.data + "/TaxonomyTables"
  // ch_phylotrees   = params.data + "/Phylotrees"
  // ch_countries_h3 = params.data + "/Countries_H3"

  // Subset data
  subset_data(ch_occ, ch_spkeys, ch_poly)

  // Estimate diversity
  estimate_diversity(
    subset_data.out.occ_counts_long,
    ch_tree
  )


// On completion
workflow.onComplete {
    println "Pipeline completed at : $workflow.complete"
    println "Duration              : ${workflow.duration}"
    println "Execution status      : ${workflow.success ? 'All done!' : 'Failed' }"
}

// On error
workflow.onError {
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

