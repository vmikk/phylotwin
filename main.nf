#!/usr/bin/env nextflow


println( "Running PhyloTwin diversity estimation pipeline\n" )

// Path to the directory with built-in phylogenetic trees
def treesDir = "${projectDir}/data/Phylotrees"

// Include Biodiverse sub-workflow
include { BIODIVERSE } from "./subworkflows/biodiverse.nf"

// Directory for publishing outputs
OUTDIR = params.userid ? params.outdir + "/" + params.userid : params.outdir
OUTDIR_1_SUB = OUTDIR + "/01.Occurrence_subset"
OUTDIR_2_DIV = OUTDIR + "/02.Diversity_estimates"
OUTDIR_3_VIZ = OUTDIR + "/03.Visualization"

// Prepare data subset for diversity estimation (based on spatial and taxonomy filters)
process subset_data {

    label "container_phylotwin"

    publishDir OUTDIR_1_SUB, mode: "${params.publish_dir_mode}", overwrite: true

    input:
      path occurrences
      path specieskeys
      path polygon

    output:
      path "aggregated_counts.parquet", emit: occ_counts_long
      path "aggregated_counts.csv",     emit: occ_counts_csv, optional: true
      path "dataset_keys.tsv",          emit: dataset_keys

    script:
    def spkeysArg = specieskeys.name != 'no_specieskeys' ? "--specieskeys ${specieskeys}" : ''
    def polyArg = polygon.name != 'no_polygon' ? "--polygon ${polygon}" : ''
    def csvArg  = params.bd_indices ? "--csv TRUE" : ""
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
      ${csvArg} \
      --threads     ${task.cpus}

    echo "..Done"
    """
}



// Estimate diversity (mostly fast indices)
process estimate_diversity {

    label "container_phylotwin"

    // Code `99` is used to indicate that no species occurrences were found with the specified filters
    errorStrategy { task.exitStatus == 99 ? 'ignore' : 'retry' }

    publishDir OUTDIR_2_DIV, mode: "${params.publish_dir_mode}", overwrite: true

    input:
      path table
      path tree

    output:
      path "diversity_estimates.txt",        emit: txt,        optional: true
      path "diversity_estimates.qs",         emit: qs,         optional: true
      path "diversity_estimates.gpkg",       emit: geopackage, optional: true
      path "diversity_estimates.geojson",    emit: geojson,    optional: true

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

    exit_code=\$?

    # ## Compress GeoJSON file
    # if [ -f diversity_estimates.geojson ]; then
    #   gzip -4 diversity_estimates.geojson
    # fi

    echo "..Done"

    ## Exit with the same code as the R script
    if [ \$exit_code -ne 0 ]; then
      exit \$exit_code
    fi

    """
}


// Plot diversity indices (interactive map - Leaflet-based choropleth)
process viz_leaflet {

    label "container_phylotwin"

    publishDir OUTDIR_3_VIZ, mode: "${params.publish_dir_mode}", overwrite: true
    tag "${index}"

    input:
      path divestimates  // table with diversity estimates
      val index          // index name

    output:
      path "Choropleth_*.html", emit: choropleth
      path "Leaflet_*.qs",      emit: leaflet, optional: true

    script:
    """
    Visualization_Leaflet.R \
      --divestimates    ${divestimates} \
      --variable        ${index} \
      --shortid         TRUE \
      --antimeridianfix TRUE \
      --redundancy      ${params.redundancy}

      # --variable2
      # --saveqs
    """
}



workflow {

  // Channels
  ch_occ  = Channel.fromPath(params.occ, type: 'dir', checkIfExists: true)
  ch_tree = Channel.value(params.tree)
    .map { tree ->
        def treePath = "${treesDir}/${tree}"
        if (!file(treePath).exists()) {
            throw new IllegalArgumentException("Specified tree file '${tree}' does not exist in ${treesDir}")
        }
        return file(treePath)
    }

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

  // Biodiverse-based subworkflow
  if (params.bd_indices) {

    BIODIVERSE(
      subset_data.out.occ_counts_csv,
      ch_tree
    )

  } // end of Biodiverse subworkflow

  // Plot diversity indices into HTML files
  if (params.noviz == false) {

      // Prepare channel with diversity index names for plotting
      // ch_vizindices = Channel.of("Richness", "PD", "SES.PD")    // for debugging
      ch_vizindices = Channel.from(params.viz.split(','))

      // Replicate diversity estimates for each index
      ch_divests = ch_vizindices.combine(estimate_diversity.out.qs).map { it[1] }

      // Plot diversity indices
      viz_leaflet(
        ch_divests,
        ch_vizindices
      )

  } // end of noviz

}



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

