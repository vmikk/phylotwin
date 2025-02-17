#!/usr/bin/env nextflow


ver = "version " + workflow.manifest.version
if (workflow.commitId) { ver += " revision " + workflow.commitId.substring(0, 7) }
println( "Running PhyloTwin diversity estimation pipeline, ${ver}\n" )

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

    containerOptions "-v ${projectDir}:${projectDir}"

    input:
      path occurrences
      path specieslist
      path polygon
      val  tree

    output:
      path "aggregated_counts.parquet", emit: occ_counts_long
      path "aggregated_counts.csv",     emit: occ_counts_csv, optional: true
      path "dataset_keys.tsv",          emit: dataset_keys
      path "phylogenetic_tree.nex",     emit: phylogenetic_tree
      // path "species_selected.txt.gz",  emit: species_selected
      // path "h3_cells.txt.gz",          emit: h3_cells

    script:
    def spnamesArg = specieslist.name != 'no_specieslist' ? "--specieslist ${specieslist}" : ''
    def polyArg = polygon.name != 'no_polygon' ? "--polygon ${polygon}" : ''
    def csvArg  = params.bd_indices ? "--csv TRUE" : ""
    def duckdbArg = (workflow.containerEngine == 'docker' || workflow.containerEngine == 'singularity') ? '--duckdb_extdir "/usr/local/bin/duckdb_ext"' : ''
    """
    echo -e "Subsetting data\n"

    subset_data.R \
      --inpdir      ${occurrences} \
      --tree        ${tree} \
      --phylum      ${params.phylum} \
      --class       ${params.classs} \
      --order       ${params.order} \
      --family      ${params.family} \
      --genus       ${params.genus} \
      ${spnamesArg} \
      --resolution  ${params.resolution} \
      --country     ${params.country} \
      --latmin      ${params.latmin} \
      --latmax      ${params.latmax} \
      --lonmin      ${params.lonmin} \
      --lonmax      ${params.lonmax} \
      ${polyArg} \
      --data        ${params.data} \
      ${csvArg} \
      ${duckdbArg} \
      --threads     ${task.cpus}

    echo "..Done"
    """
}



// Estimate diversity (mostly fast indices)
process estimate_diversity {

    label "container_phylotwin"

    // Code `99` is used to indicate that no species occurrences were found with the specified filters
    errorStrategy { task.exitStatus == 99 ? 'ignore' : 'retry' }

    input:
      path table
      path tree

    output:
      path "estimate_diversity_results.qs", emit: qs, optional: true

    script:
    """
    echo -e "Estimating diversity\n"

    echo "Species abundances: " ${table}
    echo "Phylogenetic tree: "  ${tree}

    estimate_diversity.R \
      --input       ${table} \
      --tree        ${tree} \
      --output      estimate_diversity_results \
      --div         ${params.div} \
      --threads     ${task.cpus}

    exit_code=\$?

    # ## Compress GeoJSON file
    # if [ -f estimate_diversity_results.geojson ]; then
    #   gzip -4 estimate_diversity_results.geojson
    # fi

    echo "..Done"

    ## Exit with the same code as the R script
    if [ \$exit_code -ne 0 ]; then
      exit \$exit_code
    fi

    """
}


// Combine Biodiverse results with `estimate_diversity` results
// Prepare spatial data for visualization
process aggregate_spatial_results {

    label "container_phylotwin"

    publishDir OUTDIR_2_DIV, mode: "${params.publish_dir_mode}", overwrite: true

    input:
      path estimate_diversity_results
      path biodiverse_results

    output:
      path "diversity_estimates.geojson",    emit: geojson,    optional: true
      path "diversity_estimates.gpkg",       emit: geopackage, optional: true
      path "diversity_estimates.txt",        emit: txt,        optional: true

    script:
    def bdArg = biodiverse_results.name != 'no_biodiverse' ? "${biodiverse_results}" : 'NA'
    """
    combine_results.R \
      --estdiv     ${estimate_diversity_results} \
      --biodiv     ${bdArg} \
      --resolution ${params.resolution} \
      --output     diversity_estimates
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
  ch_spnames = params.specieslist ? Channel.fromPath(params.specieslist) : Channel.fromPath("no_specieslist", checkIfExists: false)

  // Pipeline data
  ch_data = Channel.fromPath(params.data, type: 'dir', checkIfExists: true)
  // ch_taxatables   = params.data + "/TaxonomyTables"
  // ch_phylotrees   = params.data + "/Phylotrees"
  // ch_countries_h3 = params.data + "/Countries_H3"

  // Subset data
  subset_data(ch_occ, ch_spnames, ch_poly, ch_tree)

  // Estimate diversity [should be always enabled, as we need Record summaries and Redundancy estimates]
  estimate_diversity(
    subset_data.out.occ_counts_long,
    subset_data.out.phylogenetic_tree
  )

  // Biodiverse-based subworkflow
  if (params.bd_indices) {

    BIODIVERSE(
      subset_data.out.occ_counts_csv,
      subset_data.out.phylogenetic_tree
    )

    // Combine Biodiverse results with `estimate_diversity` results
    aggregate_spatial_results(
      estimate_diversity.out.qs,
      BIODIVERSE.out.bdres
    )

    // end of Biodiverse subworkflow
  } else {
    // If Biodiverse is not enabled, use `estimate_diversity` results

    aggregate_spatial_results(
      estimate_diversity.out.qs,
      Channel.fromPath("no_biodiverse", checkIfExists: false))

  }

  // Plot diversity indices into HTML files
  if (params.noviz == false) {

      // Prepare channel with diversity index names for plotting
      // ch_vizindices = Channel.of("Richness", "PD", "SES.PD")    // for debugging
      ch_vizindices = Channel.from(params.viz.split(','))

      // Replicate diversity estimates for each index
      ch_divests = ch_vizindices.combine(aggregate_spatial_results.out.txt).map { it[1] }

      // Plot diversity indices
      viz_leaflet(
        ch_divests,
        ch_vizindices
      )

  } // end of noviz

}



// On completion
workflow.onComplete {
    println "Completed at      : ${workflow.complete}"
    println "Duration          : ${workflow.duration}"
    println "Success           : ${workflow.success}"
    println "Execution status  : ${workflow.success ? 'All done!' : 'Failed' }"
}

// On error
workflow.onError {
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

