
// Subworkflow for Biodiverse-based diversity estimation

// How many randomization iterations should be per Biodiverse thread?
if( params.bd_threads > 1 ) {
    randomization_chunks = (0..<params.bd_threads)
    iterations_per_thread = params.rnd / params.bd_threads
    iterations_per_thread = Math.ceil(iterations_per_thread).toInteger()
}
else {
    randomization_chunks = [ 1 ]
    iterations_per_thread = params.rnd
}
biodiverse_args = "function=" + params.bd_randname + " max_iters=" + iterations_per_thread


// Create Biodiverse input files
process prep_biodiv {

    label "container_biodiverse"

    // cpus 1

    input:
      path occurrences    // Occurrence data, CSV
      path tree           // Phylogenetic tree, Newick

    output:
      path "occ.bds",          emit: BDS
      path "tree.bts",         emit: BTS
      path "occ_analysed.bds", emit: BDA
      path "occ.bds.csv",      emit: BDOBS

    script:
    """
    ## Prepare Biodiverse input file
    ## NB! column numbers are zero-based here
    ## Latitude = Y, Longitude = X
    
    ## Input occurrence data format:
    # H3, Latitude, Longitude, specieskey, total_records
    #  0,    1,         2,         3,           4

    echo -e "\n\n---- Preparing occurrence data ----\n\n"
    Biodiverse_00_create_bds.pl \
      --csv_file ${occurrences} \
      --out_file "occ.bds" \
      --label_column_number     '3' \
      --sampcount_column_number '4' \
      --group_column_number_x   '2' \
      --group_column_number_y   '1' \
      --cell_size_x '0' \
      --cell_size_y '0'

    ## Prepare the tree for Biodiverse
    echo -e "\n\n---- Preparing phylogenetic tree ----\n\n"
    Biodiverse_00_create_bts.pl \
      --input_tree_file ${tree} \
      --out_file "tree.bts"

    ## Run the analyses
    echo -e "\n\n---- Running Biodiverse analyses ----\n\n"
    Biodiverse_02_biodiverse_analyses.pl \
      --input_bds_file "occ.bds" \
      --input_bts_file "tree.bts" \
      --calcs ${params.bd_indices}

    """
}



// Estimate phylogenetic diversity with Biodiverse
process phylodiv {

    label "container_biodiverse"

   // cpus 1

    input:
      path(BDA)
      val(chunkid)

    output:
      path "Biodiv_randomized_${chunkid}.bds", emit: BDArand

    script:
    """
    Biodiverse_03_run_randomisation.pl \
      --basedata ${BDA} \
      --bd_name  ${BDA} \
      --out_file "Biodiv_randomized.bds" \
      --rand_name 'rand' \
      --iterations ${iterations_per_thread} \
      --args ${biodiverse_args} \
      seed=${chunkid}

    ## Add chunk ID into the file name
    mv "Biodiv_randomized.bds" "Biodiv_randomized_${chunkid}.bds"

    """
}



// Prepare a shapefile to spatially constrain randomizations
process prep_shapefile {

    label "container_phylotwin"

    // cpus 1

    input:
      path(occurrences)
      path(polygons)

    output:
      path "shapefile*", emit: shapefile

    script:
    """
    Biodiverse_12_Prepare_SpatialConstraints.R \
      --input ${occurrences} \
      --randconstrain ${polygons} \
      --threads ${task.cpus} \
      --output "shapefile"
    """
}


// Estimate phylogenetic diversity using spatially-constrained randomizations
process phylodiv_constrianed {

    label "container_biodiverse"

   // cpus 1

    input:
      path(BDA)
      path(polygons)
      val(chunkid)

    output:
      path "Biodiv_randomized_${chunkid}.bds", emit: BDArand

    script:
    """
    Biodiverse_03_run_randomisation.pl \
      --basedata ${BDA} \
      --bd_name  ${BDA} \
      --out_file "Biodiv_randomized.bds" \
      --rand_name 'rand' \
      --iterations ${iterations_per_thread} \
      --args ${biodiverse_args} \
      seed=${chunkid} \
      spatial_conditions_for_subset='sp_points_in_same_poly_shape (file => "shapefile.shp")'

    ## Add chunk ID into the file name
    mv "Biodiv_randomized.bds" "Biodiv_randomized_${chunkid}.bds"

    """
}


// Create a file with paths to all chunks with randomization results
process rand_filelist {

    // container image is required for Cloud only
    // label "container_phylotwin"

    input:
      path(randfiles)

      //// To avoid name collisions if analysis was done for multiple datasets
      // path(randfiles, stageAs: "?/*")

    output:
      path "randomization_results.txt", emit: RND

    shell:
    $/
    echo "${randfiles}" \
      | sed -z 's/, /\n/g; s/^\[//; s/\]//' \
      > randomization_results.txt
    /$
}
// Groovy allows an alternative syntax for string definitions
// which uses the $ as escape character in place of \ character.
// These strings are delimited with an opening $/ and and a closing /$



// Aggregate the randomization results - with Biodiverse
process aggregate_rnds_biodiv {

    label "container_biodiverse"

   // cpus 1

    input:
      path RND
      path BDArand

    output:
      path "Biodiverse.bds", emit: Biodiv

    script:
    """
    Biodiverse_05_reintegrate_basedatas_post_rand.pl \
      --glob ${RND} \
      --output_prefix Biodiverse
    """
}


// Export Biodiverse results into CSV
process div_to_csv {

    label "container_biodiverse"

    // cpus 1

    input:
      path Biodiv

    output:
      path "bd_out/*.csv", emit: csvs

    script:
    """
    Biodiverse_04_load_bds_and_export_results.pl \
      --input_bds_file ${Biodiv} \
      --output_csv_prefix 'RND'

    mkdir -p bd_out
    mv RND_*.csv bd_out/

    """
}


// Merge Biodiverse results into a single table
process merge_biodiverse_results {

    label "container_phylotwin"

    // cpus 1

    input:
      path inps

    output:
      path "Biodiverse_results.txt", emit: bdres

    script:
    """
    biodiverse_combine_tabs.R \
      --inpdir . \
      --prefix RND \
      --resolution ${params.resolution} \
      --output Biodiverse_results.txt

    """
}



workflow BIODIVERSE {
    take:
        occurrences
        tree

    main:

      // Prepare Biodiverse input files
      prep_biodiv(occurrences, tree)

      // Channel of randomization chunks
      rnd_ch = Channel.fromList( randomization_chunks )


    // If no spatial constraints (for randomization) are provided
    if(params.randconstrain == null){

      // Perform unconstrained randomizations
      phylodiv(prep_biodiv.out.BDA, rnd_ch)

      // Prepare a file with paths of `phylodiv` output (multiple chunks of randomizations)
      // and create a new channel from it
      rand_filelist(phylodiv.out.BDArand.collect())

      // Aggregate randomization results (with Biodiverse script)
      aggregate_rnds_biodiv(
          rand_filelist.out.RND,
          phylodiv.out.BDArand.collect())

    } else {

      // If spatial constraints are provided, split dataset in parts (for each polygon)

      // A channel with spatial polygons
      polygons = file(params.randconstrain)

      // Prepare a shapefile with polygons
      prep_shapefile(merge_occ.out.occurrences, polygons)

      // Run spatially-constrained randomizations
      phylodiv_constrianed(
        prep_biodiv.out.BDA,
        prep_shapefile.out.shapefile,
        rnd_ch)

      // Collect `phylodiv_constrianed` output (multiple chunks of randomizations)
      rand_filelist(phylodiv_constrianed.out.BDArand.collect())

      // Aggregate randomization results (with Biodiverse script)
      aggregate_rnds_biodiv(
          rand_filelist.out.RND,
          phylodiv_constrianed.out.BDArand.collect())

    } // end of randomizations


    // // Split occurrences by polygons and run randomizations independently
    // polygons = file(params.randconstrain)
    //
    // // Split dataset
    // split_by_polygons(merge_occ.out.occurrences, polygons)
    //
    // // Channel with spatially-constrained datasets
    // ch_spatconstr = split_by_polygons.out.occsplit.flatten()
    // ch_spatconstr.view()
    //
    // // Prepare data for Biodiverse
    // prep_biodiv(ch_spatconstr, merge_occ.out.tree)
    //
    // // Channel with the number of randomization chunks
    // rnd_ch_tmp = Channel.fromList( randomization_chunks )
    //
    // // Apply randomization chunks for each spatially-contrained dataset (cartesian product)
    // rnd_ch = prep_biodiv.out.BDA.combine(rnd_ch_tmp)
    // rnd_ch.view()


    // Output results as CSV
    div_to_csv(aggregate_rnds_biodiv.out.Biodiv)

    emit:
        bdres = merge_biodiverse_results.out.bdres

} // end of sub-workflow

