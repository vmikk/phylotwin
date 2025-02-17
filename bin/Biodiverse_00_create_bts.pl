#!/usr/bin/env perl

## Script to prepare phylogenetic tree for Biodiverse

## Parameters:
# --input_tree_file  = Phylogenetic tree in Newick or Nexusformat (can be gzipped)
# --out_file         = Output file (bts)

## Based on `create_bts.pl` (4c70164, 2022-09-22)
# https://github.com/vmikk/PhyloNext/blob/main/bin/00_create_bts.pl


use strict;
use warnings;
use Carp;        #  warnings and dropouts
use File::Spec;  #  for the cat_file sub
use English qw ( -no_match_vars );
use IO::Zlib;    # for gzip support

use Biodiverse::BaseData;
use Biodiverse::ElementProperties;  #  for remaps
use Biodiverse::ReadNexus;
use Biodiverse::Tree;

use Getopt::Long::Descriptive;

$| = 1;

my ($opt, $usage) = describe_options(
  '%c <arguments>',
  [ 'input_tree_file=s',   'The input tree file in Newick or Nexus format', { required => 1 } ],
  [ 'out_file=s',  'The output biodiverse tree file .bts', { required => 1 }],
  [],
  [ 'help',       "print usage message and exit" ],
);

 
if ($opt->help) {
    print($usage->text);
    exit;
}

print "Preparing phylogenetic tree for Biodiverse\n";

my $tree_inp_file = $opt->input_tree_file;
my $tree_out_file = $opt->out_file;

# Read the input file (auto-detecting gzip compression)
my $is_gzipped = $tree_inp_file =~ /\.gz$/i;
my $base_filename = $is_gzipped ? $tree_inp_file =~ s/\.gz$//r : $tree_inp_file;
my $is_nexus = $base_filename =~ /\.(nex|nexus)$/i;
my $is_newick = $base_filename =~ /\.(tre|nwk|newick)$/i;

die "Unrecognized file extension. Expected .nex, .nexus for Nexus format or .tre, .nwk, .newick for Newick format (case insensitive, can be .gz compressed)\n"
    unless $is_nexus || $is_newick;

# Create a ReadNexus object
my $read_nex = Biodiverse::ReadNexus->new();

if ($is_nexus) {
    # For Nexus files
    if ($is_gzipped) {
        my $content = '';
        my $zh = IO::Zlib->new($tree_inp_file, "rb") 
            or die "Cannot open gzipped file $tree_inp_file: $!";
        while (my $line = $zh->getline()) {
            $content .= $line;
        }
        $zh->close;
        
        # Write to temporary file
        my $temp_file = $tree_inp_file . "_temp";
        open my $temp_fh, '>', $temp_file or die "Cannot create temporary file: $!";
        print $temp_fh $content;
        close $temp_fh;
        
        # Import from temporary file
        $read_nex->import_data(
            file => $temp_file,
            use_element_properties => 0,
        );
        
        # Clean up
        unlink $temp_file;
    } else {
        $read_nex->import_data(
            file => $tree_inp_file,
            use_element_properties => 0,
        );
    }
} else {
    # For Newick files
    my $content = '';
    if ($is_gzipped) {
        my $zh = IO::Zlib->new($tree_inp_file, "rb") 
            or die "Cannot open gzipped file $tree_inp_file: $!";
        while (my $line = $zh->getline()) {
            $content .= $line;
        }
        $zh->close;
    } else {
        open my $fh, '<', $tree_inp_file 
            or die "Cannot open file $tree_inp_file: $!";
        while (my $line = <$fh>) {
            $content .= $line;
        }
        close $fh;
    }
    
    # Import the Newick tree data from the content
    my $success = $read_nex->import_newick(
        data => $content
    );
}

# Get the tree object
my $tree_array = $read_nex->get_tree_array;
my $tree = $tree_array->[0];  # Get first tree if multiple trees exist

$tree->save (filename => $tree_out_file);

print "Tree saved to $tree_out_file\n";
