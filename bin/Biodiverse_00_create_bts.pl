#!/usr/bin/env perl

## Script to prepare phylogenetic tree for Biodiverse

## Parameters:
# --input_tree_file  = Phylogenetic tree in Newick format (can be gzipped)
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
  [ 'input_tree_file=s',   'The input tree file in Newick format', { required => 1 } ],
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
my $newick_content;
if ($tree_inp_file =~ /\.gz$/) {
    my $zh = IO::Zlib->new($tree_inp_file, "rb") 
        or die "Cannot open gzipped file $tree_inp_file: $!";
    local $/;  # enable slurp mode
    $newick_content = <$zh>;
    $zh->close;
} else {
    open my $fh, '<', $tree_inp_file 
        or die "Cannot open file $tree_inp_file: $!";
    local $/;  # enable slurp mode
    $newick_content = <$fh>;
    close $fh;
}

# Create a ReadNexus object which can also handle Newick format
my $read_nex = Biodiverse::ReadNexus->new();

# Import the Newick tree data from the content
my $success = $read_nex->import_newick(
    data => $newick_content
);

# Get the tree object
my $tree_array = $read_nex->get_tree_array;
my $tree = $tree_array->[0];  # Get first tree if multiple trees exist

$tree->save (filename => $tree_out_file);

print "Tree saved to $tree_out_file\n";
