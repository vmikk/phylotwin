#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Biodiverse::Tree;

# Check if a file argument was provided
die "Usage: $0 <tree_file.bts>\n" unless @ARGV == 1;

my $tree_file = $ARGV[0];

# Check if file exists
die "Error: File '$tree_file' does not exist\n" unless -e $tree_file;

my $tree = Biodiverse::Tree->new(file => $tree_file);

# Print the entire tree structure
print Dumper($tree);

# Print node names
print "\nTree nodes:\n";
my $nodes = $tree->get_node_refs;
foreach my $node (@$nodes) {
    print $node->get_name . "\n";
}

