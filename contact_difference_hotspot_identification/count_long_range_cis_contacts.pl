#!/usr/bin/perl -w
# written by Lingbin
# Count long-range cis contacts for each genomic window.
# Usage: perl count_long_range_cis_contacts.pl <input_file> <output_prefix>
# Input: uncompressed, headerless, tab-delimited cooler dump --join output:
#   chr1, start1, end1, chr2, start2, end2, count
# Supply raw integer counts with each matrix pixel represented once.
# Keep cis pixels with abs(start1 - start2) strictly greater than each threshold:
#   25,000 / 50,000 / 75,000 / 100,000 bp, measured between bin starts.
# Add each retained pixel count to both endpoint windows.
# Outputs: <output_prefix>.filter25kb.position, .filter50kb.position,
#          .filter75kb.position and .filter100kb.position.
# Each output has four columns without a header: chr, start, end, count.
# Only windows encountered in retained pixels are written; coordinates are unchanged.
# This script can process either specific or total contact counts.
# To match call_contact_difference_hotspots.R, use an output prefix ending in:
#   <sample>.SelfMapping.<hap>.HPRC.25kb.Unique.position  (specific counts)
#   <sample>.SelfMapping.<hap>.HPRC.25kb.Total.position   (total counts)
# Place each output in its corresponding counts root under <sample>/.
# Specific and total inputs must use the same sample/haplotype, bins and filters.
# Example after creating /path/to/specific_contact_counts/HG00000/:
# perl /path/to/count_long_range_cis_contacts.pl /path/to/HG00000.SelfMapping.Hap1.HPRC.25kb.Unique.bed /path/to/specific_contact_counts/HG00000/HG00000.SelfMapping.Hap1.HPRC.25kb.Unique.position
use strict;
use warnings;

if (@ARGV != 2) { die "Usage: $0 <input_file> <output_prefix>\n"; }
my ($input_file, $prefix) = @ARGV;

my @thresholds = (25000, 50000, 75000, 100000);
my @labels     = ("25kb", "50kb", "75kb", "100kb");

my @counts = ({}, {}, {}, {});           # One count hash per threshold.

open my $fh, '<', $input_file or die "Could not open '$input_file': $!\n";
while (<$fh>) {
    chomp;
    my @f = split("\t");
    next if @f < 7;
    my ($chr1,$s1,$e1, $chr2,$s2,$e2, $count) = @f[0,1,2,3,4,5,6];

    # Keep contacts on the same chromosome.
    next if ($chr1 ne $chr2);

    # Identify both endpoint windows.
    my $pos1 = "$chr1\t$s1\t$e1";
    my $pos2 = "$chr2\t$s2\t$e2";
    my $dist = abs($s1 - $s2);            # Distance between bin starts.

    for my $i (0 .. $#thresholds) {
        # Retain distances strictly greater than the threshold.
        next if ($dist <= $thresholds[$i]);

        # Add the contact count to both endpoints.
        $counts[$i]{$pos1} += $count;
        if ($pos1 ne $pos2) {            # Count a distinct second endpoint.
            $counts[$i]{$pos2} += $count;
        }
    }
}
close $fh;

for my $i (0 .. $#thresholds) {
    my $outfile = "${prefix}.filter${labels[$i]}.position";
    open my $out, '>', $outfile or die "Cannot write $outfile: $!\n";
    foreach my $position (sort keys %{$counts[$i]}) {
        print $out "$position\t$counts[$i]{$position}\n";
    }
    close $out;
}
