#!/usr/bin/perl -w 
# written by Lingbin
# 只保留 cis(同染色体)接触,删除 trans,并过滤近距离 cis。
# 输出4种距离阈值(过滤相邻 1/2/3/4 个窗口 = 25/50/75/100 kb)的每窗口接触数。
# 用法: perl this.pl <input_file> <output_prefix>
use strict;
use warnings;

if (@ARGV != 2) { die "Usage: $0 <input_file> <output_prefix>\n"; }
my ($input_file, $prefix) = @ARGV;

my @thresholds = (25000, 50000, 75000, 100000);
my @labels     = ("25kb", "50kb", "75kb", "100kb");

my @counts = ({}, {}, {}, {});           # 每种阈值一个哈希

open my $fh, '<', $input_file or die "Could not open '$input_file': $!\n";
while (<$fh>) {
    chomp;
    my @f = split("\t");
    next if @f < 7;
    my ($chr1,$s1,$e1, $chr2,$s2,$e2, $count) = @f[0,1,2,3,4,5,6];

    # 删除 trans:不同染色体的接触直接跳过
    next if ($chr1 ne $chr2);

    # 到这里都是 cis(同染色体)
    my $pos1 = "$chr1\t$s1\t$e1";
    my $pos2 = "$chr2\t$s2\t$e2";
    my $dist = abs($s1 - $s2);            # 起点距离

    for my $i (0 .. $#thresholds) {
        # 过滤近距离 cis:距离 <= 阈值 -> 跳过
        next if ($dist <= $thresholds[$i]);

        # 保留:累加两端
        $counts[$i]{$pos1} += $count;
        if ($pos1 ne $pos2) {            # 距离过滤后 pos1 必然 != pos2,这里是保险
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
