#!/usr/bin/perl

while (<>) {
    ($weight,$onode,$priprob,$nodeinfo,$rnode,$trainseq)=split;
    print "$nodeinfo $rnode $trainseq 0 0\n";
}
