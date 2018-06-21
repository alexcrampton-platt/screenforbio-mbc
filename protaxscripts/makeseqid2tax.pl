#!/usr/bin/perl

$level=shift;

# use full length names for taxa
while (<>) {
    ($id,$nimi) = split;
    @a=split(/,/,$nimi);
    if (scalar(@a) >= $level) {
	$tax=join(',',@a[0 .. ($level-1)]);
	print "$id\t$tax\n";
    }
}
