#!/usr/bin/perl

while (<>) {
    if (!(/^#/)) {
	($score, $name1, $name2, $alnSize2) = (split)[0,1,6,8];
	$f = $score/$alnSize2;
	$key = "$name2\t$name1";
	if (!exists($dat{$key}) or ($f > $dat{$key})) {
	    $dat{$key} = $f;
	}
    }
}

foreach $i (sort keys %dat) {
    $f = sprintf("%.3f",$dat{$i});
    print "$i\t$f\n";
}
