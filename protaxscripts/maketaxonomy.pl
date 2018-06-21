#!/usr/bin/perl

print "0\t0\t0\troot\n";
$count=1;

while (<>) {
    ($k,$p,$c,$o,$f,$g,$s)=split;
    # fill NA with dummy names
    if ($s ne "NA") {
	$specific=$s;
	if ($g eq "NA") {
	    $g = "dummy_genus_$specific";
	}
	else {
	    $specific = $g;
	}
	if ($f eq "NA") {
	    $f = "dummy_family_$specific";
	}
	else {
	    $specific = $f;
	}	
	if ($o eq "NA") {
	    $o = "dummy_order_$specific";
	}
    }
    if (!exists($nid{$o})) {
	$nid{$o}=$count;
	print "$nid{$o}\t0\t1\t$o\n";
	$count++;
    }
    $cur = "$o,$f";
    if (!exists($nid{$cur})) {
	$nid{$cur}=$count;
	print "$nid{$cur}\t$nid{$o}\t2\t$cur\n";
	$count++;
    }
    $parent = $cur;
    $cur = "$o,$f,$g";
    if (!exists($nid{$cur})) {
	$nid{$cur}=$count;
	print "$nid{$cur}\t$nid{$parent}\t3\t$cur\n";
	$count++;
    }
    $parent = $cur;
    $cur = "$o,$f,$g,$s";
    if (!exists($nid{$cur})) {
	$nid{$cur}=$count;
	print "$nid{$cur}\t$nid{$parent}\t4\t$cur\n";
	$count++;
    }
}
