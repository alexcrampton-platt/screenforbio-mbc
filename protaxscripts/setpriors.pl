#!/usr/bin/perl

$usage= "setpriors.pl proportion seleced_species taxonomy_file";

$proportion=shift;
if (($proportion <= 0) or ($proportion >= 1)) {
    die "ERROR: proportion should be between 0..1 (was $proportion).\nusage: $usage\n";
}

$file=shift;
open(FD,$file);
while (<FD>) {
    ($nid)=split;
    $select{$nid}=1;
}
close(FD);

$num_selected=0;
$maxlevel=0;
%num_nodes = ();
while (<>) {
    s/\s+$//;
    ($id,$pid,$level,$name,$prior) = split; 
    next if ($id eq ""); # skip empty lines
    $num_selected++ if (exists($select{$id}));
    $dat{$id} = "$id\t$pid\t$level\t$name";
    push(@nids,$id);
    push(@{$node2children{$pid}},$id) if ($id != $pid);
    $node2level{$id} = $level;
    $maxlevel = $level if ($level > $maxlevel);
    $num_nodes{$level}++;
} 

# equal prior for species nodes
# if 90% for selected, 10% for other nodes (if $proportion = 0.9)

$p1 = $proportion/$num_selected;
$p0 = (1-$proportion)/($num_nodes{$maxlevel}-$num_selected);

foreach $nid (sort {$b <=> $a} keys %node2level) {
    if ($node2level{$nid} == $maxlevel) {
	if (exists($select{$nid})) {
	    $nodeprob{$nid} = $p1;
	}
	else {
	    $nodeprob{$nid} = $p0;
	}
    }
    else {
	# sum probabilities of child nodes
	$sumprob=0;
	foreach $cnid (@{$node2children{$nid}}) {
	    $sumprob += $nodeprob{$cnid};
	}
	$nodeprob{$nid} = $sumprob;
    }
}

foreach $nid (@nids) {
    $f = sprintf("%.6g",$nodeprob{$nid});
    print "$dat{$nid}\t$f\n";
}
