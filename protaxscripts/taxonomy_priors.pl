#!/usr/bin/perl

$maxlevel=0;
%num_nodes = ();
while (<>) {
    s/\s+$//;
    ($id,$pid,$level,$name,$prior) = split; 
    next if ($id eq ""); # skip empty lines
    $dat{$id} = $_;
    push(@nids,$id);
    push(@{$node2children{$pid}},$id) if ($id != $pid);
    $node2level{$id} = $level;
    $maxlevel = $level if ($level > $maxlevel);
    $num_nodes{$level}++;
} 

# equal prior for species nodes

foreach $nid (sort {$b <=> $a} keys %node2level) {
    if ($node2level{$nid} == $maxlevel) {
	$nodeprob{$nid} = 1/$num_nodes{$maxlevel};
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
    $f = sprintf("%.10f",$nodeprob{$nid});
    print "$dat{$nid}\t$f\n";
}
