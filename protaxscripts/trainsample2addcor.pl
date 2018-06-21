#!/usr/bin/perl

if (scalar(@ARGV) < 3) {
    die "usage: trainsample2addcor.pl classification_file seq2taxname_file taxonomy_file\n";
}

$infile=shift;
$seq2taxfile=shift;
$taxonomyfile=shift;

open(FD,$taxonomyfile);
while (<FD>) {
    ($id,$pid,$level,$name)=split;
    $node2tname{$id}=$name;
    $tname2node{$name}=$id;
}
close(FD);

open(FD,$seq2taxfile);
while (<FD>) {
    ($foo1,$foo2,$seqid,$tname)=split;
    push(@seqids,$seqid);
    push(@cortnames,$tname);
}
close(FD);

open(FD,$infile);
$i=0;
while (<FD>) {
    ($foo1,$foo2,$seqid,$node,$logprob)=split;
    if ($node eq "") {
	$prob = 0;
    }
    else {
	$prob = exp($logprob);
    }

    if ($seqid ne $seqids[$i]) {
	die "two files not in syncrony ($seqids[$i] vs $seqid), line $i+1.\n";
    }
    # node may be e.g 4 or 4,unk
    if ($node eq "") {
	$tname = "";
    }
    else {
	$foo=$node;
	$foo =~ s/,unk//;
	if (!exists($node2tname{$foo})) {
	    die "ERROR: cannot find taxname for node ($foo,$node), seqid '$seqid'.\n";
	}
	$tname = $node2tname{$foo};
	if ($node =~ /unk/) {
	    $tname = $tname . ',unk';
	}
    }
    
    if ($tname eq $cortnames[$i]) {
	$cor=1;
    }
    else {
	$cor=0;
    }
    if ($node eq "") {
	$node = NA;
	$tname = NA;
    }
    
    print "$seqid\t$prob\t$cor\t$node\t$tname\tNA\t$cortnames[$i]\n";
    $i++;
}
close(FD);
