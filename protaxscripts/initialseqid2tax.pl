#!/usr/bin/perl

$taxfile=shift;
open(FD,$taxfile);
while (<FD>) {
    ($nid,$pid,$level,$name)=split;
    if ($level == 4) {
	($o,$f,$g,$s)=split(/,/,$name);
	$k= $g . "," . $s;
	$short2full{$k} = $name;
    }
}
close(FD);

# all fasta id names contain 3 parts genus_species_numid
while (<>) {
    if (/^>/) {
	($id)=split;
	$id =~s/^>//;
	($g,$s,$numid)=split(/_/,$id);
	$k= $g . "," . $s;
	if (!exists($short2full{$k})) {
	    die "ERROR: cannot find full taxon name for '$id'.\n";
	}
	else {
	    print "$id\t$short2full{$k}\n";
	}
    }
}
