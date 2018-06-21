#!/usr/bin/perl

if (scalar(@ARGV < 2)) {
    die "usage: add_taxonomy_info.pl taxonomy_file classification_result_file\n";
}

$taxfile=shift;
$cfile=shift;

open(FD,$taxfile) or die "ERROR: cannot read taxonomy file '$taxfile'. $!\n";
while (<FD>) {
    ($id,$pid,$level,$name)=split;
    $taxinfo{$id} = "$level\t$name";
}
close(FD);

open(FD,$cfile) or die "ERROR: cannot read classification output file '$cfile'. $!\n";
while (<FD>) {
    ($seqid,$taxid,$p)=split;
    ($id,$unk) = split(/,/,$taxid);
    if ($unk eq 'unk') {
	($level,$taxname)=split(/\t/,$taxinfo{$id});
	$level++;
	print "$seqid\t$taxid\t$p\t$level\t$taxname,unk\n";
    }
    else {
	print "$seqid\t$taxid\t$p\t$taxinfo{$id}\n";
    }
}
close(FD);
