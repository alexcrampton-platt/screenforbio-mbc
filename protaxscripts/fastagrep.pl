#!/usr/bin/perl

$file=shift;
open(FD,$file);
while (<FD>) {
    ($id)=split;
    $ota{$id}=1;
}
close(FD);

$ok=0;
while (<>) {
    if (/^>/) {
	($id)=split;
	$id =~ s/^>//;
	if (exists($ota{$id})) {
	    $ok=1;
	}
	else {
	    $ok=0;
	}
    }

    print if ($ok);
}
