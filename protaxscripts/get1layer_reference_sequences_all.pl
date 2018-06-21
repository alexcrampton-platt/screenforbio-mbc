#!/usr/bin/perl

use warnings;
use strict;

sub readTaxonomyFile {
    # taxonomy file format: id(integer) pid(integer) level(integer) name(string)
    
    my ($file) = @_;
    my %node2parent = ();
    my %node2children = ();
    my %node2level = ();    
    my %name2node = ();    
    my ($id, $pid, $level, $name);
    my $maxid=0;
    open(FD,$file) or die "ERROR (readTaxonomyFile): cannot read file '$file'. $!\n";
    while (<FD>) {
	($id,$pid,$level,$name) = split;
	next if ($id eq ""); # skip empty lines
	$node2parent{$id} = $pid;
	push(@{$node2children{$pid}},$id) if ($id != $pid);
	$node2level{$id} = $level;
	$name2node{$name} = $id;
	$maxid = $id if ($id > $maxid);
    } 
    close(FD);
    return (\%node2parent,\%node2children,\%node2level, \%name2node, $maxid);
}

sub readSequence2TaxonomyFile {
    my ($file, $taxname2node) = @_;
    my %seqid2taxonomy = ();
    my %node2seqids = ();

    my ($seqid,$taxname);
    open(FD, $file) or die "ERROR (readSequence2TaxonomyFile): cannot read file '$file'. $!\n";
    while (<FD>) {
	($seqid,$taxname) = split;
	next if ($seqid eq ""); # skip empty lines
	die "ERROR (readSequence2TaxonomyFile): cannot find taxonomy node for seqid: '$seqid' taxname: '$taxname'.\n" if (!exists($taxname2node->{$taxname}));
	$seqid2taxonomy{$seqid} = $taxname2node->{$taxname};
	push(@{$node2seqids{ $taxname2node->{$taxname} }}, $seqid);
    }
    close(FD);

    return (\%seqid2taxonomy,\%node2seqids);
}

sub printTaxonomy {
    my ($n2pref, $n2cref, $n2lref, $n2seqs, $unknoderef) = @_;
#    foreach my $i (sort {$n2lref->{$a} <=> $n2lref->{$b}} keys %$n2lref) {
    foreach my $i (sort {$n2lref->{$a} <=> $n2lref->{$b} || $a <=> $b} keys %$n2lref) {
	print "LEVEL $n2lref->{$i}: '$i' parent: $n2pref->{$i}";
	if (exists($unknoderef->{$i})) {
	    print ", this is UNK node";
	}
	else {
	    if (exists($n2cref->{$i})) {
		print " child nodes: @{$n2cref->{$i}}";
	    }
	    else {
		print " child nodes: -";
	    }
	    if (exists($n2seqs->{$i})) {
		print ", seqids: @{$n2seqs->{$i}}";
	    }    
	}
	print "\n";
    }
    return (0);
}

sub arrayPair {
    my ($aref,$wref) = @_;
    my @a = @$aref;
    my @w = @$wref;
    return (\@a,\@w);
}

sub allRseqs {
    my ($n2lref,$n2cref,$n2seqref, $debug_mode) = @_;
    my %node2rseq = ();
    my $rootnode = 0;
    my ($node, $cnode, $count);
    my (@rid_array,@weight_array,$wtext);

    # take rseqs both from direct seqs and from children's rseqs, proceed from leaves towards the root
    # weighting is 1/num_rseqs_of_cnode -> equal weighting of cnodes independent of differences in num_rseqs
    foreach $node (sort {$n2lref->{$b} <=> $n2lref->{$a}} keys %$n2lref) {
	next if ($node == $rootnode);	

	@rid_array = ();
	@weight_array = ();

	if (exists($n2cref->{$node})) {
	    foreach $cnode (@{$n2cref->{$node}}) {
		if (exists($node2rseq{$cnode})) {
		    $count = scalar(@{$node2rseq{$cnode}->[0]});
		    print "level $n2lref->{$node} node '$node' cnode '$cnode' $count child rseqs @{$node2rseq{$cnode}->[0]}\n" if ($debug_mode > 1);
		    push(@rid_array, @{$node2rseq{$cnode}->[0]});
		    $wtext = sprintf("%.3g",1/$count);
		    push(@weight_array, ($wtext) x $count);
		}
	    }
	}

	if (exists($n2seqref->{$node})) {
	    # take rseqs labeled exactly to this node
	    $count = scalar(@{$n2seqref->{$node}});
	    print "level $n2lref->{$node} node '$node' $count direct rseqs @{$n2seqref->{$node}}\n" if ($debug_mode > 1);
	    push(@rid_array, @{$n2seqref->{$node}});
	    $wtext = sprintf("%.3g",1/$count);
	    push(@weight_array, ($wtext) x $count);	    
	}

	if (scalar(@rid_array) > 0) {
	    @{$node2rseq{$node}} = arrayPair(\@rid_array, \@weight_array);	
	}
    }
    return (\%node2rseq);
}

sub writeRseqsOneLayer {
    my ($n2rid, $rseq_separator, $outfile, $n2lref, $target_layer) = @_;

    my ($node, $s, $w, @rseqs,@weights, $i);
    open (FD, ">$outfile") or die "ERROR (writeRseqs): cannot write to file '$outfile'. $!\n";
    foreach $node (sort {$a <=> $b} keys %$n2rid) {
	if ($n2lref->{$node} == $target_layer) {
	    $s = join($rseq_separator, @{$n2rid->{$node}->[0]});
	    $w = join($rseq_separator, @{$n2rid->{$node}->[1]});
	    print FD "$node\t$s\t$w\n";
	}
    }
    close (FD);
    
    return (0);
}

my $usage = "usage: get1layer_reference_sequences_all.pl layer taxonomyfile seqid2taxonomy_file outfile (debug_mode)\n";

if (scalar(@ARGV) < 4) {
    die $usage;
}

my $target_layer = shift;
my $taxonomyfile=shift;
my $seqid2taxidfile=shift;
my $outfile = shift;
die $usage if (!defined($outfile));
my $debug_mode = shift;
$debug_mode = 0 if (!defined($debug_mode));

my ($n2pref,$n2cref,$n2lref,$tname2nref)=readTaxonomyFile($taxonomyfile);
my ($s2tref, $n2sidref) = readSequence2TaxonomyFile($seqid2taxidfile,$tname2nref);
printTaxonomy($n2pref, $n2cref, $n2lref, $n2sidref, {}) if ($debug_mode);

my ($n2ridref) = allRseqs($n2lref,$n2cref,$n2sidref, $debug_mode);
writeRseqsOneLayer($n2ridref, ',', $outfile, $n2lref, $target_layer);
