#!/usr/bin/perl

use warnings;
use strict;

my $MISSING_IS_ZERO=1;

sub readTaxonomyFile {
    # taxonomy file format: id(integer) pid(integer) level(integer) name(string) prior(float)
    
    my ($file) = @_;
    my %node2parent = ();
    my %node2children = ();
    my %node2level = ();    
    my %name2node = ();    
    my ($id, $pid, $level, $name, $prior);
    my $maxid=0;
    my %node2prior = ();
    open(FD,$file) or die "ERROR (readTaxonomyFile): cannot read file '$file'. $!\n";
    while (<FD>) {
	($id,$pid,$level,$name,$prior) = split;	
	next if ($id eq ""); # skip empty lines
	$node2parent{$id} = $pid;
	push(@{$node2children{$pid}},$id) if ($id != $pid);
	$node2level{$id} = $level;
	$name2node{$name} = $id;
	$maxid = $id if ($id > $maxid);
	die "ERROR (readTaxonomyFile): cannot find prior (5th column) file '$file'. $!\n" if ($prior eq "");
	$node2prior{$id} = $prior;
    } 
    close(FD);
    return (\%node2parent,\%node2children,\%node2level, \%name2node, $maxid, \%node2prior);
}

sub addMissingTaxonomy {
    my ($n2pref,$n2cref,$n2lref, $maxid) = @_;

    # add one extra child for each node (missing 'unknown' taxonomy branch) except leaf nodes

    my %unknodes = ();
    my ($node, $newnode);
    foreach $node (keys %$n2cref) {
	$maxid++;
	$newnode = $maxid;
	$n2lref->{$newnode} = $n2lref->{$node} + 1;
	$n2pref->{$newnode} = $node;
	$unknodes{$newnode} = 1;
    }
    foreach $newnode (keys %unknodes) {
	$node = $n2pref->{$newnode};
	push(@{$n2cref->{$node}},$newnode);
    }
    return (\%unknodes);
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

sub min {
    return ((sort {$a <=> $b} @_)[0]);
}

sub childIndex {
    my ($node, $n2pref, $n2cref) = @_;

    # find index of $node in its parents child array, return -1 if not found

    my $pnode = $n2pref->{$node};
    my $ind = 0;
    my $found = -1;
    foreach my $i (@{$n2cref->{$pnode}}) {
	$found=$ind, last if ($node eq $i);
	$ind++;
    }
    
    return ($found);
}

sub simMeanMax {
    my ($qseqid, $tseqidref, $tseqweightref, $seqsimref, $seq_dont_use_ref) = @_;
    my $sim=0;
    my $max=0;
    my $best2 = 0;
    foreach my $tseqid (@$tseqidref) {
	next if (exists($seq_dont_use_ref->{$tseqid}));
	if (exists($seqsimref->{$qseqid}->{$tseqid})) {
	    $sim = $seqsimref->{$qseqid}->{$tseqid};
	    if ($sim > $max) {
		$best2 = $max;
		$max = $sim;
	    }
	    elsif ($sim > $best2) {
		$best2 = $sim;
	    }
	}
	elsif (exists($seqsimref->{$tseqid}->{$qseqid})) {
	    $sim = $seqsimref->{$tseqid}->{$qseqid};
	    if ($sim > $max) {
		$best2 = $max;
		$max = $sim;
	    }
	    elsif ($sim > $best2) {
		$best2 = $sim;
	    }
	}
    }
    return ($best2, $max);
}

sub createX {
    my ($node, $n2cref, $n2rseqref, $n2rseqweightref, $seq_dont_use_ref, $qseqid, $remove_i, $unknoderef, $seqsimref, $n2prior, $debug_mode) = @_;

    my ($i,$rseqid,$cnode,$sumprob);

    my %cnode_count = ();
    $sumprob=0;
    for ($i=0; $i<(scalar(@{$n2cref->{$node}})); $i++) {
	$cnode = $n2cref->{$node}->[$i];
	if (($i != $remove_i) and (!exists($unknoderef->{$cnode}))) {
	    $sumprob += $n2prior->{$cnode};
	}
	next if (($i == $remove_i) || exists($unknoderef->{$cnode}) || (!exists($n2rseqref->{$cnode})));
	foreach $rseqid (@{$n2rseqref->{$cnode}}) {
	    $cnode_count{$i}++ if (!exists($seq_dont_use_ref->{$rseqid}));		
	}
    }
        
    # store all matrix X elements in one array row by row
    

    my @X = (0, 0, 0, 0);
    # unk prior is 0.026: x/(x+sum) = p -> x = p/(1-p) * sum, 0.026/(1-0.026) = 0.0267
    my @xweight = (sprintf("%.6g",0.0267*$sumprob));

    my ($mean,$max, $seqcov, $wref);
    my $use_weighting;
    for ($i=0; $i<(scalar(@{$n2cref->{$node}})); $i++) {
	$cnode = $n2cref->{$node}->[$i];
	next if (($i == $remove_i) || exists($unknoderef->{$cnode}));
	push(@xweight, $n2prior->{$cnode});

	if (exists($cnode_count{$i})) {
	    if (exists($n2rseqweightref->{$cnode})) {
		$wref = $n2rseqweightref->{$cnode};
	    }
	    else {
		$wref = [];
	    }

	    if ($debug_mode > 1) {
		my $j=0;
		my $winfo = "rseq weight: NONE";
		$use_weighting = 0;
		if (exists($n2rseqweightref->{$cnode})) {
		    $use_weighting = 1 if (scalar(@{$n2rseqref->{$cnode}}) == scalar(@{$n2rseqweightref->{$cnode}}));
		}

		foreach my $sid (@{$n2rseqref->{$cnode}}) {
		    $winfo = "rseq weight: $n2rseqweightref->{$cnode}->[$j]" if ($use_weighting);
		    $j++;
		    if (exists($seqsimref->{$qseqid}->{$sid})) {
			print " cnode '$cnode' seqsim('$qseqid','$sid') =  $seqsimref->{$qseqid}->{$sid} $winfo\n";
		    }
		    elsif (exists($seqsimref->{$sid}->{$qseqid})) {
			print " cnode '$cnode' seqsim('$qseqid','$sid') =  $seqsimref->{$sid}->{$qseqid} (symmetric) $winfo\n";
		    }
		    else {
			if (exists($seq_dont_use_ref->{$sid})) {
			    print " cnode '$cnode' seqsim('$qseqid','$sid') =  flagged dont_use rseq ($winfo)\n";
			}
			else {
			    print " cnode '$cnode' seqsim('$qseqid','$sid') =  doesn't exist ($winfo)\n";
			}
		    }
		}
	    } # end debug_mode > 1
    
	    # calculate mean and max similarity between query seq and child node's refseqs
	    
	    ($mean, $max) = simMeanMax($qseqid, $n2rseqref->{$cnode}, $wref, $seqsimref, $seq_dont_use_ref);
	    #$seqcov = log($cnode_count{$i});
	    if ($cnode_count{$i} > 1) {
		$seqcov = 1;
	    }
	    else {
		$seqcov = 0;
	    }
	    push(@X, (0, 1, sprintf("%.6g",$max), sprintf("%.6g",$max - $mean)));
	}
	else {
	    # known taxon without refseqs
	    push(@X, (1, 0, 0, 0));
	}
    }
    return (\@X, \@xweight);
}


sub readTaxonomyRseqIdFile {
    my ($filename, $rseq_separator) = @_;
    my %node2rseq = ();
    my %node2rseqweight = ();
    
    open(FD, $filename) or die "ERROR (readTaxonomyRseqIdFile): cannot read file '$filename'. $!\n";
    while (<FD>) {
	my ($node, $rseqs, $weights) = split;
	@{$node2rseq{$node}} = split(/$rseq_separator/, $rseqs);
	@{$node2rseqweight{$node}} = split(/$rseq_separator/, $weights) if ($weights ne "");
    }
    close(FD);
    
    return (\%node2rseq, \%node2rseqweight);
}

sub addRseqs2directSeqs {
    my ($noderef,$n2rseq, $n2seqref, $s2node) = @_;
    my ($node,$rseqid);
    foreach $node (@$noderef) {
	if (exists($n2rseq->{$node})) {
	    @{$n2seqref->{$node}} = () if (!exists($n2seqref->{$node}));
	    foreach $rseqid (@{$n2rseq->{$node}}) {
		# add rseq if it's not node's direct seq (i.e. already in node's sequence list)
		push (@{$n2seqref->{$node}}, $rseqid) if (!exists($s2node->{$rseqid}) or ($s2node->{$rseqid} ne $node));
	    }
	}
    }
    return (0);
}

sub readSeqSimilarityFile {
    my ($file) = @_;
    my %sim = ();
    my ($sid1,$sid2,$seqsim);

    open(FD,$file) or die "ERROR (readSeqSimilarityFile): cannot read file '$file'. $!\n";
    while (<FD>) {
	($sid1,$sid2,$seqsim) = split;
	$sim{$sid1}{$sid2} = $seqsim;
    }
    close(FD);

    return (\%sim);
}

if (scalar(@ARGV) < 6) {
    die "usage: create1layer_xdata.pl trainsample_file taxonomyfile seqid2taxonomy_file taxonomy2rseqid_file seqsimfile outfile (debug_mode 0,1,2)\n";
}

my $trainsamplefile = shift;
my $taxonomyfile=shift;
my $seqid2taxfile = shift;
my $tax2ridfile=shift;
my $seqsimfile=shift;
my $outfile = shift;
my $debug_mode = shift;
$debug_mode = 0 if (!defined($debug_mode));

print "MISSING SIMILARITIES treated as ZEROs: $MISSING_IS_ZERO\n";

my ($n2pref,$n2cref,$n2lref,$tname2nref,$maxid,$n2prior) = readTaxonomyFile($taxonomyfile);
my $unknoderef = addMissingTaxonomy($n2pref,$n2cref,$n2lref,$maxid);
my ($s2tref, $n2sidref) = readSequence2TaxonomyFile($seqid2taxfile,$tname2nref);
my ($n2ridref, $n2rweightref) = readTaxonomyRseqIdFile($tax2ridfile, ',');
my @allnodes = keys %$n2pref;
addRseqs2directSeqs(\@allnodes,$n2ridref, $n2sidref, $s2tref);

my $seqsimref = readSeqSimilarityFile($seqsimfile);
my @allxdat;
my @allxweight; # weights for X matrix rows (branches may have user-defined prior probabilities)

open(FD,"$trainsamplefile") or die "ERROR: cannot read trainsample_file '$trainsamplefile'. $!\n";
my $i=0;
my ($seqid,$yi, $remove_i,%seq_dont_use);
my ($weight,$onode,$priprob,$nodeinfo,$nodeclass,$onerefseq,$node,$trainseq);
while (<FD>) {
    ($weight,$onode,$priprob,$nodeinfo,$node,$trainseq) = split;
    $i++;
    print "$i): onode '$onode' prob $priprob" if ($debug_mode);
    $remove_i = -1;
    %seq_dont_use = ();
    ($nodeclass,$onerefseq) = split(/,/,$nodeinfo);

    if ($nodeclass eq "nseq2") {
	# node belongs to known taxonomy with refseqs, +1 because of unk row in X
	$yi = 1 + childIndex($node, $n2pref, $n2cref);
	print " #seqs >=2, trainseq '$trainseq'\n" if ($debug_mode);

	$seq_dont_use{$trainseq} = 1;
    }
    elsif ($nodeclass eq "nseq1") {
	# node represents known species with only 1 refseq, +1 because of unk row in X
	$yi = 1 + childIndex($node, $n2pref, $n2cref);
	print " #seqs ==1, rnode '$node' trainseq '$trainseq', refseq '$onerefseq'\n" if ($debug_mode);
	
	$seq_dont_use{$trainseq} = 1;
	foreach $seqid (@{$n2ridref->{$node}}) {
	    $seq_dont_use{$seqid} = 1 if ($seqid ne $onerefseq);
	}
    }
    elsif ($nodeclass eq "nseq0") {
	# node represents known taxonomy without refseqs, +1 because of unk row in X
	$yi = 1 + childIndex($node, $n2pref, $n2cref);
	print " #seqs ==0 in taxonomy, rnode '$node' trainseq '$trainseq'\n" if ($debug_mode);

	$seq_dont_use{$trainseq} = 1; 
	foreach $seqid (@{$n2ridref->{$node}}) {
	    $seq_dont_use{$seqid} = 1;
	}
    }
    elsif ($nodeclass eq "unk") {
        # node represents missing taxonomy branch, 0 is unk row index in X
	$yi = 0;
	$remove_i = childIndex($node, $n2pref, $n2cref);
	print " #seqs ==0 UNK, rnode '$node' trainseq '$trainseq'\n" if ($debug_mode);

	$seq_dont_use{$trainseq} = 1;
	foreach $seqid (@{$n2ridref->{$node}}) {
	    $seq_dont_use{$seqid} = 1;
	}
    }
    else {
	die "ERROR: unknown nodeclass '$nodeclass':\n$_\n";
    }

    # for each node along path from present node to root, create X matrix based on sequence similarity
    # get sequence similarity between training seq and each node's reference seqs

    my $pnode = $n2pref->{$node};
    my @xdat=($weight, $priprob);
    my ($x, $xwref) = createX($pnode, $n2cref, $n2ridref, $n2rweightref, \%seq_dont_use, $trainseq, $remove_i, $unknoderef, $seqsimref, $n2prior, $debug_mode);
#    push(@xdat, [$n2lref->{$pnode}, $yi, @$x]);
    push(@xdat, [0, $yi, @$x]);
    my @xweight = ([@$xwref]);

    # root node == 0
#    while ($pnode != 0) {
#	$node = $pnode;
#	$pnode = $n2pref->{$node};
#	$yi = 1 + childIndex($node, $n2pref, $n2cref);
#	$x = createX($pnode, $n2cref, $n2ridref, $n2rweightref, \%seq_dont_use, $trainseq, -1, $unknoderef, $seqsimref, $debug_mode);
#	push(@xdat, [$n2lref->{$pnode}, $yi, @$x]);
#    }

    if ($debug_mode > 1) {
	print "XDAT: weight: $xdat[0] priprob: $xdat[1] X-matrix rows:\n";
	for (my $j=2; $j<scalar(@xdat); $j++) {
	    print "   @{$xdat[$j]}\n";
	}
    }
    
    push(@allxdat, [@xdat]);
    push(@allxweight, [@xweight]);
}
close(FD);

open(FD,">$outfile") or die "ERROR: cannot write to '$outfile'. $!\n";

# offset: array start index in output file (0 or 1)
my $offset = 1;
my $betadim = 4;
my $x;
for (my $i=0; $i<scalar(@allxdat); $i++) {
    my @xdat = @{$allxdat[$i]};
    my @xweight = @{$allxweight[$i]};
    my $weight = shift(@xdat);
    my $priprob = shift(@xdat);
    $priprob = sprintf("%.10f",$priprob);
    my $nlevels = scalar(@xdat);
    my @xout = ();
    for (my $j=0; $j<$nlevels; $j++) {
	# -2: first two elements for prob and yi
	my $nrows = (scalar(@{$xdat[$j]}) - 2)/$betadim;
	my $n = scalar(@{$xdat[$j]}) -1;
	$x = join(',', ($xdat[$j]->[0] + $offset,$xdat[$j]->[1] + $offset,$nrows,$betadim,@{$xweight[$j]},@{$xdat[$j]}[2..$n]));
	push(@xout, $x);
    }
    $x = join(';', @xout);
    print FD "$weight\t$priprob\t$x\n";
}
close(FD);
