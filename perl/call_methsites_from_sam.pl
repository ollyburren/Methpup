#!/usr/bin/perl

use strict;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/basename/;

######################
## ARGUMENT HANDLING #
######################

my $USAGE=<<EOL;
perl $0 -[s]amfile -[m]eth_file -[o]utdir {-v[erbose] -r[elaxed] -[a]lign_fasta}

	samfile - path to input SAM formatted file (output from bowtie2)
	align_fasta - path to FASTA formatted sequences used to run bowtie step - NOT CURRENTLY IMPLEMENTED.
	meth_file - path to tab delim file for locations of meth site for each insert.
	relaxed - boolean, if set to 0 then reads w/o C/T at methsite will be output
	outdir - path to write out results
	
Original algorithm Xin Yang, converted to Perl by Olly Burren
EOL

my($sam,$out,$meth,$help,$verbose,$fasta,$relaxed);

GetOptions (
	'samfile|s=s' => \$sam,
	'outdir|o=s' => \$out,
	'methfile|m=s' => \$meth,
	'align_fasta|a=s'=>\$fasta,
	'help|h' => \$help,
	'relaxed|r'=> \$relaxed,
	'verbose|v' => \$verbose
	);

my $ERROR_FLAG;

if($help){
	print STDERR $USAGE."\n";
	exit(1);
}

if(!$sam){
	print STDERR "[ERROR] -[s]amfile parameter required\n";
	$ERROR_FLAG++;
}

if(!$out){
	print STDERR "[ERROR] -[o]utfile parameter required\n";
	$ERROR_FLAG++;
}

if(!$meth){
	print STDERR "[ERROR] -[m]eth_file parameter required\n";
	$ERROR_FLAG++;
}

unless(-e $sam){
	print STDERR "[ERROR] cannot locate samfile $sam\n";
	$ERROR_FLAG++;
}

unless(-e $meth){
	print STDERR "[ERROR] cannot locate meth_file $meth\n";
	$ERROR_FLAG++;
}


#unless(-e $fasta){
#	print STDERR "[ERROR] cannot locate align_fasta $fasta\n";
#	$ERROR_FLAG++;
#}

unless(-d $out){
	print STDERR "[ERROR] cannot locate $out for writing\n";
}

if($ERROR_FLAG){
	print STDERR $USAGE."\n";
	exit(1);
}

#########
## MAIN #
#########
(my $stub = basename($sam)) =~s/\.sam$//;

my $meths = &read_meth_sites($meth);
#my $sequences = &read_align_fasta($fasta,$meths);


open(SAM,$sam) || die "Cannot open SAM file for reading $sam\n";
my %RHASH;
while(<SAM>){
	chomp;
	my $res = &process_indel($_);
	my ($g,$lm,$s,$l)=@$res;
	$RHASH{$g}{$s}{$lm}++;
}

my %CALLS;
my %EFF;
foreach my $g(keys %RHASH){
	foreach my $s(keys %{$RHASH{$g}}){
		foreach my $lm(keys %{$RHASH{$g}{$s}}){
			my $count = $RHASH{$g}{$s}{$lm};
			my $cm = &call_meth($g,$lm,$s,$count,$meths,$relaxed);
			## NOTE: Add this in once I have a source for information
			## to calculate BS efficiency
			## my $eff = &bs_eff($g,$lm,$s,$count,$meths,$relaxed);
			$CALLS{$g}{$cm->[1]}+=$cm->[2];
		}
	}
}

## Here we output to format expected by Methplot R library
## https://github.com/XinYang6699/Methplot
## Actual format returned is more efficient perhaps recommend 
## Methplot supports this

open(MP,">$out/$stub.methplot") || die "Cannot open $out/$stub.methplot for writing\n";
foreach my $g(keys %CALLS){
	foreach my $s(keys %{$CALLS{$g}}){
		for(my $i;$i<$CALLS{$g}{$s};$i++){
			last if $s eq 'NOCALL';
			## note that we do not filter results 
			## that are not C or T
			print MP join("\t",$g,$s)."\n";
		}
	}
}
close(MP);


## output bp calls at each position format favoured by Dan

my $c_to_t = &c_to_t(\%CALLS,$meths);

foreach my $k(keys%{$c_to_t}){
	my @header='gene';
	my @value=$k;
	foreach my $r(@{$c_to_t->{$k}}){
		push @header,$r->{name};
		push @value,$r->{count};
	}
	open(CTOT,">$out/${stub}.$k.tab") || die "Cannot open  >$out/${stub}.$k.tab for writing\n";
	print CTOT join("\t",@header)."\n".join("\t",@value)."\n";
	close(CTOT)
}

## some messing around to generate PWM's for
## use is seqLogo's

##my $foobar = &generate_pwm(\%CALLS);

##foreach my $k(keys %$foobar){
##	print "##$k\n";
##	print $foobar->{$k}."\n";
#}

## print $foobar;


###############
# SUBROUTINES #
###############

sub read_align_fasta{
	my ($file,$ms) = @_;
	open(FASTA,"$file") || die "Cannot open $file\n";
	my %results;
	my ($gene,$seq);
	while(<FASTA>){
		next if/^#/ | /^$/;
		chomp;
		##allow spaces but trim off end
		if(/^>(.*)\s*$/){
			$gene=$1;
		}else{
			(my $tmp = $_) =~s/\s//g;
			$results{$gene}.=$tmp;
		}
	}
	my %return;
	foreach my $g(keys %{$ms}){
		## count Y's or C's not paired with a G
		next unless $results{$g};
		## get 
		my $ct =()= $results{$g} =~ /[YC][^G]/g;
		## remove the number that are methylated
		#$ct= $ct - scalar(@{$ms->{$g}});
		$return{$g}=$ct;
	}
	return \%return;
}

		
		

#output pwm format for use by seqLogo

sub generate_pwm{
	my ($rhash)=@_;
	my %r=%$rhash;
	my %return;
	foreach my $g(keys %r){
		my %fs;
		foreach my $s(keys %{$r{$g}}){
			my $total = $r{$g}{$s};
			next if $s eq 'NOCALL';
			my @seq = split('',$s);
			for(my $i=0;$i<@seq;$i++){
				#$return{$g}{$i}{$seq[$i]}+=$total;
				$fs{$i}{$seq[$i]}+=$total;
				$fs{$i}{TOTAL}+=$total;
			}
		}
		my @pwm;
		my $string="\n";
		foreach my $bp(qw/A C G T/){
			my @posum;
			for(my $j=0;$j<scalar(keys %fs);$j++){
				push @posum,$fs{$j}{$bp}/$fs{$j}{TOTAL};
			}
			$string.=join("\t",@posum)."\n";
		}
		
		$return{$g}=$string;
	}
	return \%return;
}

## produces output that Dan is familiar with and has asked for

sub c_to_t{
	my ($rhash,$meths)=@_;
	my %r=%$rhash;
	my %ret;
	my %return;
	foreach my $g(keys %r){
		foreach my $s(keys %{$r{$g}}){
			my $total = $r{$g}{$s};
			if ($s eq 'NOCALL'){
				$ret{$g}{'NOCALL'}+=$total;
				next;
			}elsif($s =~ /[^CT]/){
				$ret{$g}{'BADCALL'}+=$total;
				next;
			}
			## counts of c to t
			my $c =()= $s =~ /C/g;
			my $t =()= $s =~ /T/g;
			$ret{$g}{"C${c}T${t}"}+=$total;
		}
	}
	foreach my $g (keys(%ret)){
		my %gt = %{$ret{$g}};
		my @keys = grep{/^C/}keys %gt;
		my @sortedk = sort{
			$a =~ /C([0-9]+)T([0-9]+)/;
			my ($ac,$at) = $a=~/C([0-9]+)T([0-9]+)/;
			my ($bc,$bt) = $b=~/C([0-9]+)T([0-9]+)/;
			$ac <=> $bc || $bt <=> $at;
		}@keys;
		## get the max number of C's and max number of T's
		$sortedk[-1]=~/C([0-9]+)/g;
		my $maxC=$1;
		$sortedk[0]=~/C[0-9]+T([0-9]+)/g;
		my $maxT=$1;
		my $counter = $maxC < $maxT ? $maxT : $maxC;
		#print $sortedk[-1]."\t".$maxC."\t".$maxT."\n";
		for(my $j=0;$j<=$counter;$j++){
			my $k = "C${j}T".($counter-$j);
			push @{$return{$g}},{name=>$k,count=>$gt{$k}||0};
		}
		#foreach my $k(@sortedk){
		#	(my $seq = $gt{$k})=~s/N0//;
		#	push @{$return{$g}},{name=>$k,count=>$gt{$k}};
		#}
		push @{$return{$g}},{name=>'NOCALL',count=>$gt{NOCALL}||0};
		push @{$return{$g}},{name=>'BADCALL',count=>$gt{BADCALL}||0};
	}
	return \%return;
}
## read in tab delim file with coordinates for methylations sites.

sub read_meth_sites{
	my $meth=shift;
	open(METH,$meth) || die "Cannot open meth_file $meth\n";
	my %METHS;
	while(<METH>){
		chomp;
		my @vals = grep{/\d+/}split("\t",$_);
		$METHS{$vals[0]}=[@vals[1..$#vals]];
	}
	close(METH);
	return \%METHS;
}

## call methylation sites for a particular sequence

sub call_meth{
	my($g,$lm,$s,$count,$msites,$relaxed)=@_;
	my @seq=split("",$s);
	my @c = @{$msites->{$g}};
	my $ret;
	if(($c[0]- $lm) > 0 && ($c[-1]-$lm)<=length($s)){
		my @res;
		foreach my $si(@c){
			push @res,$seq[$si-$lm]
		}
		$ret = join("",@res);
		$ret = 'NOCALL' if $ret!~/^[CT]+$/ && ! $relaxed; 
	}else{
		$ret = 'NOCALL';
	}
	return [$g,$ret,$count];
}

## preprocess SAM files so so that insertions are removed and deletions
## are added back and softclips are removed.

sub process_indel{
	my $samline=shift;
	my @out;
	my @vals = split("\t",$_);
	push @out,@vals[2,3];
	my (@numb,@patb);
	push @numb,$1 while($vals[5]=~/(\d+)/g);
	push @patb,$1 while($vals[5]=~/(\D+)/g);
	##perfect match
	if(@patb==0){
		push @out,$vals[9];
	}else{
		my @seq=split("",$vals[9]);
		@seq=@seq[0..($#seq-$numb[-1])]if $patb[-1] eq 'S'; ##REPLACE
		## compatibility with Xin's code possible bug
		## if the REPLACE line is replaced with the following then
		## get the same output as Xin. It seems her code removes
		## an extra base replaces with the last base constittutively 
		## for 3' soft clipped alignments. It shouldn't make any difference
		## to the downstream analysis (it's impossible to have a methylation site 
		## thats callable at the end of an alignment).
		
		## my $lastseq=$seq[-1]; 
		## @seq=@seq[0..($#seq-$numb[-1]-1)]if $patb[-1] eq 'S';
		## push @seq,$lastseq if $patb[-1] eq 'S';
		
		#get postions for deletions
		my @nd;
		for(my $i=0;$i<@patb;$i++){
			push @nd,$i if $patb[$i] eq 'D'
		}
		my @posb=@numb;
		if(@nd>0){
			$posb[$_]=0 foreach @nd; 
		}
		@posb=cumsum(@posb);
		#check if there is an insertion
		my @ni;
		for(my $i=0;$i<@patb;$i++){
			push @ni,$i if $patb[$i] eq 'I'
		}
		if(@ni>0){
			$seq[$posb[$_]]='I' foreach @ni; 
		}
		#check if there is a deletion
		if(@nd>0){
			@nd =reverse @nd;
			foreach my $i(@nd){
				splice @seq,$posb[$i],0,'D'x $numb[$i]
			}
		}
		@seq=splice(@seq,$posb[0]) if $patb[0] eq 'S';  
		push @out, join('',grep{$_ ne 'I'}@seq);
	}
	return \@out;
}

## compute the cumulative sum over an array.

sub cumsum{
	my @array=@_;
	my @tmp=$array[0];
	for(my $i=1;$i<@array;$i++){
		my $prev=$tmp[$i-1];
		push @tmp,$prev+$array[$i];
	}
	return @tmp;
}
			
			
	
	


