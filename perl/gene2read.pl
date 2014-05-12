#!/usr/bin/perl

use strict;
use Getopt::Long;
use Data::Dumper;
use File::Find;
use File::Basename;

my %STEPS=(
	pipeline => 'Demultiplex',
	trim => 'Trimmomatic',
	cutadapt => 'Cutadapt',
	flash => 'FLASH');

# suppress output if # of readcounts is below this threshold
my $MIN_READCOUNT=12;
# position in FA to start to generate unique sequence to lookup
my $START_POS=8;
# length of unique sequence lookup
my $SLENGTH=13;

######################
## ARGUMENT HANDLING #
######################


my $USAGE=<<EOL;
perl $0 -[d]ir -[f]asta > outfile

	dir - top level directory containing well based Methpup analysis.
	fasta - path to FASTA formatted design file for experiment.
	
Outputs a tab delimited file of results to STDOUT
EOL

my ($fasta,$dir,$help);
GetOptions (
	'dir|d=s' => \$dir,
	'fasta|f=s' => \$fasta,
	'help|h'=>\$help);

my $ERROR_FLAG;

if($help){
	print STDERR $USAGE."\n";
	exit(1);
}

if(!$dir){
	print STDERR "[ERROR] -[d]ir parameter required\n";
	$ERROR_FLAG++;
}

if(!$fasta){
	print STDERR "[ERROR] -[f]asta parameter required\n";
	$ERROR_FLAG++;
}

unless(-d $dir){
	print STDERR "[ERROR] cannot locate dir $dir\n";
	$ERROR_FLAG++;
}

unless(-e $fasta){
	print STDERR "[ERROR] cannot locate fasta $fasta\n";
	$ERROR_FLAG++;
}

if($ERROR_FLAG){
	print STDERR $USAGE."\n";
	exit(1);
}



#my $FASTA_FILE='/home/oliver/GIT_REPOS/Methpup/NGS12_reference.fa';
#my $dir  =     '/stats/oliver/NGS12_BISSEQ_LATEST/pipeline/5_A6';



my $greps = getIndex($fasta,$START_POS,$SLENGTH);
my $id = basename($dir);
                                         
my %RESULTS;
find(sub{
		next unless /(\.fwd\.fq\.gz$)|(\.fwd\.fq\.gz\.extendedFrags)/;
		next if /sam.gz$/;
		next if -d  $File::Find::name;    
		my $tdir = $STEPS{basename(dirname($File::Find::name))} || 'Demultiplex';       
		foreach my $s(keys %$greps){
			my $count = `zcat $File::Find::name | grep -c $s`;
			chomp($count);
			next if $count <$MIN_READCOUNT;
			my $g = $greps->{$s};
			$RESULTS{$tdir}{$g}=$count;
		}
},$dir);

## what is in the sam file - follows this recipe
my $samfile = "$dir/bowtie2/$id.fwd.fq.gz.extendedFrags.sam.gz";
foreach my $s(keys %$greps){
	my $g = $greps->{$s};
	my $count = `zcat $samfile | grep -v '^\@'  | grep -c $g `;
	chomp($count);
		next if $count <$MIN_READCOUNT;
	$RESULTS{'Bowtie'}{$g}=$count;
}

foreach my $step(qw/Demultiplex Trimmomatic Cutadapt FLASH Bowtie/){
	foreach my $g(keys %{$RESULTS{$step}}){
		print join("\t",$g,$step,$RESULTS{$step}{$g},$id)."\n";
	}
}
	

sub getIndex{
	my($fa,$start,$minlen)=@_;
	open(FASTA,"$fa") || die "Cannot open $fa\n";
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
	close(FASTA);
	my %return;                            
	foreach my $g(keys %results){
		my $seq = substr($results{$g},$start,$minlen);
		    
		if($return{$seq}){
			die "COLLISION $g\n";
		}else{
			my @match;
			foreach my $sg(keys %results){
				next if $g eq $sg;
				push @match,$sg if $results{$sg}=~m/$seq/;
			}
			#my @match = grep{return /$seq/}values %results;
			unless($seq=~m/^[AGCT]+$/){
				print STDERR "$g $seq overlaps methylation site truncating\n";
				$seq=~s/(^[AGCT]+)[^AGCT]+/\1/;
			}
			if(@match!=0){
				print STDERR "$g $seq matches ".join(",",@match)." product(s)\n";
			}
			$return{$seq}=$g;
		}
	}
	
	## we need to check here that the sequence doesn't match 
	## any other products
	my %MATCH;
	return \%return;
}
