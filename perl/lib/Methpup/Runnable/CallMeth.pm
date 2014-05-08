package Methpup::Runnable::CallMeth;

use strict;
use base('Methpup::RunnableI');
use Data::Dumper;
use File::Basename;

#tests <= (
#	test.ctrl => [ ds1,ds2] ..
#)

my $SUBDIR_NAME = 'methylation_sites';

use constant _inputs_expected => {
	meth_site_file=>'file',
	genomic_fa_file=>'file',
	samfiles=>'hashref',
	bowtie_sub=>'string',
};

use constant _defaults => {
};

	
sub run{
	my $self=shift;
	my %samfiles = %{$self->inputs->{samfiles}};
	
	foreach my $s(keys(%samfiles)){
		$self->verbose("Processing $s");
		my $file = $samfiles{$s};
		my $dname = dirname($file);
		my $subdir = $self->inputs->{bowtie_sub};
		$dname =~ s/$subdir$/$SUBDIR_NAME/;
		(my $fpattern = basename($file)) =~s/\.fastq.gz$//;
		File::Path::remove_tree($dname,{keep_root => 1,result=> \my $del_dirs});
		Macd::Utils::mkpath_wrapper($dname);
		#forward;
		my @params = "-s $file";
		push @params,"-o $dname";
		push @params,"-m ".$self->inputs->{meth_site_file};
		push @params,"-g ".$self->inputs->{genomic_fa_file};
		push @params,"-i $s";
		my $cmd = join(" ",'perl',$self->binary,@params);
		$self->verbose($cmd);
		$self->step->add_job(Macd::Step::Job->new(command=>$cmd));
	}
	if($self->step->execute()){
		$self->step->wait_on_complete();
	}
	return 1;
}

sub filelist{
	my $self=shift;
	my $hidewarning=shift;
	#if($self->outputs->{filelist}){
	#	return $self->outputs->{filelist};
	#}
	my %samfiles = %{$self->inputs->{samfiles}};
	my $subdir = $self->inputs->{bowtie_sub};
	my %ret;
	foreach my $s(keys(%samfiles)){
		#(my $fpattern = basename($samfiles{$s})) =~s/\.sam.gz$/*/;
		(my $fpattern = basename($samfiles{$s})) =~s/^([^\.]+)\..*$/\1*/;
		#$fpattern.='.extendedFrags.fastq.gz';
		(my $dname=dirname($samfiles{$s}))=~s/$subdir$/$SUBDIR_NAME/;
		next unless -d $dname;
		opendir(DIR,$dname) || die "Cannot open $dname\n";
		my @files = grep{/$fpattern/}readdir(DIR);
		foreach my $f(@files){
			if(-e $dname."/".$f){
				push @{$ret{$s}},$dname.'/'.$f;
			}elsif(-d $self->log_dir){
				$self->msg("[WARNING] missing expected file ".$f." - skipping further analysis") unless $hidewarning;
			}
		}
	}
	$self->outputs->{filelist}=\%ret;
	return \%ret;
}

sub subdir{
	return $SUBDIR_NAME;
}


sub _skip_message{
	my $self=shift;
	return "[".ref($self)."]: The output path ".$self->inputs->{out_dir}."/".$SUBDIR_NAME."/ already exists and contains files. Would you like to skip this step ";
}

sub run_conditions{
	#return 0;
	my $self=shift;
	my @k = keys(%{$self->filelist('hidewarning')});
	my $flag=0;
	$flag =1 if @k > 0;
	return $flag;
}
1;



