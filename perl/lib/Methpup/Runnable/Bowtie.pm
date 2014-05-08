package Methpup::Runnable::Bowtie;

use strict;
use base('Methpup::RunnableI');
use Data::Dumper;
use File::Basename;

#tests <= (
#	test.ctrl => [ ds1,ds2] ..
#)

my $SUBDIR_NAME = 'bowtie2';

use constant _inputs_expected => {
	ref_seq_file=>'file',
	pe_file_list=>'hashref',
	flash_sub=>'string',
	phred=>'number',
};

use constant _defaults => {
	#ops=>'--np 0 --no-unal --local --no-head'
	ops=>'--np 0 --no-unal --local'
};
	
sub run{
	my $self=shift;
	## first we build bowtie indexes
	if(!-d $self->inputs->{index_dir}){
		Macd::Utils::mkpath_wrapper($self->inputs->{index_dir});
	}
	my $cmd = dirname($self->binary)."/bowtie2-build ".$self->inputs->{ref_seq_file}.' '.$self->inputs->{index_dir};
	$self->verbose($cmd);
	$self->step->add_job(Macd::Step::Job->new(command=>$cmd));
	if($self->step->execute()){
		$self->step->wait_on_complete();
	}
	sleep(3);
	$cmd='cp '.$self->inputs->{ref_seq_file}.' '.$self->inputs->{index_dir};
	`$cmd`;
	#return 1;
	my %pe_list = %{$self->inputs->{pe_file_list}};
	foreach my $pef(keys(%pe_list)){
		$self->verbose("Bowtie2 processing $pef");
		my $file = $pe_list{$pef};
		my $dname = dirname($file);
		my $subdir = $self->inputs->{flash_sub};
		$dname =~ s/$subdir$/$SUBDIR_NAME/;
		(my $fpattern = basename($file)) =~s/\.fastq.gz$//;
		File::Path::remove_tree($dname,{keep_root => 1,result=> \my $del_dirs});
		Macd::Utils::mkpath_wrapper($dname);
		#forward;
		my @params = '--phred'.$self->inputs->{phred};
		push @params, $self->inputs->{ops};
		push @params, '-x '.$self->inputs->{index_dir};
		push @params, '-U '.$file;
		push @params, "-S  $dname/$fpattern.sam ; gzip $dname/$fpattern.sam";
		my $cmd = join(" ",$self->binary,@params);
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
	my %pe_list = %{$self->inputs->{pe_file_list}};
	my $subdir = $self->inputs->{flash_sub};
	my %ret;
	foreach my $pef(keys(%pe_list)){
		(my $fpattern = basename($pe_list{$pef})) =~s/\.fastq.gz$/.sam.gz/;
		#$fpattern.='.extendedFrags.fastq.gz';
		(my $dname=dirname($pe_list{$pef}))=~s/$subdir$/$SUBDIR_NAME/;
		$fpattern="$dname/$fpattern";
		if(-e $fpattern){
			$ret{$pef}=$fpattern;
		}elsif(-d $self->log_dir){
			$self->msg("[WARNING] missing expected file ".$fpattern." - skipping further analysis") unless $hidewarning;
		}## else we are silent
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



