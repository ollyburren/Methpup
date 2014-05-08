package Methpup::Runnable::FLASH;

use strict;
use base('Methpup::RunnableI');
use Data::Dumper;
use File::Basename;

#tests <= (
#	test.ctrl => [ ds1,ds2] ..
#)

my $SUBDIR_NAME = 'flash';

use constant _inputs_expected => {
	maxol=>'number',
	minol=>'number',
	pe_file_list=>'hashref',
	cutadapt_sub=>'string',
};

use constant _defaults => {
};
	
sub run{
	my $self=shift;
	my %pe_list = %{$self->inputs->{pe_file_list}};
	foreach my $pef(keys(%pe_list)){
		$self->verbose("FLASH processing $pef");
		my @files = @{$pe_list{$pef}};
		my $dname = dirname($files[1]);
		my $subdir = $self->inputs->{cutadapt_sub};
		$dname =~ s/$subdir$/$SUBDIR_NAME/;
		(my $fpattern = basename($files[0])) =~s/\.1\.fq.gz$//;
		File::Path::remove_tree($dname,{keep_root => 1,result=> \my $del_dirs});
		Macd::Utils::mkpath_wrapper($dname);
		#forward;
		my @params = '-M '.$self->inputs->{maxol};
		push @params, '-m '.$self->inputs->{minol};
		push @params, '-z';
		push @params, '-o '.$fpattern;
		push @params, '-d '.$dname;
		push @params, join(" ",@files);
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
	my $subdir = $self->inputs->{cutadapt_sub};
	my %ret;
	foreach my $pef(keys(%pe_list)){
		(my $fpattern = basename($pe_list{$pef}->[0])) =~s/\.1\.fq.gz$//;
		$fpattern.='.extendedFrags.fastq.gz';
		(my $dname=dirname($pe_list{$pef}->[0]))=~s/$subdir$/$SUBDIR_NAME/;
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



