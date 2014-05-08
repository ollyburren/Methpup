package Methpup::Runnable::Cutadapt;

use strict;
use base('Methpup::RunnableI');
use Data::Dumper;
use File::Basename;

#tests <= (
#	test.ctrl => [ ds1,ds2] ..
#)

my $SUBDIR_NAME = 'cutadapt';

use constant _inputs_expected => {
	forseq=>'string',
	revseq=>'string',
	linformin=>'number',
	linrevmin=>'number',
	pe_file_list=>'hashref',
	trimmo_sub=>'string',
};

use constant _defaults => {
};
	
sub run{
	my $self=shift;
	my %pe_list = %{$self->inputs->{pe_file_list}};
	foreach my $pef(keys(%pe_list)){
		$self->verbose("Cutadapt processing $pef");
		my @files = @{$pe_list{$pef}};
		my $dname = dirname($files[1]);
		my $trimsubdir = $self->inputs->{trimmo_sub};
		$dname =~ s/$trimsubdir/$SUBDIR_NAME/;
		File::Path::remove_tree($dname,{keep_root => 1,result=> \my $del_dirs});
		Macd::Utils::mkpath_wrapper($dname);
		#forward;
		my @params = '-a '.$self->inputs->{forseq};
		push @params, '-O '.$self->inputs->{linformin};
		push @params, $files[0];
		push @params, '| gzip -';
		push @params, " > $dname/".basename($files[0]);
		my $cmd = join(" ",$self->binary,@params);
		$self->verbose($cmd);
		$self->step->add_job(Macd::Step::Job->new(command=>$cmd));
		#reverse
		@params = '-a '.$self->inputs->{revseq};
		push @params,  '-O '.$self->inputs->{linrevmin};
		push @params, $files[1];
		push @params, '| gzip -';
		push @params, " > $dname/".basename($files[1]);
		$cmd = join(" ",$self->binary,@params);
		$self->verbose($cmd);
		$self->step->add_job(Macd::Step::Job->new(command=>$cmd));
	}
	if($self->step->execute()){
		$self->step->wait_on_complete();
	}
	return 1;
}

sub cutadapt_filelist{
	my $self=shift;
	my $hidewarning=shift;
	#if($self->outputs->{cutadapt_filelist}){
	#	return $self->outputs->{cutadapt_filelist};
	#}
	my %pe_list = %{$self->inputs->{pe_file_list}};
	my $trimsubdir = $self->inputs->{trimmo_sub};
	my %ret;
	foreach my $pef(keys(%pe_list)){
		my @files = map{
											my $dname=dirname($_);
											$dname=~s/$trimsubdir$/$SUBDIR_NAME/;
											$dname.'/'.basename($_);
								}@{$pe_list{$pef}};
		if(scalar(grep{-e $_}@files)==2){
			$ret{$pef}=\@files;
		}elsif(-d $self->log_dir){
			$self->msg("[WARNING] missing expected file ".dirname($files[1])." - skipping further analysis") unless $hidewarning;
		}## else we are silent
	}
	$self->outputs->{cutadapt_filelist}=\%ret;
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
	my @k = keys(%{$self->cutadapt_filelist('hidewarning')});
	my $flag=0;
	$flag =1 if @k > 0;
	return $flag;
}
1;



