package Methpup::Runnable::Reports;

use strict;
use base('Methpup::RunnableI');
use Data::Dumper;
use File::Basename;
use IO::Zlib;

#tests <= (
#	test.ctrl => [ ds1,ds2] ..
#)

my $SUBDIR_NAME = 'report';

use constant _inputs_expected => {
	cm_files=>'hashref',
	callmeth_sub=>'string'
};
	
sub run{
	my $self=shift;
	my @dirs;
	foreach my $k(keys %{$self->inputs->{cm_files}}){
		my $d = dirname($self->inputs->{cm_files}->{$k}->[0]);
		my $subdir = $self->inputs->{callmeth_sub};
		$d =~s/$subdir$//;
		my $dname = $d.'/'.$SUBDIR_NAME.'/';
		$self->verbose("Reporting on $d\n");
		File::Path::remove_tree($dname,{keep_root => 1,result=> \my $del_dirs});
		Macd::Utils::mkpath_wrapper($dname);
		my @params = $d;
		my $cmd ='perl '.$self->binary." $d > $dname/counts.tab";
		$self->verbose($cmd);
		$self->step->add_job(Macd::Step::Job->new(command=>$cmd));
	}
	if($self->step->execute()){
		$self->step->wait_on_complete();
	}
	return 1;
}



sub subdir{
	return $SUBDIR_NAME;
}

sub trim_filelist{
	my $self=shift;
	my $hidewarning=shift;
	#if($self->outputs->{trim_filelist}){
	#	return $self->outputs->{trim_filelist};
	#}
	my %cm_files = %{$self->inputs->{cm_files}};
	my $subdir = $self->inputs->{callmeth_sub};
	my %ret;
	foreach my $k(keys(%cm_files)){
		my $d = dirname($self->inputs->{cm_files}->{$k}->[0]);
		#$fpattern.='.extendedFrags.fastq.gz';
		(my $dname=$d) =~ s/$subdir$/$SUBDIR_NAME/;
		next unless -d $dname;
		opendir(DIR,$dname) || die "Cannot open $dname\n";
		my @files = grep{/counts\.tab/}readdir(DIR);
		foreach my $f(@files){
			if(-e $dname."/".$f){
				push @{$ret{$k}},$dname.'/'.$f;
			}elsif(-d $self->log_dir){
				$self->msg("[WARNING] missing expected file ".$f." - skipping further analysis") unless $hidewarning;
			}
		}
	}
	$self->outputs->{filelist}=\%ret;
	return \%ret;
}


sub _skip_message{
	my $self=shift;
	return "[".ref($self)."]: The output path ".$self->inputs->{out_dir}."/".$SUBDIR_NAME."/ already exists and contains files. Would you like to skip this step ";
}

sub run_conditions{
	#return 0;
	my $self=shift;
	my @k = keys(%{$self->trim_filelist('hidewarning')});
	my $flag=0;
	$flag =1 if @k > 0;
	return $flag;
}
1;



