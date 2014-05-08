package Methpup::Runnable::Reports;

use strict;
use base('Methpup::RunnableI');
use Data::Dumper;
use File::Basename;

#tests <= (
#	test.ctrl => [ ds1,ds2] ..
#)

my $JAVA_JAR = 'java -jar';
my $SUBDIR_NAME = 'trim';

use constant _inputs_expected => {
	working_dir=>'file',
	steps=>'array'
};

	
sub run{
	my $self=shift;
	my %pe_list = %{$self->inputs->{pe_file_list}};
	foreach my $pef(keys(%pe_list)){
		$self->verbose("Trimmomatic processing $pef");
		my @files = @{$pe_list{$pef}};
		my $dname = dirname($files[1])."/$SUBDIR_NAME/";
		if(-d $dname){
			#if($self->skip_ask){
			#	$self->msg("Skipping $pef");
			#}
		}else{
			File::Path::remove_tree($dname,{keep_root => 1,result=> \my $del_dirs});
			Macd::Utils::mkpath_wrapper($dname);
		}
		my @params = 'PE';
		push @params, "-phred".$self->inputs->{phred};
		push @params, "-trimlog ${dname}${pef}_trimmomatic.log";
		push @params, join(" ",@files);
		push @params, join(" ",map{s/\.fq\.gz$//;"${dname}$_.fq.gz ${dname}${_}.unpaired.fq.gz"}map{basename($_)}@files);
		push @params, 'HEADCROP:'.$self->inputs->{linker_length}.' '.$self->inputs->{ops};
		my $cmd = join(" ",$JAVA_JAR,$self->binary,@params);
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
	my %pe_list = %{$self->inputs->{pe_file_list}};
	my %ret;
	foreach my $pef(keys(%pe_list)){
		my @files = map{dirname($_)."/trim/".basename($_)}@{$pe_list{$pef}};
		if(scalar(grep{-e $_}@files)==2){
			$ret{$pef}=\@files;
		}elsif(-d $self->log_dir){
			$self->msg("[WARNING] missing expected file ".dirname($files[1])." - skipping further analysis") unless $hidewarning;
		}## else we are silent
	}
	$self->outputs->{trim_filelist}=\%ret;
	return \%ret;
}


sub _skip_message{
	my $self=shift;
	return "[".ref($self)."]: The output path ".$self->inputs->{out_dir}." already exists and contains files. Would you like to skip this step ?";
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



