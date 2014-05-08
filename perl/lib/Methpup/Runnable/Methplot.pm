package Methpup::Runnable::Methplot;

use strict;
use base('Methpup::RunnableI');
use Data::Dumper;
use File::Basename;
use File::Find;

#tests <= (
#	test.ctrl => [ ds1,ds2] ..
#)

my $SUBDIR_NAME = 'cutadapt';

use constant _inputs_expected => {
	script=>'file',
	in_dir=>'file',
	out_dir=>'file'
};

use constant _defaults => {
};
	
sub run{
	my $self=shift;
	my $odir = $self->inputs->{out_dir};
	if(-d $odir){
		#File::Path::remove_tree($self->inputs->{out_dir},{keep_root => 1,result=> \my $del_dirs});
	}else{
		Macd::Utils::mkpath_wrapper($odir);
	}
	## prepare plots for all files with extenstion
	find(sub{
			next unless m/\.methplot$/;
			next if m/eff/; ## don't plot bs files at this stage.
			my @params = $self->inputs->{script};
			push @params,"mp.file=\\'$File::Find::name\\'";
			push @params ,"out.dir=\\'$odir\\'";
			my $cmd = join(" ",$self->binary,@params);
			$self->verbose($cmd);
			$self->step->add_job(Macd::Step::Job->new(command=>$cmd));
	},$self->inputs->{in_dir});
	if($self->step->execute()){
		$self->step->wait_on_complete();
	}
	return 1;
}

sub _skip_message{
	my $self=shift;
	return "[".ref($self)."]: The output path ".$self->inputs->{out_dir}."/".$SUBDIR_NAME."/ already exists and contains files. Would you like to skip this step ";
}

sub run_conditions{
	#return 0;
	my $self=shift;
	my $flag=0;
	find(sub{
			if(/^methplot/){
				$flag++;
			}
	},$self->inputs->{out_dir});
	return $flag;
}
1;



