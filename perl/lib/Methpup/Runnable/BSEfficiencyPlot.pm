package Methpup::Runnable::BSEfficiencyPlot;

use strict;
use base('Methpup::RunnableI');
use Data::Dumper;
use File::Basename;
use File::Find;

#tests <= (
#	test.ctrl => [ ds1,ds2] ..
#)


use constant _inputs_expected => {
	in_dir=>'file',
	out_dir=>'file',
	file_pattern=>'string'
};

use constant _defaults => {
};

#efficiency.methplot

sub run{
	my $self=shift;
	my $idir = $self->inputs->{in_dir};
	my $rdir = $self->inputs->{out_dir};
	my $fpattern = $self->inputs->{file_pattern};
	if(-d $rdir){
		$self->debug("Would remove $rdir");
		#File::Path::remove_tree($rdir,{keep_root => 1,result=> \my $del_dirs});
	}
	Macd::Utils::mkpath_wrapper($rdir);
	my $fcmd = "find $idir -name '$fpattern' | xargs cat > ${rdir}bs.efficiency.tab";
	$self->verbose($fcmd);
	`$fcmd`;
	## in.file out.dir
	my @params = $self->inputs->{script};
	push @params,"in.file=\\'${rdir}bs.efficiency.tab\\'";
	push @params ,"out.dir=\\'$rdir\\'";
	my $cmd = join(" ",$self->binary,@params);
	$self->verbose($cmd);
	$self->step->add_job(Macd::Step::Job->new(command=>$cmd));
	if($self->step->execute()){
		$self->step->wait_on_complete();
	}
	return 1;
}


sub _skip_message{
	my $self=shift;
	return "[".ref($self)."]: This step has already been run. Would you like to skip it ";
}

sub run_conditions{
	my $self=shift;
	return 1 if -e  $self->inputs->{out_dir}."/bs.efficiency.pdf";
	return 0;
}
1;



