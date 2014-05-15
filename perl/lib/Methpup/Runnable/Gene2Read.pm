package Methpup::Runnable::Gene2Read;

use strict;
use base('Methpup::RunnableI');
use Data::Dumper;
use File::Basename;
use File::Find;

#tests <= (
#	test.ctrl => [ ds1,ds2] ..
#)

my $SUBDIR_NAME = 'gene2read';

use constant _inputs_expected => {
	in_dir=>'file',
	fasta_file=>'file',
};

use constant _defaults => {
};

	
sub run{
	my $self=shift;
	my $idir = $self->inputs->{in_dir};
	opendir(DIR,$idir) || die "Cannot open $idir\n";
	my @dirs = map{"$idir/$_"}grep{-d "$idir/$_" && $_ !~/^[\.]+/}readdir(DIR);
	foreach my $d(@dirs){
		my $rdir = "$d/$SUBDIR_NAME/";
		if(-d $rdir){
			$self->verbose("Would remove $rdir");
			File::Path::remove_tree($rdir,{keep_root => 1,result=> \my $del_dirs});
		}
		Macd::Utils::mkpath_wrapper($rdir);
		my $id = basename($d);
		my $ofile = "${rdir}$id.gene2read.tab";
		my @params = "-d $d";
		push @params, "-f ".$self->inputs->{fasta_file};
		push @params, "> $ofile";
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
	my @return;
	find(sub{
			next unless /gene2read.tab$/;
			push @return, $File::Find::name;
	},$self->inputs->{in_dir});
	return @return;
}

sub subdir{
	return $SUBDIR_NAME;
}


sub _skip_message{
	my $self=shift;
	return "[".ref($self)."]: This step has already been run. Would you like to skip it ";
}

sub run_conditions{
	my $self=shift;
	my @k = $self->filelist;
	my $flag=0;
	$flag =1 if @k > 0;
	return $flag;
}
1;



