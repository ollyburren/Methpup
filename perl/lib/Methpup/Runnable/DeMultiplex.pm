package Methpup::Runnable::DeMultiplex;

use strict;
use base('Methpup::RunnableI');
use Data::Dumper;
use File::Basename;

#tests <= (
#	test.ctrl => [ ds1,ds2] ..
#)

use constant _inputs_expected => {
	forward_file => 'file',
	reverse_file => 'file',
	tag_file => 'file',
	out_dir => 'file',
	t=>'number',
	F=>'string',
};

use constant _defaults => {
	t=>30,
	F=>'ILM1.8',
};
	
sub run{
	my $self=shift;
	my $odir = $self->inputs->{out_dir};
	if(-d $odir){
		File::Path::remove_tree($self->inputs->{out_dir},{keep_root => 1,result=> \my $del_dirs});
	}else{
		Macd::Utils::mkpath_wrapper($odir);
	}
	##unlicensed version does not support gzipped files
	$self->msg("Unzipping files for use with ".$self->binary);
	my $f = $self->gunzip($self->{inputs}->{forward_file});
	my $r = $self->gunzip($self->{inputs}->{reverse_file});
	$self->step->add_job(Macd::Step::Job->new(command=>$f->[0])) unless -e $f->[1];
	$self->step->add_job(Macd::Step::Job->new(command=>$r->[0])) unless -e $r->[1];
	if($self->step->execute()){
		$self->step->wait_on_complete();
	}
	my %params;
	$params{-b}=$self->{inputs}->{tag_file};
	$params{-f}=join(" ",$f->[1],$r->[1]);
	$params{-t}=$self->{inputs}->{t};
	$params{-d}=$odir;
	$params{-F}=$self->{inputs}->{F};
	my $cmd = join(" ", $self->{binary},map{"$_ $params{$_}"}keys %params,'--GZIP');
	$self->msg("Demultiplexing with ".$self->binary);
	$self->verbose($cmd);
	my $job = Macd::Step::Job->new(command=>$cmd);
	$self->step->add_job($job);
	if($self->step->execute()){
		$self->step->wait_on_complete();
	}
	$self->msg("Finished demultiplexing with ".$self->binary);
	sleep(10);
	## next step is rename all the dirs in the target dir with wells
	my $lu = $self->_get_dlu();
	opendir (DIR, $odir) or die $!;
	#my @ofile=(basename($self->inputs->{forward_file}),basename($self->inputs->{reverse_file}));
	foreach my $d(grep{-d "$odir/$_"}readdir(DIR)){
		next if $d =~/^\./;
		if(my $dest = $lu->{$d}){
			delete $lu->{$d};
			`mv $odir/$d $odir/$dest`;
			$self->verbose("mv $odir/$d $odir/$dest");
			## next we move individual files to somewhere more sensible
			my $fwd_file = basename($self->{inputs}->{forward_file});
			my $rev_file = basename($self->{inputs}->{reverse_file});
			my $fwd_des_file = "$dest.fwd.fq.gz";
			my $rev_des_file = "$dest.rev.fq.gz";
			`mv $odir/$dest/$fwd_file $odir/$dest/$fwd_des_file`;
			$self->verbose("mv $odir/$dest/$fwd_file $odir/$dest/$fwd_des_file");
			`mv $odir/$dest/$rev_file $odir/$dest/$rev_des_file`;
			$self->verbose("mv $odir/$dest/$rev_file $odir/$dest/$rev_des_file");
		}else{
			$self->msg("Barcode $d not found in ".$self->inputs->{tag_file});
		}
	}
	open(MISS,">$odir/missing.tab") || die "Cannot open $odir/missing.tab\n";
	foreach my $k(keys %$lu){
		print MISS "$k ".$lu->{$k}."\n";
	}
	close(MISS);
	## finally clear up by deleting unzipped source files
	$self->verbose("Removing ".join(", ",$f->[1],$r->[1]));
	unlink($f->[1]);
	unlink($r->[1]);
	#$self->outputs->{demux_filelist}=\@demuxed;
	return 1;
}

## this routine returns a hashref of all the things demultiplexed.

sub demux_filelist{
	my $self=shift;
	my $hidewarning=shift;
	#if($self->outputs->{demux_filelist}){
	#	return $self->outputs->{demux_filelist};
	#}
	my $odir = $self->inputs->{out_dir};
	#my $ffile = basename($self->inputs->{forward_file});
	#my $rfile = basename($self->inputs->{reverse_file});
	my $lu = $self->_get_dlu();
	my %ret;
	foreach my $v(values(%$lu)){
		#my @fs=("$odir/$v/$ffile","$odir/$v/$rfile");
		my @fs=("$odir/$v/$v.fwd.fq.gz","$odir/$v/$v.rev.fq.gz");
		foreach my $tf(@fs){
			$self->msg("[WARNING] missing expected file $tf - skipping barcode") if(!-e $tf && !$hidewarning);
			## we skip this 
			next;
		}
		$ret{$v}=\@fs;
		#last;
	}
	delete $ret{'NOT_CALLED'};
	$self->outputs->{demux_filelist}=\%ret;
	return \%ret;
}

sub _get_dlu{
	my $self=shift;
	my %ret;
	my $tfile = $self->inputs->{tag_file};
	open(IN,$tfile) || die "Cannot open $tfile\n";
	while(<IN>){
		chomp;
		next if /^(Distance|Format)\s/;
		my ($dest,$bc)=split(/\s+/,$_);
		$ret{$bc}=$dest;
	}
	$ret{NC}='NOT_CALLED';
	close(IN);
	return \%ret;
}

sub _skip_message{
	my $self=shift;
	return "[".ref($self)."]: The output path ".$self->inputs->{out_dir}."/pipeline/ already exists and contains files. Would you like to skip this step ";
}

sub run_conditions{
	#return 0;
	my $self=shift;
	my $odir = $self->inputs->{out_dir};
	my $flag=0;
	if(-d $odir){
		opendir (DIR, $odir) or die $!;
		if(scalar(grep{-d "$odir/$_" && $_ !~/^[\.]+/}readdir(DIR))>0){
			$flag=1;
		}
		closedir(DIR);
	}
	return $flag;
}

sub gunzip{
	my $self = shift;
	my $file = shift;
	(my $ofile = basename($file))=~s/\.gz$//;
	my $out =  $self->inputs->{out_dir}."/$ofile";
	my $cmd = 'gunzip -c '.$file." > $out";
	return [$cmd,$out];
}
1;



