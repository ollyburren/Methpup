#!/usr/bin/perl

## Name: Methpup.pl
## Function: Run pipeline to call methylation status on multiplexed, pulled down PE illumina sequencing data
## Author: Xin Yang (bash), Olly Burren (Perl conversion)
## Date: 07/04/2014

#########
#PRAGMA #
#########

use strict;
use FindBin qw($Bin) ;
use lib ("$Bin/perl/lib","$Bin/../macd/lib");
use Getopt::Long;
use Config::IniFiles;
use Macd;
use Methpup::Runnable;
use Data::Dumper;

############
#CONSTANTS #
############

use constant PHRED=>33;
use constant LINFORMIN=>5;
use constant LINREVMIN=>5;
use constant MINLEN=>10;
use constant MAXLEN=>222;

##########
#GLOBALS #
##########

my $DEP_FILE = "$Bin/dependencies.tab";
my %DEP_LU;
if(! -e $DEP_FILE){
	print "[ERROR] Cannot find $DEP_FILE - please check installation\n";
	exit(1);
}else{
	open(DEP,"$DEP_FILE") || die "[ERROR] Cannot read $DEP_FILE please check permissions\n";
	while(<DEP>){
		chomp;
		next if /(^parameter)|(^\s*$)/;
		my ($p,$n,$v,$url)=split("\t",$_);
		$DEP_LU{$p}={name=>$n,url=>$url,version=>$v};
	}
	close(DEP);
}

my $ERROR_FLAG=0;

my %DEFAULTS = (
	phred=>PHRED,
	formin=>LINFORMIN,
	revmin=>LINREVMIN,
	minol=>MINLEN,
	maxol=>MAXLEN);
	

#print Dumper(\%DEP_LU);
	


#######################################
#OPTION HANDLING AND INI FILE PARSING #
#######################################

my $USAGE=<<EOL;
perl $0 --[i]ni_file {-[v]erbose -[h]elp)

	ini_file - path to a suitable inifile  http://github.com/ollyburren/Methpup/ini/default.ini) for a template
	verbose - enable verbose output to STDERR - currently not implemented
	help - print this message.

For help and suggestions please contact Olly Burren (http://github.com/ollyburren)
EOL

my ($ini_file,$verbose,$help,$ERROR_FLAG);
GetOptions (
	'ini_file|i=s' => \$ini_file,
	'help|h' => \$help,
	'verbose|v' => \$verbose
	);

if($help){
	print STDERR $USAGE."\n";
	exit(1);
}
	
if(!-e $ini_file){
	print STDERR "Cannot find ini_file:$ini_file. See http://github.com/ollyburren/ini/default.ini for details\n";
	exit(1);
}

if($verbose){
	print STDERR "Verbose flag not yet implemented.\n";
}

##################
## CONFIGURATION #
##################

my $cfg = Config::IniFiles->new( -file => $ini_file );

###################
## SETUP PROGRAMS #
###################

my %PROGRAMS;
foreach my $p($cfg->Parameters('PROGRAMS')){
	my %pro = %{$DEP_LU{$p}};
	if(my $set = $cfg->val('PROGRAMS',$p)){
		if(-e $set){
			$PROGRAMS{$p}=$set;
		}else{
			
			print "[ERROR] Can't find location for $pro{name}: $set, $pro{version}. Try $pro{url}\n";
		}
	}else{
		print "[ERROR] Please install: $p($pro{name}-$pro{version}) - $pro{url}. Then update [PROGRAMS] section of ini\n";
	}
}

########################
## SETUP PROJECT FILES #
########################

my %PROJECT;
##grab file based parameters 
foreach my $p($cfg->Parameters('PROJECT files')){
	my $set = $cfg->val('PROJECT files',$p);
	if(! -e $set){
		print STDERR "[ERROR] Cannot access path for $p:$set\n";
		$ERROR_FLAG=1;
	}else{
		$PROJECT{$p}=$cfg->val('PROJECT files',$p);
	}
}

###########################
## SETUP PROJECT SETTINGS #
###########################

foreach  my $p($cfg->Parameters('PROJECT settings')){
	if(my $set = $cfg->val('PROJECT settings',$p)){
		$PROJECT{$p}=$set
	}elsif($set = $DEFAULTS{$p}){
		print "$p: using default $set\n";
		$PROJECT{$p}=$set
	}else{
		print STDERR "[ERROR] no value set for [PROJECT settings]$p\n";
		$ERROR_FLAG=1;
	}
}
my $MACD_CONF_FILE = $PROJECT{mac_d_conf_file};
my $BASE_DIR = $PROJECT{base_dir};
my $PROJECT = $PROJECT{name};
my $PROJECT_DIR = "$BASE_DIR/$PROJECT";
my $QUEUE_ENGINE = $PROJECT{queue_engine};
my $DRIVER = Macd::GRIDDriverFactory->instantiate($QUEUE_ENGINE,inifile => $MACD_CONF_FILE);
my $F_FASTQ = $PROJECT{fastq_input_fwd};
my $R_FASTQ = $PROJECT{fastq_input_rev};
my $TAGFILE = $PROJECT{novobarcode_index_file};

##############################
## SET UP EACH GENE SETTINGS #
##############################
#my %GENES;
#foreach my $sect(grep{/GENE/}$cfg->Sections()){
#	(my $gname = $sect)=~s/GENE (.*)/\1/;
#	unless($GENES{$gname}{forseq} = $cfg->val($sect,'forseq')){
#		print STDERR "[ERROR] Need a forseq for linker for $sect\n";
#	}
#	unless($GENES{$gname}{revseq} = $cfg->val($sect,'revseq')){
#		print STDERR "[ERROR] Need a revseq for linker for $sect\n";
#	}
#	unless($GENES{$gname}{linker_length} = $cfg->val($sect,'linker_length')){
#		print STDERR "[ERROR] Need a linker_length for linker for $sect\n";
#	}
#}

#print Dumper(\%PROGRAMS)."\n";
#print Dumper(\%PROJECT)."\n";
#print Dumper(\%GENES)."\n";

################
## DEMULTIPLEX #
################

my $odir = "$PROJECT_DIR/pipeline/";

my $demux = Methpup::Runnable::DeMultiplex->new(
		log_dir=>"$PROJECT_DIR/log/demultiplex/",
		binary=>$PROGRAMS{novobc_bin},
		debug_flag=>1,
		#env_settings=>$RLIB_ENV,                         
		macd_driver=>$DRIVER,
		inputs=>{
			forward_file=>$F_FASTQ,
			reverse_file=>$R_FASTQ,
			tag_file=>$TAGFILE,
			#output_dir=>$odir
			
		},
		outputs=>{
			out_dir=>$odir
		}
)->run_step;


################
## TRIMMOMATIC #
################

my $trimmo = Methpup::Runnable::Trimmomatic->new(
	log_dir=>"$PROJECT_DIR/log/trimmomatic/",
	binary=>$PROGRAMS{trimmo_jar},
	debug_flag=>1,          
	macd_driver=>$DRIVER,
	inputs=>{
		linker_length=>$PROJECT{linker_length},### SORT THIS OUT LATER
		phred=>$PROJECT{phred},
		pe_file_list=>$demux->demux_filelist('hidewarning'),
	},
	outputs=>{
		out_dir=>$odir
	})->run_step;
#die Dumper($trimmo);
#############
## CUTADAPT #
#############


my $cutadapt = Methpup::Runnable::Cutadapt->new(
	log_dir=>"$PROJECT_DIR/log/cutadapt/",
	binary=>$PROGRAMS{cutadapt_bin},
	debug_flag=>1,          
	macd_driver=>$DRIVER,
	inputs=>{
		forseq=>$PROJECT{forseq},
		revseq=>$PROJECT{revseq},
		linformin=>$PROJECT{formin},        
		linrevmin=>$PROJECT{revmin},
		trimmo_sub=>$trimmo->subdir,
		pe_file_list=>$trimmo->trim_filelist('hidewarning')
	},
	outputs=>{
		out_dir=>$odir
	})->run_step;
#die(Dumper($cutadapt));
##########
## FLASH #
##########

my $flash = Methpup::Runnable::FLASH->new(
	log_dir=>"$PROJECT_DIR/log/flash/",
	binary=>$PROGRAMS{flash_bin},
	debug_flag=>1,          
	macd_driver=>$DRIVER,
	inputs=>{
		maxol=>$PROJECT{maxol},
		minol=>$PROJECT{minol},
		cutadapt_sub=>$cutadapt->subdir,
		pe_file_list=>$cutadapt->cutadapt_filelist('hidewarning')
	},
	outputs=>{
		out_dir=>$odir
	})->run_step;

###########
## BOWTIE #
###########

my $bowtie = Methpup::Runnable::Bowtie->new(
	log_dir=>"$PROJECT_DIR/log/bowtie/",
	binary=>$PROGRAMS{bowtie_bin},
	debug_flag=>1,          
	macd_driver=>$DRIVER,
	inputs=>{
		index_dir=>"$PROJECT_DIR/bowtie2_index/",
		ref_seq_file=>$PROJECT{'ref_seq_file'},
		phred=>$PROJECT{phred},
		flash_sub=>$flash->subdir,
		pe_file_list=>$flash->filelist('hidewarning')
	},
	outputs=>{
		out_dir=>$odir
	})->run_step;

sleep(5);


###########################
## CALL METHYLATION SITES #
###########################

my $cm = Methpup::Runnable::CallMeth->new(
	log_dir=>"$PROJECT_DIR/log/callmeth/",
	binary=>"$Bin/perl/call_methsites_from_sam.pl",
	debug_flag=>1,          
	macd_driver=>$DRIVER,
	inputs=>{
		meth_site_file=>$PROJECT{meth_site_file},
		bowtie_sub=>$bowtie->subdir,
		samfiles=>$bowtie->filelist('hidewarning')
	},
	outputs=>{
		out_dir=>$odir
	})->run_step;