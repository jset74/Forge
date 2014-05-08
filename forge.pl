#!/usr/bin/env perl

=head1 NAME

forge.pl - Functional element Overlap analysis of the Results of GWAS Experiments.

=head1 SYNOPSIS

forge.pl options (-f file) (-snps snplist)

=head1 DESCRIPTION

Analyse a set of SNPs for their overlap with DNase 1 hotspots compared to matched background SNPs. Identifies enrichment in DHS by tissue and plots graphs and table to display. Arbitrarily a minumum of 20 SNPs is required.  Note that if no SNPs are given the script will run on Pulmonary funciotn GWAS as an example output.

Several outputs are made.

A straight base R graphics pdf chart of the data.

A polychart (https://github.com/Polychart/polychart2) interactive javascript graphic using rCharts (http://ramnathv.github.io/rCharts/).

A dimple (http://dimplejs.org) d3 interactive graphic using rCharts.

A table using the Datatables (https://datatables.net) plug-in for the jQuery Javascript library, again accessed through rCharts.

In each of the graphics the colouring should be consistent. Blue (Z < 2.58), light red or pink (2.58 =< Z < 3.39), red or dark red (Z >= 3.39 ) for the 99% and 99.9% cIs. Or whatever other thresholds are specified.

Forge functions, plotting options and stats are provided by Forge::Forge, Forge::Plot and Forge::Stats modules.

=head1 OPTIONS

=over

=item B<data>

Data set to analyse. Either ENCODE data ('encode') or Roadmap Epigenome data ('erc'). erc by default.

=item B<peaks>

Use peaks instead of hotspots. Peaks are more stringent DNase1 peaks calls representing DNase hypersensitive sites, rather than hotspots which are regions of generalised DNAs1 sensitivity or open chromatin. Default is to use hotspots.

=item B<bkgd>

Specify whether the background matches should be picked from general set of arrays used in GWAS ('gwas') or from the Illumina_HumanOmni2.5 ('omni'). General GWAS arrays include

Affy_GeneChip_100K_Array
Affy_GeneChip_500K_Array
Affy_SNP6
HumanCNV370-Quadv3
HumanHap300v2
HumanHap550v3.0
Illumina_Cardio_Metabo
Illumina_Human1M-duoV3
Illumina_Human660W-quad

Defaults to 'gwas'. In both cases SNPs have to be on the arrays AND in the 1000 genomes phase 1 integrated call data set at phase1/analysis_results/integrated_call_sets.

=item B<label>

Supply a label that you want to use for the plotting titles, and filenames.

=item B<f>

Supply the name of a file containing a list of SNPs. Format must be given by the -format flag. If not supplied the analysis is performed either on snps provided as rsids in a comma separated list through the snps option or on a set of data from a gwas study on Pulmonary_function (http://www.ncbi.nlm.nih.gov/pubmed/21946350, http://www.ncbi.nlm.nih.gov/pubmed/20010835 and http://www.ncbi.nlm.nih.gov/pubmed/20010834). Note that 20 SNPs are required at a minimum.

=item B<snps>

Can provide the snps as rsids in a comma separated list.

=item B<min_snps>

Specify the minimum number of SNPs to be allowed. Default is 5 now are using binomial test.

=item B<thresh>

Alter the default binomial p value thresholds. Give a comma separate list of two e.g. 2.58,3.39 for the defaults

=item B<format>

if f is specified, specify the file format as follow:

rsid = list of snps as rsids each on a separate line. Optionally can add other fields after the rsid which are ignored, unless the pvalue filter is specified, in which case Forge assumes that thesecond field is the minus log10 pvalue

bed  = File given is a bed file of locations (chr\tbeg\tend) aka Personal Genome SNP format.  bed format should be 0 based and the chromosome should be given as chrN. Hoever will also accept chomosomes as just N (ensembl) and 1-based format where beg and end are the same

vcf = File given is a vcf file.

tabix = File contains SNPs in tabix format.

ian = 1-based chr\tbeg\tend\trsid\tpval\tminuslog10pval

=item B<filter>

Set a filter on the SNPs based on the -log10 pvalue.  This works for files in the 'ian' or 'rsid' format. Give a value as the lower threshols and only SNPs with -log10 pvalues >= to the threshold will be analysed. Defaiult is no filtering.

=item B<bkgrd>

Output background stats for investigation.

=item B<reps>

THe number of background matching sets to pick and analyse. Default 1000.

=item B<ld>

Apply filter for SNPs in LD at eithr r2 >= 0.8 ("high LD"), or r2 >= 0.1 ("independent SNPs"). Specify ld 0.8, or ld 0.1. Default is to filter at r2 >= 0.8.  With ld filter specified, forge will report SNPs removed due to LD with another SNP in the list and will randomly pick one for each LD block.

To turn off LD fitlering specify -nold

=item B<nold>

Turn off LD filtering.

=item B<noplot>

Just make the data file, don't plot.

=item B<help|h|?>

Print a brief help message and exits.

=item B<man|m>

Print this perldoc and exit.

=back

=head1 LICENCE

forge.pl Functional analysis of GWAS SNPs

Copyright (C) 2013  EMBL - European Bioinformatics Institute

This program is free software: you can redistribute it and/or modify it under the terms of
the GNU General Public License as published by the Free Software Foundation, either version 3
of the License, or (at your option) any later version. This program is distributed in the hope
that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. Neither
the institution name nor the name forge.pl can be used to endorse or promote products derived from
this software without prior written permission. For written permission, please contact
dunham@ebi.ac.uk. Products derived from this software may not be called forge.pl nor may forge.pl
appear in their names without prior written permission of the developers. You should have received
a copy of the GNU General Public License along with this program.  If not, see http://www.gnu.org/licenses/.

=head1 AUTHOR

Ian Dunham, EBI

=head1 CONTACT

Ian Dunham <dunham@ebi.ac.uk>

=cut

use strict;
use 5.010;
use warnings;


use Config::IniFiles;
use Cwd;
use Getopt::Long;
use File::Basename;
use Pod::Usage;

use lib "lib";

my ($bkgd, $data, $peaks, $label, $file, $format, $min_snps, $bkgrdstat, $noplot, $reps, $help, $man, $thresh, $ld, $nold, $filter, @snplist);
my ($datadir, $dsn, $rlibs, $user, $pass, $output);


GetOptions (
    'data=s'     => \$data,
    'peaks'      => \$peaks,
    'bkgrd'      => \$bkgrdstat,
    'label=s'    => \$label,
    'f=s'        => \$file,
    'format=s'   => \$format,
    'snps=s'     => \@snplist,
    'min_snps=i' => \$min_snps,
    'noplot'     => \$noplot,
    'reps=i'     => \$reps,
    'thresh=s'   => \$thresh,
    'ld=f'       => \$ld,
    'nold'       => \$nold,
    'filter=f'   => \$filter,
    'help|h|?'   => \$help,
    'man|m'      => \$man,
    'datadir=s'  => \$datadir, # e.g /usr/local/forge/
    'dsn=s'      => \$dsn, # e.g dbi:SQLite:/usr/local/forge/forge.db
    'rlibs=s'     => \$rlibs, # e.g /usr/local/libs/R
    'out=s'      => \$output, # e.g /tmp
);

pod2usage(1) if ($help);
pod2usage(-verbose => 2) if ($man);

my $dirname = dirname(__FILE__);
my $cfg = Config::IniFiles->new( -file => "$dirname/forge.ini" );

$datadir ||= $cfg->val('Data', 'datadir');
$dsn ||= $cfg->val('Data', 'dsn');
$user ||= $cfg->val('Data', 'user');
$pass ||= $cfg->val('Data', 'pass');

$rlibs ||= $cfg->val('Files', 'rlibs');

$output ||= $cfg->val('Files', 'output') || getcwd;

use Bio::Analysis::Forge;

my $forge = Bio::Analysis::Forge->new(
    {
        'output'     => $output,
        'r_libs'     => $rlibs,
        'datadir'    => $datadir, 
	'dsn'        => $dsn,
	'user'       => $user,
	'pass'       => $pass,
	'noplot'     => $noplot,
	'label'      => $label,
        'data'       => $data,
        'peaks'      => $peaks,
        'bkgrdstat'  => $bkgrdstat,
	'reps'       => $reps,
	'thresh'     => $thresh,
	'min_snps'   => $min_snps,	
	'nold'       => $nold,
	'ld'         => $ld,
    });

$forge or die "System error: Failed to initialise Forge analysis.";

# snps need to come either from a file or a list
my @snps;

# A series of data file formats to accept.
if (defined $file){
    @snps = @{$forge->parse_input($file, $format) || []};
}
elsif (@snplist){
    @snps = split(/,/,join(',',@snplist));
}
else{
# Test SNPs from gwascatalog_21_03_2012  Pulmonary_function.snps.bed
# If no options are given it will run on the default set of SNPs
    warn "No SNPs given, so running for example on Pulmonary function set from the GWAS catalogue.";
    @snps = qw(rs2865531 rs2395730 rs12914385 rs11168048 rs1529672 rs357394 rs13147758 rs3769124 rs2647044 rs12504628 rs1541374 rs2869967 rs1928168 rs3094548 rs3867498 rs9978142 rs4762767 rs6903823 rs11172113 rs9310995 rs2571445 rs2070600 rs11727189 rs3734729 rs2906966 rs1036429 rs16909898 rs3995090 rs12477314 rs2544527 rs2284746 rs993925 rs2277027 rs1344555 rs1455782 rs2855812 rs2838815 rs11001819 rs12716852 rs2798641 rs4129267 rs7068966 rs12899618 rs153916 rs1551943 rs730532 rs1980057 rs3820928 rs2036527 rs10516526 rs2857595 rs3817928 rs310558 rs808225 rs12447804);
}

if (my $jobid = $forge->run(\@snps)) {
    warn "Done . Results are in $output/$jobid \n";
} else {
    warn "Error : ",$forge->error, "\n";
}

