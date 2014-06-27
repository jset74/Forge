package Bio::Analysis::Forge;
use strict;

use DBI;
use Sort::Naturally;
use Storable;
use IO::Handle;
use File::Path qw(make_path);
use Data::Dumper;
use Scalar::Util qw(looks_like_number);
use Forge::Stats;
use Forge::Plot;
use Forge::Forge; # ld_filter . need to use all functions from plot and forge

our $VERSION = '1.1';
# percentile bins for the bkgrd calculations. This is hard coded so there are enough SNPs to choose from, but could later be altered.
my $per = 10;
# number of sets to analyse for bkgrd. Again currently hardcoded to 1000
my $reps = 1000;

# internal
sub  tissues :lvalue { $_[0]->{'_tissues'}; }
sub  cells :lvalue { $_[0]->{'_cells'}; }
sub  jobid :lvalue { $_[0]->{'_jobid'}; }
sub  error :lvalue { $_[0]->{'_error'}; }
sub  sh :lvalue { $_[0]->{'_sh'}; }
sub  dbh :lvalue { $_[0]->{'_dbh'}; }
sub  logfile :lvalue { $_[0]->{'_logfile'}; }

# global params that can be passed to the constructor
sub  dsn :lvalue { $_[0]->{'dsn'}; }
sub  datadir :lvalue { $_[0]->{'datadir'}; }
sub  user :lvalue { $_[0]->{'user'}; }
sub  pass :lvalue { $_[0]->{'pass'}; }
sub  output :lvalue { $_[0]->{'output'}; }
sub  analysis :lvalue { $_[0]->{'data'}; }
sub  bkgd :lvalue { $_[0]->{'bkgd'}; }
sub  label :lvalue { $_[0]->{'label'}; }
sub  noplot :lvalue { $_[0]->{'noplot'}; }
sub  r_libs :lvalue { $_[0]->{'r_libs'}; }

sub  t1 :lvalue { $_[0]->{'tmin'}; }
sub  t2 :lvalue { $_[0]->{'tmax'}; }
sub  repetitions :lvalue { $_[0]->{'repetitions'}; }
sub  bkgrdstat :lvalue { $_[0]->{'bkgrdstat'}; }

sub r2 :lvalue { $_[0]->{'r2'}; }
sub nold :lvalue { $_[0]->{'nold'}; }
sub ld :lvalue { $_[0]->{'ld'}; }

# Constructor - requires location of the forge data as DSN string and directory with forge settings
sub new {
    my ($class, $params) = @_;

    my @chars = ("A".."Z", "a".."z", "0".."9");
    my $str = '';
    $str .= $chars[rand @chars] for 1..4;
    $params->{_jobid} = "forge/$str/".time;


    $params->{min_snps} ||= 5;
    $params->{repetitions} ||= 100;

    $params->{data} ||= 'erc';
    $params->{bkgd} ||= 'gwas';
    $params->{user} ||= '';
    $params->{pass} ||= '';
    $params->{tmin} ||= 0.01;
    $params->{tmax} ||= 0.05;

# Set r2 LD thresholds
    unless (defined $params->{nold}){
	my $ld = $params->{ld} || 0.8;
	unless ($ld == 0.1 || $ld == 0.8){
	    warn "You have specified LD filtering, but given an invalid value $ld. the format is ld 0.0, ld 0.8, or ld 0.1";
	    return undef;
	}
	my $r2;
	($r2 = $ld) =~ s /\.//;
	$params->{r2} = "r".$r2;
    }

    unless (looks_like_number($params->{tmin}) && looks_like_number($params->{tmax})){
	warn "You must specify numerical thersholds";
	return undef;
    }

    my $dsn = $params->{dsn};

    unless ($dsn) {
	warn "ERROR: Need the Forge data connection string, e.g. dbi:SQLite:/usr/local/forge/forge.db";
	return undef;
    }

    unless ($params->{datadir}) {
	warn "ERROR: Need a location of Forge data parameters ( omni.snp_params.10 )";
	return undef;
    }

    unless (-e $params->{datadir}) {
	warn "ERROR: Forge data folder does not exist";
	return undef;
    }


    my $dbh = DBI->connect($dsn, $params->{user}, $params->{pass});
    unless ($dbh) {
	warn "ERROR: System error: $DBI::errstr ";
	return;
    }
    $params->{_dbh} = $dbh;

    return bless $params;
}

# function to submit a job in a new thread 
# for deployment in a website
sub submit {
    my ($self, $snps) = @_;
    my $pid = fork();
    if (not defined $pid) {
	warn "FORK ERROR: Resources not avilable.";
    } elsif ($pid == 0) {
	$self->run($snps);
	if ($self->error) {
	    $self->debug($self->error);
	}
     	exit(0);
    } else {
	#    waitpid($pid,0);
    }  
    return $self->jobid;
}

# the main function that does all the work
sub run {
    my ($self, $snps) = @_;

    my $t = time;
# set up folder where all results and logs will go    
    my $folder = $self->output.'/'.$self->jobid;
    make_path($folder);
    unless (-e $folder) {
	$self->error = "Failed to create $folder : $!";
	return;
    }
    my $sfh;
# this is where will print the job status 
# so the website can check on the progress
    my $sfilename = "$folder/STATUS";
    if (open $sfh, ">$sfilename") { 
	$sfh->autoflush(1);
	$self->sh = $sfh;
	print $sfh "LOAD#0\n";
    } else {
	$self->error =  "ERROR: Cannot open $sfilename: $!"; 
	return;
    }

# all info from the analysis will go into log
    $self->logfile = "$folder/LOG";


# Remove redundancy in the input

    my %nonredundant;
    foreach my $snp (@$snps){
	$nonredundant{$snp}++;
    }
    foreach my $snp (keys %nonredundant){
	if ($nonredundant{$snp} > 1) {
	    warn "$snp is present " . $nonredundant{$snp} . " times in the input. Analysing only once."
	}
    }

    my @snps = keys %nonredundant;

# Perform ld filter unless -nold is specified.

    my ($ld_excluded, $output, $input);
    unless (defined $self->nold) {
	$input = scalar @snps;
	($ld_excluded, @snps) = ld_filter(\@snps, $self->r2, $self->dbh);
	$output = scalar @snps;
    }


    eval {
	$self->get_cells();

	my $rows = $self->get_bits(\@snps);

# unpack the bitstrings and store the overlaps by cell.
	my $test = $self->process_bits($rows);

	if ($self->bkgrdstat){
	    open my $bfh, ">", "$folder/bkgrd.stats" or die "cannot open $folder/bkgrd.stats";
	    my (@maf, @tss, @gc);
	    foreach my $rsid (keys %{$$test{'SNPS'}}){
		my ($maf, $tss, $gc) = split "\t", $$test{'SNPS'}{$rsid}{'PARAMS'};
		push @maf, $maf;
		push @tss, $tss;
		push @gc, $gc;
	    }
	    print $bfh join("\t", "test", "maf", @maf);
	    print $bfh join("\t", "test", "tss", @tss);
	    print $bfh join("\t", "test", "gc", @gc);
	}
    

# Identify SNPs that weren't found and warn about them.
	my @missing;
	foreach my $rsid (@$snps){
	    if (defined $self->ld) {
		next if exists $$ld_excluded{$rsid}; # if the snps are not in the 1000 genomes set they will not be found by the LD filter, so we have to recheck here.
	    }
	    unless (exists $$test{'SNPS'}{$rsid}){
		push @missing, $rsid;
	    }
	}

# only pick background snps matching snps that had bitstrings originally.
	my @foundsnps = keys %{$$test{'SNPS'}};
	$self->debug("A total of " . scalar @foundsnps . " SNPs have been analysed.\n" ) ;
	
	if (scalar @missing > 0) {
	    $self->debug("\nThe following " . scalar @missing . " SNPs have not been analysed because they were not found in the 1000 genomes phase 1 integrated call data\n");
	    $self->debug(join("\n", @missing) . "\n");
	}

	if (defined $self->ld) {
	    if ($output < $input) {
		self->debug("For ".$self->label.", $input SNPs provided, " . scalar @snps . " retained, " . scalar @missing . " not analysed, "  . scalar(keys %$ld_excluded) . " LD filtered at ".$self->ld.".");
	    }
	}

	
# identify the gc, maf and tss, and then make bkgrd picks
	my $picks = $self->match(\%$test);

# for bgrd set need to get distribution of counts instead
# make a hash of data -> cell -> bkgrd-Set -> overlap counts
	my %bkgrd; #this hash is going to store the bkgrd overlaps

# Get the bits for the background sets and process
	my $num = scalar(keys %$picks);
	my $ic = 0;

# Get the bits for the background sets and process
	my $backsnps;

	foreach my $bkgrd (keys %{$picks}){
	    $rows = $self->get_bits(\@{$$picks{$bkgrd}});
	    $backsnps += scalar @$rows; #$backsnps is the total number of background SNPs analysed

	    unless (scalar @$rows == scalar @foundsnps){
		$self->debug("Background " . $bkgrd . " only " . scalar @$rows . " SNPs out of " . scalar @foundsnps . "\n");
	    }
	    
	    my $result = $self->process_bits($rows);
	    foreach my $cell (keys %{$$result{'CELLS'}}){
		push @{$bkgrd{$cell}}, $$result{'CELLS'}{$cell}{'COUNT'}; # accumulate the overlap counts by cell
	    }
	    $self->status( sprintf("RUNNING#%d\n", $ic++ / $num * 100) ) ;


	    if ($self->bkgrdstat){
		open my $bfh, ">", "$folder/bkgrd.stats" or die "cannot open $folder/bkgrd.stats";
		my (@maf, @tss, @gc);
		foreach my $rsid (keys %{$$test{'SNPS'}}){
		    my ($maf, $tss, $gc) = split "\t", $$test{'SNPS'}{$rsid}{'PARAMS'};
		    push @maf, $maf;
		    push @tss, $tss;
		    push @gc, $gc;
		}
		print $bfh join("\t", "test", "maf", @maf);
		print $bfh join("\t", "test", "tss", @tss);
		print $bfh join("\t", "test", "gc", @gc);
	    }
	}

	$self->dbh->disconnect();
	$self->status( "RESULTS#0\n");

#Having got the test overlaps and the bkgd overlaps now calculate Zscores and output the table to be read into R for plotting.
	my $filename = $folder."/chart.tsv";

	open my $ofh, ">", "$filename";
	unless ($ofh) {
	    $self->error = "ERROR: Cannot open $filename: $!"; 
	    return;
	}
	
	print $ofh join("\t", "Zscore", "Cell", "Tissue", "File", "SNPs", "Number", "Accession") ."\n";
	my $n =1;
	my $tissues = $self->tissues;
	my $cells = $self->cells;

	my %tissuecount;
	foreach my $cell (keys %$tissues){
	    my $tissue = $$tissues{$cell}{'tissue'};
	    $tissuecount{$tissue}++;
	}

	my $tissuecount = scalar keys %tissuecount;
	my $t1 = $self->t1/$tissuecount; # bonferroni correction by number of tissues
	my $t2 = $self->t2/$tissuecount;
	
	$self->t1 = -log10($t1);
	$self->t2 = -log10($t2);
	


	my $pos = 0;
	my $snpcount = scalar @foundsnps;

	open my $bfh, ">", "$folder/background.tsv" or die "Cannot open $folder/background.tsv";

	
	foreach my $cell (sort {ncmp($$tissues{$a}{'tissue'},$$tissues{$b}{'tissue'}) || ncmp($a,$b)} @$cells){ # sort by the tissues alphabetically (from $tissues hash values)
    # ultimately want a data frame of names(results)<-c("Zscore", "Cell", "Tissue", "File", "SNPs")
	    print $bfh join("\t", @{$bkgrd{$cell}});

	    my $teststat = $$test{'CELLS'}{$cell}{'COUNT'}; #number of overlaps for the test SNPs

    # binomial pvalue, probability of success is derived from the background overlaps over the tests for this cell
    # $backsnps is the total number of background SNPs analysed
    # $tests is the number of overlaps found over all the background tests
	    my $tests;
	    foreach (@{$bkgrd{$cell}}){
		$tests+= $_;
	    }
	    my $p = sprintf("%.6f", $tests/$backsnps);

    # binomial probability for $teststat or more hits out of $snpcount snps
    # sum the binomial for each k out of n above $teststat
	    my $pbinom;
	    foreach my $k ($teststat .. $snpcount){
		$pbinom += binomial($k, $snpcount, $p);
	    }
	    if ($pbinom >1) {
		$pbinom = 1;
	    }
	    $pbinom = -log10($pbinom);
    # Z score calculation
	    my $mean = mean(@{$bkgrd{$cell}});
	    my $sd = std(@{$bkgrd{$cell}});
	    my $zscore;
	    if ($sd == 0){
		$zscore = "NA";
	    }
	    else{
		$zscore = sprintf("%.3f", ($teststat-$mean)/$sd);
	    }
	    if ($pbinom >= $self->t2){
		$pos++;
	    }
	    my $snp_string = "";
	    $snp_string = join(",", @{$$test{'CELLS'}{$cell}{'SNPS'}}) if defined $$test{'CELLS'}{$cell}{'SNPS'}; # This gives the list of overlapping SNPs for use in the tooltips. If there are \
a lot of them this can be a little useless
    my ($shortcell, undef) = split('\|', $cell); # undo the concatenation from earlier to deal with identical cell names.
	    print $ofh join("\t", $zscore, $pbinom, $shortcell, $$tissues{$cell}{'tissue'}, $$tissues{$cell}{'file'}, $snp_string, $n, $$tissues{$cell}{'acc'}) . "\n";
	    $n++;
	}


	close $ofh;

	my $cellcount = scalar @$cells;
	my $fdr = $self->fdr($pos, $snpcount, $cellcount);
	$self->debug("$pos positive lines at FDR = $fdr at Z >= 3.39\n");


	unless (defined $self->noplot){
	    $self->label ||= $self->jobid;
	    #Plotting and table routines
	    $self->status( "RESULTS#10\n" );
	    $self->Chart($filename); # basic pdf plot
	    $self->status( "RESULTS#20\n" );
	    $self->dChart($filename); # rCharts Dimple chart
	    $self->status("RESULTS#60\n");
	    $self->table($filename); # Datatables chart
	}
	$self->status( "RESULTS#100\n" );
    };

   
    $self->status ( "COMPLETE#100\n" );

    if ($!) {
	$self->error("ERROR:$!");
    }
	
    $sfh->close();    
    $self->debug(sprintf("Time taken: %ds \n",  time - $t));

    return $self->jobid;
}


### Subroutines ###

sub match{
    my $self = shift;
    # identifies the bins that each of the snps lies in, and then picks 100 matching SNPs.

    my $snps = shift; # ahash that contians the found snps

    my $bkgd = $self->bkgd;
    my $datadir = $self->datadir;

    my ($bins, $params, %bins, %params);
    # load up the stored hashes that contain the bins of snps by gc, maf, and tss distance. There is one for each of the bkgd set possibilities(gwas or omni).
    # These are precalculated according to the parameters that are hard coded above.
    # the hash to load is defined by the bkgd option - defualts to 'gwas'
    if ($bkgd =~ /omni/i){
        $bins = $datadir . "/omni.snp_bins.$per";
        $params = $datadir . "/omni.snp_params.$per";
    }
    else{
        $bins = $datadir . "/snp_bins.$per";
        $params = $datadir . "/snp_params.$per";
    }

    if (-e $bins && -e $params){
        %bins = %{ retrieve($bins) };
        %params = %{ retrieve($params)};
    }
    my (%picks);

    my $num2 = keys %{$$snps{'SNPS'}};
    my $iic = 0;

    foreach my $rs (keys %{$$snps{'SNPS'}}){
	$self->status (sprintf("LOAD#%d\n", $iic++ / $num2 * 80 + 20));
        srand;
        my ($maf, $tss, $gc) = split("\t", join("\t", $$snps{'SNPS'}{$rs}{'PARAMS'}));
        #$rs is the test snp, $rsid os the matched snp.
        my ($i, $j, $k) = assign ($gc, $tss, $maf, \%params);

        my $range = scalar @{$bins{$i}{$j}{$k}};
        for (my $n = 1; $n <= $reps; $n++) {
            my ($snp_string, $rsid);
            while (1){
                my $pick = int(rand($range));
                $snp_string = ${$bins{$i}{$j}{$k}}[$pick]; #pick the $pick'th element in the array as the chosen SNP "
                (undef, undef, undef, $rsid) = split /\t/,  $snp_string;
                last unless $rsid eq $rs; # must not pick the test snp itself.
            }
            push @{$picks{$n}}, $rsid; # each $n array is a set of snps matching the test set/
        }
    }
    return \%picks;
}

sub process_bits{
    # Processes the bitstrings to get a count of overlaps for each cell type.
    my ($self, $rows, $analysis) = @_;
    my $cells = $self->cells;
    my $analysis ||= $self->analysis;

    my %test;
    foreach my $row (@{$rows}){
        my ($location, $rsid, $sum, $bit, $maf, $tss, $gc);
        if ($analysis eq "erc"){
            ($location, $rsid, undef, undef, $bit, $sum, $maf, $tss, $gc) = @$row;
        }
       else{
            ($location, $rsid, $bit, $sum, undef, undef, $maf, $tss, $gc) = @$row;
       }

        $test{'SNPS'}{$rsid}{'SUM'} = $sum;
        $test{'SNPS'}{$rsid}{'PARAMS'} = join("\t", $maf, $tss, $gc);
        my @bits = split "", $bit;
        my $index = 0;
        foreach my $cell (@$cells){
            $test{'CELLS'}{$cell}{'COUNT'} += $bits[$index];
            push @{$test{'CELLS'}{$cell}{'SNPS'}}, $rsid if $bits[$index] == 1;
            $index++;
        }
    }
    return \%test;
}

#get the bitstrings for an array of snps from the forge db
sub get_bits {
    my $self = shift;
    my ($snps) = @_;
    my $args = "\'".join ("','", @$snps)."\'";
    my $sql = "SELECT * FROM bits WHERE rsid IN ($args)";
    return $self->dbh->selectall_arrayref($sql);
}

sub fetch_rsid{
    #gets the rsid for a SNP where a location is given
    my ($self, $loc) = @_;
    my $rsid = $self->dbh->selectrow_array("SELECT rsid FROM bits WHERE location = '$loc'") || return "no RSID match for $loc";
    return $rsid;
}

# get the cell list array and the hash that connects the cells and tissues
sub get_cells{
    my $self = shift;
    my $analysis = shift || $self->analysis;
    my $dbh = $self->dbh;

    # read the correct cell list based on data (erc -encode). Also gets the tissue names for the cells.

    my $table = join('_', "cells", $analysis);
    my $sth = $dbh->prepare("SELECT shortcell,tissue,file,acc FROM $table");
    $sth->execute();
    my $ver = $sth->fetchall_arrayref();
    $sth->finish();

    my ($cells, $tissues, $acc);
    foreach my $row (@$ver){
        my $cell = shift @$row;
        my $tissue = shift @$row;
        my $file = shift @$row;
        my $acc = shift @$row;
        $cell = "$cell|$file"; # Sometimes the same cell is used twice, with a differnt file so need to record separately (e.g. WI-38).
        push @$cells, $cell;
        $$tissues{$cell}{'tissue'} = $tissue; # this is the hash that is used to connect cells and tissues and ultimately provide the sorting order
        $$tissues{$cell}{'file'} = $file;
        $$tissues{$cell}{'acc'} = $acc;
    }
    #print Dumper %$tissues;
    $self->cells = $cells;
    $self->tissues = $tissues;

    return 1;
}


sub Chart{
    # This is the original code using standard R plot to generate a static pdf.
    my ($self, $filename) = @_;

    my $Rdir = $self->output.'/'.$self->jobid;
    my $tissues = $self->tissues;
    my $cells = $self->cells;
    my $label = $self->label;
    my $t1 = $self->t1;
    my $t2 = $self->t2;

    my $chart = "chart.pdf";
    my $rfile = $Rdir. "/chart.R";

    # make plot, first calculate where dividing lines are:
    my (@lines, @label_pos, @labels, @tissue_txt);
    my $n =1;
    my $last = '0';
    foreach my $cell (sort {ncmp($$tissues{$a}{'tissue'},$$tissues{$b}{'tissue'}) || ncmp($a,$b)} @$cells) {
        my $tissue = $$tissues{$cell}{'tissue'};
        push @tissue_txt, $tissue;
        unless ($tissue eq $last){
            push @lines, $n-0.5;
            push @labels, $tissue;
        }
        $n++;
        $last = $tissue;
    }

    my $length = scalar @tissue_txt;
    my $index =0;
    foreach my $value (@lines) {
        if (defined $lines[$index+1]){
            my $pos = $value + ($lines[$index+1]-$value)/2;
            push @label_pos, $pos;
        }
        else {
            my $pos = $value + ($length-$value)/2;
            push @label_pos, $pos;
        }
        $index++;
    }

    open my $rfh, ">", "$rfile";
    if (my $rlib = $self->r_libs) {
        print $rfh ".libPaths(\"$rlib\")\n";
    }

    print $rfh "setwd(\"$Rdir\")
results<-read.table(\"$filename\",header=TRUE,sep=\"\t\")
results\$Class<-cut(results\$Zscore, breaks =c(min(results\$Zscore), $t1, $t2, max(results\$Zscore)), labels=FALSE, include.lowest=TRUE) # 99.9 and 99% CIs 1, 2, 3
pdf(\"$chart\", width=22.4, height=7)
palette(c(\"steelblue3\",\"pink2\",\"red\"))
ymin1 = min(results\$Zscore, na.rm=TRUE)*1.1
ymax1 = max(results\$Zscore, na.rm=TRUE)*1.1
ymax = max(c(abs(ymin1),ymax1))
ymin = -ymax
par(mar=c(9,4,3,1)+0.1)
plot(results\$Zscore,ylab=\"Z score\",xlab=\"\",main=\"Proportion of SNPs, DNase1 sites (probably TF sites) which are present in cell lines for $label\",ylim=c(ymin,ymax), las=2, las=2, pch=19,col=results\$Class, xaxt='n')
axis(1, seq(1,length(results\$Cell)),labels=results\$Cell, las=2, cex.axis=0.7)
mtext(1,text=\"Cell\",line=7,cex=1.2)
#abline(h=-$t1, col=\"lightpink1\") # Z score of 2.58 = 99 % probability
abline(h=$t1, col=\"lightpink1\")
#abline(h=-$t2, col=\"lightpink1\", lty=2)
abline(h=$t2, col=\"lightpink1\", lty=2)
text(c(-1),$t1+0.2,c(\"Z = $t1\"),col=\"lightpink1\",adj=1,cex=0.8)
#text(c(-1),-$t1+0.16,c(\"1%\"),col=\"lightpink1\",adj=1,cex=0.8)
text(c(-1),$t2+0.2,c(\"Z = $t2\"),col=\"lightpink1\",adj=1,cex=0.8)
#text(c(-1),-$t2+0.16,c(\"0.1%\"),col=\"lightpink1\",adj=1,cex=0.8)
palette(\"default\")\n";

    foreach my $pos (@lines){
        print $rfh "lines(c($pos,$pos),c(-22,22),lty=6,col=\"#00000070\")\n" unless $pos == 0.5;
    }
    $index = 0;
    foreach my $tissue (@labels){
        print $rfh "text(c(" . $label_pos[$index] . "),ymax,c(\"" . $tissue . "\"),col=\"burlywood3\",adj=1,srt=90,cex=0.8)\n";
        $index++;
    }
    print $rfh "dev.off()\n";
#run the R code
    system "R --no-save --quiet --slave < $rfile";
}


sub rChart{
    # Makes a polcharts plot : note X axis labelling is problematical
    my ($self, $filename) = @_;
    my $chart = "rchart.htm";
    my $Rdir = $self->output.'/'.$self->jobid;
    my $rfile = "$Rdir/rChart.R";
    my $label = $self->label;
    my $data = $self->analysis;
    (my $lab = $label) =~ s/^\w/_/g;
    $lab .= ".".$self->analysis;
    my $t1 = $self->t1;
    my $t2 = $self->t2;

    open my $rcfh, ">", "$rfile";
    if (my $rlib = $self->r_libs) {
        print $rcfh ".libPaths(\"$rlib\")\n";
    }

    print $rcfh "setwd(\"$Rdir\")
results<-read.table(\"$filename\", header = TRUE, sep=\"\t\")
results\$Colour<- 0 + (results\$Zscore < $t2) + (results\$Zscore < $t1)  # 99.9 and 99% CIs
require(rCharts)
r1 <- rPlot(Zscore ~ Cell, data=results, color=\"bin(Colour, 0.25)\", type=\"point\", tooltip = \"function(item){ return (item.Zscore + '\\\\n' + item.Cell + '\\\\n' + item.Tissue + '\\\\n' + item.File + '\\\\n' + item.SNPs + '\\\\n' + item.Accession + '\\\\n')}\")
#r1\$guides(color=list(scale = list(type = \'gradient\', lower = \'\#CCC\', upper = \'\#000\'))) # optional code to make a grey scale
r1\$addParams(width = 2000, height=600, title=\"$label overlaps with $data DHS\")
ymin1 = min(results\$Zscore, na.rm=TRUE)*1.2
ymax1 = max(results\$Zscore, na.rm=TRUE)*1.2
ymax = max(c(abs(ymin1),ymax1))
ymin = -ymax
r1\$guides(x = list(numticks = length(unique(results\$Cell)), levels=results\$Cell), y = list(min = ymin, max = ymax))
r1\$save('$chart', cdn = F)
##r1\$show() #makes a temp file\n";

    system "R --no-save --quiet --slave < $rfile";
}

sub dChart{
    # Make simple interactive chart.
    my ($self, $filename) = @_;
    my $chart = "dchart.htm";
    my $Rdir = $self->output.'/'.$self->jobid;
    my $rfile = "$Rdir/dChart.R";
    my $label = $self->label;
    my $data = $self->analysis;
    (my $lab = $label) =~ s/^\w/_/g;
    $lab .= ".".$self->analysis;
    my $t1 = $self->t1;
    my $t2 = $self->t2;


    open my $rcfh, ">", "$rfile";
    if (my $rlib = $self->r_libs) {
        print $rcfh ".libPaths(\"$rlib\")\n";
    }

    print $rcfh "setwd(\"$Rdir\")
results<-read.table(\"$filename\", header = TRUE, sep=\"\t\")
results\$Class<-cut(results\$Zscore, breaks =c(min(results\$Zscore), $t1, $t2, max(results\$Zscore)), labels=FALSE, include.lowest=TRUE) # 99.9 and 99% CIs 1, 2, 3
require(rCharts)
d1 <- dPlot(
  y = \"Zscore\",
  x = c(\"Cell\", \"Tissue\", \"SNPs\", \"Number\", \"Accession\"),
  groups = \"Class\",
  data = results,
  type = \"bubble\",
  width = 2000,
  height = 1500,
  bounds = list(x=90,y=50,height=600,width=1850),
  id = \"chart.$lab\"
)\n";
    if ($data =~ /erc/){
	print $rcfh "d1\$xAxis( type = \"addCategoryAxis\", grouporderRule = \"Tissue\", orderRule = \"Cell\")\n";
    } else {
	print $rcfh "d1\$xAxis( type = \"addCategoryAxis\", grouporderRule = \"Tissue\", orderRule = \"Number\")\n";
    }

    print $rcfh "d1\$yAxis( type = \"addMeasureAxis\" )
d1\$colorAxis(
   type = \"addColorAxis\",
   colorSeries = \"Class\",
   palette = c(\"lightblue\",\"pink\",\"red\") )
d1\$addParams(title=\"$label overlaps with $data DHS\")
d1\$save('$chart', cdn = F)\n";

    system "R --no-save --quiet --slave < $rfile";
}

sub table{
    # Make Datatables table
    my ($self, $filename) = @_;

    my $chart = "table.htm";
    my $Rdir = $self->output.'/'.$self->jobid;
    my $t1 = $self->t1;
    my $t2 = $self->t2;


    my $rfile = "$Rdir/table.R";
    open my $rcfh, ">", "$rfile";
    if (my $rlib = $self->r_libs) {
        print $rcfh ".libPaths(\"$rlib\")\n";
    }

    print $rcfh "setwd(\"$Rdir\")
    data<-read.table(\"$filename\", header = TRUE, sep=\"\t\")
    results<-data.frame(data\$Cell, data\$Tissue, data\$Accession, data\$Zscore, data\$SNPs)
    names(results)<-c(\"Cell\", \"Tissue\", \"Accession\", \"Zscore\", \"SNPs\")
    require(rCharts)
    dt <- dTable(
      results,
      sScrollY= \"600\",
      bPaginate= F,
      sScrollX= \"100%\",
      sScrollXInner= \"110%\"
    )
    dt\$save('$chart', cdn = F)";
    system "R --no-save --quiet --slave < $rfile";
}

sub debug {
    my ($self,$msg) = @_;
    warn $msg;
    if ($self->logfile) {
	if (open F, ">>".$self->logfile) {
	    print F $msg;
	    close F;
	}
    }
}

sub status {
    my ($self,$msg) = @_;
    warn $msg;
    my $sh = $self->sh;
    print $sh $msg;
}

sub parse_input {
    my ($self, $file, $format) = @_;

    my @snps = ();
    die "Error: please specify file format\n";
    open my $fh, "<", $file or die "cannot open file $file : $!";
    if ($format =~ /rsid/){
        while (<$fh>){
            chomp;
            my @rsid = split /\:/, $_;
            my $rsid = pop @rsid; # take the last one for want of a better idea.
            push @snps, $rsid;
        }
    }
    elsif ($format =~ /ian/){
        while (<$fh>){
            my ($chr, $beg, $end, $rsid, undef) = split "\t", $_;
            my @rsid = split /\:/, $rsid; # to deal with multiple rsids
            $rsid = pop @rsid; # take the last one for want of a better idea. Can't take all as they are the same thing.
            push @snps, $rsid;
        }
    }
    elsif ($format =~ /vcf/){
        while (<$fh>){
            next if /^#/;
            my ($chr, $beg, $rsid) = split "\t", $_;
            unless ($chr =~ /^chr/){
                $chr = "chr". $chr;
            }
            if ($rsid =~/^rs\d+/){
                push @snps, $rsid;
            }
            else {
                my $loc = "$chr:$beg-$beg";
                #get the rsid from the db
                $rsid = $self->fetch_rsid($loc);
                push @snps, $rsid if defined $rsid;
            }
        }
    }
    elsif ($format =~ /bed|tabix/){
        while (<$fh>){
            my $loc;
            if ($format =~/bed/){
                next if /^track/;
                my ($chr, $beg, $end) = split "\t", $_;
                unless ($chr =~ /^chr/){
                    $chr = "chr". $chr;
                }
                $loc = "$chr:$end-$end";
            }
            elsif ($format =~ /tabix/){
                chomp;
                $loc = $_;
            }
            #get the $rsid from the db
            my $rsid = $self->fetch_rsid($loc);
            push @snps, $rsid if defined $rsid;
        }
    }

    return \@snps;
}

1;

