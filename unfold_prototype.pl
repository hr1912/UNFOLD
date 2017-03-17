#! /usr/bin/env perl
use warnings; 
use strict;
use Getopt::Long;
use File::Basename qw(fileparse);
use threads;
use Data::Dumper;

=head
---------------------------------------------------------------------

Author: Ruan Hang <ruanhang@novogene.cn>
---------------------------------------------------------------------
=cut

#----------------------------
# initalizing 
#----------------------------

if (@ARGV <2) {
	die "Usage: $0 -exp use_exprience_value(default yes) -soap soap_exe -R R_exe -reads reads_lst -ref reference_index -threads mapping_threads -para max_parallelled_running_samples -out output_dir\n";
}

my ($exp_value, $soap_exe, $R_exe, $reads_lst, $ref_idx, $threads, $para, $out) = ('','','','','', 0, 0,'');

GetOptions ('exp=s' => \$exp_value,
			'soap=s' => \$soap_exe,
			'R=s' => \$R_exe,
			'reads=s' => \$reads_lst,
			'ref=s' => \$ref_idx,
			'threads=i' => \$threads,
			'para=i' => \$para,
			'out=s' => \$out);

if ($exp_value eq 'no') {
	$exp_value = 0;
} else {
	$exp_value = 1;
}

$threads = 4 if (!$threads);
$para = 4 if (!$para);
$out = 'pipe_out' if (!$out);

my $node_processor_num = `grep 'processor' /proc/cpuinfo | sort -u | wc -l`;
$node_processor_num =~ s/\s+//g;
print "Running node has $node_processor_num processors\n";

$soap_exe = "/WORK/SD/ruanhang/bin/soap" if (!$soap_exe);
$R_exe = "/PROJ/SD/share/.0/bin/R" if (!$R_exe);

`mkdir $out`;`mkdir $out/{soap,images,stats,json}`;

#----------------------------
# reads mapping
#----------------------------

open READS, "$reads_lst" || die "Can not open $reads_lst, exiting ...\n";

my (%sample, %sample_details);

my $sample_counter=0;

while (<READS>) {
	chomp;
	$sample_counter ++;
	&queue($_);
}

foreach my $thr (threads->list()) {
	my %out = $thr->join();
	if (%out) {
		$sample{$out{'readsName'}}=\%out; 				
	}
}

print "All reads mapped and info stored...\n";		

close READS;

sub queue {
	my $rd = $_[0];
	
	while (threads->list() >= $para) {
		foreach my $thr (threads->list(threads::joinable)) {
			my %out = $thr->join();
			if (%out) {
				$sample{$out{'readsName'}}=\%out;
			}
		} 
		sleep(5);
	}
	
	my ($thr) = threads->create(\&map,$rd);
} 

sub map {
	my $full_path = $_[0];
	my ($rd,$dir)=fileparse($full_path);
	$rd =~ s/\s+$//g;
	
	if (-e "$out/soap/$rd.se") {
		#print "Sample $sample_counter $rd mapped ... skipping mapping ...\n";
	} else {
		#print "Mapping Sample $sample_counter $rd ...\n";
		`$soap_exe -p $threads -r 0 -a $full_path -D $ref_idx -o $out/soap/$rd.se 2>$out/soap/$rd.soap.log`;
	}
		
	return &stats("$out/soap/$rd");
}

#----------------------------
# bins info loading
#----------------------------

my $lowess_fp = "/PROJ/SD/PrenatalDiagnosis/GC_bias_test/files/bin.lowess";
my $bin_fp = "/PROJ/SD/PrenatalDiagnosis/GC_bias_test/files/bin.stat";

my %bin;
=head
-----------------------------
bin hash structure:
Key: bin name "chr.num"
Value: 	[0] N%
		[1] GC%
-----------------------------
=cut

open BIN, "$bin_fp" || die "Can not open $bin_fp\n";

while (<BIN>) {

	my ($name, $pos, $n, $gc, $read_c) = (split)[0..4];
	next if ($n != 0.0 || !$read_c);
	$bin{$name} = [$gc,$read_c];

}

close BIN;

my $median = &getMedian(\%bin);

print "median:$median\n";

sub getMedian {
	my %hash = %{$_[0]};	
	my $size = scalar(keys %hash);
	my $counter = 0;
	
	foreach (sort {$hash{$b}->[1] <=> $hash{$a}->[1]} keys %hash) {
		return $hash{$_}->[1] if ($counter >= ($size + 1)/2);
		$counter ++;
	}
}

my %corr_factor;

open LOWESS, "$lowess_fp" || die "Can not open $lowess_fp\n";

while (<LOWESS>) {
	my @tab = split;
	$tab[1] = sprintf("%.1f",$tab[1]);
	$corr_factor{$tab[1]} = $median / 50 / $tab[2] ;
}

close LOWESS;

#----------------------------
# GC correction
#----------------------------

foreach my $sample_id(sort keys %sample) {
	foreach my $chr(keys %{$sample{$sample_id}}) {
		next if ($chr =~/reads|Per/);
		foreach my $bin(keys %{$sample{$sample_id}{$chr}}) {
			my $num;
			if ($bin =~ /orig_(\d+)/) {
				$num = $1;
			}
			my $name = "$chr.$num";
			my $gc;
			if (exists $bin{$name}) {
				$gc = $bin{$name}[0];
			} else {
				#print STDERR "$name does not exists...\n";
				next;
			}
			#print "$name\t$gc\n" if (!exists $corr_factor{$gc}); 	
			$sample{$sample_id}{$chr}{'readsPerOrig'} += $sample{$sample_id}{$chr}{$bin};
			$sample{$sample_id}{$chr}{"corr_$num"} = $sample{$sample_id}{$chr}{$bin} * $corr_factor{$gc};
			$sample{$sample_id}{$chr}{'readsPerCorr'} += $sample{$sample_id}{$chr}{"corr_$num"};				
		}
		$sample{$sample_id}{$chr}{'readsPerOrig'} /= $sample{$sample_id}{'readsNum'};
		$sample{$sample_id}{$chr}{'readsPerCorr'} /= $sample{$sample_id}{'readsNum'};
	}
}

#print Dumper(\%sample);

#----------------------------
# calculate z score 
#----------------------------

my (%sum, %mean, %sd, %z_score);
my (%sum_ir, %mean_ir, %sd_ir, %z_score_ir);
my (%sum_gc, %mean_gc, %sd_gc, %z_score_gc);
my %ir = ("chr21" => "chr14", "chr13" => "chr4", "chr18" => "chr8", "chrX" => "chr7", "chrY" => "chr7");

my (%exp, %exp_ir, %exp_gc);

$exp{'chr21'}{'mean'} = 0.01298466; $exp{'chr21'}{'sd'} = 0.00011756;
$exp{'chr13'}{'mean'} = 0.03451643; $exp{'chr13'}{'sd'} = 0.00048865;
$exp{'chr18'}{'mean'} = 0.02879179; $exp{'chr18'}{'sd'} = 0.00021333;
$exp{'chrX'}{'mean'} = 0.03987724; $exp{'chrX'}{'sd'} = 0.00130241;
$exp{'chrY'}{'mean'} = 0.00019243; $exp{'chrY'}{'sd'} = 0.00011629;

$exp_ir{'chr21'}{'mean'} = 0.40753555; $exp_ir{'chr21'}{'sd'} = 0.00408000;
$exp_ir{'chr13'}{'mean'} = 0.54804976; $exp_ir{'chr13'}{'sd'} = 0.00319680;
$exp_ir{'chr18'}{'mean'} = 0.56654396; $exp_ir{'chr18'}{'sd'} = 0.00350219;
$exp_ir{'chrX'}{'mean'} = 0.73482015; $exp_ir{'chrX'}{'sd'} = 0.02539289;
$exp_ir{'chrY'}{'mean'} = 0.00354226; $exp_ir{'chrY'}{'sd'} = 0.00213251;

$exp_gc{'chr21'}{'mean'} = 0.01258996; $exp_gc{'chr21'}{'sd'} = 0.00013175;
$exp_gc{'chr13'}{'mean'} = 0.03527885; $exp_gc{'chr13'}{'sd'} = 0.00054287;
$exp_gc{'chr18'}{'mean'} = 0.02882543; $exp_gc{'chr18'}{'sd'} = 0.00024029;
$exp_gc{'chrX'}{'mean'} = 0.04018914; $exp_gc{'chrX'}{'sd'} = 0.00132073;
$exp_gc{'chrY'}{'mean'} = 0.00019285; $exp_gc{'chrY'}{'sd'} = 0.00011806;


###########################
my $chrY_cut = 0.00015;
my $zscore_cut = 3.00;
###########################

foreach my $sample_id(sort keys %sample) {
	foreach my $chr(keys %{$sample{$sample_id}}) {
		next if ($chr =~ /reads|Per/);
		$sum{$chr} += $sample{$sample_id}{$chr}{'readsPerOrig'};
		$sd{$chr} += $sample{$sample_id}{$chr}{'readsPerOrig'} * $sample{$sample_id}{$chr}{'readsPerOrig'};
		$sum_gc{$chr} += $sample{$sample_id}{$chr}{'readsPerCorr'};
		$sd_gc{$chr} += $sample{$sample_id}{$chr}{'readsPerCorr'} * $sample{$sample_id}{$chr}{'readsPerCorr'};
	}
	
	foreach my $chr(qw !chr21 chr13 chr18 chrX chrY!) { 
		$sum_ir{$chr} += $sample{$sample_id}{$chr}{'readsPerOrig'} / $sample{$sample_id}{$ir{$chr}}{'readsPerOrig'};
		$sd_ir{$chr} += $sample{$sample_id}{$chr}{'readsPerOrig'} * $sample{$sample_id}{$chr}{'readsPerOrig'} / $sample{$sample_id}{$ir{$chr}}{'readsPerOrig'} / $sample{$sample_id}{$ir{$chr}}{'readsPerOrig'};
	}
}

foreach my $chr(keys %sum) {
	$mean{$chr} = $sum{$chr} / $sample_counter;
	$sd{$chr} = $sd{$chr} / $sample_counter - $mean{$chr} * $mean{$chr};
	$sd{$chr} = sqrt($sd{$chr});
}
 
foreach my $chr(keys %sum_gc) {
	$mean_gc{$chr} = $sum_gc{$chr} / $sample_counter;
	$sd_gc{$chr} = $sd_gc{$chr} / $sample_counter - $mean_gc{$chr} * $mean_gc{$chr};
	$sd_gc{$chr} = sqrt($sd_gc{$chr});
} 

foreach my $chr(qw !chr21 chr13 chr18 chrX chrY!) {
	$mean_ir{$chr} = $sum_ir{$chr} / $sample_counter;
	$sd_ir{$chr} = $sd_ir{$chr} / $sample_counter - $mean_ir{$chr} * $mean_ir{$chr};
	$sd_ir{$chr} = 0 if ($sd_ir{$chr} < 0);
	$sd_ir{$chr} = sqrt($sd_ir{$chr});
}


if ($exp_value) {
	
foreach my $chr(qw !chr21 chr13 chr18 chrX chrY!) {
	$mean_ir{$chr} = $exp_ir{$chr}{'mean'};
	$sd_ir{$chr} = $exp_ir{$chr}{'sd'};
	$mean{$chr} = $exp{$chr}{'mean'};
	$sd{$chr} = $exp{$chr}{'sd'};
	$mean_gc{$chr}	= $exp_gc{$chr}{'mean'};
	$sd_gc{$chr} = $exp_gc{$chr}{'sd'};
}
	  
}
	
foreach my $sample_id(sort keys %sample) {
	foreach my $chr(keys %{$sample{$sample_id}}) {
		next if ($chr =~ /reads|Per/);
		if ($sd{$chr}) {
			$z_score{$sample_id}{$chr} = ($sample{$sample_id}{$chr}{'readsPerOrig'} - $mean{$chr}) / $sd{$chr};
		} else {
			$z_score{$sample_id}{$chr} = 0;
		}
		if ($sd_gc{$chr}) {
			$z_score_gc{$sample_id}{$chr} = ($sample{$sample_id}{$chr}{'readsPerCorr'} - $mean_gc{$chr}) / $sd_gc{$chr};
		} else {
			$z_score_gc{$sample_id}{$chr} = 0; 	
		}
	}

	foreach my $chr(qw !chr21 chr13 chr18 chrX chrY!) {
		$z_score_ir{$sample_id}{$chr} = ($sample{$sample_id}{$chr}{'readsPerOrig'}/$sample{$sample_id}{$ir{$chr}}{'readsPerOrig'} 
			- $mean_ir{$chr} ) / $sd_ir{$chr};
	}

	#print "$sample{$sample_id}{'chrY'}\n";

	if ($sample{$sample_id}{'chrY'}{'readsPerCorr'}  < $chrY_cut && $sample{$sample_id}{'chrY'}{'readsPerOrig'} < $chrY_cut) {
		$sample_details{$sample_id}{'sex'} = 'female';
	} else {
		$sample_details{$sample_id}{'sex'} = 'male';
	}

	foreach my $chr(qw !chr21 chr13 chr18!) {
		if ($z_score{$sample_id}{$chr} > $zscore_cut && $z_score_ir{$sample_id}{$chr} > $zscore_cut && $z_score_gc{$sample_id}{$chr} > $zscore_cut) {
			$sample_details{$sample_id}{$chr} = 'high';
		} else {
			$sample_details{$sample_id}{$chr} = 'low';
		}
	}
} 

#####################################################

$sample_counter=0;
%sum=%sum_gc=%sum_ir=();
%sd=%sd_gc=%sd_ir=();

foreach my $sample_id (sort keys %sample) {
	
	next if ($sample_details{$sample_id}{'chr21'} eq 'high' 
		|| $sample_details{$sample_id}{'chr13'} eq 'high'
		|| $sample_details{$sample_id}{'chr18'} eq 'high');

	$sample_counter ++;

	foreach my $chr(keys %{$sample{$sample_id}}) {
		next if ($chr =~ /reads|Per/);
		$sum{$chr} += $sample{$sample_id}{$chr}{'readsPerOrig'};
		$sd{$chr} += $sample{$sample_id}{$chr}{'readsPerOrig'} * $sample{$sample_id}{$chr}{'readsPerOrig'};
		$sum_gc{$chr} += $sample{$sample_id}{$chr}{'readsPerCorr'};
		$sd_gc{$chr} += $sample{$sample_id}{$chr}{'readsPerCorr'} * $sample{$sample_id}{$chr}{'readsPerCorr'};
	}
	
	foreach my $chr(qw !chr21 chr13 chr18 chrX chrY!) { 
		$sum_ir{$chr} += $sample{$sample_id}{$chr}{'readsPerOrig'} / $sample{$sample_id}{$ir{$chr}}{'readsPerOrig'};
		$sd_ir{$chr} += $sample{$sample_id}{$chr}{'readsPerOrig'} * $sample{$sample_id}{$chr}{'readsPerOrig'} / $sample{$sample_id}{$ir{$chr}}{'readsPerOrig'} / $sample{$sample_id}{$ir{$chr}}{'readsPerOrig'};
	}
}

print "Stats based on $sample_counter samples\n";

foreach my $chr(keys %sum) {
	$mean{$chr} = $sum{$chr} / $sample_counter;
	$sd{$chr} = $sd{$chr} / $sample_counter - $mean{$chr} * $mean{$chr};
	$sd{$chr} = 0 if ($sd{$chr} < 0);
	$sd{$chr} = sqrt($sd{$chr});
}
 
foreach my $chr(keys %sum_gc) {
	$mean_gc{$chr} = $sum_gc{$chr} / $sample_counter;
	$sd_gc{$chr} = $sd_gc{$chr} / $sample_counter - $mean_gc{$chr} * $mean_gc{$chr};
	$sd_gc{$chr} = 0 if ($sd_gc{$chr} < 0);
	$sd_gc{$chr} = sqrt($sd_gc{$chr});
} 

foreach my $chr(qw !chr21 chr13 chr18 chrX chrY!) {
	$mean_ir{$chr} = $sum_ir{$chr} / $sample_counter;
	$sd_ir{$chr} = $sd_ir{$chr} / $sample_counter - $mean_ir{$chr} * $mean_ir{$chr};
	$sd_ir{$chr} = 0 if ($sd_ir{$chr} < 0);
	$sd_ir{$chr} = sqrt($sd_ir{$chr});
}

######################################################


print "-------------------------------------------\n";

foreach my $chr(qw !chr21 chr13 chr18 chrX chrY!) {
	printf "%s\tmean\tsd\n", $chr;
	printf "%s\t%.8f\t%.8f\n",'orig',$mean{$chr},$sd{$chr};
	printf "%s\t%.8f\t%.8f\n",'ir',$mean_ir{$chr},$sd_ir{$chr};
	printf "%s\t%.8f\t%.8f\n\n",'gc',$mean_gc{$chr},$sd_gc{$chr};
}

print "--------------------------------------------\n";

#----------------------------
# output csv
#----------------------------

open CSV, ">$out/stats/stats.csv" || die "Can not creat $out/stats/stats.csv, exiting ...\n";

print CSV "Sample,chr21.z.score,chr21.z.score.ir,chr21.z.score.gc,chr21.risk,chr13.z.score,chr13.z.score.ir,chr13.z.score.gc,chr13.risk,chr18.z.score,chr18.z.score.ir,chr18.z.score.gc,chr18.risk,chrX,chrX.gc,chrX.z.score,chrX.z.score.ir,chrX.zcore.gc,chrY,chrY.gc,chrY.z.score,chrY.z.score.ir,chrY.z.score.gz,sex\n";

foreach my $sample_id(sort keys %sample) {
	printf CSV "%s,%.3f,", $sample_id,$z_score{$sample_id}{'chr21'};
	printf CSV "%.3f,%.3f,%s,", $z_score_ir{$sample_id}{'chr21'},$z_score_gc{$sample_id}{'chr21'},$sample_details{$sample_id}{'chr21'};
	printf CSV "%.3f,", $z_score{$sample_id}{'chr13'};
	printf CSV "%.3f,%.3f,%s,", $z_score_ir{$sample_id}{'chr13'},$z_score_gc{$sample_id}{'chr13'},$sample_details{$sample_id}{'chr13'};
	printf CSV "%.3f,", $z_score{$sample_id}{'chr18'};
	printf CSV "%.3f,%.3f,%s,", $z_score_ir{$sample_id}{'chr18'},$z_score_gc{$sample_id}{'chr18'},$sample_details{$sample_id}{'chr18'};
	printf CSV "%.3f,%.3f,%.3f,", $sample{$sample_id}{'chrX'}{'readsPerOrig'}*100 ,$sample{$sample_id}{'chrX'}{'readsPerCorr'}*100, $z_score{$sample_id}{'chrX'};
	printf CSV "%.3f,%.3f,", $z_score_ir{$sample_id}{'chrX'},$z_score_gc{$sample_id}{'chrX'};
	printf CSV "%.3f,%.3f,%.3f,", $sample{$sample_id}{'chrY'}{'readsPerOrig'}*100, $sample{$sample_id}{'chrY'}{'readsPerCorr'}*100, $z_score{$sample_id}{'chrY'};
	printf CSV "%.3f,%.3f,%s\n", $z_score_ir{$sample_id}{'chrY'},$z_score_gc{$sample_id}{'chrY'},$sample_details{$sample_id}{'sex'};
}

close CSV;

open CSV, ">$out/results.csv" || die "Can not creat $out/results.csv, exiting ...\n";

print CSV "Sample,mapping.rate,identical.rate,identical.num,chr21.risk,chr21.z.score,chr21.z.score.ir,chr21.z.score.gc,chr13.risk,chr13.z.score,chr13.z.score.ir,chr13.z.score.gc,chr18.risk,chr18.z.score,chr18.z.score.ir,chr18.z.score.gc,chrY.percent,chrY.percent.gc,sex\n";

foreach my $sample_id(sort keys %sample) {
	printf CSV "%s,%.2f,%.2f,%d,", $sample_id, $sample{$sample_id}{'mapPer'}, $sample{$sample_id}{'identicalPer'}, $sample{$sample_id}{'readsNum'}; 
	printf CSV "%s,%.3f,%.3f,%.3f,", $sample_details{$sample_id}{'chr21'}, $z_score{$sample_id}{'chr21'}, $z_score_ir{$sample_id}{'chr21'}, $z_score_gc{$sample_id}{'chr21'};
	printf CSV "%s,%.3f,%.3f,%.3f,", $sample_details{$sample_id}{'chr13'}, $z_score{$sample_id}{'chr13'}, $z_score_ir{$sample_id}{'chr13'}, $z_score_gc{$sample_id}{'chr13'};
	printf CSV "%s,%.3f,%.3f,%.3f,", $sample_details{$sample_id}{'chr18'}, $z_score{$sample_id}{'chr18'}, $z_score_ir{$sample_id}{'chr18'}, $z_score_gc{$sample_id}{'chr18'};
	printf CSV "%.3f,%.3f,%s\n", $sample{$sample_id}{'chrY'}{'readsPerOrig'}*100, $sample{$sample_id}{'chrY'}{'readsPerCorr'}*100, $sample_details{$sample_id}{'sex'};
}

close CSV;
#----------------------------
# plotting
#----------------------------

foreach my $chr(qw !chr21 chr13	chr18!) {

my $plot_module=<<"END";
library(ggplot2)
png(type="cairo", "$out/images/$chr.zscore.png")
a<-read.csv("$out/stats/stats.csv")
p<-ggplot(a,aes(x="Original",y=$chr.z.score)) + geom_point(aes(color=$chr.risk),shape=2,size=3) + geom_hline(color="darkblue",linetype=3,yintercept=3) + labs(list(title="$chr z-score distribution",x="Samples", y="Z-score"));
p
dev.off()
END
 
open R, "|$R_exe --vanilla --slave" or die $!;
print R $plot_module;
close R;

$plot_module=<<"END";
library(ggplot2)
png(type="cairo", "$out/images/$chr.zscore_ir.png")
a<-read.csv("$out/stats/stats.csv")
p<-ggplot(a,aes(x="Internal Reference correted",y=$chr.z.score.ir)) + geom_point(aes(color=$chr.risk),shape=2,size=3) + geom_hline(color="darkblue",linetype=3,yintercept=3) + labs(list(title="$chr z-score distribution",x="Samples", y="Z-score_ir"));
p
dev.off()
END
 
open R, "|$R_exe --vanilla --slave" or die $!;
print R $plot_module;
close R;

$plot_module=<<"END";
library(ggplot2)
png(type="cairo", "$out/images/$chr.zscore_gc.png")
a<-read.csv("$out/stats/stats.csv")
p<-ggplot(a,aes(x="GC corrected",y=$chr.z.score.gc)) + geom_point(aes(color=$chr.risk),shape=2,size=3) + geom_hline(color="darkblue",linetype=3,yintercept=3) + labs(list(title="$chr z-score distribution",x="Samples", y="Z-score_gc"));
p
dev.off()
END

open R, "|$R_exe --vanilla --slave" or die $!;
print R $plot_module;
close R;

$plot_module=<<"END";
library(ggplot2)
png(type="cairo", "$out/images/$chr.zscore_all.png")
a<-read.csv("$out/stats/stats.csv")
p<-ggplot(a,aes(x="Original",y=$chr.z.score)) + geom_point(aes(color=$chr.risk),shape=2,size=3) + geom_point(aes(x="Internal Reference corrected",y=$chr.z.score.ir,color=$chr.risk),shape=2,size=3) + geom_point(aes(x="GC corrected",y=$chr.z.score.gc,color=$chr.risk),shape=2,size=3)+ geom_hline(color="darkblue",linetype=3,yintercept=3) + labs(list(title="$chr z-score distribution",x="Samples", y="Z-score"));
p
dev.off()
END

open R, "|$R_exe --vanilla --slave" or die $!;
print R $plot_module;
close R;
}

my $plot_module=<<"END";
library(ggplot2)
png(type="cairo", "$out/images/chrY.per.png")
a<-read.csv("$out/stats/stats.csv")
p<-ggplot(a,aes(x="Original",y=chrY))  + geom_point(aes(color=sex),shape=2,size=3) + geom_hline(color="darkblue",linetype=3,yintercept=0.015)+ labs(list(title="chrY mapping percent distribution",x="Samples", y="chrY %"));
p
dev.off()
END

open R, "|$R_exe --vanilla --slave" or die $!;
print R $plot_module;
close R;

$plot_module=<<"END";
library(ggplot2)
png(type="cairo", "$out/images/chrY.per_gc.png")
a<-read.csv("$out/stats/stats.csv")
p<-ggplot(a,aes(x="GC corrected",y=chrY.gc))  + geom_point(aes(color=sex),shape=2,size=3) + geom_hline(color="darkblue",linetype=3,yintercept=0.015)+ labs(list(title="chrY mapping percent distribution",x="Samples", y="chrY % GC Corrected"));
p
dev.off()
END

open R, "|$R_exe --vanilla --slave" or die $!;
print R $plot_module;
close R;

$plot_module=<<"END";
library(ggplot2)
png(type="cairo", "$out/images/chrY.per_all.png")
a<-read.csv("$out/stats/stats.csv")
p<-ggplot(a,aes(x="Original",y=chrY)) + geom_point(aes(color=sex),shape=2,size=3) + geom_point(aes(x="GC corrected",y=chrY.gc,color=sex),shape=2,size=3) + geom_hline(color="darkblue",linetype=3,yintercept=0.015)+ labs(list(title="chrY mapping percent distribution",x="Samples", y="chrY %"));
p
dev.off()
END

open R, "|$R_exe --vanilla --slave" or die $!;
print R $plot_module;
close R;

#----------------------------
# output json
#----------------------------

open JSON, ">$out/json/results.json" or die "Can not create $out/json/results.json, exiting ...\n";

print JSON "{\"$out/json/results.json\":[\n";

my $temp_c = 0;

foreach my $sample_id(sort keys %sample) {
my $each=<<"END";
{"Sample ID":"$sample_id",
"Sex":"$sample_details{$sample_id}{'sex'}",
"Chr21 Trisomy Risk":"$sample_details{$sample_id}{'chr21'}",
"Chr13 Trisomy Risk":"$sample_details{$sample_id}{'chr13'}",
"Chr18 Trisomy Risk":"$sample_details{$sample_id}{'chr18'}"
END

print JSON $each;
$temp_c ++;
if ($temp_c == $sample_counter) {
	print JSON "}\n";
} else {
	print JSON "},\n";
}

}

print JSON "]\n}\n";

close JSON;

#----------------------------
# partial subroutines
#----------------------------

sub stats {
	my $prefix = $_[0];
	my ($mapping_log, $mapping_out) = ("$prefix.soap.log", "$prefix.se");
	my ($rd_c, $mapped_c, $valid_c) = (0,0,0);
	my $mapping_time = 0;

	my %stat;

	$stat{'readsName'} = $prefix;
	if (-e $mapping_log) { 

	$rd_c = `grep "Total Reads" $mapping_log`;
	$rd_c =~ s/\D+//g;

	$mapped_c = `grep ^Alignment $mapping_log`;
	$mapped_c =~ s/^\D+//g;
	$mapped_c =~ s/\(\S+\)//g;
	
	$mapping_time = `grep "Elapsed Time" $mapping_log`;
	$mapping_time =~ s/^\D+//g;
	$mapping_time =~ s/\s+$//g;
	
	} else {
		print STDERR "$mapping_log does not exist ...\n";
	}

	open SE, "$mapping_out" || die "Can not open $mapping_out, exiting ...\n";

	while(<SE>) {
		my ($ref, $pos,$mismatch) = (split)[7..9];
		next if ($ref eq 'chrM');
		next if ($mismatch);
		$valid_c ++;
	
		my $bin_c = int($pos/50000) + 1;
		my $name= "$ref.$bin_c";
	
		$stat{'readsNum'} ++;	
		$stat{$ref}{"orig_$bin_c"} ++;
		
	}
	
	if ($rd_c) {
		printf "%s #reads %d #mapped %d %.2f%% #valid %d %.2f%% mapping time %.2f\n", 
			$prefix, $rd_c, $mapped_c, $mapped_c * 100 / $rd_c, $valid_c, $valid_c * 100 / $rd_c, $mapping_time;
	}

	$stat{'mapPer'} = $mapped_c*100/$rd_c;
	$stat{'identicalPer'} = $valid_c*100/$rd_c;

	#print Dumper(\%stat);

	return %stat;
}
	
