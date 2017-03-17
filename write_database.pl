#! /PROJ/SD/share/software/ActivePerl-5.16/bin/perl
use warnings;
use strict;
use DBI;
use Encode;

=head
---------------------------------------------------------------------

Author: Ruan Hang <ruanhang@novogene.cn>
---------------------------------------------------------------------
=cut

my $fp_csv = shift || die "$0 [pre_report.csv] [allchrom.csv]\n";
my $fp_all_chrom = shift || die "$0 [pre_report.csv] [allchrom.csv]\n";

open CSV, "$fp_csv" || die $!;
open ALL, "$fp_all_chrom" || die $!;

my %hash;
my $date = `date '+%Y-%m-%d %T'`;

while (<CSV>) {
	chomp;
	my @tab = split /\,/,$_;
	$hash{"$tab[0]"}{"chr13_zscore"} = $tab[1];
	$hash{"$tab[0]"}{"chr13_result"} = encode("gbk",decode("utf8",$tab[2] eq "+" ? '高危' : '未见明显异常'));
	if ($tab[1] > 3.0 || $tab[1] < -3.0) {
		$hash{"$tab[0]"}{"chr13_flag"} = "!";
	}
	$hash{"$tab[0]"}{"chr18_zscore"} = $tab[3];
	$hash{"$tab[0]"}{"chr18_result"} = encode("gbk",decode("utf8",$tab[4] eq "+" ? '高危' : '未见明显异常'));
	if ($tab[3] > 3.0 || $tab[3] < -3.0) {
		$hash{"$tab[0]"}{"chr18_flag"} = "!";
	}
	$hash{"$tab[0]"}{"chr21_zscore"} = $tab[5];
	$hash{"$tab[0]"}{"chr21_result"} = encode("gbk",decode("utf8",$tab[6] eq "+" ? '高危' : '未见明显异常'));
	if ($tab[5] > 3.0 || $tab[5] < -3.0) {
		$hash{"$tab[0]"}{"chr21_flag"} = "!";
	} 
}

close CSV;

print STDERR "pre_report file loaded\n";

my ($all_zscore_low,$all_zscore_high) = (-4,4);

<ALL>;
while (<ALL>) {
	chomp;
	my @tab = split /\,/,$_;
	my $libID = $tab[1];
	my @chr_zscore;
	my ($ab_chr, $ab_val) = ('','');
	for (my $i = 1; $i <= 24; $i ++) {
		next if ($i == 13 || $i == 18 || $i == 21);
		$chr_zscore[$i] = $tab[$i*3+1];
		if ($chr_zscore[$i] < $all_zscore_low || $chr_zscore[$i] > $all_zscore_high) {
			$i = 'X' if ($i==23);
			$i = 'Y' if ($i==24); 		
			$ab_chr.="chr$i,";
			$ab_val.="$chr_zscore[$i],";
		}
	}
	
	if ($ab_chr) {	
		$ab_chr =~ s/,$//;
		$ab_val = s/,$//;
		$hash{"$libID"}{"chrxx_name"} = $ab_chr;
		$hash{"$libID"}{"chrxx_zscore"} = $ab_val;
	}

}
	
my $Orcl_Sid="novoorcl";
my $Orcl_Host="192.168.7.200";
my $Orcl_User="novo_cq";
my $Orcl_Pass="novo_cq";

$ENV{ORACLE_SID}=$Orcl_Sid;
$ENV{NLS_LANG}='American_America.ZHS16GBK';

my $Orcl_DB = DBI->connect("DBI:Oracle:host=$Orcl_Host;sid=$Orcl_Sid", $Orcl_User, $Orcl_Pass)
	|| die ($DBI::errstr . "\n");

print STDERR "Oracle database connected\n";

my $obj_txt2_str = encode("gbk",decode("utf8",'合格'))." $date";

foreach my $sample_code(keys %hash) {
	my ($Update, $tmp);

	$Update = "UPDATE K_WORKFLOW_OBJ SET OBJ_TXT2 = \'$obj_txt2_str\' WHERE OBJ_TXT1 = \'$sample_code\' AND PPOINT_ID = '038654b3-f03b-4d77-9360-b84fa2376bb8'";
	$tmp = $Orcl_DB->prepare($Update);
	$tmp->execute();

	if (exists $hash{$sample_code}{chrxx_name}) {
		$Update = "UPDATE CQ_REPORT_DOC SET CHRXX_NAME = \'$hash{$sample_code}{chrxx_name}\' WHERE SAMPLE_CODE = \'$sample_code\'";
		$tmp = $Orcl_DB->prepare($Update);
		$tmp->execute();
	}

	if (exists $hash{$sample_code}{chrxx_zscore}) {
		$Update = "UPDATE CQ_REPORT_DOC SET CHRXX_ZVALUE = \'$hash{$sample_code}{chrxx_zscore}\' WHERE SAMPLE_CODE = \'$sample_code\'";
		$tmp = $Orcl_DB->prepare($Update);
		$tmp->execute();
	}

	$Update = "UPDATE CQ_REPORT_DOC SET CHR13_ZVALUE = \'$hash{$sample_code}{chr13_zscore}\' WHERE SAMPLE_CODE = \'$sample_code\'";
	$tmp = $Orcl_DB->prepare($Update);
	$tmp->execute();

	$Update = "UPDATE CQ_REPORT_DOC SET CHR13_RESULT = \'$hash{$sample_code}{chr13_result}\' WHERE SAMPLE_CODE = \'$sample_code\'";
	$tmp = $Orcl_DB->prepare($Update);
	$tmp->execute();

	if (exists $hash{$sample_code}{"chr13_flag"}) {
		$Update = "UPDATE CQ_REPORT_DOC SET CHR13_WARING = \'$hash{$sample_code}{chr13_flag}\' WHERE SAMPLE_CODE = \'$sample_code\'";
		$tmp = $Orcl_DB->prepare($Update);
		$tmp->execute();
	}

	$Update = "UPDATE CQ_REPORT_DOC SET CHR18_ZVALUE = \'$hash{$sample_code}{chr18_zscore}\' WHERE SAMPLE_CODE = \'$sample_code\'";
	$tmp = $Orcl_DB->prepare($Update);
	$tmp->execute();

	$Update = "UPDATE CQ_REPORT_DOC SET CHR18_RESULT = \'$hash{$sample_code}{chr18_result}\' WHERE SAMPLE_CODE = \'$sample_code\'";
	$tmp = $Orcl_DB->prepare($Update);
	$tmp->execute();

	if (exists $hash{$sample_code}{"chr18_flag"}) {
		$Update = "UPDATE CQ_REPORT_DOC SET CHR18_WARING = \'$hash{$sample_code}{chr18_flag}\' WHERE SAMPLE_CODE = \'$sample_code\'";
		$tmp = $Orcl_DB->prepare($Update);
		$tmp->execute();
	}
			
	$Update = "UPDATE CQ_REPORT_DOC SET CHR21_ZVALUE = \'$hash{$sample_code}{chr21_zscore}\' WHERE SAMPLE_CODE = \'$sample_code\'";
	$tmp = $Orcl_DB->prepare($Update);
	$tmp->execute();

	$Update = "UPDATE CQ_REPORT_DOC SET CHR21_RESULT = \'$hash{$sample_code}{chr21_result}\' WHERE SAMPLE_CODE = \'$sample_code\'";
	$tmp = $Orcl_DB->prepare($Update);
	$tmp->execute();

	if (exists $hash{$sample_code}{"chr21_flag"}) {
		$Update = "UPDATE CQ_REPORT_DOC SET CHR21_WARING = \'$hash{$sample_code}{chr21_flag}\' WHERE SAMPLE_CODE = \'$sample_code\'";
		$tmp = $Orcl_DB->prepare($Update);
		$tmp->execute();
	}

	print STDERR "sample $sample_code proceed\n";

}

print STDERR "All sample proceed at $date\n";

my $Table = "SELECT * FROM CQ_REPORT_DOC";
my $sth= $Orcl_DB->prepare($Table);
$sth->execute();

for (my $i=0;$i<$sth->{NUM_OF_FIELDS};$i++) {
	printf "%s\t", $sth->{NAME}[$i];
}

print "\n";

while (my @row = $sth->fetchrow_array() ) {
	foreach (@row) {
		$_ = "_" if !defined($_);
		print "%s\t", $_;
	}
	print "\n";
}

End {
	$Orcl_DB->disconnect if (defined ($Orcl_DB));
}
