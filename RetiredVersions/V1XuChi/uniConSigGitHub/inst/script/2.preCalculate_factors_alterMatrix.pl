#!/usr/bin/perl

my $jaccardcutoff=1;
#usage: perl $0 trainingCellLines conceptFile
use strict;
use warnings;
#my ($cycleN)=@ARGV;
#my $cycleN=1;
#for(my $cycleN=1;$cycleN<=20;$cycleN++){
my $cl_info="Homo_sapiens.gene_info";
my $conAddr="Conceptgene_20150526.db";
#my $file=$ARGV[0];
#my $inputAddr="./training_cellLine/$ARGV[0]";

open (CL,"$cl_info") or die "OpenError: $cl_info, $!\n";
my %idAll; my %id2name;
while(<CL>){
	next if(/^\#/);
	$_=~s/[\n\r]//g;
	my @data=split(/\t/,$_);
	$idAll{$data[1]}=1;
	$id2name{$data[1]}=$data[2];
}
close CL;
#open (INPUT,"$inputAddr") or die "OpenError: $inputAddr, $!\n";
#my %input;my %geneCount;
#
#while(<INPUT>){
#	chomp;
#	$_=~s/[\n\r]//g;
#	my @data=split(/\t/);
#	$input{"$data[0]\t$data[1]"}{$data[2]}=1;
#}
#
#close INPUT;

my $date=&getTime();
print "$date->{hour}\:$date->{minute} $date->{second}s $date->{date}\t";
print "Loading input list finished\n";

open (DB,"$conAddr") or die "OpenError: Conceptgene_test.db $!\n";
my %DBconcept_id;my %DBid_concept; my $n=0;
while(<DB>){
	if($n==0){
		$n++;
		next;
	}
	$_=~s/[\n\r]//g;
	my @data=split(/\t/);
	if((exists $data[2]) and ($data[2]=~/\d+/)){
		$DBconcept_id{"$data[0]\^$data[1]"}{$data[2]}=1;
		$DBid_concept{$data[2]}{"$data[0]\^$data[1]"}=1;
	}
}
close DB;
$date=&getTime();
print "$date->{hour}\:$date->{minute} $date->{second}s $date->{date}\t";
print "Loading concept database finished\n";

#print SIG "InputCat\tInput\tInputSize(n=)\tConceptCat\tConcept\tConceptSize(n=)\tOverlap(n=)\tTotal(n=)\tJaccardIndex\n";
my %excluded;
foreach my $concept (keys %DBconcept_id){
	if (keys %{$DBconcept_id{$concept}} < 5){
		delete $DBconcept_id{$concept};
		$excluded{$concept}=1;
	}
}
$date=&getTime();
print "$date->{hour}\:$date->{minute} $date->{second}s $date->{date}\t";

print "Begin Loading Jaccard score matrix\n";
my %conceptscore; my %score; my %scoreanno;
my %jaccardMatrix;
my $conceptFile=$conAddr;
$conceptFile=~s/\.db$//;
open (JM,"./conceptMatrix_rm1overlap_chi.xls") or die "OpenError: conceptMatrix_rm1overlap_chi.xls, $!\n";
my $line=0; my @dbTitle;
while(<JM>){
	chomp;
	$_=~s/[\n\r]//g;
	if($line==0){
		@dbTitle=split(/\t/);
		shift @dbTitle; #the first element is ""(null)
		$line++;
		next;
	}
	my @data=split(/\t/);
	my $key=shift @data;
	for (my $i=0;$i<=$#dbTitle;$i++){
		if (not exists $jaccardMatrix{"$dbTitle[$i]\+$key"}){ #load half of the matrix to mem, to save mem usage
			$jaccardMatrix{"$key\+$dbTitle[$i]"}=$data[$i]; #load the matrix to %jaccardMatrix
		}
	}
}
close JM;

#foreach my $concept (keys %DBconcept_id){
#	foreach my $genelist (keys %input){
#		my %common=();
#		my %total=();
#		foreach my $geneID (keys %{$DBconcept_id{$concept}}){
#			#print "$geneID\t$genelist\t"; print sort {$b<=>$a} keys %{$input{$genelist}}; print "\n";
#			$common{$geneID}=1 if exists $input{$genelist}{$geneID};
#			$total{$geneID}=1;
#		}
#		foreach my $geneID (keys %{$input{$genelist}}){
#			$total{$geneID}=1;
#		}
#		my $c=keys %common;
#		my $t=keys %total;
#		$jaccardMatrix{$concept}{"genelist"}{"weight"}{$genelist}=$c/$t;#print "$concept\t$genelist\t$c/$t\t$jaccardMatrix{$concept}{genelist}{weight}{$genelist}\n";
#		$jaccardMatrix{$concept}{"genelist"}{"weightdel"}{$genelist}=($c-1)/$t;
#	}
#}
$date=&getTime();
print "$date->{hour}\:$date->{minute} $date->{second}s $date->{date}\t";
print "Loading Jaccard score matrix finished\nBegin calculation of factors\n";

#for my $cutoffTmp (5..5){
	#my $cutoff=0.1;
	#for my $cycleD (1..1){
		#my $cycleD=1;
		#my $dupGeneList="/ihome/xiaosongwang/xuc13/randDup/cutoff$cutoff/$cutoff\_toBeDup$cycleD\_Conceptgene_20150526.db_cleaned";
		#system("mkdir /ihome/xiaosongwang/xuc13/randDup/uniConsigV3/cutoff$cutoff/gl$cycleD")
		#system("mkdir /ihome/xiaosongwang/xuc13/randDup/uniConsigV4/cutoff$cutoff/gl$cycleD");
		#print "Calculation is using $inputAddr as input geneList\nDuplicating $dupGeneList database\n";
		#open (DUP,"$dupGeneList") or die "OpenError: $dupGeneList, $!\n";
		#my %dup;
		#print "Dealing with Dup$cycleD\n";
		#while(<DUP>){
		#	chomp;
		#	$_=~s/[\n\r]//g;
		#	my @data=split(/\t/,$_);
		#	$dup{"$data[0]\^$data[1]"}=1;
		#}
		#calculate the consig score for each geneID
			my %copyDBid_concept=%DBid_concept;
		#foreach my $genelist (sort keys %input){
			open(OUT,">./preCal_factors_rm1overlap_$conceptFile\.xls") or die "OpenError: ./preCal_factors_$conceptFile\.xls, $!\n";
			print OUT "ID\tECN\tfactors\n";
			foreach my $geneID (sort keys %idAll){
				#my $sumV7=0;
				my $penV7=0;
				#my $m=0;
				my %factors;
				foreach my $concept (keys %{$DBid_concept{$geneID}}){
					#print "$geneID\t$concept\n";
					#my @tmp=keys $jaccardMatrix{$concept}{"genelist"}{"weight"}; print "@tmp\n";
					#next if ((not exists $jaccardMatrix{$concept}{"genelist"}{"weight"}{$genelist}) or ($jaccardMatrix{$concept}{"genelist"}{"weight"}{$genelist}>$jaccardcutoff));
					#$m++;
					if (exists $excluded{$concept}){
						next;
					}
					my $consumV7=0;
					#if($geneID==1240121){
					#	print "$geneID\t$concept\n";
					#}
					foreach my $copyConcept (keys %{$copyDBid_concept{$geneID}}){
						if(exists $jaccardMatrix{"$concept\+$copyConcept"}){
							if($jaccardMatrix{"$concept\+$copyConcept"}>0.05){
								$consumV7=$consumV7+$jaccardMatrix{"$concept\+$copyConcept"};
							}
						}elsif(exists $jaccardMatrix{"$copyConcept\+$concept"}){
							if($jaccardMatrix{"$copyConcept\+$concept"}>0.05){
								$consumV7=$consumV7+$jaccardMatrix{"$copyConcept\+$concept"};
							}
						}
					}
					$factors{$concept}=$consumV7;
					#if (exists $input{$genelist}{$geneID}){
					#	$sumV7=$sumV7+$jaccardMatrix{$concept}{"genelist"}{"weightdel"}{$genelist}/$consumV7;
					#	my @concept=split(/\^/,$concept);
					#	$conceptscore{$genelist}{$geneID}{$concept[1]}=$jaccardMatrix{$concept}{"genelist"}{"weightdel"}{$genelist};
					#}else{
					#	$sumV7=$sumV7+$jaccardMatrix{$concept}{"genelist"}{"weight"}{$genelist}/$consumV7;
					#	my @concept=split(/\^/,$concept);
					#	$conceptscore{$genelist}{$geneID}{$concept[1]}=$jaccardMatrix{$concept}{"genelist"}{"weight"}{$genelist};
					#}
					if($consumV7==0){
						print "$geneID\t$concept\n";
					}else{
						$penV7=$penV7+1/$consumV7;
					}
				}
				print OUT "$geneID\t$penV7";
				foreach my $concept (sort keys %factors){
					print OUT "\t$concept\_$factors{$concept}";
				}
				print OUT "\n";
				#if ($m==0){
				#	$score{$genelist}{$geneID}{V7}=0;
				#	$scoreanno{$genelist}{$geneID}{V7}="0\t0\t";
				#}else{
				#	$score{$genelist}{$geneID}{V7}=$sumV7/sqrt($penV7);
				#	$scoreanno{$genelist}{$geneID}{V7}="$sumV7\t$penV7\t";
				#}
				#print "---------------------------------------------------\n";
				#print "$geneID\t$genelist\t$score{$genelist}{$geneID}\n";
			}
			close OUT;
			#my %maxconsig; 
			#foreach my $geneid (sort {$score{$genelist}{$b}{V7}<=>$score{$genelist}{$a}{V7}} keys %{$score{$genelist}}){
			#	$maxconsig{$genelist}=$score{$genelist}{$geneid}{V7};
			#	print "$maxconsig{$genelist}\n";
			#last;
			#}
			#$date=&getTime();
			#print "$date->{hour}\:$date->{minute} $date->{second}s $date->{date}\t";
			#print "$genelist consig score calculation finished\nBegin outputing results\n";
			#my @nameData=split(/\t/,$genelist);
			#open (V7,">./uniConSig/curated_noDup_$nameData[0]\_$nameData[1]\_$conceptFile\.xls") or die "OpenError: ./uniConSig/curated_noDup_$nameData[0]\_$nameData[1]\_$conceptFile\.xls, $!\n";
			#print V7 "genelist\t\tgeneSym\tGeneID\tsumECWV7\tECNV7\tUniConSigScoreV7\n";
			#foreach my $geneid (sort {$score{$genelist}{$b}{V7}<=>$score{$genelist}{$a}{V7}} keys %{$score{$genelist}}){
			#	$input{$genelist}{$geneid}="" if (not exists $input{$genelist}{$geneid});
			#	if($geneid == 1299078){
			#		print "$id2name{$geneid}\t$geneid\t$scoreanno{$genelist}{$geneid}{V7}$score{$genelist}{$geneid}{V7}\n";
			#	}
			#	print V7 "$genelist\t$id2name{$geneid}\t$geneid\t$scoreanno{$genelist}{$geneid}{V7}$score{$genelist}{$geneid}{V7}\n";
			#}
			$date=&getTime();
			print "$date->{hour}\:$date->{minute} $date->{second}s $date->{date}\t";
			print "Pre-calculation finished!\n";
			#close RE;close SIG;
			#close V7;
		#}
	#}
#}

sub getTime{
	my $time=shift||time();
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime($time);
	$mon++;
	$sec=($sec<10)?"0$sec":$sec;
	$min=($min<10)?"0$min":$min;
	$hour=($hour<10)?"0$hour":$hour;
	$mday=($mday<10)?"0$mday":$mday;
	$mon=($mon<10)?"0$mon":$mon;
	$year+=1900;
	my $weekday=('Sun','Mon','Tue','Wed','Thu','Fri','Sat')[$wday];
	return {
		'second'=>$sec,
		'minute'=>$min,
		'hour'=>$hour,
		'day'=>$mday,
		'month'=>$mon,
		'year'=>$year,
		'weekNo'=>$wday,
		'wday'=>$weekday,
		'yday'=>$yday,
		'date'=>"$year$mon$mday"
	};
}
