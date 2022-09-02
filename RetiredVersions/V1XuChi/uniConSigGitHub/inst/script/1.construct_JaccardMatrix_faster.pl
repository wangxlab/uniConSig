#!/usr/bin/perl
#calculate the Jaccard score matrix, the matrix will be like:
#				concept1	concept2	concept3	...	
#		concept1	J11		J12			J13			J1i		
#		concept2	J21		J22			J23			J2i		
#		concept3	J31		J32			J33			J3i	
#		...
#my $conAddr=.\Conceptgene_20150525.db;

open (DB,"./Conceptgene_20150526.db") or die "OpenError: Conceptgene_20150526.db $!\n";
open (OUT,">conceptMatrix_rm1overlap_chi.xls") or die "OpenError: conceptMatrix_chi.xls, $!\n";

my %DBconcept_id;my %DBid_concept;

while(<DB>){
	chomp;
	$_=~s/[\n\r]//g;
	my @data=split(/\t/);
	if((exists $data[2]) and ($data[2]=~/\d+/)){
		$DBconcept_id{"$data[0]\^$data[1]"}{$data[2]}=1;
		$DBid_concept{$data[2]}{"$data[0]\t$data[1]"}=1;
	}
}
close DB;
$date=&getTime();
print "$date->{hour}\:$date->{minute} $date->{second}s $date->{date}\t";
print "Loading concept database finished, begin filtering\n";

foreach my $concept (keys %DBconcept_id){
	delete $DBconcept_id{$concept} if (keys %{$DBconcept_id{$concept}} < 5);
}
$date=&getTime();
print "$date->{hour}\:$date->{minute} $date->{second}s $date->{date}\t";
print "Filtering finished, begin construction\n";

my %copyDBconcept_id=%DBconcept_id; my %conceptscore; my %score; my %scoreanno;
my %jaccardMatrix;
my $titleConcept=join("\t", sort keys %DBconcept_id);
print OUT "\t$titleConcept\n";

foreach my $concept (sort keys %DBconcept_id){
	print OUT "$concept";
	foreach my $copyConcept (sort keys %copyDBconcept_id){
		if(exists $jaccardMatrix{$copyConcept}{$concept}){
			print OUT "\t$jaccardMatrix{$copyConcept}{$concept}";
			next;
		}
		my %common=();
		my %total=();
		foreach my $geneID (keys %{$DBconcept_id{$concept}}){
			$common{$geneID}=1 if exists $copyDBconcept_id{$copyConcept}{$geneID};
			$total{$geneID}=1;
		}
		foreach my $geneID (keys %{$copyDBconcept_id{$copyConcept}}){
			$total{$geneID}=1;
		}
		my $c=keys %common;
		my $t=keys %total;
		if($c==1){
			$c=0;
		}
		my $score=$c/$t;
		$jaccardMatrix{$concept}{$copyConcept}=$c/$t;
		print OUT "\t$score";
	}
	print OUT "\n";
}
$date=&getTime();
print "$date->{hour}\:$date->{minute} $date->{second}s $date->{date}\t";
print "Construction finished\n";
close OUT;
close DB;

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
