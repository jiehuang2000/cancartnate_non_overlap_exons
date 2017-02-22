%LET path = /folders/myshortcuts/BIOS_511_SAS/SAS_files_practice/Catherine/hg19_bed;
/*
PROC IMPORT DATAFILE = "&path/refGene_hg19_112716_3primeUTRExons_normChr_zinbaFormat.bed"
			OUT = UTR3exon
			DBMS = dlm
			REPLACE;
			DELIMITER = '090D'x;
			GETNAMES= No;
			GUESSINGROWS= 1000;  */
		
DATA DavisLab.UTR3exon;
	INFILE "&path/refGene_hg19_112716_3primeUTRExons_normChr_zinbaFormat.bed" DLM='09'x TRUNCOVER;
	input name :$60. chr :$5. start :10. stop :10. strand $1. sigval $1.;
	run;

DATA DavisLab.UTR5exon;
	INFILE "&path/refGene_hg19_112716_5primeUTRExons_normChr_zinbaFormat.bed" DLM='09'x TRUNCOVER;
	input name :$60. chr :$5. start :10. stop :10. strand $1. sigval $1.;
	run;
	
DATA DavisLab.CodingExon;
	INFILE "&path/refGene_hg19_112716_codingExons_normChr_removeDup_zinbaFormat.bed" DLM='09'x TRUNCOVER;
	input name :$60. chr :$5. start :10. stop :10. strand $1. sigval $1.;
	run;
	
DATA DavisLab.allexon;
	SET DavisLab.UTR3exon DavisLab.UTR5exon DavisLab.CodingExon;
	LENGTH ID $20 TYPE $5 EXON_ORDER 4;
	NMchr = SCAN(name,1,'_');
	ID = SCAN(name, 2, '_');
	TYPE = SCAN(name, 3, '_');
	If type = 'utr5' then type = 'a_utr5';
	EXON_ORDER = SCAN(name, 4, '_');
	IF NMchr = "NM";
	DROP NMchr;
	RUN;
	
PROC SORT DATA=DavisLab.allexon;
	*BY ID TYPE EXON_ORDER;
	BY CHR START STOP;
	RUN;

PROC MEANS DATA=DAVISLAB.ALLEXON MAX;
    VAR EXON_ORDER;
	CLASS TYPE;
	RUN;
	
* FIND OUT ALL THE DUPLICATED EXONS;	
DATA allexon2;
	SET DAVISLAB.ALLEXON;
	RETAIN REP 0;
	CHRP = LAG(CHR); STARTP = LAG(START); STOPP = LAG(STOP); REPP = LAG(REP);
	IF CHR = CHRP AND START = STARTP AND STOP = STOPP THEN REP = 1;
	ELSE REP=0;
	DROP CHRP -- REPP;
	RUN;

DATA allexon3;
	merge allexon2 (firstobs=1) allexon2 (firstobs=2 rename= (name=namep chr=chrp start=startp stop=stopp strand=strandp rep=repp));
	rep2 = max(rep, repp);
	KEEP NAME -- EXON_ORDER REP2;
	RENAME REP2 = repeated_exon;
	run;

PROC SORT DATA=allexon3;
	BY ID TYPE EXON_ORDER;
	RUN;	



*ALL GENE(isoforms) and the number of repeated exon in the gene;
DATA NON_OVERLAP_GENE;
	SET ALLEXON3;
	BY ID;
	RETAIN COUNT_REP 0;
	IF FIRST.ID THEN COUNT_REP = 0;
	COUNT_REP = COUNT_REP +  REPEATED_EXON;
	IF LAST.ID THEN OUTPUT;
	KEEP ID COUNT_REP;
	RUN;

*We only want genes with 0 repeated exons;

PROC SORT DATA=DavisLab.allexon; BY ID TYPE EXON_ORDER; RUN;
PROC SORT DATA=NON_OVERLAP_GENE; BY ID; RUN;

DATA DAVISLAB.genes_with_non_overlap_exon;
	MERGE DAVISLAB.ALLEXON NON_OVERLAP_GENE;
	BY ID;
	IF COUNT_REP = 0;
	SIGVAL = 0;
	KEEP NAME -- SIGVAL;
	RUN;


proc export data=DAVISLAB.GENES_WITH_NON_OVERLAP_EXON 
   outfile="&PATH/GENES_WITH_NON_OVERLAP_EXON.bed" 
   dbms=tab;
run;



data temp; set davislab.genes_with_non_overlap_exon; diff=stop-start; keep diff;run;
proc sort; by descending diff;run;
proc print data=temp(obs=50); run;

* the longest exon is 22562 bp;

data coord;
	INFILE "/folders/myshortcuts/BIOS_511_SAS/SAS_files_practice/Catherine/hg19_bed/WT1.coord" 
	       DLM='09'x TRUNCOVER firstobs=2;
	input name2 :$60. chr :$5. start :10. stop :10. strand :$1. sigval $8. seq $32767.;
	len = countw(seq,'09'x);
	diff = stop - start + 1;
	run;	
proc sort data=coord; by chr start stop; run;
proc sort data=DAVISLAB.genes_with_non_overlap_exon; by chr start stop; run;

data davislab.wt1coord2;
	merge DAVISLAB.genes_with_non_overlap_exon (in=in1) coord (in=in2);
	by chr start stop;
	if in2;
	run;







	
	
	