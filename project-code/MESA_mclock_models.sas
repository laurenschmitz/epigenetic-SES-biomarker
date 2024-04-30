

/*import MESA phenotypes created in MESA phenos.sas code*/
proc import
	out=mesa
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\data\MESA_clock_phenos_v2.csv'
	dbms=csv
	replace;
run;

data mesa;
	set mesa;

	if momschl2=1 then meducq=4;
	else if momschl2=2 then meducq=3;
	else if momschl2=3 then meducq=2;
	else if momschl2=4 then meducq=1;
	else if momschl2 in (5 6) then meducq=0;

	if dadschl2=1 then feducq=4;
	else if dadschl2=2 then feducq=3;
	else if dadschl2=3 then feducq=2;
	else if dadschl2=4 then feducq=1;
	else if dadschl2 in (5 6) then feducq=0;

	peducq=mean(meducq, feducq);
	cSES=peducq;
run;

proc export
	data=mesa
	outfile='S:\MESA\Methylation\Lauren_DNAmAge\data\MESA_clock_phenos_v3.csv'
	dbms=csv
	replace;
run;



/*import SES methylation biomarkers created in HRS_ElasticNet_XValid.R code*/
proc import
	out=mclock
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\results\HRS_Enet\MESA_SESindex_mclock_450K_raw.csv'
	dbms=csv
	replace;
run;

proc import
	out=mclock2
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\results\HRS_Enet\MESA_cSES_mclock_450K_raw.csv'
	dbms=csv
	replace;
run;


/*import main exam6 data*/
proc import
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\data\SHARe_Exam6Main.txt'
	dbms=tab
	out=exam6
	replace;
	guessingrows=600;
run;


/*import events data*/
proc import
	out=events
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\data\SHARe_EventsThruYear2019_DS.txt'
	dbms=tab
	replace;
run;


/*import dementia events data*/
proc import
	out=demevents
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\data\SHARe_ThruYear2018Events_Dementia.txt'
	dbms=tab
	replace;
run;


/*import APOE data*/
proc import
	out=apoe(keep=pin apoe_imputed e4 rename=(pin=subject_id))
	datafile='S:\MESA\KT ApoE Social\APOE_MESA.csv'
	dbms=csv
	replace;
run;

proc sort;
	by subject_id;
run;




/*import MPOA data*/
proc import
	datafile='S:\MESA\Methylation\ELA_DNAmAge\Data\MESA_DunedinPACE_new_all_CpG.csv'
	dbms=csv
	out=mpoa
	replace;
run;






/*import ID mapping files*/
proc import
	datafile='R:\MESA\MESA_EPI_PCA\ID_mapping.csv'
	dbms=csv
	out=xwalk (rename=(dbgapid=sidno phenoid=idno))
	replace;
run;

proc sort;
	by idno;
run;

proc import
	datafile='R:\MESA\Methylation_081413\MESA_Epi_METH_idmap.txt'
	dbms=tab
	out=list
	replace;
run;

proc import
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\data\MESA_Epi_METH_idmap_SIDNO.txt'
	dbms=tab
	out=list2
	replace;
run;

proc sort data=list out=test;
	by chip_meth pos_meth age1c gender site bcell;
run;

proc sort data=list2;
	by chip_meth pos_meth age1c gender site bcell;
run;

data xwalk2;
	merge test list2(in=a);
	by chip_meth pos_meth age1c gender site bcell;
	if a;
	keep idno sidno;
run;

proc sort;
	by idno;
run;




/*import PCs*/
data pcs;
	infile 'R:\MESA\MESA_EPI_PCA\Combined.txt';
	input idno PC1-PC50;
	keep idno pc1-pc4;
run;

proc sort;
	by idno;
run;





/*merge data with IDNO*/
data mesa2;
	merge xwalk(rename=(sidno=subject_id)) xwalk2(rename=(sidno=subject_id)) mesa(in=a) mclock mclock2 mpoa pcs;
	by idno;
	if a;
run;

proc sort;
	by subject_id;
run;






/*merge data by SIDNO*/
data mesa3;
	merge mesa2(in=a) exam6(in=b) events demevents apoe /*minda(where=(subject_id^=.))*/;
	by subject_id;
	if a;

	if mi=1 or ang=1 or chf=1 or chda=1 or cvda=1 then heart=1;
	else if mi=0 and ang=0 and chf=0 and chda=0 and cvda=0 then heart=0;

	if strk=1 or tia=1 then stroke=1;
	else if strk=0 and tia=0 then stroke=0;

	if emphys6=1 then copd=1;
	else if emphys6=0 then copd=0;
/******decided not to use count of conditions in MESA*********
	if arthritism6=1 then arthr=1;
	else if arthritism6=0 then arthr=0;

	array psy[10] ptsdm6 bipolm6 schizm6 depressm6 anxietym6 ocdm6 autismm6 adhdm6 dyslexm6 psyothm6;
	do i=1 to 10;
		if psy[i]=9 then psy[i]=.;
	end;

	npsych=n(ptsdm6, bipolm6, schizm6, depressm6, anxietym6, ocdm6, autismm6, adhdm6, dyslexm6, psyothm6);
	if ptsdm6=1 or bipolm6=1 or schizm6=1 or depressm6=1 or anxietym6=1 or ocdm6=1 or autismm6=1 or adhdm6=1 or dyslexm6=1 or psyothm6=1 then psych=1;
	else if npsych>7 then psych=0;

	nconditions=n(htn6c, dm036t, heart, stroke, copd, arthr, psych);
	if nconditions>4 then conditions=sum(htn6c, dm036t, heart, stroke, copd, arthr, psych);
*/
run;




proc freq data=mesa3;
	table icddementia_dh dth genhlth6 htn6c dm036t heart stroke copd /*arthr psych npsych conditions nconditions*/;* / missing;
run;

proc means data=mesa3;
	var genhlth6;
run;




/*standardize and take means of conditions to create MSI variables*/
proc stdize data=mesa3 out=mesa4 oprefix sprefix=z;
	var htn6c dm036t heart stroke sesindex_mclock cses_mclock sesindex cses;
run;

data mesa4;
	set mesa4;

	if n(zhtn6c, zdm036t, zheart, zstroke)=4 then MSI=mean(zhtn6c, zdm036t, zheart, zstroke);
run;

proc stdize data=mesa4 out=mesa4 oprefix sprefix=z;
	var MSI;
run;



proc means;
	var zMSI;
run;





/*This macro will do all health outcomes at once*/
%macro sesmodels(pred, ses, model) / minoperator;
%let out=/*conditions*/ zMSI genhlth6 dth icddementia_dh;
%do i=1 %to 4;
	%let outi=%scan(&out, &i);
	%if &i in (1 2) %then %do;
		proc glm data=mesa4 plots=diagnostics;
			class race gender smoke / ref=first;
			model &outi=&pred age5c race gender smoke &ses / clparm solution;
			ods output ParameterEstimates=estimates;
		run;

		data estimates&i;
			set estimates;
			if parameter="&pred";
			drop biased--tvalue;
			rename estimate=estimate&model probt=probt&model lowercl=lowercl&model uppercl=uppercl&model;
		run;
	%end;	
/*
	%else %if &i=1 %then %do;
		proc genmod data=mesa4;
			class race gender smoke / ref=first param=ref;
			model &outi=&pred age5c race gender smoke &ses / dist=poisson;
			ods output ParameterEstimates=estimates;
		run;

		data estimates&i;
			retain dependent parameter IRR;
			set estimates;
			if parameter="&pred";
			dependent="&outi";
			IRR=exp(estimate);
			rename IRR=estimate&model probchisq=probt&model;
			keep dependent parameter IRR probchisq;
		run;
	%end;
*/
	%else %if &i=3 %then %do;
		proc logistic data=mesa4 desc;
			class race gender smoke / ref=first param=ref;
			model &outi=&pred age5c race gender smoke &ses / expb rsq;
			ods output ParameterEstimates=estimates oddsratios=ors;
			ods output rsquare=rsquare(keep=label1 cvalue1 cvalue2);
		run;

		data ors;
			set ors;
			if effect="&pred";
			keep lowercl uppercl; 
			rename lowercl=lowercl&model uppercl=uppercl&model;
		run;

		data estimates&i;
			retain dependent variable expest;
			set estimates;
			if variable="&pred";
			dependent="&outi";
			rename expest=estimate&model probchisq=probt&model variable=parameter;
			keep dependent variable expest probchisq;
		run;

		data estimates&i;
			merge estimates&i ors;
		run;
	%end;

	%else %if &i=4 %then %do;
		proc logistic data=mesa4 desc;
			class race gender smoke e4 / ref=first param=ref;
			model &outi=&pred age5c race gender smoke e4 &ses / expb rsq;
			ods output ParameterEstimates=estimates oddsratios=ors;
			ods output rsquare=rsquare(keep=label1 cvalue1 cvalue2);
		run;

		data ors;
			set ors;
			if effect="&pred";
			keep lowercl uppercl; 
			rename lowercl=lowercl&model uppercl=uppercl&model;
		run;

		data estimates&i;
			retain dependent variable expest;
			set estimates;
			if variable="&pred";
			dependent="&outi";
			rename expest=estimate&model probchisq=probt&model variable=parameter;
			keep dependent variable expest probchisq;
		run;

		data estimates&i;
			merge estimates&i ors;
		run;
	%end;

	data estimatesall&model;
		set estimates1 estimates2 estimates3 estimates4 /*estimates5*/;
	run;
%end;
%mend;

%sesmodels(zSESindex_mclock,,1)
%sesmodels(zSESindex,,1star)
%sesmodels(zSESindex_mclock,SESindex,2)
%sesmodels(zSESindex_mclock,dnamgrimage,3a)
%sesmodels(zSESindex_mclock,SESindex dnamgrimage,4a)
%sesmodels(zSESindex_mclock,dunedinpace,3b)
%sesmodels(zSESindex_mclock,SESindex dunedinpace,4b)
%sesmodels(zSESindex_mclock,cses,5)

/*run appropriate models above before this step*/
data estimatesall_SES;
	merge estimatesall1star estimatesall1 estimatesall2 estimatesall3a estimatesall4a estimatesall3b estimatesall4b estimatesall5;
run;

ods csv;
proc print data=estimatesall_SES noobs;
run;
ods csv close;



%sesmodels(zcSES_mclock,,1)
%sesmodels(zcSES,,1star)
%sesmodels(zcSES_mclock,cSES,2)
%sesmodels(zcSES_mclock,dnamgrimage,3a)
%sesmodels(zcSES_mclock,cSES dnamgrimage,4a)
%sesmodels(zcSES_mclock,dunedinpace,3b)
%sesmodels(zcSES_mclock,cSES dunedinpace,4b)
%sesmodels(zcSES_mclock,sesindex,5)

/*run appropriate models above before this step*/
data estimatesall_cSES;
	merge estimatesall1star estimatesall1 estimatesall2 estimatesall3a estimatesall4a estimatesall3b estimatesall4b estimatesall5;
run;


ods csv;
proc print data=estimatesall_cSES noobs;
run;
ods csv close;












/*macro for finding freqs of categorical variables for Table of descriptors*/
%macro freqs(var);
proc freq data=mesa4;
	where zsesindex_mclock^=.;
	table &var / out=freqstot;
run;

data freqs&var;
	retain Variable;
	length Category $20.;
	set freqstot;

	Category=vvalue(&var);
	Variable="&var";

	pct=put(percent, 6.2);
	stat=cat(count, ' (', strip(pct), ')'); 
	keep variable category stat;
run;
%mend;


%freqs(gender)
%freqs(race)
%freqs(smoke)
%freqs(educ3)
%freqs(income)
%freqs(wealth)
%freqs(meducq)
%freqs(feducq)
%freqs(dth)
%freqs(icddementia_dh)





/*macro for finding means of continuous variables for Table1*/
%macro means(var);
proc means data=mesa4 mean std;
	where zsesindex_mclock^=.;
	var &var;
	ods output summary=means;
run;

data means&var;
	set means;
	Variable="&var";
	mean=put(&var._mean, 20.2);
	std=put(&var._stddev, 20.2);
	stat=cat(strip(mean), ' (', strip(std), ')');
	keep variable stat;
run;
%mend;


%means(age5c)
%means(zoccprest2)
%means(a_factor_ana)
%means(A_SOCENV_CEB2)
%means(SESindex)
%means(dnamgrimage)
%means(dunedinpace)
%means(cSES)
%means(zMSI)
%means(genhlth6)




/*constructing Table 1 of descriptors*/
data table1;
	retain variable category;
	length variable $12.;
	set meansage5c freqsgender freqsrace freqssmoke freqseduc3 freqsincome meanszoccprest2 meansa_factor_ana meansA_SOCENV_CEB2 meanssesindex meansdnamgrimage meansdunedinpace;
run;


ods rtf;
proc print data=table1 noobs;
run;
ods rtf close;







/*look at predictiveness of the methylation clock on SES index*/
ods pdf;
proc glm data=mesa4;
	class race gender;
	model SESindex=SESindex_mclock age5c gender /*race bcell tcell nkcell neutro pc1 pc2 pc3 pc4 smoke*/ / solution;
*	ods output modelanova=ma overallanova=oa;
run;

proc glm data=mesa4;
	class race gender;
	model SESindex=age5c gender /*race bcell tcell nkcell neutro pc1 pc2 pc3 pc4 smoke*/ / solution;
run;
ods pdf close;




/*look at predictiveness of the methylation clock on cSES index*/
ods pdf;
proc glm data=mesa4;
	class race gender;
	model cSES=cSES_mclock age5c gender race bcell tcell nkcell neutro pc1 pc2 pc3 pc4 smoke/**/ / solution;
*	ods output modelanova=ma overallanova=oa;
run;

proc glm data=mesa4;
	class race gender;
	model cSES=age5c gender race bcell tcell nkcell neutro pc1 pc2 pc3 pc4 smoke/**/ / solution;
run;
ods pdf close;






proc sql;
select 
    ma.Dependent,
    ma.source, 
    ma.SS / oa.SS as type3_PartialRsquare format=percentn7.2,
    ma.probF
from 
    ma, oa
where oa.source="Corrected Total";
quit;





/*clock correlations*/
option orientation=landscape;
ods csv;
proc corr data=mesa4 nomiss nosimple;
	var SESindex_mclock cSES_mclock SESindex cses dnamgrimage dunedinpace;
run;
ods csv close;

