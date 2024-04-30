
libname library 'S:\MESA\dbGaP2013July_v11\31735\PhenoGenotypeFiles';



proc format;
	value educ
		0='NO SCHOOLING'
		1='GRADES 1-8'
		2='GRADES 9-11'
		3='COMPLETED HIGH SCHOOL/GED'
		4='SOME COLLEGE BUT NO DEGREE'
		5='TECHNICAL SCHOOL CERTIFICATE'
		6='ASSOCIATE DEGREE'
		7='BACHELOR DEGREE'
		8='GRADUATE OR PROFESSIONAL SCHOOL'
		.='MISSING';

	value peduc
		1='NO SCHOOLING'
		2='SOME SCHOOLING BUT DID NOT COMPLETE HIGH SCHOOL'
		3='HIGH SCHOOL DEGREE'
		4='SOME COLLEGE BUT NO DEGREE'
		5='COLLEGE DEGREE'
		6='GRADUATE OR PROFESSIONAL SCHOOL'
		999='MISSING';

	value yn
		0=' No'
		1='Yes'
		999=' zMissing';

	value smoke 
		0='NEVER SMOKED'
		1='FOMER SMOKER - QUIT >1 YEAR AGO'
		2='FOMER SMOKER - QUIT <1 YEAR AGO'
		3='CURRENT SMOKER'
		4='DO NOT KNOW'
		9='DO NOT KNOW'
		.='MISSING';

	value educx
		0=' No degree'
		1=' HS Degree/GED'
		2='Any College Certificate/Degree'
		999=' zMissing';

	value educxx
		5='5: No degree'
		4='4: HS Degree/GED'
		3='3: Some College - No 4-year Degree'
		2='2: 4-year degree'
		1='1: graduate or professional school'
		999='Missing';

	value smokex
		0='Never'
		1=' Former (>1 Year)'
		2=' Current'
		999=' Missing';

	value obese
		0='Normal: BMI<=25'
		1=' Overweight: 25<BMI<30'
		2=' Obese: BMI>=30'
		999='Missing';

	value drink
		0='Non-Drinker'
		1=' <=2 Drinks/Day'
		2=' >2 Drinks/Day'
		999=' Missing';

	value exer
		0=' Low: <600 MET-min/wk'
		1=' Moderate: 600-2999 MET-min/wk'
		2='High: >=3000 MET-min/wk'
		999=' zMissing';

	value inc
		1='1: <$20,000'
		2='2: $20,000 - 39,999'
		3='3: $40,000 - 74,999'
		4='4: $75,000+'
		999='Missing';

	value incx
		5='5: <$20,000'
		4='4: $20,000 - 34,999'
		3='3: $35,000 - 49,999'
		2='2: $50,000 - 74,999'
		1='1: $75,000+'
		999='Missing';

	value tert
		1='High'
		2='Medium'
		3='Low';

	value agecat
		1='<65'
		2='65+';

run;






/*Import ID maps and lists*/

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
	by sidno;
run;




proc import
	datafile='R:\MESA\MESA_EPI_PCA\ID_mapping.csv'
	dbms=csv
	out=xwalk (rename=(dbgapid=sidno phenoid=idno))
	replace;
run;

proc sort;
	by sidno;
run;





/*Import occupation data*/

proc import
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\data\MESA_Occupation_Data.txt'
	dbms=tab
	out=occ
	replace;
run;


proc import
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\data\final_MESA_w2012ops.csv'
	dbms=csv
	out=occ2(keep=idno occprest_overall)
	replace;
run;

proc sort;
	by idno;
run;









/*Import main exam data*/

proc import
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\data\SHARe_Exam1Main.txt'
	dbms=tab
	out=exam1x
	replace;
	guessingrows=600;
run;

proc import
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\data\SHARe_Exam2Main.txt'
	dbms=tab
	out=exam2x
	replace;
	guessingrows=600;
run;

proc import
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\data\SHARe_Exam3Main.txt'
	dbms=tab
	out=exam3x
	replace;
	guessingrows=600;
run;

proc import
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\data\SHARe_Exam4Main.txt'
	dbms=tab
	out=exam4x
	replace;
	guessingrows=600;
run;

proc import
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\data\SHARe_Exam5Main.txt'
	dbms=tab
	out=exam5x
	replace;
	guessingrows=600;
run;





/*Import neighborhood and census data*/

proc import
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\data\SHARe_AncilMesaNeighborSurveyScales.txt'
	dbms=tab
	out=socenv
	replace;
	guessingrows=600;
run;

proc import
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\data\SHARe_AncilMesaNeighborCensusTract\SHARe_AncilMesaNeighborCensusTract.txt'
	dbms=tab
	out=census
	replace;
	guessingrows=600;
run;






/*Merge crosswalks with main exam and social data*/

data mesax;	
	merge xwalk xwalk2
		exam1x(keep=sidno evsmk1 cursmk1 cigsday1 bmi1c pamvcm1c curjob1 educ1 age1c gender1 race1c site1c)
		exam2x(keep=sidno smkstat2 cigsday2 bmi2c pamvcm2c age2c empstat2 curjob2 momschl2 dadschl2)
		exam3x(keep=sidno smkstat3 cigsday3 bmi3c pamvcm3c age3c hometyp3 owncar3 ownland3 invest3 empstat3 curjob3)
		exam4x(keep=sidno smkstat4 cigsday4 bmi4c age4c empstat4 curjob4)
		exam5x(keep=sidno smkstat5 cigsday5 bmi5c pamvcm5c exercm5c age5c curalc5--liqwk5 marital5 empstat5 curjob5 income5 site5c)
		census(keep=sidno a_f1_pc2 a_factor_ana exam where=(exam=5))
		socenv(keep=sidno a_socenv_ceb exam where=(exam=5));
	by sidno;
run;

proc sort;
	by idno;
run;








/*Merge on occupation data*/

data mesa2x;
	merge list(in=a) mesax occ(keep=idno socatcd1-socatcd3) occ2;
	by idno;
	if a;
run;







/*Import clock data*/

proc import
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\data\DNAm_Age\DNAm_Clock\MESA_DNAm_Age1264.csv'
	dbms=csv
	out=clocks
	replace;
run;


proc import
	datafile='S:\MESA/Methylation/Lauren_DNAmAge/data/DNAm_Age/DNAm_Clock/MESA_PoAm_Age1264.csv'
	dbms=csv
	out=clocks2(rename=(idno=idnox))
	replace;
run;

data clocks2;
	set clocks2;
	idno=input(substr(idnox, 2), 10.);
	drop idnox;
run;


proc import
	datafile='S:\MESA/Methylation/Lauren_DNAmAge/data/DNAm_Age/MESA_new_DNAmGrimAge_Exam5Age_everyone.csv'
	dbms=csv
	out=clocks3
	replace;
run;







/*Merge on clock data and create phenotypes*/

data mesa3;
	merge mesa2x clocks(keep=idno Hannum_DNAm--DNAmTIMP1) clocks2 clocks3;
	by idno;

	if educ1 in (0 1 2) then educ=0;
	else if educ1 in (3 4) then educ=1;
	else if educ1 in (5 6 7 8) then educ=2;

	if educ in (0 1) then educ2=1;
	else if educ=2 then educ2=0;

	if educ1 in (0 1 2) then educ3=5;
	else if educ1=3 then educ3=4;
	else if educ1 in (4 5 6) then educ3=3;
	else if educ1=7 then educ3=2;
	else if educ1=8 then educ3=1;

	if educ1=0 then edyears=0;
	else if educ1=1 then edyears=4.5;
	else if educ1=2 then edyears=10;
	else if educ1=3 then edyears=12;
	else if educ1=4 then edyears=13;
	else if educ1 in (5 6) then edyears=14;
	else if educ1=7 then edyears=16;
	else if educ1=8 then edyears=18;
	else edyears=.; 

	edyears2=18-edyears;

	peduc1=max(momschl2, dadschl2);

	if peduc1 in (1 2) then peduc=0;
	else if peduc1 in (3 4) then peduc=1;
	else if peduc1 in (5 6) then peduc=2;

	if peduc1=1 then pedyears=0;
	else if peduc1=2 then pedyears=6;
	else if peduc1=3 then pedyears=12;
	else if peduc1=4 then pedyears=13;
	else if peduc1=5 then pedyears=16;
	else if peduc1=6 then pedyears=18;
	else pedyears=.;

	if momschl2=1 then momedyears=0;
	else if momschl2=2 then momedyears=10;
	else if momschl2=3 then momedyears=12;
	else if momschl2=4 then momedyears=14;
	else if momschl2=5 then momedyears=16;
	else if momschl2=6 then momedyears=18;
	else momedyears=0;

	if dadschl2=1 then dadedyears=0;
	else if dadschl2=2 then dadedyears=10;
	else if dadschl2=3 then dadedyears=12;
	else if dadschl2=4 then dadedyears=14;
	else if dadschl2=5 then dadedyears=16;
	else if dadschl2=6 then dadedyears=18;
	else dadedyears=0;

	if smkstat5=0 then smoke=0;
	else if smkstat5=1 then smoke=1;
	else if smkstat5 in (2 3) then smoke=2;

	if 0<bmi5c<=25 then obesity=0;
	else if 25<bmi5c<30 then obesity=1;
	else if bmi5c>=30 then obesity=2;

	if curalc5=0 then drinks=0;
	else drinks=sum(of rwinewk5--liqwk5);

	if drinks=0 then drinkscat=0;
	else if 0<drinks<=14 then drinkscat=1;
	else if drinks>14 then drinkscat=2;

	if 0<=exercm5c<600 then exercise=0;
	else if 600<=exercm5c<3000 then exercise=1;
	else if exercm5c>=3000 then exercise=2;

	if substr(socatcd1, 1, 1)='1' then white_collar=1;
	else if substr(socatcd1, 1, 1) in ('2' '3' '4' '5' '6') then white_collar=0;
	else if (substr(socatcd1, 1, 1)='7' or socatcd1='') and (substr(socatcd2, 1, 1)='1' or substr(socatcd3, 1, 1)='1') then white_collar=1;
	else if (substr(socatcd1, 1, 1)='7' or socatcd1='') and (substr(socatcd2, 1, 1) in ('2' '3' '4' '5' '6') or substr(socatcd3, 1, 1) in ('2' '3' '4' '5' '6')) then white_collar=0;

	if marital5=1 then married=1;
	else if marital5 in (2 3 4 5) then married=0;

	if hometyp3 in (2 3) then homeown=1;
	else if hometyp3 in (1 4) then homeown=0;

	if owncar3 in (1 2) then carown=1;
	else if owncar3=0 then carown=0;

	if ownland3 in (1 2) then landown=1;
	else if ownland3=0 then landown=0;

	if (homeown=1 or landown=1) and carown=1 and invest3=1 then wealth=1;
	else if (homeown=1 or landown=1) and (carown=1 or invest3=1) then wealth=2;
	else if (homeown=1 or landown=1) and carown=0 and invest3=0 then wealth=3;
	else if homeown=0 and landown=0 and carown=1 and invest3=1 then wealth=3;
	else if homeown=0 and landown=0 and (carown=1 or invest3=1) then wealth=4;
	else if homeown=0 and landown=0 and carown=0 and invest3=0 then wealth=5;

	wealth2=sum(homeown, landown, carown, invest3);

	if income5 in (1 2 3 4 5) then income=5;
	else if income5 in (6 7 8) then income=4;
	else if income5 in (9 10) then income=3;
	else if income5=11 then income=2;
	else if income5 in (12 13 14 15) then income=1;

	if income5=1 then incomec=2500;
	else if income5=2 then incomec=6500;
	else if income5=3 then incomec=10000;
	else if income5=4 then incomec=14000;
	else if income5=5 then incomec=18000;
	else if income5=6 then incomec=22500;
	else if income5=7 then incomec=27500;
	else if income5=8 then incomec=32500;
	else if income5=9 then incomec=37500;
	else if income5=10 then incomec=45000;
	else if income5=11 then incomec=62500;
	else if income5=12 then incomec=87500;
	else if income5=13 then incomec=112500;
	else if income5=14 then incomec=137500;
	else if income5=15 then incomec=150000;


	A_SOCENV_CEB2=-1*A_SOCENV_CEB;
	A_F1_PC2r=-1*A_F1_PC2;

	zA_F1_PC2=(A_F1_PC2+0.3167194)/1.1091279;
	zSOCENV2=(A_SOCENV_CEB2-0.4805792)/2.6959082;
	occprest2=-1*occprest_overall;
	zoccprest2=-1*(occprest_overall-45.7414145)/13.4912960;

	if -4<zA_F1_PC2<=-0.4816959 then zA_F1_PC2q=1;
	else if -0.4816959<zA_F1_PC2<=0.1058026 then zA_F1_PC2q=2;
	else if 0.1058026<zA_F1_PC2<=0.3640933 then zA_F1_PC2q=3;
	else if 0.3640933<zA_F1_PC2<=0.8223676 then zA_F1_PC2q=4;
	else if zA_F1_PC2>0.8223676 then zA_F1_PC2q=5;

	if -18<a_factor_ana<=-3.5758540 then a_factor_anaq=5;
	else if -3.5758540<a_factor_ana<=-1.0182591 then a_factor_anaq=4;
	else if -1.0182591<a_factor_ana<=1.0854188 then a_factor_anaq=3;
	else if 1.0854188<a_factor_ana<=3.5741545 then a_factor_anaq=2;
	else if a_factor_ana>3.5741545 then a_factor_anaq=1;

	if -4<zSOCENV2<=-0.9525260 then zSOCENV2q=1;
	else if -0.9525260<zSOCENV2<=-0.3929517 then zSOCENV2q=2;
	else if -0.3929517<zSOCENV2<=0.2065530 then zSOCENV2q=3;
	else if 0.2065530<zSOCENV2<=0.9268341 then zSOCENV2q=4;
	else if zSOCENV2>0.9268341 then zSOCENV2q=5;

	if -4<zoccprest2<=-1.0362835 then zoccprest2q=1;
	else if -1.0362835<zoccprest2<=-0.1721240 then zoccprest2q=2;
	else if -0.1721240<zoccprest2<=0.3022588 then zoccprest2q=3;
	else if 0.3022588<zoccprest2<=0.9567827 then zoccprest2q=4;
	else if zoccprest2>0.9567827 then zoccprest2q=5;

	curjob=curjob5;
	array job[4] curjob1-curjob4;
	array emp[4] empstat2-empstat5;
	do i=1 to 4;
		if curjob=. and emp[5-i]=0 then curjob=job[5-i];
	end;

	SESn=n(educ3, income, wealth, a_factor_anaq, zSOCENV2q, zoccprest2q);

	if SESn>=3 then SESindex=sum(educ3, income, wealth, a_factor_anaq, zSOCENV2q, zoccprest2q)/SESn;

	if 0<sesindex<=2.3 then ses=1;
	else if 2.3<sesindex<3.2 then ses=2;
	else if sesindex>=3.2 then ses=3;

	if 0<age5c<65 then agecat=1;
	else if age5c>=65 then agecat=2;

	format educ3 educxx. income incx.;

	drop sidno curalc5--liqwk5;
run;



/*
proc means data=mesa3;
	var A_F1_PC2 A_SOCENV_CEB2 occprest_overall occprest2;
run;

proc freq data=mesa3;
	table educ1 educ3 income income5 wealth;
run;

proc means data=mesa3 min p20 p40 p60 p80 max;
	var zoccprest2 zSOCENV2 zA_F1_PC2 a_factor_ana;
run;

proc freq data=mesa3;
	table zoccprest2q zSOCENV2q zA_F1_PC2q curjob*income;
run;

proc means data=mesa3;
	var SESindex;
run;

proc freq data=mesa3;
	table SESindex ses sesn;
run;
*/






/*Create factor scores*/

ods pdf;
proc factor data=mesa3 nfactors=1 score out=mesa4;
	var edyears2 income wealth zA_F1_PC2 zSOCENV2 zoccprest2;
run;
ods pdf close;

proc factor data=mesa4 nfactors=1 score out=mesa4 prefix=Factorx;
	var edyears2 income wealth zoccprest2;
run;


/*
proc means data=mesa4;
	var factor1 factorx1;
run;

proc corr data=mesa4;
	var sesindex factor1 factorx1;
run;

proc freq data=mesa4;
	table factor1 factorx1;
run;
*/

proc freq data=mesa4;
	table sespc1 sespc1x;
run;



/*Create tertiles for factor scores*/

data mesa4;
	set mesa4;

	if -3<factor1<=-0.514376443 then sespc1=1;
	else if -0.514376443<factor1<=0.34564250787 then sespc1=2;
	else if factor1>0.3456425078 then sespc1=3;

	if -3<factorx1<=-0.51839736 then sespc1x=1;
	else if -0.51839736<factorx1<=0.3346353108 then sespc1x=2;
	else if factorx1>0.3346353108 then sespc1x=3;

	format ses sespc1 sespc1x tert.;
run;


	







/*Save raw data*/

proc export
	data=mesa4
	outfile='S:\MESA\Methylation\Lauren_DNAmAge\data\MESA_clock_phenos_v2.csv'
	dbms=csv
	replace;
run;







/*Add labels and formats*/

data MESA_clock_phenos;
	set mesa4;

	format educ1 educ. momschl2 dadschl2 peduc. evsmk1 cursmk1 white_collar married homeown carown landown invest3 yn. smkstat2-smkstat5 smoke. educ peduc educx. smoke smokex. obesity obese. drinkscat drink. exercise exer. agecat agecat.;

	label 
		idno='ID Number'
		chip_meth='Methylation Chip Number'
		pos_meth='Methylation Chip Position'
		age1c='Age at Exam 1'
		bcell='WBC Proportions: B cells'
		tcell='WBC Proportions: T cells'
		nkcell='WBC Proportions: Natural Killer cells'
		neutro='WBC Proportions: Neutrophils'
		bmi1c='BMI at Exam 1'
		educ1='Education: Highest Level Completed'
		evsmk1='Smoked at least 100 Cigarettes in Lifetime'
		cursmk1='Smoked in Last 30 Days'
		cigsday1='Average # Cigarettes Smoked Per Day'
		pamvcm1c='Moderate and Vigorous Physical Activity Total (MET-min/wk)'
		age2c='Age at Exam 2'
		bmi2c='BMI at Exam 2'
		smkstat2='Smoking Status'
		cigsday2='Average # Cigarettes Smoked Per Day'
		pamvcm2c='Moderate and Vigorous Physical Activity Total (MET-min/wk)'
		age3c='Age at Exam 3'
		bmi3c='BMI at Exam 3'
		smkstat3='Smoking Status'
		cigsday3='Average # Cigarettes Smoked Per Day'
		pamvcm3c='Moderate and Vigorous Physical Activity Total (MET-min/wk)'
		age4c='Age at Exam 4'
		bmi4c='BMI at Exam 4'
		smkstat4='Smoking Status'
		cigsday4='Average # Cigarettes Smoked Per Day'
		age5c='Age at Exam 5'
		bmi5c='BMI at Exam 5'
		smkstat5='Smoking Status'
		cigsday5='Average # Cigarettes Smoked Per Day'
		pamvcm5c='Moderate and Vigorous Physical Activity Total (MET-min/wk)'
		exercm5c='Total Intentional Exercise (MET-min/wk)'
		A_F1_PC2='Neighborhood Socioeconomic Disadvantage Score - Average'
		A_SOCENV_CEB='Neighborhood Social Environment Score - Average'
		A_SOCENV_CEB2='Neighborhood Social Environment Score - Average - Reverse Coded'
		zA_F1_PC2='Standardized Neighborhood Socioeconomic Disadvantage Score - Average'
		zSOCENV2='Standardized Neighborhood Social Environment Score - Average - Reverse Coded'
		Hannum_DNAm='Hannum epigenetic age using Morgan’s code'
		Lin_99_CpG='Lin epigenetic age using Morgan’s code'
		Yang_385_CpG='Yang epigenetic age using Morgan’s code'
		Levine_DNAm_PhenoAge='Levine epigenetic age using Morgan’s code'
		Horvath_391CPG='Horvath epigenetic age using Morgan’s code'
		Horvath_DNAm='Horvath epigenetic age from online calculator'
		PredictedGender='Predicted gender from online calculator'
		DNAmGrimAge='Grim epigenetic age from online calculator by using advanced analysis option'
		DNAmADM='Grim epigenetic age component - ADM'
		DNAmB2M='Grim epigenetic age component - B2M'
		DNAmCystatinC='Grim epigenetic age component - Cystatin C'
		DNAmGDF15='Grim epigenetic age component - GDF15'
		DNAmLeptin='Grim epigenetic age component - Leptin'
		DNAmPACKYRS='Grim epigenetic age component - Pack Years'
		DNAmPAI1='Grim epigenetic age component - PAI1'
		DNAmTIMP1='Grim epigenetic age component - TIMP1'
		educ='Education: Highest Level Completed'
		educ2='Education: Less than College Degree'
		edyears='Education: Years of Education'
		peduc='Parental Education: Highest Level Completed'
		pedyears='Parental Education: Years of Education'
		smoke='Smoking Status'
		obesity='Obesity Category'
		drinks='Total Drinks per Week'
		drinkscat='Total Drinks per Week Category'
		exercise='Total Intentional Exercise Category'
		white_collar='Management/Professional Occupation Indicator'
		occprest_overall='Occupational Prestige Score'
		zoccprest2='Standardized Occupational Prestige Score - Reverse Coded'
		socatcd1='Exam 1: BLS Standard occupation category code'
		socatcd2='Exam 2: BLS Standard occupation category code'
		socatcd3='Exam 3: BLS Standard occupation category code';
run;






/*
proc freq data=MESA_clock_phenos;
	table site;
run;


proc corr data=MESA_clock_phenos;
	var Hannum_DNAm Lin_99_CpG Yang_385_CpG Levine_DNAm_PhenoAge Horvath_391CPG Horvath_DNAm DNAmGrimAge;
	with age5c;
run;



ods pdf;
proc glm data=MESA_clock_phenos;
	class site;
	model A_F1_PC2=site / solution;
	means site;
run;

proc glm data=MESA_clock_phenos;
	class site;
	model A_SOCENV_CEB=site / solution;
	means site;
run;

proc glm data=MESA_clock_phenos;
	class site;
	model occprest_overall=site / solution;
	means site;
run;

proc freq data=MESA_clock_phenos;
	table (drinkscat educ exercise obesity peduc smoke)*site / chisq;
run;
ods pdf close;



proc sort data=MESA_clock_phenos out=bysite;
	by site;
run;

ods pdf;
proc reg data=bysite;
	ods select fitplot;
	var A_SOCENV_CEB A_F1_PC2;
	model A_F1_PC2=A_SOCENV_CEB;
run;

proc reg data=bysite;
	ods select fitplot;
	by site;
	var A_SOCENV_CEB A_F1_PC2;
	model A_F1_PC2=A_SOCENV_CEB;
run;

proc sgplot data=bysite;	
	reg x=A_SOCENV_CEB y=A_F1_PC2 / group=site;
run;
ods pdf close;

proc means data=MESA_clock_phenos;
	class race;
	var zSOCENV2;
run;
*/








/*Create age accelaration variables*/

proc reg data=MESA_clock_phenos plots=none;
	model Hannum_DNAm Lin_99_CpG Yang_385_CpG Levine_DNAm_PhenoAge Horvath_391CPG Horvath_DNAm DNAmGrimAge PoAm_age DNAmADM DNAmB2M DNAmCystatinC DNAmGDF15 DNAmLeptin DNAmPACKYRS DNAmPAI1 DNAmTIMP1=age5c;
	output out=MESA_clock_phenos2 r=Hannum_DNAm_AA Lin_99_CpG_AA Yang_385_CpG_AA DNAm_PhenoAge_AA Horvath_391_AA Horvath_DNAm_AA DNAmGrimAge_AA PoAm_age_AA DNAmADM_AA DNAmB2M_AA DNAmCystatinC_AA DNAmGDF15_AA DNAmLeptin_AA DNAmPACKYRS_AA DNAmPAI1_AA DNAmTIMP1_AA; 
run;
quit;

/*
proc reg data=MESA_clock_phenos plots=none;
	where race='WHITE, CAUCASIAN';
	model Hannum_DNAm Lin_99_CpG Yang_385_CpG Levine_DNAm_PhenoAge Horvath_391CPG Horvath_DNAm DNAmGrimAge DNAmADM DNAmB2M DNAmCystatinC DNAmGDF15 DNAmLeptin DNAmPACKYRS DNAmPAI1 DNAmTIMP1=age5c;
	output out=MESA_clock_phenos2ea r=Hannum_DNAm_AAx Lin_99_CpG_AAx Yang_385_CpG_AAx DNAm_PhenoAge_AAx Horvath_391_AAx Horvath_DNAm_AAx DNAmGrimAge_AAx DNAmADM_AA DNAmB2M_AAx DNAmCystatinC_AAx DNAmGDF15_AAx DNAmLeptin_AAx DNAmPACKYRS_AAx DNAmPAI1_AAx DNAmTIMP1_AAx; 
run;
quit;

proc reg data=MESA_clock_phenos plots=none;
	where race='BLACK, AFRICAN-AME';
	model Hannum_DNAm Lin_99_CpG Yang_385_CpG Levine_DNAm_PhenoAge Horvath_391CPG Horvath_DNAm DNAmGrimAge DNAmADM DNAmB2M DNAmCystatinC DNAmGDF15 DNAmLeptin DNAmPACKYRS DNAmPAI1 DNAmTIMP1=age5c;
	output out=MESA_clock_phenos2aa r=Hannum_DNAm_AAx Lin_99_CpG_AAx Yang_385_CpG_AAx DNAm_PhenoAge_AAx Horvath_391_AAx Horvath_DNAm_AAx DNAmGrimAge_AAx DNAmADM_AA DNAmB2M_AAx DNAmCystatinC_AAx DNAmGDF15_AAx DNAmLeptin_AAx DNAmPACKYRS_AAx DNAmPAI1_AAx DNAmTIMP1_AAx; 
run;
quit;

proc reg data=MESA_clock_phenos plots=none;
	where race='HISPANIC';
	model Hannum_DNAm Lin_99_CpG Yang_385_CpG Levine_DNAm_PhenoAge Horvath_391CPG Horvath_DNAm DNAmGrimAge DNAmADM DNAmB2M DNAmCystatinC DNAmGDF15 DNAmLeptin DNAmPACKYRS DNAmPAI1 DNAmTIMP1=age5c;
	output out=MESA_clock_phenos2ha r=Hannum_DNAm_AAx Lin_99_CpG_AAx Yang_385_CpG_AAx DNAm_PhenoAge_AAx Horvath_391_AAx Horvath_DNAm_AAx DNAmGrimAge_AAx DNAmADM_AA DNAmB2M_AAx DNAmCystatinC_AAx DNAmGDF15_AAx DNAmLeptin_AAx DNAmPACKYRS_AAx DNAmPAI1_AAx DNAmTIMP1_AAx;
run;
quit;

data MESA_clock_phenos2;
	merge MESA_clock_phenos2 
		MESA_clock_phenos2ea(keep=idno Hannum_DNAm_AAx Lin_99_CpG_AAx Yang_385_CpG_AAx DNAm_PhenoAge_AAx Horvath_391_AAx Horvath_DNAm_AAx DNAmGrimAge_AAx DNAmADM_AA DNAmB2M_AAx DNAmCystatinC_AAx DNAmGDF15_AAx DNAmLeptin_AAx DNAmPACKYRS_AAx DNAmPAI1_AAx DNAmTIMP1_AAx)
		MESA_clock_phenos2aa(keep=idno Hannum_DNAm_AAx Lin_99_CpG_AAx Yang_385_CpG_AAx DNAm_PhenoAge_AAx Horvath_391_AAx Horvath_DNAm_AAx DNAmGrimAge_AAx DNAmADM_AA DNAmB2M_AAx DNAmCystatinC_AAx DNAmGDF15_AAx DNAmLeptin_AAx DNAmPACKYRS_AAx DNAmPAI1_AAx DNAmTIMP1_AAx)
		MESA_clock_phenos2ha(keep=idno Hannum_DNAm_AAx Lin_99_CpG_AAx Yang_385_CpG_AAx DNAm_PhenoAge_AAx Horvath_391_AAx Horvath_DNAm_AAx DNAmGrimAge_AAx DNAmADM_AA DNAmB2M_AAx DNAmCystatinC_AAx DNAmGDF15_AAx DNAmLeptin_AAx DNAmPACKYRS_AAx DNAmPAI1_AAx DNAmTIMP1_AAx);
	by idno;
run;
*/









/*Standardize clocks and create missing categories*/

proc means data=MESA_clock_phenos2;
	where age5c^=.;
	var Hannum_DNAm_AA Lin_99_CpG_AA Yang_385_CpG_AA DNAm_PhenoAge_AA Horvath_391_AA Horvath_DNAm_AA DNAmGrimAge_AA
		Hannum_DNAm Lin_99_CpG Yang_385_CpG Levine_DNAm_PhenoAge Horvath_391CPG Horvath_DNAm DNAmGrimAge PoAm_age;
run;

/*
proc means data=MESA_clock_phenos2;
	var Hannum_DNAm_AAx Lin_99_CpG_AAx Yang_385_CpG_AAx DNAm_PhenoAge_AAx Horvath_391_AAx Horvath_DNAm_AAx DNAmGrimAge_AAx;
run;
*/

data MESA_clock_phenos3;
	set MESA_clock_phenos2;
	if age5c^=.;

	array cats[*] educ educ3 peduc momschl2 dadschl2 smoke obesity drinkscat exercise married homeown carown landown invest3 income wealth;
	array miss[*] educmiss educ3miss peducmiss momschl2miss dadschl2miss smokemiss obesitymiss drinkscatmiss exercisemiss marriedmiss homeownmiss carownmiss landownmiss invest3miss incomemiss wealthmiss;
	do i=1 to dim(cats);
		if cats[i]=. then cats[i]=999;

		if cats[i]=999 then miss[i]=1;
		else miss[i]=0;
	end;

	if abs(Hannum_DNAm_AA)>5*4.6315609 then Hannum_DNAm_AA=.;
	if abs(Lin_99_CpG_AA)>5*6.4890115 then Lin_99_CpG_AA=.;
	if abs(Yang_385_CpG_AA)>5*0.0057556 then Yang_385_CpG_AA=.;
	if abs(DNAm_PhenoAge_AA)>5*5.9074863 then DNAm_PhenoAge_AA=.;
	if abs(Horvath_391_AA)>5*3.3859589 then Horvath_391_AA=.;
	if abs(Horvath_DNAm_AA)>5*5.3733591 then Horvath_DNAm_AA=.;
	if abs(DNAmGrimAge_AA)>5*4.1334598 then DNAmGrimAge_AA=.;

	if abs(Hannum_DNAm-73.0160058)>5*8.5472361 then Hannum_DNAm=.;
	if abs(Lin_99_CpG-82.5411814)>5*10.5781521 then Lin_99_CpG=.;
	if abs(Yang_385_CpG-0.0894892)>5*0.0059374 then Yang_385_CpG=.;
	if abs(Levine_DNAm_PhenoAge-71.7649907)>5*9.2759005 then Levine_DNAm_PhenoAge=.;
	if abs(Horvath_391CPG-69.1860598)>5*7.5468500 then Horvath_391CPG=.;
	if abs(Horvath_DNAm-63.9096833)>5*8.6272638 then Horvath_DNAm=.;
	if abs(DNAmGrimAge-79.3998827)>5*8.0972077 then DNAmGrimAge=.;
	if abs(PoAm_age-1.0596210)>5*0.0753813 then PoAm_age=.;
run;





/*
proc freq data=MESA_clock_phenos3;
	table educmiss peducmiss smokemiss obesitymiss drinkscatmiss exercisemiss marriedmiss homeownmiss carownmiss landownmiss invest3miss incomemiss;
run;

proc means data=MESA_clock_phenos3;
	var zoccprest2 zA_F1_PC2 zSOCENV2;
run;

ods pdf;
proc means data=MESA_clock_phenos3 N mean std median min max;
	var age5c Horvath_DNAm Horvath_391CPG Hannum_DNAm Levine_DNAm_PhenoAge Lin_99_CpG Yang_385_CpG DNAmGrimAge PoAm_age 
		sesindex occprest2 A_F1_PC2 A_SOCENV_CEB2 zoccprest2 zA_F1_PC2 zSOCENV2;
run;

proc freq data=MESA_clock_phenos3;
	table gender race married educ3 income wealth smoke drinkscat obesity site;
run;
ods pdf close;


ods pdf;
proc corr data=MESA_clock_phenos3;
	var age5c Horvath_DNAm Horvath_391CPG Hannum_DNAm Levine_DNAm_PhenoAge Lin_99_CpG Yang_385_CpG DNAmGrimAge PoAm_age
		Horvath_DNAm_AA Horvath_391_AA Hannum_DNAm_AA DNAm_PhenoAge_AA Lin_99_CpG_AA Yang_385_CpG_AA DNAmGrimAge_AA PoAm_age_AA;
run;
ods pdf close;

ods pdf;
proc means;
	where momschl2miss=0;
	var momedyears;
run;

proc means;
	where dadschl2miss=0;
	var dadedyears;
run;

proc freq;
	table momschl2miss dadschl2miss;
run;
ods pdf close;
*/









/*Macro for running regressions and creating tables*/

%macro mix(out, pred);
proc mixed data=MESA_clock_phenos3;
	class /*&pred(ref=first)*/ gender race married chip_meth pos_meth;
	model &out=&pred age5c gender race married bcell tcell nkcell neutro / solution;
	random chip_meth pos_meth;
	ods output SolutionF=solutions1; 
run;

proc mixed data=MESA_clock_phenos3;
	class /*&pred(ref=first)*/ gender race married chip_meth pos_meth smoke drinkscat obesity;
	model &out=&pred age5c gender race married smoke drinkscat obesity bcell tcell nkcell neutro / solution;
	random chip_meth pos_meth;
	ods output SolutionF=solutions2; 
run;

proc mixed data=MESA_clock_phenos3;
	class /*&pred(ref=first)*/ gender race married chip_meth pos_meth smoke drinkscat obesity site(ref=first);
	model &out=&pred age5c gender race married smoke drinkscat obesity site bcell tcell nkcell neutro / solution;
	random chip_meth pos_meth;
	ods output SolutionF=solutions3; 
run;

data solutions1;
	retain clock Effect Category;
	length clock $16. Category $30.;
	set solutions1;
	clock="&out";

	array preds[*] /*&pred*/ married;
	array preds2[*] $ gender race;
	Category='';
	do i=1 to dim(preds);
		if strip(Category) in ('' '_') then Category=vvalue(preds[i]);
	end; 

	do i=1 to dim(preds2);
		if strip(Category) in ('' '_') then Category=vvalue(preds2[i]);
	end;

	if Effect not in ('bcell' 'tcell' 'nkcell' 'neutro');
	keep clock effect category estimate stderr probt;
run;

proc sort;
	by clock effect category;
run;

data solutions2;
	retain clock Effect Category;
	length clock $16. Category $30.;
	set solutions2;
	clock="&out";

	array preds[*] /*&pred*/ married smoke drinkscat obesity;
	array preds2[*] $ gender race;
	Category='';
	do i=1 to dim(preds);
		if strip(Category) in ('' '_') then Category=vvalue(preds[i]);
	end; 

	do i=1 to dim(preds2);
		if strip(Category) in ('' '_') then Category=vvalue(preds2[i]);
	end;

	if Effect not in ('bcell' 'tcell' 'nkcell' 'neutro');
	keep clock effect category estimate stderr probt;
	rename estimate=estimate2 stderr=stderr2 probt=probt2;
run;

proc sort;
	by clock effect category;
run;

data solutions3;
	retain clock Effect Category;
	length clock $16. Category $30.;
	set solutions3;
	clock="&out";

	array preds[*] /*&pred*/ married smoke drinkscat obesity;
	array preds2[*] $ gender race site;
	Category='';
	do i=1 to dim(preds);
		if strip(Category) in ('' '_') then Category=vvalue(preds[i]);
	end; 

	do i=1 to dim(preds2);
		if strip(Category) in ('' '_') then Category=vvalue(preds2[i]);
	end;

	if Effect not in ('bcell' 'tcell' 'nkcell' 'neutro');
	keep clock effect category estimate stderr probt;
	rename estimate=estimate3 stderr=stderr3 probt=probt3;
run;

data blank;
	set solutions3;
	keep clock;
run;

proc sort nodupkey;
	by clock;
run;

data solutions3;
	set blank solutions3;
	n=_N_;
run;

proc sort;
	by clock effect category;
run;

data all;
	merge solutions1 solutions2 solutions3;
	by clock effect category;
run;

proc sort;
	by n;
run;

data &out&pred;
	set all;
	drop n;
run;
%mend;


%mix(DNAmGrimAge, sesindex)
%mix(Levine_DNAm_PhenoAge, sesindex)
%mix(Hannum_DNAm, sesindex)
%mix(Horvath_391CPG, sesindex)
%mix(Horvath_DNAm, sesindex)
%mix(Lin_99_CpG, sesindex)
%mix(Yang_385_CpG, sesindex)
%mix(PoAm_age, sesindex)

/*
%mix(DNAmGrimAge_AA, sespc1)
%mix(DNAm_PhenoAge_AA, sespc1)
%mix(Hannum_DNAm_AA, sespc1)
%mix(Horvath_391_AA, sespc1)
%mix(Horvath_DNAm_AA, sespc1)
%mix(Lin_99_CpG_AA, sespc1)
%mix(Yang_385_CpG_AA, sespc1)
%mix(PoAm_age, sespc1)

%mix(DNAmGrimAge_AA, sespc1x)
%mix(DNAm_PhenoAge_AA, sespc1x)
%mix(Hannum_DNAm_AA, sespc1x)
%mix(Horvath_391_AA, sespc1x)
%mix(Horvath_DNAm_AA, sespc1x)
%mix(Lin_99_CpG_AA, sespc1x)
%mix(Yang_385_CpG_AA, sespc1x)
%mix(PoAm_age, sespc1x)
*/




/*Create tables*/

data sestable1;
	set DNAmGrimAgesesindex Levine_DNAm_PhenoAgesesindex Hannum_DNAmsesindex Horvath_391CPGsesindex Horvath_DNAmsesindex Lin_99_CpGsesindex Yang_385_CpGsesindex PoAm_agesesindex;
run;


/*
data sespc1table1;
	set DNAmGrimAge_AAsespc1 DNAm_PhenoAge_AAsespc1 Hannum_DNAm_AAsespc1 Horvath_391_AAsespc1 Horvath_DNAm_AAsespc1 Lin_99_CpG_AAsespc1 Yang_385_CpG_AAsespc1 PoAm_agesespc1;
run;

data sespc1xtable1;
	set DNAmGrimAge_AAsespc1x DNAm_PhenoAge_AAsespc1x Hannum_DNAm_AAsespc1x Horvath_391_AAsespc1x Horvath_DNAm_AAsespc1x Lin_99_CpG_AAsespc1x Yang_385_CpG_AAsespc1x PoAm_agesespc1x;
run;
*/



ods rtf;
proc print data=sestable1(drop=clock) noobs;
run;
ods rtf close;








/*Individual SES components*/

%mix(DNAmGrimAge, zoccprest2)
%mix(DNAmGrimAge, za_f1_pc2)
%mix(DNAmGrimAge, zsocenv2)
%mix(PoAm_age, zoccprest2)
%mix(PoAm_age, za_f1_pc2)
%mix(PoAm_age, zsocenv2)


*uncomment &pred in mix() for these;
%mix(DNAmGrimAge, educ3)
%mix(DNAmGrimAge, income)
%mix(DNAmGrimAge, wealth)
%mix(PoAm_age, educ3)
%mix(PoAm_age, income)
%mix(PoAm_age, wealth)


data occtable;
	set DNAmGrimAgezoccprest2 PoAm_agezoccprest2;
run;

data pc2table;
	set DNAmGrimAgeza_f1_pc2 PoAm_ageza_f1_pc2;
run;

data socenvtable;
	set DNAmGrimAgezsocenv2 PoAm_agezsocenv2;
run;

data eductable;
	set DNAmGrimAgeeduc3 PoAm_ageeduc3;
run;

data inctable;
	set DNAmGrimAgeincome PoAm_ageincome;
run;

data wealthtable;
	set DNAmGrimAgewealth PoAm_agewealth;
run;


ods csv;
proc print data=occtable noobs;
run;
ods csv close;

ods csv;
proc print data=pc2table noobs;
run;
ods csv close;

ods csv;
proc print data=socenvtable noobs;
run;
ods csv close;

ods csv;
proc print data=eductable noobs;
run;
ods csv close;

ods csv;
proc print data=inctable noobs;
run;
ods csv close;

ods csv;
proc print data=wealthtable noobs;
run;
ods csv close;








/*Macro for stratified regression analyses and table creation*/

%macro mix2(out, pred, strata);
proc sort data=MESA_clock_phenos3(rename=(Levine_DNAm_PhenoAge=DNAm_PhenoAge)) out=sort;
	by &strata;
run;

%if &strata=gender %then %let cov=race;
%else %if &strata=race %then %let cov=gender;
%else %if &strata=agecat %then %let cov=gender race;

proc mixed data=sort;
	by &strata;
	class /*&pred(ref=first)*/ gender race married chip_meth pos_meth;
	model &out=&pred age5c &cov married bcell tcell nkcell neutro / solution;
	random chip_meth pos_meth;
	ods output SolutionF=solutions1; 
run;

proc mixed data=sort;
	by &strata;
	class /*&pred(ref=first)*/ gender race married chip_meth pos_meth smoke drinkscat obesity;
	model &out=&pred age5c &cov married smoke drinkscat obesity bcell tcell nkcell neutro / solution;
	random chip_meth pos_meth;
	ods output SolutionF=solutions2; 
run;

proc mixed data=sort;
	by &strata;
	class /*&pred(ref=first)*/ gender race married chip_meth pos_meth smoke drinkscat obesity site;
	model &out=&pred age5c &cov married smoke drinkscat obesity site bcell tcell nkcell neutro / solution;
	random chip_meth pos_meth;
	ods output SolutionF=solutions3; 
run;

data solutions1;
	retain clock strata Effect Category;
	length clock $16. Category $30.;
	set solutions1;
	clock="&out";
	strata=put(&strata, 12.);

	array preds[*] /*&pred*/ married;
	array preds2[*] $ &cov;
	Category='';
	do i=1 to dim(preds);
		if strip(Category) in ('' '_') then Category=vvalue(preds[i]);
	end; 

	do i=1 to dim(preds2);
		if strip(Category) in ('' '_') then Category=vvalue(preds2[i]);
	end;

	if Effect not in ('bcell' 'tcell' 'nkcell' 'neutro');
	keep clock strata effect category estimate stderr probt;
run;

proc sort;
	by clock strata effect category;
run;

data solutions2;
	retain clock strata Effect Category;
	length clock $16. Category $30.;
	set solutions2;
	clock="&out";
	strata=put(&strata, 12.);

	array preds[*] /*&pred*/ married smoke drinkscat obesity;
	array preds2[*] $ &cov;
	Category='';
	do i=1 to dim(preds);
		if strip(Category) in ('' '_') then Category=vvalue(preds[i]);
	end; 

	do i=1 to dim(preds2);
		if strip(Category) in ('' '_') then Category=vvalue(preds2[i]);
	end;

	if Effect not in ('bcell' 'tcell' 'nkcell' 'neutro');
	keep clock strata effect category estimate stderr probt;
	rename estimate=estimate2 stderr=stderr2 probt=probt2;
run;

proc sort;
	by clock strata effect category;
run;

data solutions3;
	retain clock strata Effect Category;
	length clock $16. Category $30.;
	set solutions3;
	clock="&out";
	strata=put(&strata, 12.);

	array preds[*] /*&pred*/ married smoke drinkscat obesity;
	array preds2[*] $ &cov site;
	Category='';
	do i=1 to dim(preds);
		if strip(Category) in ('' '_') then Category=vvalue(preds[i]);
	end; 

	do i=1 to dim(preds2);
		if strip(Category) in ('' '_') then Category=vvalue(preds2[i]);
	end;

	if Effect not in ('bcell' 'tcell' 'nkcell' 'neutro');
	keep clock strata effect category estimate stderr probt;
	rename estimate=estimate3 stderr=stderr3 probt=probt3;
run;

data blank;
	set solutions3;
	keep clock strata;
run;

proc sort nodupkey;
	by clock strata;
run;

data solutions3;
	set blank solutions3;
run;

proc sort;
	by clock strata;
run;

data solutions3;
	set solutions3;
	n=_N_;
run;

proc sort;
	by clock strata effect category;
run;

data all;
	merge solutions1 solutions2 solutions3;
	by clock strata effect category;
run;

proc sort;
	by n;
run;

data &out&pred&strata;
	set all;
	drop n;
run;
%mend;



%mix2(DNAmGrimAge, sesindex, race)
%mix2(DNAm_PhenoAge, sesindex, race)
%mix2(Hannum_DNAm, sesindex, race)
%mix2(Horvath_391CPG, sesindex, race)
%mix2(Horvath_DNAm, sesindex, race)
%mix2(Lin_99_CpG, sesindex, race)
%mix2(Yang_385_CpG, sesindex, race)
%mix2(PoAm_age, sesindex, race)

%mix2(DNAmGrimAge, sesindex, gender)
%mix2(DNAm_PhenoAge, sesindex, gender)
%mix2(Hannum_DNAm, sesindex, gender)
%mix2(Horvath_391CPG, sesindex, gender)
%mix2(Horvath_DNAm, sesindex, gender)
%mix2(Lin_99_CpG, sesindex, gender)
%mix2(Yang_385_CpG, sesindex, gender)
%mix2(PoAm_age, sesindex, gender)

%mix2(DNAmGrimAge, sesindex, agecat)
%mix2(DNAm_PhenoAge, sesindex, agecat)
%mix2(Hannum_DNAm, sesindex, agecat)
%mix2(Horvath_391CPG, sesindex, agecat)
%mix2(Horvath_DNAm, sesindex, agecat)
%mix2(Lin_99_CpG, sesindex, agecat)
%mix2(Yang_385_CpG, sesindex, agecat)
%mix2(PoAm_age, sesindex, agecat)



data sestable1strat;
	retain clock strata;
	set DNAmGrimAgesesindex DNAmGrimAgesesindexrace DNAmGrimAgesesindexgender DNAmGrimAgesesindexagecat
		Levine_DNAm_PhenoAgesesindex DNAm_PhenoAgesesindexrace DNAm_PhenoAgesesindexgender DNAm_PhenoAgesesindexagecat 
		Hannum_DNAmsesindex Hannum_DNAmsesindexrace Hannum_DNAmsesindexgender Hannum_DNAmsesindexagecat 
		Horvath_391CPGsesindex Horvath_391CPGsesindexrace Horvath_391CPGsesindexgender Horvath_391CPGsesindexagecat 
		Horvath_DNAmsesindex Horvath_DNAmsesindexrace Horvath_DNAmsesindexgender Horvath_DNAmsesindexagecat 
		Lin_99_CpGsesindex Lin_99_CpGsesindexrace Lin_99_CpGsesindexgender Lin_99_CpGsesindexagecat 
		Yang_385_CpGsesindex Yang_385_CpGsesindexrace Yang_385_CpGsesindexgender Yang_385_CpGsesindexagecat 
		PoAm_agesesindex PoAm_agesesindexrace PoAm_agesesindexgender PoAm_agesesindexagecat;

	if effect='SESindex' or (strata='' and effect='');
run;





ods rtf;
proc print data=sestable1strat(drop=clock effect) noobs;
run;
ods rtf close;



















%macro bsmix(out, pred);
proc mixed data=MESA_clock_phenos3;
	class gender race married chip_meth pos_meth smoke drinkscat obesity;
	model &out=&pred age5c gender race married smoke drinkscat obesity bcell tcell nkcell neutro / solution;
	random chip_meth pos_meth;
	ods output SolutionF=solutions1; 
run;

proc mixed data=MESA_clock_phenos3;
	class gender race married chip_meth pos_meth smoke drinkscat obesity;
	model &out=&pred age5c gender race married smoke drinkscat obesity bcell tcell nkcell neutro momedyears momschl2miss / solution;
	random chip_meth pos_meth;
	ods output SolutionF=solutions2; 
run;

proc mixed data=MESA_clock_phenos3;
	class gender race married chip_meth pos_meth smoke drinkscat obesity;
	model &out=&pred age5c gender race married smoke drinkscat obesity bcell tcell nkcell neutro momedyears momschl2miss dadedyears dadschl2miss / solution;
	random chip_meth pos_meth;
	ods output SolutionF=solutions3; 
run;

data solutions1;
	retain clock Effect Category;
	length clock $16. Category $30.;
	set solutions1;
	clock="&out";

	array preds[*] married smoke drinkscat obesity;
	array preds2[*] $ gender race;
	Category='';
	do i=1 to dim(preds);
		if strip(Category) in ('' '_') then Category=vvalue(preds[i]);
	end; 

	do i=1 to dim(preds2);
		if strip(Category) in ('' '_') then Category=vvalue(preds2[i]);
	end;

	if Effect not in ('bcell' 'tcell' 'nkcell' 'neutro');
	keep clock effect category estimate stderr probt;
run;

proc sort;
	by clock effect category;
run;

data solutions2;
	retain clock Effect Category;
	length clock $16. Category $30.;
	set solutions2;
	clock="&out";

	array preds[*] married smoke drinkscat obesity;
	array preds2[*] $ gender race;
	Category='';
	do i=1 to dim(preds);
		if strip(Category) in ('' '_') then Category=vvalue(preds[i]);
	end; 

	do i=1 to dim(preds2);
		if strip(Category) in ('' '_') then Category=vvalue(preds2[i]);
	end;

	if Effect not in ('bcell' 'tcell' 'nkcell' 'neutro');
	keep clock effect category estimate stderr probt;
	rename estimate=estimate2 stderr=stderr2 probt=probt2;
run;

proc sort;
	by clock effect category;
run;

data solutions3;
	retain clock Effect Category;
	length clock $16. Category $30.;
	set solutions3;
	clock="&out";

	array preds[*] married smoke drinkscat obesity;
	array preds2[*] $ gender race;
	Category='';
	do i=1 to dim(preds);
		if strip(Category) in ('' '_') then Category=vvalue(preds[i]);
	end; 

	do i=1 to dim(preds2);
		if strip(Category) in ('' '_') then Category=vvalue(preds2[i]);
	end;

	if Effect not in ('bcell' 'tcell' 'nkcell' 'neutro');
	keep clock effect category estimate stderr probt;
	rename estimate=estimate3 stderr=stderr3 probt=probt3;
run;

proc sort;
	by clock effect category;
run;

data &out&pred;
	merge solutions1 solutions2 solutions3;
	by clock effect category;
run;
%mend;


%bsmix(DNAmGrimAge, sesindex)
%bsmix(Levine_DNAm_PhenoAge, sesindex)
%bsmix(Hannum_DNAm, sesindex)
%bsmix(Horvath_391CPG, sesindex)
%bsmix(Horvath_DNAm, sesindex)
%bsmix(Lin_99_CpG, sesindex)
%bsmix(Yang_385_CpG, sesindex)
%bsmix(PoAm_age, sesindex)




/*Create tables*/

data sestable1;
	set DNAmGrimAgesesindex Levine_DNAm_PhenoAgesesindex Hannum_DNAmsesindex Horvath_391CPGsesindex Horvath_DNAmsesindex Lin_99_CpGsesindex Yang_385_CpGsesindex PoAm_agesesindex;
run;


ods csv;
proc print data=sestable1 noobs;
run;
ods csv close;

