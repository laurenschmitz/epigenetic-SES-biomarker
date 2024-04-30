
libname library 'H:\HRS\RAND\randhrs1992_2020v1_SAS';
libname track 'H:\HRS\trk2018v2a';
libname hrs98 'H:\HRS\1998\DATA\SAS';
libname hrs00 'H:\HRS\2000\DATA\SAS';
libname hrs02 'H:\HRS\2002\DATA\SAS';
libname hrs04 'H:\HRS\2004\DATA\SAS';
libname hrs06 'H:\HRS\2006\DATA\SAS';
libname hrs08 'H:\HRS\2008\DATA\SAS';
libname hrs10 'H:\HRS\2010\DATA\SAS';
libname hrs12 'H:\HRS\2012\DATA\SAS';
libname hrs14 'H:\HRS\2014\DATA\SAS';
libname hrs16 'H:\HRS\2016\DATA\SAS';
libname hrs20 'H:\HRS\2020\DATA\SAS';
libname vbs 'H:\HRS\VBS';
libname hcap 'H:\HRS\HCAP\HC16\HC16summary';
libname langa 'H:\HRS\LangaWeir2020';
libname apoe 'H:\HRS\APOE_serotonin_release 2021';
libname xwalk 'R:\Health_Retirement_Study\Genotype_Jan2015_Processing\Preliminary_GWAS\Phenotype';
libname factors 'S:\Social_Epigenomics';
libname acs 'R:\Health_Retirement_Study\HRS-CDR-Census-v2-SAS\Tract';
libname geo 'R:\Health_Retirement_Study\GeoCode\Detail\HRSxDetail2018v2\data\SAS'; 
libname edimp 'S:\Health_Retirement_Study\Episodic_Memory_R03\From Colter\Pheno_March17_Jessica';
libname occ 'R:\Health_Retirement_Study\Industry Occupation Data\2016\ioxwave16\data';
libname cses 'S:\MESA\Methylation\Lauren_DNAmAge\data\cSESmeasures';
libname clocks 'R:\Health_Retirement_Study\Epigenetic Clocks\Clocks';
libname data 'S:\MESA\Methylation\Lauren_DNAmAge\data';



proc format;
	value race
		1='EA'
		2='AA'
		3='HA'
		4='Other';

	value smoke
		0='0:Never'
		1='1:Former'
		2='2:Current';
run;



/*import xwalk and methylation sample files*/
proc import
	datafile='R:\Health_Retirement_Study\Methylation_Apr_2021\transfer2\xwalk.csv'
	out=methy(rename=(hhid=chhid pn=cpn id=sample))
	dbms=csv
	replace;
run;

data methy;
	length hhid $6. PN $3.;
	set methy;

	PN=cat('0', put(cpn,2.));

	if length(strip(chhid))=5 then hhid=cat('0', put(chhid,5.));
	else if length(strip(chhid))=6 then hhid=put(chhid,6.);

	drop cpn chhid;
run;

proc sort;
	by FID;
run;



/*import cell proportions*/
proc import 
	datafile='R:\Health_Retirement_Study\Methylation_Apr_2021\transfer2\DNAm_cellpropest_HRSn4018.csv'
	dbms=csv
	out=cells
	replace;
run;

proc sort;
	by sample;
run;



/*import technical variables*/
proc import 
	datafile='R:\Health_Retirement_Study\Methylation_Apr_2021\transfer04.13.2021\finalPD.csv'
	dbms=csv
	out=plate(drop=PN)
	replace;
run;



/*combine info files above*/
data info1;
	merge methy plate;
	by FID;
run;

proc sort;
	by sample;
run;

data info;
	merge info1 cells;
	by sample;
run;

proc sort;
	by hhid pn;
run;





/*Cleans up each wave of Leave-Behind data for the social environment variables.  These are then inputs to the HRS_socenv.R program which creates the final scales and saves them to HRS_socenv_06_16.csv which is imported below*/ 
/*
data lb06;
	merge track.Trk2018tr_r(keep=hhid pn gender kage kinsamp) geo.Hrsxgeo18v8b_r(where=(year=2006) keep=hhid pn year linkcen2010 stfips10) hrs06.H06lb_r(keep=hhid pn klb021a--klb021h);
	by hhid pn;
	if kinsamp=1;

	klb021b=8-klb021b;
	klb021d=8-klb021d;
	klb021h=8-klb021h;

	hhidpn=hhid||pn;

	m1=.;
	m2=.;
	m3=.;
run;

proc transpose data=lb06 out=so_co(rename=(col1=so_co));
	var klb021a klb021c klb021e klb021g;
	by hhidpn linkcen2010 gender kage stfips10;
run;

proc transpose data=lb06 out=aes_qual(rename=(col1=aes_qual));
	var klb021b klb021f klb021h m1;
	by hhidpn;
run;

proc transpose data=lb06 out=safe(rename=(col1=safe));
	var klb021d m1 m2 m3;
	by hhidpn;
run;

data lb06x;
	merge so_co aes_qual(keep=aes_qual) safe(keep=safe);
	if so_co^=. or aes_qual^=. or safe^=.;
run;
 
proc export
	data=lb06x
	outfile='S:\MESA\Methylation\Lauren_DNAmAge\data\lb06.csv'
	dbms=csv
	replace;
run;





data lb08;
	merge track.Trk2018tr_r(keep=hhid pn gender lage linsamp) geo.Hrsxgeo18v8b_r(where=(year=2008) keep=hhid pn year linkcen2010 stfips10) hrs08.H08lb_r(keep=hhid pn llb021a--llb021h);
	by hhid pn;
	if linsamp=1;

	hhidpn=hhid||pn;

	m1=.;
	m2=.;
	m3=.;
run;

proc transpose data=lb08 out=so_co(rename=(col1=so_co));
	var llb021a llb021c llb021e llb021g;
	by hhidpn linkcen2010 gender lage stfips10;
run;

proc transpose data=lb08 out=aes_qual(rename=(col1=aes_qual));
	var llb021b llb021f llb021h m1;
	by hhidpn;
run;

proc transpose data=lb08 out=safe(rename=(col1=safe));
	var llb021d m1 m2 m3;
	by hhidpn;
run;

data lb08x;
	merge so_co aes_qual(keep=aes_qual) safe(keep=safe);
	if so_co^=. or aes_qual^=. or safe^=.;
run;
 
proc export
	data=lb08x
	outfile='S:\MESA\Methylation\Lauren_DNAmAge\data\lb08.csv'
	dbms=csv
	replace;
run;





data lb10;
	merge track.Trk2018tr_r(keep=hhid pn gender mage minsamp) geo.Hrsxgeo18v8b_r(where=(year=2010) keep=hhid pn year linkcen2010 stfips10) hrs10.H10lb_r(keep=hhid pn mlb021a--mlb021h);
	by hhid pn;
	if minsamp=1;

	hhidpn=hhid||pn;

	m1=.;
	m2=.;
	m3=.;
run;

proc transpose data=lb10 out=so_co(rename=(col1=so_co));
	var mlb021a mlb021c mlb021e mlb021g;
	by hhidpn linkcen2010 gender mage stfips10;
run;

proc transpose data=lb10 out=aes_qual(rename=(col1=aes_qual));
	var mlb021b mlb021f mlb021h m1;
	by hhidpn;
run;

proc transpose data=lb10 out=safe(rename=(col1=safe));
	var mlb021d m1 m2 m3;
	by hhidpn;
run;

data lb10x;
	merge so_co aes_qual(keep=aes_qual) safe(keep=safe);
	if so_co^=. or aes_qual^=. or safe^=.;
run;
 
proc export
	data=lb10x
	outfile='S:\MESA\Methylation\Lauren_DNAmAge\data\lb10.csv'
	dbms=csv
	replace;
run;






data lb12;
	merge track.Trk2018tr_r(keep=hhid pn gender nage ninsamp) geo.Hrsxgeo18v8b_r(where=(year=2012) keep=hhid pn year linkcen2010 stfips10) hrs12.H12lb_r(keep=hhid pn nlb021a--nlb021h);
	by hhid pn;
	if ninsamp=1;

	hhidpn=hhid||pn;

	m1=.;
	m2=.;
	m3=.;
run;

proc transpose data=lb12 out=so_co(rename=(col1=so_co));
	var nlb021a nlb021c nlb021e nlb021g;
	by hhidpn linkcen2010 gender nage stfips10;
run;

proc transpose data=lb12 out=aes_qual(rename=(col1=aes_qual));
	var nlb021b nlb021f nlb021h m1;
	by hhidpn;
run;

proc transpose data=lb12 out=safe(rename=(col1=safe));
	var nlb021d m1 m2 m3;
	by hhidpn;
run;

data lb12x;
	merge so_co aes_qual(keep=aes_qual) safe(keep=safe);
	if so_co^=. or aes_qual^=. or safe^=.;
run;
 
proc export
	data=lb12x
	outfile='S:\MESA\Methylation\Lauren_DNAmAge\data\lb12.csv'
	dbms=csv
	replace;
run;





data lb14;
	merge track.Trk2018tr_r(keep=hhid pn gender oage oinsamp) geo.Hrsxgeo18v8b_r(where=(year=2014) keep=hhid pn year linkcen2010 stfips10) hrs14.H14lb_r(keep=hhid pn olb020a--olb020h);
	by hhid pn;
	if oinsamp=1;

	hhidpn=hhid||pn;

	m1=.;
	m2=.;
	m3=.;
run;

proc transpose data=lb14 out=so_co(rename=(col1=so_co));
	var olb020a olb020c olb020e olb020g;
	by hhidpn linkcen2010 gender oage stfips10;
run;

proc transpose data=lb14 out=aes_qual(rename=(col1=aes_qual));
	var olb020b olb020f olb020h m1;
	by hhidpn;
run;

proc transpose data=lb14 out=safe(rename=(col1=safe));
	var olb020d m1 m2 m3;
	by hhidpn;
run;

data lb14x;
	merge so_co aes_qual(keep=aes_qual) safe(keep=safe);
	if so_co^=. or aes_qual^=. or safe^=.;
run;
 
proc export
	data=lb14x
	outfile='S:\MESA\Methylation\Lauren_DNAmAge\data\lb14.csv'
	dbms=csv
	replace;
run;





data lb16;
	merge track.Trk2018tr_r(keep=hhid pn gender page pinsamp) geo.Hrsxgeo18v8b_r(where=(year=2016) keep=hhid pn year linkcen2010 stfips10) hrs16.H16lb_r(keep=hhid pn plb020a--plb020h);
	by hhid pn;
	if pinsamp=1;

	hhidpn=hhid||pn;

	m1=.;
	m2=.;
	m3=.;
run;

proc transpose data=lb16 out=so_co(rename=(col1=so_co));
	var plb020a plb020c plb020e plb020g;
	by hhidpn linkcen2010 gender page stfips10;
run;

proc transpose data=lb16 out=aes_qual(rename=(col1=aes_qual));
	var plb020b plb020f plb020h m1;
	by hhidpn;
run;

proc transpose data=lb16 out=safe(rename=(col1=safe));
	var plb020d m1 m2 m3;
	by hhidpn;
run;

data lb16x;
	merge so_co aes_qual(keep=aes_qual) safe(keep=safe);
	if so_co^=. or aes_qual^=. or safe^=.;
run;
 
proc export
	data=lb16x
	outfile='S:\MESA\Methylation\Lauren_DNAmAge\data\lb16.csv'
	dbms=csv
	replace;
run;
*/



/*import social environment scales created above and in separate R program*/
proc import
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\data\HRS_socenv_06_16.csv'
	out=socenv
	dbms=csv
	replace;
run;

proc sort;
	by hhid pn;
run;





/*import census variables including factor ana*/
proc sort data=geo.Hrsxgeo18v8b_r out=geo(where=(year=2016));
	by LINKCEN2010;
run;

data census;
	merge geo(keep=hhid pn LINKCEN2010 in=a) acs.Tract10_pct_201216acs(keep=LINKCEN2010 ppv1216 mhinc1216) factors.Factorscore_1216(rename=(GEOID10=LINKCEN2010));
	by LINKCEN2010;
	if a;
run;

proc sort;
	by hhid pn;
run;





/*create occupational prestige variable across all waves.  This is saved to data.HRS_occ_prestige*/
/*
proc import
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\data\Occ_Prestige_Data_1980-2010_for_Scott.csv'
	out=occp
	dbms=csv
	replace;
run;

proc sort data=occp out=occp80 nodupkey;
	by occ_1980;
run;

proc sort data=occp out=occp00 nodupkey;
	by occ_2000;
run;

proc sort data=occp out=occp10 nodupkey;
	by occ_2010;
run;

data ioxwave16;
	set occ.ioxwave16;

	occ_1992=max(v2720, v2821, v3407, v3432, v4302, v4402, v4437);
	occ_1994=max(w3609, w4312b, w7007);
	occ_1996=max(e2666, e2732, e3134);
	occ_1998=max(f3187, f3255, f3652);
	occ_2000=max(g3436, g3505, g3962, g4339);
	occ_2002=max(hj168, hk010, hm025b);
	occ_2004=max(jj062_00, jj168_00, jk010_00);
	if occ_2004 in (. 999 9999) then occ_2004=max(JMW201_1_00, JMW201_2_00, JMW201_3_00, JMW201_4_00);
	occ_2006=max(KJ168, KK010);
	if occ_2006 in (. 999 9999) then occ_2006=max(KMW201_1, KMW201_2, KMW201_3, KMW201_4);
	occ_2008=max(LJ062, LJ168, LK010);
	if occ_2008 in (. 999 9999) then occ_2008=max(LMW201_1, LMW201_2, LMW201_3, LMW201_4);
	occ_2010=max(MJ062_2010, MJ168_2010, MK010_2010);
	if occ_2010 in (. 999 9999) then occ_2010=max(MMW201_1_2010, MMW201_2_2010, MMW201_3_2010, MMW201_4_2010);
	occ_2012=max(NJ062, NJ168, NK010);
	if occ_2012 in (. 999 9999) then occ_2012=max(NMW201_1, NMW201_2, NMW201_3, NMW201_4);
	occ_2014=max(OJ168, OK010);
	if occ_2014 in (. 999 9999) then occ_2014=max(OMW201_1, OMW201_2, OMW201_3, OMW201_4);
	occ_2016=max(PJ168, PK010, PMW201_1, PMW201_2, PMW201_3, PMW201_4);
	if occ_2016 in (. 999 9999) then occ_2016=max(PMW201_1, PMW201_2, PMW201_3, PMW201_4);
run;

proc sort data=ioxwave16;
	by occ_1992;
run;

data occ;
	merge ioxwave16(in=a) occp80(keep=occ_1980 occ_prest_1980 rename=(occ_1980=occ_1992 occ_prest_1980=occ_prest_1992));
	by occ_1992;
	if a;
	if occ_1992=. then occ_prest_1992=.;
run;

proc sort;
	by occ_1994;
run;
 
data occ;
	merge occ(in=a) occp80(keep=occ_1980 occ_prest_1980 rename=(occ_1980=occ_1994 occ_prest_1980=occ_prest_1994));
	by occ_1994;
	if a;
	if occ_1994=. then occ_prest_1994=.;
run;

proc sort;
	by occ_1996;
run;
 
data occ;
	merge occ(in=a) occp80(keep=occ_1980 occ_prest_1980 rename=(occ_1980=occ_1996 occ_prest_1980=occ_prest_1996));
	by occ_1996;
	if a;
	if occ_1996=. then occ_prest_1996=.;
run;

proc sort;
	by occ_1998;
run;
 
data occ;
	merge occ(in=a) occp80(keep=occ_1980 occ_prest_1980 rename=(occ_1980=occ_1998 occ_prest_1980=occ_prest_1998));
	by occ_1998;
	if a;
	if occ_1998=. then occ_prest_1998=.;
run;

proc sort;
	by occ_2000;
run;
 
data occ;
	merge occ(in=a) occp80(keep=occ_1980 occ_prest_1980 rename=(occ_1980=occ_2000 occ_prest_1980=occ_prest_2000));
	by occ_2000;
	if a;
	if occ_2000=. then occ_prest_2000=.;
run;

proc sort;
	by occ_2002;
run;
 
data occ;
	merge occ(in=a) occp80(keep=occ_1980 occ_prest_1980 rename=(occ_1980=occ_2002 occ_prest_1980=occ_prest_2002));
	by occ_2002;
	if a;
	if occ_2002=. then occ_prest_2002=.;
run;

proc sort data=occ;
	by occ_2004;
run;
 
data occ;
	merge occ(in=a) occp00(keep=occ_2000 occ_prest_2000 rename=(occ_2000=occ_2004 occ_prest_2000=occ_prest_2004));
	by occ_2004;
	if a;
	if occ_2004=. then occ_prest_2004=.;
run;

proc sort data=occ;
	by occ_2006;
run;
 
data occ;
	merge occ(in=a) occp00(keep=occ_2000 occ_prest_2000 rename=(occ_2000=occ_2006 occ_prest_2000=occ_prest_2006));
	by occ_2006;
	if a;
	if occ_2006=. then occ_prest_2006=.;
run;

proc sort data=occ;
	by occ_2008;
run;
 
data occ;
	merge occ(in=a) occp00(keep=occ_2000 occ_prest_2000 rename=(occ_2000=occ_2008 occ_prest_2000=occ_prest_2008));
	by occ_2008;
	if a;
	if occ_2008=. then occ_prest_2008=.;
run;

proc sort data=occ;
	by occ_2010;
run;
 
data occ;
	merge occ(in=a) occp10(keep=occ_2010 occ_prest_2010);
	by occ_2010;
	if a;
run;

proc sort data=occ;
	by occ_2012;
run;
 
data occ;
	merge occ(in=a) occp10(keep=occ_2010 occ_prest_2010 rename=(occ_2010=occ_2012 occ_prest_2010=occ_prest_2012));
	by occ_2012;
	if a;
run;

proc sort data=occ;
	by occ_2014;
run;
 
data occ;
	merge occ(in=a) occp10(keep=occ_2010 occ_prest_2010 rename=(occ_2010=occ_2014 occ_prest_2010=occ_prest_2014));
	by occ_2014;
	if a;
run;

proc sort data=occ;
	by occ_2016;
run;
 
data occ;
	merge occ(in=a) occp10(keep=occ_2010 occ_prest_2010 rename=(occ_2010=occ_2016 occ_prest_2010=occ_prest_2016));
	by occ_2016;
	if a;

	occ_prest=mean(of occ_prest_1992--occ_prest_2016);
run;

proc standard data=occ mean=0 std=1 out=occz;
	var occ_prest;
run;

data occz;
	set occz;
	occ_prest=-1*occ_prest;
run;

data data.HRS_occ_prestige;
	set occz;
	keep hhid pn occ_1992--occ_prest;
run;

proc sort;
	by hhid pn;
run;
*/









/*combine datasets and create variables*/
data phenos;
	merge info(in=a)
		track.trk2018tr_r(keep=hhid pn page firstiw) 
		library.Randhrs1992_2020v1(keep=hhid pn raeduc raedyrs rameduc rafeduc ragender raracem rahispan 
			R13VGACTX R13MDACTX R13LTACTX R13DRINK R13DRINKD R13DRINKN R13SMOKEN R13SMOKEV r13bmi 
			h1atotb h2atotb h4atotb h5atotb h6atotb h7atotb h8atotb h9atotb h10atotb h11atotb h12atotb h13atotb 
			h1itot h2itot h3itot h4itot h5itot h6itot h7itot h8itot h9itot h10itot h11itot h12itot h13itot h13amort 
			h13amrtb h13atran h13astck h13aira h13abond h13acd r14iwstat r13conde r14conde r12bpsys r13bpsys r14bpsys 
			r12bpdia r13bpdia r14bpdia r13shlt r14shlt r14shltc r13diabe r14diabe r13hibpe r14hibpe r13hearte r14hearte 
			r13stroke r14stroke r13cogtot)
		hrs98.h98a_r(keep=hhid pn f993--f996)
		hrs00.h00a_r(keep=hhid pn g1080--g1083)
		hrs02.h02b_r(keep=hhid pn hb020--hb023)
		hrs04.h04b_r(keep=hhid pn jb020--jb023)
		hrs06.h06b_r(keep=hhid pn kb020--kb023)
		hrs08.h08b_r(keep=hhid pn lb020--lb023)
		hrs10.h10b_r(keep=hhid pn mb020--mb023)
		hrs12.h12b_r(keep=hhid pn nb020--nb023)
		hrs14.h14b_r(keep=hhid pn ob020--ob023)
		hrs16.h16b_r(keep=hhid pn pb020--pb023)
		hrs16.h16c_r(keep=hhid pn pc116 pc117)
		hrs16.h16pr_r(keep=hhid pn pz205)
		census
		socenv
		data.hrs_occ_prestige(keep=hhid pn occ_prest)
		edimp.Hrs_cog_varmerge_all_may18(keep=hhid pn rameduc_v2 rafeduc_v2)
		cses.Cses_measures
		hcap.Hc16hp_f(keep=hhid pn r1hcapdx)
		apoe.Apoe_serotonin_release
		langa.Cogfinalimp_9520wide(keep=hhid pn cogtot27_imp2016 cogfunction2016 cogtot27_imp2018 cogfunction2018);
	by hhid pn;
	if a;

	if RAEDYRS=.M then RAEDYRS=.;

	if raeduc=5 then COLLEDUC=0;
	else if raeduc>0 then COLLEDUC=1;

	if raeduc=1 then HSEDUC=1;
	else if raeduc>1 then HSEDUC=0;

	if 0<=raedyrs<12 then educ=5;
	else if raedyrs=12 then educ=4;
	else if 12<raedyrs<16 then educ=3;
	else if raedyrs=16 then educ=2;
	else if raedyrs>16 then educ=1;

	raedyrsr = -raedyrs;

	if RAHISPAN=1 then race=3;
	else if RARACEM=1 then race=1;
	else if RARACEM=2 then race=2;
	else if RARACEM=3 then race=4;

	if rameduc_v2 in (7.5 8.5) then rameduc2=round(myrs);
	else rameduc2=rameduc_v2;

	if rafeduc_v2 in (7.5 8.5) then rafeduc2=round(fyrs);
	else rafeduc2=rafeduc_v2;

	peducyrs=max(rameduc2, rafeduc2);

	if 0<=peducyrs<12 then peduc=0;
	else if 12<=peducyrs<16 then peduc=1;
	else if peducyrs>=16 then peduc=2;

	if 0<=peducyrs<8 then peduc2=1;
	else if peducyrs>=8 then peduc2=0;

	if 0<=rameduc2<8 then meduc2=1;
	else if rameduc2>=8 then meduc2=0;

	if 0<=rafeduc2<8 then feduc2=1;
	else if rafeduc2>=8 then feduc2=0;

	if 0<=rameduc2<8 then meducq=4;
	else if 8<=rameduc2<12 then meducq=3;
	else if rameduc2=12 then meducq=2;
	else if 12<rameduc2<16 then meducq=1;
	else if rameduc2>=16 then meducq=0;

	if 0<=rafeduc2<8 then feducq=4;
	else if 8<=rafeduc2<12 then feducq=3;
	else if rafeduc2=12 then feducq=2;
	else if 12<rafeduc2<16 then feducq=1;
	else if rafeduc2>=16 then feducq=0;

*	if peducyrs in (7.5 8.5) then peducyrs=.;
*	if rameduc_v2=8.5 then meducq=.;
*	if rafeduc_v2=8.5 then feducq=.;
*	if rameduc_v2 in (7.5 8.5 .D) then rameduc_v2=.;
*	if rafeduc_v2 in (7.5 8.5 .D) then rafeduc_v2=.;
	if rameduc2=.D then rameduc2=.;
	if rafeduc2=.D then rafeduc2=.;
	
	peducq=mean(meducq, feducq);

	if pz205=5 or pc116=5 then smoke=0;
	else if (pz205=1 or pc116=1) and pc117=5 then smoke=1;
	else if (pz205=1 or pc116=1) and pc117=1 then smoke=2;

	if r13drink=0 then drinkspw=0;
	else if r13drink=1 then drinkspw=r13drinkd*r13drinkn;

	if drinkspw=0 then drinkscat=0;
	else if 0<drinkspw<=14 then drinkscat=1;
	else if drinkspw>14 then drinkscat=2;

	hhinc = mean(1.72*h1itot, 1.62*h2itot, 1.53*h3itot, 1.47*h4itot, 1.40*h5itot, 1.34*h6itot, 1.28*h7itot, 1.19*h8itot, 1.12*h9itot, 1.09*h10itot, 1.05*h11itot, 1.01*h12itot, h13itot);
	hhincr = -hhinc;
	zhhincr = (82496.92-hhinc)/108216.53;
	lnhhinc=log(hhinc+1);
	if lnhhinc=0 then lnhhincx=.;
	else lnhhincx=lnhhinc;
	lnhhincxr=-lnhhincx;

	if 0<=hhinc<27051.83 then hhincq=5;
	else if 27051.83<=hhinc<48327.82 then hhincq=4;
	else if 48327.82<=hhinc<74297.45 then hhincq=3;
	else if 74297.45<=hhinc<117465.11 then hhincq=2;
	else if hhinc>=117465.11 then hhincq=1;


	if h13amort>0 or h13amrtb>0 then homelandown=1;
	else homelandown=0;

	if h13atran>0 then carown=1;
	else carown=0;

	if h13astck>0 or h13aira>0 or h13abond>0 or h13acd>0 then invest=1;
	else invest=0;

	if homelandown=1 and carown=1 and invest=1 then wealth=1;
	else if homelandown=1 and (carown=1 or invest=1) then wealth=2;
	else if homelandown=1 and carown=0 and invest=0 then wealth=3;
	else if homelandown=0 and carown=1 and invest=1 then wealth=3;
	else if homelandown=0 and (carown=1 or invest=1) then wealth=4;
	else if homelandown=0 and carown=0 and invest=0 then wealth=5;

	wealth2 = mean(1.72*h1atotb, 1.62*h2atotb, 1.47*h4atotb, 1.40*h5atotb, 1.34*h6atotb, 1.28*h7atotb, 1.19*h8atotb, 1.12*h9atotb, 1.09*h10atotb, 1.05*h11atotb, 1.01*h12atotb, h13atotb);
	wealth2r = -wealth2;
	zwealth2r = (456760.70-wealth2)/954464.16;
	lnwealth2=log(wealth2+663710);
	if lnwealth2=0 or lnwealth2>16 then lnwealth2x=.;
	else lnwealth2x=lnwealth2;
	lnwealth2xr=-lnwealth2x;

	if -663709.00<=wealth2<27731.57 then wealth2q=5;
	else if 27731.57<=wealth2<118929.21 then wealth2q=4;
	else if 118929.21<=wealth2<286247.41 then wealth2q=3;
	else if 286247.41<=wealth2<649582.86 then wealth2q=2;
	else if wealth2>=649582.86 then wealth2q=1;

	if -5<F1_PC2_1216<=-0.7990601 then F1_PC2q=1;
	else if -0.7990601<F1_PC2_1216<=-0.2121291 then F1_PC2q=2;
	else if -0.2121291<F1_PC2_1216<=0.2142798 then F1_PC2q=3;
	else if 0.2142798<F1_PC2_1216<=0.6017173 then F1_PC2q=4;
	else if F1_PC2_1216>0.6017173 then F1_PC2q=5;

	if -12<factor_ana_1216<=-4.5189606 then factor_anaq=5;
	else if -4.5189606<factor_ana_1216<=-1.9455142 then factor_anaq=4;
	else if -1.9455142<factor_ana_1216<=0.3971477 then factor_anaq=3;
	else if 0.3971477<factor_ana_1216<=3.4494588 then factor_anaq=2;
	else if factor_ana_1216>3.4494588 then factor_anaq=1;

	factor_anar = -factor_ana_1216;

	if -6<socenv<=-2.1986680 then socenvq=1;
	else if -2.1986680<socenv<=-0.9538195 then socenvq=2;
	else if -0.9538195<socenv<=0.2688697 then socenvq=3;
	else if 0.2688697<socenv<=1.8231023 then socenvq=4;
	else if socenv>1.8231023 then socenvq=5;

	if -4<occ_prest<=-0.8206871 then occ_prestq=1;
	else if -0.8206871<occ_prest<=-0.2626447 then occ_prestq=2;
	else if -0.2626447<occ_prest<=0.3113287 then occ_prestq=3;
	else if 0.3113287<occ_prest<=0.7744443 then occ_prestq=4;
	else if occ_prest>0.7744443 then occ_prestq=5;

	SESn=n(educ, hhincq, wealth2q, factor_anaq, socenvq, occ_prestq);

	if SESn>=3 then SESindex=sum(educ, hhincq, wealth2q, factor_anaq, socenvq, occ_prestq)/SESn;

	cfsi020=max(f993, g1080, hb020, jb020, kb020, lb020, mb020, nb020, ob020, pb020);
	cfsi021=max(f994, g1081, hb021, jb021, kb021, lb021, mb021, nb021, ob021, pb021);
	cfsi022=max(f995, g1082, hb022, jb022, kb022, lb022, mb022, nb022, ob022, pb022);
	cfsi023=max(f996, g1083, hb023, jb023, kb023, lb023, mb023, nb023, ob023, pb023);

	if cfsi020 in (5 6) then cfsi1=1;
	else if cfsi020 in (1 3) then cfsi1=0;

	if cfsi021=1 then cfsi2=1;
	else if cfsi021=5 then cfsi2=0;

	if cfsi022=1 then cfsi3=1;
	else if cfsi022=5 then cfsi3=0;

	if cfsi023 in (1 6 7) then cfsi4=1;
	else if cfsi023=5 then cfsi4=0;

	cfsi=sum(of cfsi1-cfsi4);

	cSES=mean(peducq, cfsi);

	if r14iwstat=5 then dead18=1;
	else dead18=0;

	conde_change=r14conde-r13conde;

	bpsys_change=r14bpsys-r12bpsys;
	bpdia_change=r14bpdia-r12bpdia;

	if r1hcapdx in (2 3) then hcap_MCI=1;
	else if r1hcapdx=1 then hcap_MCI=0;

	if r1hcapdx=3 then hcap_dementia=1;
	else if r1hcapdx in (1 2) then hcap_dementia=0;

	if cogfunction2016 in (2 3) then langaweir16_CIND=1;
	else if cogfunction2016=1 then langaweir16_CIND=0;

	if cogfunction2016=3 then langaweir16_dementia=1;
	else if cogfunction2016 in (1 2) then langaweir16_dementia=0;

	if cogfunction2018 in (2 3) then langaweir18_CIND=1;
	else if cogfunction2018=1 then langaweir18_CIND=0;

	if cogfunction2018=3 then langaweir18_dementia=1;
	else if cogfunction2018 in (1 2) then langaweir18_dementia=0;

	if apoe in (24 34 44) then apoee4=1;
	else if apoe in (22 23 33) then apoee4=0;

	format _ALL_;
run;

/*
proc means N min p20 p40 p60 p80 max mean std;
	var hhinc wealth2;
run;
*/


/*standardize and take means of conditions to create MSI variables*/
proc stdize data=phenos out=phenos2 oprefix sprefix=z;
	var r13hibpe r14hibpe r13diabe r14diabe r13hearte r14hearte r13stroke r14stroke;
run;

data phenos2;
	set phenos2;

	MSI16=mean(zr13hibpe, zr13diabe, zr13hearte, zr13stroke);
	MSI18=mean(zr14hibpe, zr14diabe, zr14hearte, zr14stroke);
run;

proc stdize data=phenos2 out=phenos3 oprefix sprefix=z;
	var MSI16 MSI18;
run;

data phenos3;
	set phenos3;

	MSI_change=zMSI18-zMSI16;
run;



/*
proc univariate data=phenos3;
	var lnhhincx lnwealth2x;
	histogram;
run;
*/


/*factor scores*/
proc factor data=phenos3 nfactors=1 score out=phenos3 prefix=SESfactor;
	var raedyrsr lnhhincxr lnwealth2xr factor_anar socenv occ_prest;
run;

proc corr data=phenos3;
	var sesindex sesfactor1;
run;






data data.HRS_phenos;
	set phenos3;
run;


proc freq data=phenos4;
	table peducq;
run;





/*
proc factor data=phenos method=ml nfactors=1 out=factors;
	var cfsi1-cfsi4;
run;

proc means data=phenos;
	var cfsi cses_index;
run;

proc corr;
	var cfsi cses_index;
run;

ods pdf;
proc corr;
	var cfsi1--cfsi4 rameduc_v2 rafeduc_v2;
run;
ods pdf close;


proc calis method=fiml data=phenos outstat=test;
   factor
      fin ===> cfsi1-cfsi4;
run;

proc score data=factors score=test out=scores;
	var cfsi1-cfsi4;
run;
*/


/*summary tables*/
ods pdf;
proc freq;
	table rameduc_v2 rafeduc_v2 peducyrs peduc2 meduc2 feduc2 sesindex sesn cfsi cfsi1-cfsi4;
run;
ods pdf close;


proc means data=phenos N min p20 p40 p60 p80 max mean std;
	var factor_ana_1216 F1_PC2_1216 socenv occ_prest rameduc_v2 rafeduc_v2;
run;


proc univariate data=phenos noprint;
	var factor_ana_1216 F1_PC2_1216 socenv;
	histogram;
run;


proc freq;
	table meducq feducq peducq cfsi cses rameduc_v2 rafeduc_v2;
run;

proc freq;
	table r14iwstat dead18;
run;






/*save phenotype data*/
proc export
	data=phenos3
	outfile='S:\MESA\Methylation\Lauren_DNAmAge\data\HRS_methyset_phenos.csv'
	dbms=csv
	replace;
run;



data phenos3;
	set data.HRS_phenos;
run;






/*import and combine all clock data*/
proc import
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\results\HRS_Enet\HRS_SESindex_mclock.csv'
	out=mclock1(rename=(hhid=hhidx pn=pnx))
	dbms=csv
	replace;
run;

proc sort;
	by hhidx pnx;
run;

proc import
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\results\HRS_Enet\HRS_cSES_mclock.csv'
	out=mclock2(rename=(hhid=hhidx pn=pnx))
	dbms=csv
	replace;
run;

proc sort;
	by hhidx pnx;
run;

proc import
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\results\HRS_Enet\HRS_SESindex_mclock_450K.csv'
	out=mclock3(rename=(hhid=hhidx pn=pnx SESindex_mclock=SESindex_mclock450))
	dbms=csv
	replace;
run;

proc sort;
	by hhidx pnx;
run;

proc import
	datafile='S:\MESA\Methylation\Lauren_DNAmAge\results\HRS_Enet\HRS_cSES_mclock_450K.csv'
	out=mclock4(rename=(hhid=hhidx pn=pnx cSES_mclock=cSES_mclock450))
	dbms=csv
	replace;
run;

proc sort;
	by hhidx pnx;
run;

data mclocks;
	length hhid $6. PN $3.;
	merge mclock1 mclock2 mclock3 mclock4;
	by hhidx pnx;

	PN=cat('0', put(pnx,2.));

	if length(strip(hhidx))=5 then hhid=cat('0', put(hhidx,5.));
	else if length(strip(hhidx))=6 then hhid=put(hhidx,6.);

	drop pnx hhidx;
run;

proc stdize data=mclocks out=mclocks2 oprefix sprefix=z;
	var SESindex_mclock cSES_mclock SESindex_mclock450 cSES_mclock450;
run;


proc sort data=clocks.grimage_components out=grimcomps;
	by hhid pn;
run;



proc import
	datafile='S:\MESA\Methylation\ELA_DNAmAge\Data\HRS_DunedinPACE_new.csv'
	out=pace
	dbms=csv
	replace;
run;

proc sort;
	by sample;
run;





/*final phenotype data with clocks*/
data phenos4;
	merge data.hrs_phenos(in=a) clocks.epiclocka_r grimcomps(keep=hhid pn DNAmGDF_15--DNAmPACKYRS) mclocks2;
	by hhid pn;
	if a;
run;

proc sort;
	by sample;
run;

data phenos4;
	merge phenos4 pace;
	by sample;
run;

proc sort;
	by hhid pn;
run;

proc stdize data=phenos4 out=phenos4 oprefix sprefix=z;
	var SESindex cSES;
run;






data test;
	set phenos4;
	myrs_round=round(myrs);
	fyrs_round=round(fyrs);
run;

proc freq;
	where rafeduc in (7.5 8.5);
	table fyrs_round / missing;
run;

proc corr;
	where rameduc in (.D .M);
	var rameduc_v2 myrs_round;
run;

proc freq;
	where zsesindex_mclock^=.;
	table r1hcapdx r13cogtot ;*/ missing;
run;



option orientation=landscape;
ods csv;
proc corr data=phenos4 nomiss  nosimple;
	var zsesindex_mclock zcses_mclock zsesindex_mclock450 zcses_mclock450 sesindex cses HORVATH_DNAMAGE--HORVATHSKIN_DNAMAGE DNAMGRIMAGE DunedinPACE;
*	with zsesindex_mclock zcses_mclock sesindex cses HORVATH_DNAMAGE--HORVATHSKIN_DNAMAGE DNAMGRIMAGE DunedinPACE;
run;
ods csv close;

ods csv;
proc corr nomiss noprob nosimple;
	var zsesindex_mclock zcses_mclock sesindex cses DNAMGRIMAGE DNAmGDF_15--DNAmPACKYRS;
run;
ods csv close;

proc freq data=phenos4;
	table conde_change r14shltc;
run;

proc freq data=phenos4;
	where zsesindex_mclock^=.;
	table dead18*(ragender race smoke) / chisq;
run;

proc ttest data=phenos4;
	where zsesindex_mclock^=.;
	class dead18;
	var page sesindex cses zsesindex_mclock zcses_mclock;
run;

proc corr;
	var conde_change r14shltc;
run;

proc freq data=phenos4;
	where zsesindex_mclock^=.;
	table educ hhincq wealth2q factor_anaq socenvq occ_prestq;
run;
	





data test;
	set phenos4;
	sample_id=input(local_id, 6.);
run;

proc sort;
	by sample_id;
run;

proc import
	datafile='R:\Health_Retirement_Study\Genotype_Jan2015_Processing\PCA\Phase1234\Top_PCs.csv'
	out=pcs
	dbms=csv
	replace;
run;

proc sort data=pcs;
	by sample_id;
run;

data test;
	merge test(in=a) pcs(in=poo);
	by sample_id;
	if a and poo;
run;

proc sort;
	by hhid pn;
run;

data test2;
	merge test(in=a) track.trk2018tr_r(keep=hhid pn qiwwave qiwmode qiwtype qinsamp);
	by hhid pn;
	if a and sesindex_mclock^=.;
run;

proc freq;
	table qiwwave qiwmode qiwtype qinsamp r14hibpe r14diabe / missing;
run;

data phenos4;
	merge phenos4 test2(keep=hhid pn qiwtype);
	by hhid pn;

	if qiwtype=1 then in2018=1;
	else in2018=0;
run;

proc freq;
	table in2018;
run;


proc freq;
	table r1hcapdx cogfunction2016;
run;








/*This macro will do all health outcomes at once.  Can do them individually instead in next section of code*/
%macro sesmodels(pred, ses, model) / minoperator;
%let out=r14conde conde_change zMSI18 MSI_change r14shlt r14shltc dead18 langaweir18_CIND langaweir18_dementia;
%do i=1 %to 9;
	%let outi=%scan(&out, &i);
	%if &i in (3 4 5 6) %then %do;
		proc glm data=phenos4 plots=diagnostics;
			class race ragender smoke / ref=first;
			model &outi=&pred page race ragender smoke &ses / clparm solution;
			ods output ParameterEstimates=estimates;
		run;

		data estimates&i;
			set estimates;
			if parameter="&pred";
			drop biased--tvalue;
			rename estimate=estimate&model probt=probt&model lowercl=lowercl&model uppercl=uppercl&model;
		run;
	%end;	

	%else %if &i in (1 2) %then %do;
		proc genmod data=phenos4;
			class race ragender smoke / ref=first param=ref;
			model &outi=&pred page race ragender smoke &ses / dist=poisson;
			ods output ParameterEstimates=estimates;
		run;

		data estimates&i;
			retain dependent parameter IRR lowercl uppercl;
			set estimates;
			if parameter="&pred";
			dependent="&outi";
			IRR=exp(estimate);
			lowercl=exp(lowerwaldcl);
			uppercl=exp(upperwaldcl);
			rename IRR=estimate&model probchisq=probt&model lowercl=lowercl&model uppercl=uppercl&model;
			keep dependent parameter IRR probchisq lowercl uppercl;
		run;
	%end;

	%else %if &i=7 %then %do;
		proc logistic data=phenos4 desc;
			class race ragender smoke / ref=first param=ref;
			model &outi=&pred page race ragender smoke &ses / expb rsq;
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

	%else %if &i in (8 9) %then %do;
		proc logistic data=phenos4 desc;
			class race ragender smoke apoee4 / ref=first param=ref;
			model &outi=&pred page race ragender smoke apoee4 &ses / expb rsq;
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
			merge estimates ors;
			if variable="&pred";
			dependent="&outi";
			rename expest=estimate&model probchisq=probt&model variable=parameter;
			keep dependent variable expest probchisq lowercl uppercl;
		run;

		data estimates&i;
			merge estimates&i ors;
		run;
	%end;

	data estimatesall&model;
		set estimates1 /*estimates2*/ estimates3 /*estimates4*/ estimates5 /*estimates6*/ estimates7 /*estimates8*/ estimates9;
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
%sesmodels(zSESindex_mclock,cSES,5)


/*run appropriate models above before this step*/
data estimatesall_SES;
	merge estimatesall1star estimatesall1 estimatesall2 estimatesall3a estimatesall4a estimatesall3b estimatesall4b estimatesall5;
run;

%sesmodels(zcSES_mclock,,1)
%sesmodels(zcSES,,1star)
%sesmodels(zcSES_mclock,cses,2)
%sesmodels(zcSES_mclock,dnamgrimage,3a)
%sesmodels(zcSES_mclock,cses dnamgrimage,4a)
%sesmodels(zcSES_mclock,dunedinpace,3b)
%sesmodels(zcSES_mclock,cses dunedinpace,4b)
%sesmodels(zcSES_mclock,SESindex,5)

/*run appropriate models above before this step*/
data estimatesall_cSES;
	merge estimatesall1star estimatesall1 estimatesall2 estimatesall3a estimatesall4a estimatesall3b estimatesall4b estimatesall5;
run;

%sesmodels(zSESindex_mclock450,,1)
%sesmodels(zSESindex,,1star)
%sesmodels(zSESindex_mclock450,SESindex,2)
%sesmodels(zSESindex_mclock450,dnamgrimage,3a)
%sesmodels(zSESindex_mclock450,SESindex dnamgrimage,4a)
%sesmodels(zSESindex_mclock450,dunedinpace,3b)
%sesmodels(zSESindex_mclock450,SESindex dunedinpace,4b)
%sesmodels(zSESindex_mclock450,cses,5)

/*run appropriate models above before this step*/
data estimatesall_SES450;
	merge estimatesall1star estimatesall1 estimatesall2 estimatesall3a estimatesall4a estimatesall3b estimatesall4b estimatesall5;
run;

%sesmodels(zcSES_mclock450,,1)
%sesmodels(zcSES,,1star)
%sesmodels(zcSES_mclock450,cses,2)
%sesmodels(zcSES_mclock450,dnamgrimage,3a)
%sesmodels(zcSES_mclock450,cses dnamgrimage,4a)
%sesmodels(zcSES_mclock450,dunedinpace,3b)
%sesmodels(zcSES_mclock450,cses dunedinpace,4b)
%sesmodels(zcSES_mclock450,SESindex,5)

/*run appropriate models above before this step*/
data estimatesall_cSES450;
	merge estimatesall1star estimatesall1 estimatesall2 estimatesall3a estimatesall4a estimatesall3b estimatesall4b estimatesall5;
run;



ods csv;
proc print data=estimatesall_SES noobs;
run;
ods csv close;

ods csv;
proc print data=estimatesall_cSES noobs;
run;
ods csv close;

ods csv;
proc print data=estimatesall_SES450 noobs;
run;
ods csv close;

ods csv;
proc print data=estimatesall_cSES450 noobs;
run;
ods csv close;








/*individual models for aSES mclock*/

%let sesvar=;
%let sesvar=raedyrs;
%let sesvar=factor_ana_1216;
%let sesvar=socenv;
%let sesvar=SESindex;
%let sesvar=dnamgrimage;
%let sesvar=SESindex dnamgrimage;
%let sesvar=dunedinpace;
%let sesvar=SESindex dunedinpace;


ods pdf;
/*
proc glm data=phenos4 plots=diagnotics;
	class race ragender / ref=first;
	model bpsys_change=zSESindex_mclock page race ragender / solution;
run;

proc glm data=phenos4 plots=diagnotics;
	class race ragender / ref=first;
	model bpdia_change=zSESindex_mclock page race ragender / solution;
run;
*/
proc genmod data=phenos4;
	class race ragender smoke / ref=first param=ref;
	model r14conde=zSESindex_mclock page race ragender smoke &sesvar / dist=poisson;
	ods output ParameterEstimates=estimates;
run;

proc genmod data=phenos4;
	class race ragender smoke / ref=first param=ref;
	model conde_change=zSESindex_mclock page race ragender smoke &sesvar / dist=poisson;
run;

proc logistic data=phenos4 desc;
	class race ragender smoke / ref=first param=ref;
	model dead18=zSESindex_mclock page race ragender smoke &sesvar /expb rsq;
	ods output ParameterEstimates=estimates;
run;

proc GENMOD data=phenos4 desc;
	class race ragender smoke / ref=first param=ref;
	model dead18=zSESindex_mclock page race ragender smoke &sesvar / dist=bin;
	*ods output ParameterEstimates=estimates;
run;

proc glm data=phenos4 plots=diagnotics;
	class race ragender smoke / ref=first;
	model r14shlt=zSESindex_mclock page race ragender smoke &sesvar / clparm solution;
	ods output ParameterEstimates=estimates;
run;

proc glm data=phenos4 plots=diagnotics;
	class race ragender smoke / ref=first;
	model r14shltc=zSESindex_mclock page race ragender smoke &sesvar / solution;
run;

proc glm data=phenos4 plots=diagnotics;
	class race ragender smoke / ref=first;
	model zMSI18=zSESindex_mclock page race ragender smoke &sesvar / solution;
run;
quit;

proc glm data=phenos4 plots=diagnotics;
	class race ragender smoke / ref=first;
	model MSI_change=zSESindex_mclock page race ragender smoke &sesvar / solution;
run;
quit;

proc genmod data=phenos4 desc;
	class race ragender smoke apoe(ref='33') / ref=first param=ref;
	model hcap_mci=zSESindex_mclock page race ragender smoke apoe &sesvar / dist=bin;
run;

proc genmod data=phenos4 desc;
	class race ragender smoke apoe(ref='33') / ref=first param=ref;
	model hcap_dementia=zSESindex_mclock page race ragender smoke apoe &sesvar / dist=bin;
run;
ods pdf close;


/*MODEL 1*s  */
ods pdf;
proc genmod data=phenos4;
	class race ragender smoke / ref=first param=ref;
	model r14conde=SESindex page race ragender smoke / dist=poisson;
run;

proc glm data=phenos4 plots=diagnotics;
	class race ragender smoke / ref=first;
	model zMSI18=SESindex page race ragender smoke / solution;
run;
quit;

proc glm data=phenos4 plots=diagnotics;
	class race ragender smoke / ref=first;
	model r14shlt=SESindex page race ragender smoke / solution;
run;

proc logistic data=phenos4 desc;
	class race ragender smoke / ref=first param=ref;
	model dead18=SESindex page race ragender smoke / expb rsq;
	ods output ParameterEstimates=estimates oddsratios=ors;
run;

proc logistic data=phenos4 desc;
	class race ragender smoke apoe(ref='33') / ref=first param=ref;
	model langaweir18_dementia=SESindex page race ragender smoke apoe / rsq;
run;
ods pdf close;





/*individual models for cSES mclock*/

%let sesvar=;
%let sesvar=peducyrs;
%let sesvar=cfsi;
%let sesvar=cses;
%let sesvar=dnamgrimage;
%let sesvar=cses dnamgrimage;
%let sesvar=dunedinpace;
%let sesvar=cses dunedinpace;

ods pdf;
/*
proc glm data=phenos4 plots=diagnotics;
	class race ragender / ref=first;
	model bpsys_change=zcSES_mclock page race ragender / solution;
run;

proc glm data=phenos4 plots=diagnotics;
	class race ragender / ref=first;
	model bpdia_change=zcSES_mclock page race ragender / solution;
run;
*/
proc genmod data=phenos4;
	class race ragender smoke / ref=first param=ref;
	model r14conde=zcSES_mclock page race ragender smoke &sesvar / dist=poisson;
run;

proc genmod data=phenos4;
	class race ragender smoke / ref=first param=ref;
	model conde_change=zcSES_mclock page race ragender smoke &sesvar / dist=poisson;
run;

proc logistic data=phenos4 desc;
	class race ragender smoke / ref=first param=ref;
	model dead18=zcSES_mclock page race ragender smoke &sesvar / rsq;
run;

proc glm data=phenos4 plots=diagnotics;
	class race ragender smoke / ref=first;
	model r14shlt=zcSES_mclock page race ragender smoke &sesvar / solution;
run;

proc glm data=phenos4 plots=diagnotics;
	class race ragender smoke / ref=first;
	model r14shltc=zcSES_mclock page race ragender smoke &sesvar / solution;
run;

proc glm data=phenos4 plots=diagnotics;
	class race ragender smoke / ref=first;
	model zMSI18=zcSES_mclock page race ragender smoke &sesvar / solution;
run;
quit;

proc glm data=phenos4 plots=diagnotics;
	class race ragender smoke / ref=first;
	model MSI_change=zcSES_mclock page race ragender smoke &sesvar / solution;
run;
quit;
ods pdf close;


/*MODEL 1*s  */
ods pdf;
proc genmod data=phenos4;
	class race ragender smoke / ref=first param=ref;
	model r14conde=cses page race ragender smoke / dist=poisson;
run;

proc glm data=phenos4 plots=diagnotics;
	class race ragender smoke / ref=first;
	model zMSI18=cses page race ragender smoke / solution;
run;
quit;

proc glm data=phenos4 plots=diagnotics;
	class race ragender smoke / ref=first;
	model r14shlt=cses page race ragender smoke / solution;
run;

proc logistic data=phenos4 desc;
	class race ragender smoke / ref=first param=ref;
	model dead18=cses page race ragender smoke;
run;

proc logistic data=phenos4 desc;
	class race ragender smoke apoe(ref='33') / ref=first param=ref;
	model langaweir18_dementia=cses page race ragender smoke apoe;
run;
ods pdf close;







/*individual models for aSES mclock-450K*/

%let sesvar=;
%let sesvar=raedyrs;
%let sesvar=factor_ana_1216;
%let sesvar=socenv;
%let sesvar=SESindex;
%let sesvar=dnamgrimage;
%let sesvar=SESindex dnamgrimage;

proc genmod data=phenos4;
	class race ragender smoke / ref=first param=ref;
	model r14conde=zSESindex_mclock450 page race ragender smoke &sesvar / dist=poisson;
run;

proc genmod data=phenos4;
	class race ragender smoke / ref=first param=ref;
	model conde_change=zSESindex_mclock450 page race ragender smoke &sesvar / dist=poisson;
run;

proc logistic data=phenos4 desc;
	class race ragender smoke / ref=first param=ref;
	model dead18=zSESindex_mclock450 page race ragender smoke &sesvar / rsq;
run;

proc glm data=phenos4 plots=diagnotics;
	class race ragender smoke / ref=first;
	model r14shlt=zSESindex_mclock450 page race ragender smoke &sesvar / solution;
run;

proc glm data=phenos4 plots=diagnotics;
	class race ragender smoke / ref=first;
	model r14shltc=zSESindex_mclock450 page race ragender smoke &sesvar / solution;
run;

proc glm data=phenos4 plots=diagnotics;
	class race ragender smoke / ref=first;
	model zMSI18=zSESindex_mclock450 page race ragender smoke &sesvar / solution;
run;
quit;

proc glm data=phenos4 plots=diagnotics;
	class race ragender smoke / ref=first;
	model MSI_change=zSESindex_mclock450 page race ragender smoke &sesvar / solution;
run;
quit;










/*macro for finding freqs of categorical variables for Table1*/
%macro freqs(var, strata);
proc freq data=phenos4;
	where zsesindex_mclock^=.;
	table &var / out=freqstot;
run;

proc freq data=phenos4;
	where zsesindex_mclock^=.;
	table &var*&strata / chisq out=freqs outpct;
	ods output chisq=chisq;
run;

data freqs2;
	set freqstot freqs;

	if &strata=. then pct=put(percent, 6.2);
	else pct=put(pct_col, 6.2);

	NPCT=cat(count, ' (', strip(pct), ')'); 
run;

proc sort data=freqs2;
	by &var;
run;

proc transpose data=freqs2 out=freqs&var(drop=_name_);
	var NPCT;
	by &var;
run;

data chisq;
	set chisq;
	if statistic='Chi-Square';
	keep prob;
run;

data freqs&var;
	retain Variable;
	length Category $20.;
	merge freqs&var chisq;
	Category=vvalue(&var);
	Variable="&var";
	drop &var;
run;
%mend;


%freqs(race, dead18)
%freqs(ragender, dead18)
%freqs(smoke, dead18)
%freqs(colleduc, dead18)
%freqs(hseduc, dead18)
%freqs(peduc2, dead18)
%freqs(feduc2, dead18)
%freqs(meduc2, dead18)

%freqs(race, in2018)
%freqs(ragender, in2018)
%freqs(smoke, in2018)
%freqs(dead18, in2018)
%freqs(educ, in2018)
%freqs(colleduc, in2018)
%freqs(hseduc, in2018)
%freqs(peduc2, in2018)
%freqs(feduc2, in2018)
%freqs(meduc2, in2018)
%freqs(langaweir18_dementia, in2018)
%freqs(hhincq, in2018)
%freqs(meducq, in2018)
%freqs(feducq, in2018)
%freqs(wealth2q, in2018)




/*macro for finding means of continuous variables for Table1*/
%macro means(var, strata);
proc means data=phenos4 mean std;
	where zsesindex_mclock^=.;
	class &strata;
	var &var;
	types () &strata;
	ods output summary=means;
run;

data means;
	set means;
	mean=put(&var._mean, 20.2);
	std=put(&var._stddev, 20.2);
	mstd=cat(strip(mean), ' (', strip(std), ')');
run;

proc transpose data=means out=meansx(drop=_label_ _name_);
	var mstd;
run;

proc ttest data=phenos4;
	where zsesindex_mclock^=.;
	class &strata;
	var &var;
	ods output ttests=ttest;
run;

data ttest;
	set ttest;
	if method='Pooled';
	rename probt=prob;
	format probt pvalue6.4;
	keep probt;
run;

data means&var;
	retain Variable;
	merge meansx ttest;
	Variable="&var";
run;
%mend;


%means(page, dead18)
%means(factor_ana_1216, dead18)
%means(F1_PC2_1216, dead18)
%means(socenv, dead18)
%means(raedyrs, dead18)
%means(occ_prest, dead18)
%means(sesindex, dead18)
%means(peducyrs, dead18)
%means(rafeduc_v2, dead18)
%means(rameduc_v2, dead18)
%means(cfsi, dead18)
%means(cses, dead18)
%means(zsesindex_mclock, dead18)
%means(zcses_mclock, dead18)

%means(page, in2018)
%means(factor_ana_1216, in2018)
%means(raedyrs, in2018)
%means(hhinc, in2018)
%means(wealth2, in2018)
%means(occ_prest, in2018)
%means(sesindex, in2018)
%means(socenv, in2018)
%means(peducyrs, in2018)
%means(rafeduc2, in2018)
%means(rameduc2, in2018)
%means(cfsi, in2018)
%means(cses, in2018)
%means(dnamgrimage, in2018)
%means(Dunedinpace, in2018)
%means(zsesindex_mclock, in2018)
%means(zcses_mclock, in2018)
%means(r14conde, in2018)
%means(conde_change, in2018)
%means(r14shlt, in2018)
%means(r14shltc, in2018)
%means(zMSI18, in2018)
%means(MSI_change, in2018)




/*constructing Table 1 of descriptors*/
data table1;
	retain variable category;
	length variable $12.;
	set meanspage freqsragender freqsrace freqssmoke meansfactor_ana_1216 meanssocenv meansraedyrs freqscolleduc freqshseduc meansocc_prest meanssesindex meanspeducyrs freqspeduc2 meansrameduc2 freqsmeduc2 meansrafeduc2 freqsfeduc2 meanscfsi meanscses meanszsesindex_mclock meanszcses_mclock;
run;


ods rtf;
proc print data=table1 noobs;
run;
ods rtf close;






data table1a;
	retain variable category;
	length variable $12.;
	set meanspage freqsragender freqsrace freqssmoke meansfactor_ana_1216 meanssocenv meansraedyrs freqscolleduc freqshseduc meansocc_prest meanssesindex meanspeducyrs freqspeduc2 meansrameduc2 freqsmeduc2 meansrafeduc2 freqsfeduc2 meanscfsi meanscses meansdnamgrimage meansMPOA meanszsesindex_mclock meanszcses_mclock meansr14conde meanszmsi18 meansr14shlt freqsdead18 /*meansconde_change meansmsi_change meansr14shltc*/;
run;


ods rtf;
proc print data=table1a noobs;
run;
ods rtf close;




