/*#################################################################*/
/*################## Primary analysis #############################*/
/*#################################################################*/

libname lib 'W:\SUD\Data\2022\SAS data';
libname der "W:\SUD\Data\2022\derived";
libname tmplib "E:\SAS Data\Missed dose protocol_v2";


option compress=Y;

/*#####################*/
/*#######  Box 1 #######*/
/*#####################*/

proc freq data=der.Epi_hosp_14;
 table drug_start;
run; /* total 293548 episode */


data test;
 set der.Epi_hosp_14;
 by studyid_rev;
 where dt_start >= '01JAN2010'd and dt_start < '31DEC2022'd and drug_start in ("B", "M");
run;/*JM: dt_start< '31dec2022'd? */

proc freq data=test;
 table drug_start;
run; /* MHM:  total 214132 episode (M=131733, B=82399) */ 

proc sql;
    select count(distinct studyid_rev) as num_individuals
    from test
    where drug_start = "M";
quit; /* 41173 unique individuals */


proc sql;
    select count(distinct studyid_rev) as num_individuals
    from test
    where drug_start = "B";
quit; /* 34651 unique individuals for BNX */



/*#####################*/
/*#######  Box 2 #######*/
/*#####################*/

data test1_b;
set tmplib.TT_analytic_alls_b;
by studyid_rev episode;
where dt_start >= '01JAN2010'd and dt_start < '31DEC2022'd and drug_start = "B" and predrug in ("B", "") and ddose_rev ne .;
if first.episode and trial_gap <= 14 and trial_gap >= 3;
run;

proc freq data = test1_b;
table trial_cat/missing;
run; /* 14404 episode in total */


proc freq data = test1_b;
where episode = 1;
table trial_cat/missing;
run; /* 3705 first episode in total */


proc sql;
    select count(distinct studyid_rev) as num_individuals
    from test1_b;
quit; /* 10205 unique individuals */

/* BNX  82399 - 14404 = 67995 episodes were excluded;
and 34651 - 10205 = 24446 individuals were excluded */

ods html close; ods html;


data test1_m;
set tmplib.TT_analytic_alls_m;
by studyid_rev episode;
where dt_start >= '01JAN2010'd and dt_start < '31DEC2022'd and drug_start = "M" and predrug in ("M", "") and ddose_rev ne .;
if first.episode and trial_gap <= 14 and trial_gap >= 3;;
run;

proc freq data = test1_m;
table trial_cat/missing;
run; /* 15804 episode total*/


proc sql;
    select count(distinct studyid_rev) as num_individuals
    from test1_m;
quit; /* 11034 unique individuals for MET */


/* MET  131733 - 15804 = 115929 episodes were excluded;
and 41173 - 11034 = 30139 individuals were excluded */


/*Total interruption in box 2*/
data data_m;
set tmplib.TT_analytic_alls_m;
by studyid_rev episode myepisode;
where dt_start >= '01JAN2010'd and dt_start < '31DEC2022'd and drug_start = "M" and predrug in ("M", "") and ddose_rev ne .;
if trial_gap <= 14 and trial_gap >= 3;
run;

data data_b;
set tmplib.TT_analytic_alls_b;
by studyid_rev episode myepisode;
where dt_start >= '01JAN2010'd and dt_start < '31DEC2022'd and drug_start = "B" and predrug in ("B", "") and ddose_rev ne .;
if trial_gap <= 14 and trial_gap >= 3;
run;

proc freq data=data_m; table trial_cat; run; /* For MET 149995 interruptions */
proc freq data=data_b; table trial_cat; run; /* For BNX 75952 interruptions */




/*#####################*/
/*#######  Box 3 #######*/
/*#####################*/

data test2_nw_b;
 set test1_b;
 if age_eli=1 and palli_excl=0 and canc_excl=0 and preg_eli=0 and incar_grp_bl ne '2_current';
run;

proc freq data = test2_nw_b;
table trial_cat;
run; /* 11626 episode of BNX have age >= 18, no preganancy, no incarceration at time0, and no cancer or palliative care history*/

/*# Unique individuals*/
proc sql;
    select count(distinct studyid_rev) as num_individuals
    from test2_nw_b;
quit; /* 8307 unique individuals for BNX*/

/*# first episode*/
proc sql;
    select count(distinct studyid_rev) as num_individuals
    from test2_nw_b
	where episode = 1;
quit; /*2989 first episode for BNX*/


/* BNX: 14404 - 11626 = 2778 episodes were excluded */
/*10205 - 8307 = 1898 individuals were excluded */



data test2_nw_m;
set test1_m;
 if age_eli=1 and palli_excl=0 and canc_excl=0 and preg_eli=0 and incar_grp_bl ne '2_current';
run;

proc freq data = test2_nw_m;
table trial_cat;
run; /* 14172 episode of BNX have age >= 18, no preganancy, no incarceration at time0, and no cancer or palliative care history: 
      */

proc sql;
    select count(distinct studyid_rev) as num_individuals
    from test2_nw_m;
quit; /* 10053 unique individuals for MET*/



proc freq data = test2_nw_m;
where episode = 1;
table trial_cat;
run; /*2396 first episode of MET*/



/* MET  15804 - 14172 = 1632 episodes were excluded */
/* 11034 - 10053 = 981 individuals were excluded */



/*##########################################*/
/*#######  Box 4: Episodes for Trial category  #######*/
/*##########################################*/

/* BNX */
proc freq data=test2_nw_b;
 table trial_cat/missing;
run; /* Interruption of 3-5 days: 5915  
		Interruption of 6-14 days: 5711
		Total  : 11626 episodes*/


/* MET */
proc freq data=test2_nw_m;
 table trial_cat /missing;
run; /* Interruption of 3-4 days: 4844  
		Interruption of 5-14 days: 9328
		Total  : 14172 episodes*/


/*##########################################*/
/*#######  Box 5: Episodes for Treatment  #######*/
/*##########################################*/


/* BNX */
proc freq data=test2_nw_b; 
where trial_gap in (3, 4, 5);
table trt_s5/missing; run;

	 /* No dose change: 3650
		<50% reduction: 595
       >=50% reduction: 596
		     increase : 1074
		Total  		  : 5915 */


/* BNX */
proc freq data=test2_nw_b; 
where trial_cat = "6-14";
table trt_s5/missing; run;

	 /* No dose change: 2957
		<50% reduction: 776
       >=50% reduction: 826
		     increase : 1152
		Total  		  : 5711 */





/* MET */
proc freq data=test2_nw_m; 
where trial_gap in (3, 4);
table trt_s456/missing; run;

	 /* No dose change: 2871
		<50% reduction: 836
       >=50% reduction: 213
		     increase : 924
		Total  		  : 4844 */


/* MET */
proc freq data=test2_nw_m; 
where trial_cat = "5-14";
table trt_s456/missing; run;

	 /* No dose change: 5113
		<50% reduction: 1763
       >=50% reduction: 560
		     increase : 1892
		Total  		  : 4844 */



