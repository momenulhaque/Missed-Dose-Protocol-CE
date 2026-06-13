/*#################################################################*/
/*################## Primary analysis #############################*/
/*#################################################################*/

libname lib 'W:\SUD\Data\2022\SAS data';
libname der "W:\SUD\Data\2022\derived";
libname tmplib "E:\SAS Data\Missed dose protocol";
libname tmplib1 "E:\SAS Data\Missed dose protocol_v1";

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
set tmplib1.TT_analytic_alls_b;
by studyid_rev episode;
where dt_start >= '01JAN2010'd and dt_start < '31DEC2022'd and drug_start = "B" and predrug in ("B", "") and ddose_rev ne .;
if first.episode and trial_gap <= 14;
run;/*JM: the codes above have condition that the first interruption within an episode should have <=14 gap.
If you have >14 gap at first and <=14 gap at second interruption within an episode,
any of those were not selected for that episode.
what about selecting the first interruption after getting those with a gap <14 days first:

data test1_b;
set tmplib.episodes5_b(where=(trial_gap <= 14));
by studyid_rev episode;
if first.episode;
run;
If you decide to apply the condition, please modify your other sas files accordingly*/


proc freq data = test1_b;
table trial_cat/missing;
run; /* MHM: 44848 episode of BNX has interruption <= 14 days */
/*JM: you may notice your sas is ruuning slow after this step (if the freq table has many rows).
You should clear your html output window clear by running this:
ods html close; ods html; */

proc sql;
    select count(distinct studyid_rev) as num_individuals
    from test1_b;
quit; /* 23661 unique individuals */

/* BNX  82399 - 44848 = 37551 episodes were excluded;
and 34651 - 23661 = 10990 individuals were excluded */

ods html close; ods html;


data test1_m;
set tmplib1.TT_analytic_alls_m;
by studyid_rev episode;
where dt_start >= '01JAN2010'd and dt_start < '31DEC2022'd and drug_start = "M" and predrug in ("M", "") and ddose_rev ne .;
if first.episode and trial_gap <= 14;
run;
/*JM: you can get the number of episodes just by looking at the number of observations in test1_m in the log files
(i.e. NOTE: The data set WORK.test1_m has x observations and x variables.)*/

proc freq data = test1_m;
table trial_cat/missing;
run; /* MHM: 93524 episode of MET has interruption <= 14 days: 
Total 44848+93524=138372; we excluded 214132-138372=75760
BNX 44848
*/


proc sql;
    select count(distinct studyid_rev) as num_individuals
    from test1_m;
quit; /* 34651 unique individuals for MET */


/* MET  131733 - 93524 = 38209 episodes were excluded;
and 41173 - 34651 = 6522 individuals were excluded */


/*Total interruption in box 2*/
data data_m;
set tmplib1.TT_analytic_alls_m;
by studyid_rev episode myepisode;
where dt_start >= '01JAN2010'd and dt_start < '31DEC2022'd and drug_start = "M" and predrug in ("M", "") and ddose_rev ne .;
if trial_gap <= 14;
run;

data data_b;
set tmplib1.TT_analytic_all_b;
by studyid_rev episode myepisode;
where dt_start >= '01JAN2010'd and dt_start < '31DEC2022'd and drug_start = "B" and predrug in ("B", "") and ddose_rev ne .;
if trial_gap <= 14;
run;

proc freq data=data_m; table trial_cat; run; /* For MET 1535886 interruptions */
proc freq data=data_b; table trial_cat; run; /* For BNX 329459 interruptions */




/*#####################*/
/*#######  Box 3 #######*/
/*#####################*/

data test2_nw_b;
 set test1_b;
 if age_eli=1 and palli_excl=0 and canc_excl=0 and preg_eli=0 and incar_grp_bl ne '2_current';
run;

proc freq data = test2_nw_b;
table trial_cat;
run; /* MHM: 37966 episode of BNX have age >= 18, no preganancy, no incarceration at time0, and no cancer or palliative care history*/

/*# Unique individuals*/
proc sql;
    select count(distinct studyid_rev) as num_individuals
    from test2_nw_b;
quit; /*20470 unique individuals for BNX*/

/*# Incident users*/
proc sql;
    select count(distinct studyid_rev) as num_individuals
    from test2_nw_b
	where episode = 1;
quit; /*10871 incident users for BNX*/

/*Total interruption in BNX episodes */
data data_nw_b;
 set data_b;
 if age_eli=1 and palli_excl=0 and canc_excl=0 and preg_eli=0 and incar_grp_bl ne '2_current';
run;
proc freq data = data_nw_b;
table trial_cat;
run; /* 302123 interruptions */



/* BNX  44848 - 37966 = 6882 episodes were excluded */
/*23661 - 20470 = 3191 individuals were excluded 
329459-302123 = 27336 interruptions excluded */






data test2_nw_m;
set test1_m;
 if age_eli=1 and palli_excl=0 and canc_excl=0 and preg_eli=0 and incar_grp_bl ne '2_current';
run;

proc freq data = test2_nw_m;
table trial_cat;
run; /* MHM: 85679 episode of BNX have age >= 18, no preganancy, no incarceration at time0, and no cancer or palliative care history: 
      */

proc sql;
    select count(distinct studyid_rev) as num_individuals
    from test2_nw_m;
quit; /*32102 unique individuals for MET*/


/* BNX  93524 - 85679 = 7845 episodes were excluded */
/* 34651 - 32102 = 2549 individuals were excluded */

proc sql;
    select count(distinct studyid_rev) as num_individuals
    from test2_nw_m
    where episode = 1;
quit; /*15359 incident users of MET*/

/*Total interruption in MET episodes */
data data_nw_m;
 set data_m;
 if age_eli=1 and palli_excl=0 and canc_excl=0 and preg_eli=0 and incar_grp_bl ne '2_current';
run;
proc freq data = data_nw_m;
table trial_cat;
run; /* 1430697 interruptions */

/*1535886 - 1430697 = 105189 interruptions excluded */



/*##########################################*/
/*#######  Box 4: Episodes for Trial category  #######*/
/*##########################################*/

/*## Incident users */
proc freq data=test2_nw_b;
where episode = 1;
 table trial_cat/missing;
run; /* Trial 4: 9510 (1-5 days interruption)  
		Trial 5: 1361 (6-14 days interruption)
		Total  : 10871 episodes*/

/*## Prevalent new users */
proc freq data=test2_nw_b;
 table trial_cat/missing;
run; /* Trial 4: 32255 (1-5 days interruption)  
		Trial 5: 5711 (6-14 days interruption)
		Total  : 37966 episodes*/

/*## All prevalent users */
proc freq data=data_nw_b;
 table trial_cat/missing;
run; /* Trial 4: 275464 (1-5 days interruption)  
		Trial 5: 26659 (6-14 days interruption)
		Total  : 302123 episodes*/




data test2_nw_m; set test2_nw_m; 
trial_cat2 = trial_cat;
if trial_gap = 1 or trial_gap = 2 then trial_cat2 = "1-2";
if trial_gap = 3 or trial_gap = 4 then trial_cat2 = "3-4";
run;

data data_nw_m; set data_nw_m; 
trial_cat2 = trial_cat;
if trial_gap = 1 or trial_gap = 2 then trial_cat2 = "1-2";
if trial_gap = 3 or trial_gap = 4 then trial_cat2 = "3-4";
run;


/*## Incident users */
proc freq data=test2_nw_m;
where episode = 1;
 table trial_cat trial_cat2 /missing;
run; /* Trial 1: 13496 (1-3 days interruption)  
		Trial 2: 369 (4 days interruption)

        Trial 1: 12963 (1-2 days interruption)  
		Trial 2: 902 (3-4 days interruption)

		Trial 3: 1494 (5-14 days interruption)
		Total  : 15359 episodes*/

/*## Prevalent new users */
proc freq data=test2_nw_m;
 table trial_cat trial_cat2 /missing;
run; /* Trial 1: 74412 (1-3 days interruption)  
		Trial 2: 1939 (4 days interruption)

        Trial 1: 71507 (1-2 days interruption)  
		Trial 2: 4844 (3-4 days interruption)

		Trial 3: 9328 (5-14 days interruption)
		Total  : 85679 episodes*/

/*## All prevalent users */
proc freq data=data_nw_m;
 table trial_cat trial_cat2/missing;
run; 
     /* Trial 1: 1332739 (1-3 days interruption)  
		Trial 2: 21424 (4 days interruption)

        Trial 1: 1296216 (1-2 days interruption)  
		Trial 2: 57947 (3-4 days interruption)

		Trial 3: 76534 (5-14 days interruption)
		Total  : 1430697 episodes*/





/*##########################################*/
/*#######  Box 5: Episodes for Treatment  #######*/
/*##########################################*/

/*###########################*/
/*##### Trial 1: 1-3 days #############*/
/*###########################*/

/*Incident user*/
proc freq data=test2_nw_m; 
where trial_cat="1-3" and episode = 1;
table trt_p1/missing; run;

	 /* No doe change : 11133
		Dose reduction: 671
		Dose increase : 1692
		Total  		  : 13496 */

/*Prevalent new user*/
proc freq data=test2_nw_m; 
where trial_cat="1-3" ;
table trt_p1/missing; run;

	 /* No doe change : 61607
		Dose reduction: 2453
		Dose increase : 10352
		Total  		  : 74411 */
proc freq data=test2_nw_m; 
where trial_cat="1-3" ;
table trt_p1*trt_s456/missing; run;

/*All interruptions */
proc freq data=data_nw_m; 
where trial_cat="1-3" ;
table trt_p1/missing; run;
	 /* No doe change : 1204923
		Dose reduction: 40976
		Dose increase : 86840
		Total  		  : 1332739 */


/*###########################*/
/*##### Trial 1: 1-2 days #############*/
/*###########################*/

/*Incident user*/
proc freq data=test2_nw_m; 
where trial_cat2="1-2" and episode = 1;
table trt_p1/missing; run;

	 /* No dose change: 10827
		Dose reduction: 564
		Dose increase : 1572
		Total  		  : 12963 */

/*Prevalent new user*/
proc freq data=test2_nw_m; 
where trial_cat2="1-2" ;
table trt_p1/missing; run;

	 /* No doe change : 59808
		Dose reduction: 1874
		Dose increase : 9825
		Total  		  : 71507 */

/*All interruptions */
proc freq data=data_nw_m; 
where trial_cat2="1-2" ;
table trt_p1/missing; run;
	 /* No doe change : 1185221
		Dose reduction: 27712
		Dose increase : 83282
		Total  		  : 1296216 */


/*###########################*/
/*##### Trial 2: 4 #############*/
/*###########################*/
data test2_nw_m; set test2_nw_m;
trt_p2 = trt_s456;
if trt_s456 in ("1-<25%", "25-<50%") then trt_p2 = "1-<50%";
run;
data data_nw_m; set data_nw_m;
trt_p2 = trt_s456;
if trt_s456 in ("1-<25%", "25-<50%") then trt_p2 = "1-<50%";
run;

/*Incident user*/
proc freq data=test2_nw_m; 
where trial_cat="4" and episode = 1;
table trt_p2/missing; run; 
	 /* No dose change     : 189
		1-<50% reduction   : 67
		>=50% reduct/<=40mg: 16
		dose increase      : 97
		Total  		       : 369
		*/

/*Prevalent new user*/
proc freq data=test2_nw_m; 
where trial_cat="4" ;
table trt_p2/missing; run; 
	 /* No dose change     : 1072
		1-<50% reduction   : 370
		>=50% reduct/<=40mg: 100
		dose increase      : 397
		Total  		       : 1939
		*/
/*All prevalent user*/
proc freq data=data_nw_m; 
where trial_cat="4" ;
table trt_p2/missing; run; 
	 /* No dose change     : 9319
		1-<50% reduction   : 6672
		>=50% reduct/<=40mg: 3572
		dose increase      : 1861
		Total  		       : 21424
		*/

/*###########################*/
/*##### Trial 2: 3-4d #############*/
/*###########################*/

/*Incident user*/
proc freq data=test2_nw_m; 
where trial_cat2="3-4" and episode = 1;
table trt_p2/missing; run; 
	 /* No dose change     : 495
		1-<50% reduction   : 148
		>=50% reduct/<=40mg: 42
		dose increase      : 217
		Total  		       : 902
		*/

/*Prevalent new user*/
proc freq data=test2_nw_m; 
where trial_cat2="3-4" ;
table trt_p2/missing; run; 
	 /* No dose change     : 2871
		1-<50% reduction   : 836
		>=50% reduct/<=40mg: 213
		dose increase      : 924
		Total  		       : 4844
		*/
/*All prevalent user*/
proc freq data=data_nw_m; 
where trial_cat2="3-4" ;
table trt_p2/missing; run; 
	 /* No dose change     : 29021
		1-<50% reduction   : 15746
		>=50% reduct/<=40mg: 7762
		dose increase      : 5418
		Total  		       : 21424
		*/


/*###########################*/
/*##### Trial 3 #############*/
/*###########################*/

/*incident user*/
proc freq data=test2_nw_m; 
where trial_cat="5-14" and episode=1;
table trt_p2/missing; run; 
	 /* No dose change     : 688
		1-<50% reduction   : 265
		>=50% reduct/<=40mg: 107
		Dose increase 	   : 434
		Total  		       : 1494 */

/*Prevalent new user*/
proc freq data=test2_nw_m; 
where trial_cat="5-14" ;
table trt_p2/missing; run; 
	 /* No dose change     : 5113
		1-<50% reduction   : 1763
		>=50% reduct/<=40mg: 560
		Dose increase 	   : 1892
		Total  		       : 9328 */


/*All prevalent user*/
proc freq data=data_nw_m; 
where trial_cat="5-14" ;
table trt_p2/missing; run; 
	 /* No dose change     : 30936
		1-<50% reduction   : 23292
		>=50% reduct/<=40mg: 14983
		Dose increase 	   : 7323
		Total  		       : 76534 */


/*###########################*/
/*##### Trial 4 #############*/
/*###########################*/
/*Incident user*/
proc freq data=test2_nw_b; 
where trial_cat="1-5" and episode = 1;
table trt_p4s/missing; run; 
/*We can add the increase category for trial 4 also. */
     /* No doe change      : 6769
		Dose   reduction   : 1439
		Dose increase	   : 1302
		Total  		       : 9510
		*/

/*Prevalent new user*/
proc freq data=test2_nw_b; 
where trial_cat="1-5" ;
table trt_p4s/missing; run; 
/*We can add the increase category for trial 4 also. */
     /* No doe change      : 24114
		Dose   reduction   : 3374
		Dose increase	   : 4767
		Total  		       : 32255
		*/

/*All prevalent user*/
proc freq data=data_nw_b; 
where trial_cat="1-5" ;
table trt_p4s/missing; run; 
/*We can add the increase category for trial 4 also. */
     /* No doe change      : 247319
		Dose   reduction   : 14335
		Dose increase	   : 13810
		Total  		       : 275464
		*/


/*###########################*/
/*##### Trial 5 #############*/
/*###########################*/
data test2_nw_b; set test2_nw_b;
trt_p5 = trt_s5;
if trt_s5 in ("1-<25%", "25-<50%") then trt_p5 = "1-<50%";
run;

data data_nw_b; set data_nw_b;
trt_p5 = trt_s5;
if trt_s5 in ("1-<25%", "25-<50%") then trt_p5 = "1-<50%";
run;

/*Incident user*/
proc freq data=test2_nw_b; 
where trial_cat="6-14" and episode = 1;
table trt_p5/missing; run; 
     /* No dose change     : 560
		1-<50% reduction   : 219
		>=50%/<=4mg        : 269
		Dose increase 	   : 313
		Total  		       : 1361
		*/

/*Prevalent new user*/
proc freq data=test2_nw_b; 
where trial_cat="6-14";
table trt_p5/missing; run; 
     /* No dose change     : 2957
		1-<50% reduction   : 776
		>=50%/<=4mg        : 826
		Dose increase 	   : 1152
		Total  		       : 5711
		*/

/*All prevalent user*/
proc freq data=data_nw_b; 
where trial_cat="6-14";
table trt_p5/missing; run; 
     /* No dose change     : 16131
		1-<50% reduction   : 4092
		>=50%/<=4mg        : 3481
		Dose increase 	   : 2955
		Total  		       : 26659
		*/
