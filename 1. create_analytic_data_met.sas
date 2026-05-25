proc datasets kill;
quit;

options compress=y;

libname der "W:\SUD\Data\2022\derived";
libname lib "W:\SUD\Data\2022\SAS data";
libname pss "W:\SUD\Data\2022\derived\PSS";
libname tmplib "E:\SAS Data\Missed dose protocol";
libname tmplib1 "E:\SAS Data\Missed dose protocol_v1";
libname tmplib2 "E:\SAS Data\Missed dose protocol_v2";

*oat episode;
data epi;
set der.Epi_hosp_14;
by studyid_rev;
*gap between episode;
dt_end_lag=lag(dt_end);
if first.studyid_rev then dt_end_lag=.;
gap_betwn_epi= dt_start-dt_end_lag; 
gap30_eli=(gap_betwn_epi=. or gap_betwn_epi>=30);
gap90_eli=(gap_betwn_epi=. or gap_betwn_epi >= 90);
* the washout will be applied later (gap_betwn_epi=. or gap_betwn_epi>=30);
* any episodes starting with either met, bnx, srom (or any combination)=> need to resctrict;
if (find(drug_start, 'M')>0 or find(drug_start,'B')>0 or find(drug_start,'S')>0)
	and dt_start>='1jan2010'd and drug_start ne 'BUP';
*BUP indicates buprenophine monotherapy=> need to exclude;
format dt_end_lag date9.;
run;

proc freq data=epi;
table drug_start; run;*B, BS, M, MB, MBS, MS, S;

/*************************/
/*Creating my time zero*/
/************************/
proc sql;
create table epi_b as
select b.*, a.dt2, ddose, days, hosp, drug, predrug
from epi b 
left join der.days_m a 
on a.studyid_rev = b.studyid_rev 
   and a.dt2 >= b.dt_start 
   and a.dt2 <= b.dt_end
order by b.studyid_rev, b.episode, a.dt2;
quit;


data epi_b;
   set epi_b;
   by studyid_rev episode;
   lag_dt2 = lag(dt2);
   format lag_dt2 date9.;
   if first.episode then lag_dt2 = .;
   gap_days = dt2 - lag_dt2 - 1;
run;


data epi_b;
   set epi_b;
   by studyid_rev episode;
   if first.episode then myepisode = .; 
   if gap_days ne . then do;
     	myepisode + (gap_days ne 0); 
      end;
   output;
run;


data epi_b;
   set epi_b;
   by studyid_rev episode myepisode;
   if myepisode <= 0 then trial_gap = .;
   if myepisode >= 1 then do;
       if first.myepisode then trial_gap = 0;
       trial_gap + gap_days;
	   end;
    output;
run;



data epi_b;
   set epi_b;
   lag_ddose = lag(ddose);
   by studyid_rev episode myepisode;
   if myepisode > 0 and first.myepisode then do;
      ddose_prev = lag_ddose;
	  ddose_red = lag_ddose - ddose;
      ddose_rev = (ddose_red/lag_ddose)*100;
	  time0 = dt2;
      end;
	  format time0 date9.;
   if myepisode > 0 and first.myepisode;
run;


/**************************************************
	Outcome
***************************************************/
*the first date of overdose-related acute care
(only based on the primary reason [most responsible];
proc sql;
create table f_od_acute as
select w.studyid_rev, w.episode, w.myepisode, 
	min(z.dt) format=date9. as f_od_acute_dt
from epi_b as w left join der.Sud_psych_acute(
	where=((dad=1 or nacrs=1 or (msp=1 and servloc in ('I','E'))) and main_od=1)) as z
on w.studyid_rev=z.studyid_rev and w.time0 <= z.dt
group by w.studyid_rev, w.episode, w.myepisode;
quit;



*overdose death (either based on ICD code from vs or any records from coroners service);
proc sql;
create table od_death as
select w.studyid_rev, w.episode, w.myepisode, min(z.dt) format=date9. as od_death_dt
from epi_b as w left join der.Sud_psych_acute(where=(death=1 and main_od=1)) as z
on w.studyid_rev=z.studyid_rev and w.time0 <= z.dt
group by w.studyid_rev, w.episode, w.myepisode;
quit;

* all-cause acute care (first record from DAD/NACRS since OAT epi initiation);
proc sql;
create table f_hosp as
select w.studyid_rev, w.episode, w.myepisode, 
	min(z.addate_rev) format=date9. as f_hosp_dt
from epi_b as w left join der.dad_rev as z
on w.studyid_rev=z.studyid_rev and w.time0 <= z.addate_rev
group by w.studyid_rev, w.episode, w.myepisode;
quit;


proc sql;
create table f_ed as
select w.studyid_rev, w.episode, w.myepisode, 
	min(z.regdt_rev) format=date9. as f_ed_dt
from epi_b as w left join der.nacrs_rev(where=(edvisit=1)) as z
on w.studyid_rev=z.studyid_rev and w.time0 <= z.regdt_rev
group by w.studyid_rev, w.episode, w.myepisode;
quit;


data f_acute;
merge f_hosp f_ed;
by studyid_rev episode myepisode;
f_acute_dt=min(f_hosp_dt,f_ed_dt);
format f_acute_dt date9.;
run;

/*************************************************
	merge with cohort
***************************************************/
data cohort;
merge der.cohort(keep=studyid_rev f_oat_dt end_fu_dt rev_dod 
		dob gender f_oud_dt inmate_dup_ind where=(inmate_dup_ind=0))
	epi_b(in=a) f_od_acute od_death f_acute(drop=f_hosp_dt f_ed_dt)
	tmplib.prior_po(in=b);
by studyid_rev;
*years since oud diagnosis;
yr_since_oud=(time0-f_oud_dt)/365.25; 
*year of episode initiation;
yr_start=year(time0); 
*end of episode+gap;
if drug_end in ('M','MS','S') then dt_end_gap=min(end_fu_dt, dt_end+14);
else if drug_end='I' then dt_end_gap=min(end_fu_dt, dt_end+14);
else dt_end_gap=min(end_fu_dt, dt_end+14);
*po prior to OUD diagnosis;
po_preOUD=0;
if b and first_po < f_oud_dt then po_preOUD=1;
if a;
format dt_end_gap date9.;
drop inmate_dup_ind;
run;


/*************************************************
	Eligibility criteria 
***************************************************/
data canc_palli;
merge cohort(in=a) der.cohort_comorb(keep=studyid_rev first_canc_date first_palli_date);
by studyid_rev;
palli_excl=(first_palli_date ne . and first_palli_date <= time0);
canc_excl=(first_canc_date ne . and first_canc_date <= time0);
if a;
run;

*the first date of pregnacy at time0;
proc sql;
create table f_preg as
select e.studyid_rev, e.episode, e.myepisode, e.time0,
	min(p.first_hc_contact_preg_rev) format=date9. as f_preg_dt
from cohort as e left join der.peri_episode as p
on e.studyid_rev=p.studyid_rev and (
p.first_hc_contact_preg_rev <= e.time0 <= p.baby_date_of_birth_rev+90)
group by e.studyid_rev, e.episode, e.episode, e.time0;
quit;

*the first date of incarceration at time0;
proc sql;
create table f_incar as
select e.studyid_rev,e.episode, e.myepisode, e.time0,
	min(p.adm_dt_rev) format=date9. as f_incar_dt
from cohort as e left join der.inmate_epi as p
on e.studyid_rev=p.studyid_rev and (
p.adm_dt_rev <= e.time0 <= p.rel_dt_rev )
group by e.studyid_rev,e.episode, e.myepisode, e.time0;
quit;

*incarceration in the past year at time 0 (including time0);
proc sql;
create table pastyear_incar as
select w.studyid_rev, w.episode, w.myepisode,
	max(adm_dt_rev) as last_yr_adm_incar format=date9.,
	max(rel_dt_rev) as last_yr_rel_incar format=date9.
from cohort as w left join der.inmate_epi as i
on w.studyid_rev=i.studyid_rev and 
(w.time0-365 <= i.rel_dt_rev <= w.time0 or 
w.time0-365 <= i.adm_dt_rev <= w.time0 or
(i.adm_dt_rev<=w.time0-365 and w.time0<=i.rel_dt_rev))
group by w.studyid_rev, w.episode, w.myepisode;
quit;


*the first date of oop indication per year;
data oop;
set der.lha_yr(where=(year>=2010 and oop=1));
f_dt_oop=min(msp_f_dt, nacrs_f_dt,dad_f_dt, pnet_f_dt, mdy(12,31,year));
format f_dt_oop date9.;
run;


proc sql;
  create table oop_cohort as select w.studyid_rev, w.episode, w.myepisode, w.time0, w.dt_end, w.end_fu_dt,
    min(i.f_dt_oop) format=date9. as f_dt_oop from cohort as w 
  left join oop as i
    on w.studyid_rev = i.studyid_rev 
    and year(w.time0) <= year(i.f_dt_oop) 
    and year(i.f_dt_oop) <= year(w.end_fu_dt)
  group by 
    w.studyid_rev, w.episode, w.myepisode, w.time0, w.dt_end, w.end_fu_dt;
quit;


data cohort2;
merge cohort canc_palli f_preg f_incar pastyear_incar oop_cohort;
by studyid_rev episode myepisode;
*incarceration indication at time 0;
length incar_grp_bl $15.;
if last_yr_adm_incar = . and last_yr_rel_incar=. then incar_grp_bl='0_neveror>1yr';
else if last_yr_rel_incar ne . and last_yr_rel_incar < time0 then incar_grp_bl='1_pastyear';
else if last_yr_adm_incar <= time0<=last_yr_rel_incar then incar_grp_bl='2_current';
*oop;
oop_bl=0;
if f_dt_oop ne . and f_dt_oop<=time0 then oop_bl=1;
drop last_yr_adm_incar last_yr_rel_incar;
run;

proc freq; table incar_grp_bl; run;



*save cohort;
data tmplib1.TT_cohort_m;
set cohort2;
run;


/* part of SA */

data tmplib1.TT_cohort_m;
set tmplib1.TT_cohort_m;
by studyid_rev;
*end of episode+gap;
if drug_end in ('M','MS','S') then dt_end_gap30=min(end_fu_dt, dt_end+30);
else if drug_end='I' then dt_end_gap30=min(end_fu_dt, dt_end+30);
else dt_end_gap30=min(end_fu_dt, dt_end+30);
format dt_end_gap30 date9.;
if drug_end in ('M','MS','S') then dt_end_gap14=min(end_fu_dt, dt_end+14);
else if drug_end='I' then dt_end_gap14=min(end_fu_dt, dt_end+14);
else dt_end_gap14=min(end_fu_dt, dt_end+14);
format dt_end_gap14 date9.;
run;



/*###################################*/
/* Covariates creation at Time zero */
/*####################################*/


/*data tmplib.udt_days;*/
/*set lib.msp;*/
/*udt=0;*/
/*poc=0;non_oat=0;*/
/*if (fitm in (15039,15040,92518,92550,36521) or 92500<=fitm<=92514) */
/*	and fitm not in (92504,92509,92512) then udt=1;*/
/** confirmatory urine tests;*/
/**if (fitm in (92520,92521,92525,92543,92545,92546) or */
/*	92527<=fitm<=92529 or 92534<=fitm<=92539) and */
/*	fitm not in (92537) then udt=1;*/
/**if fitm = 36521 then udt=1;*/
/*if udt and servdt >= '01JAN2009'd;*/
/*if fitm in (15039,15040,36521) then poc=1;*/
/*if fitm = 15040 then non_oat=1;*/
/*run;*/

* UDT in past month ;
proc sql;
create table past_month_udt as
select w.studyid_rev, w.episode, w.myepisode, w.time0,
	max(z.poc) as poc_udt, max(1-z.poc) as lab_udt, max(z.servdt ne .) as udt
from tmplib1.TT_cohort_m as w left join tmplib.udt_days as z
on w.studyid_rev=z.studyid_rev and w.time0-30 <= z.servdt < w.time0
group by w.studyid_rev, w.episode, w.myepisode, w.time0;
quit;


proc format;
 value att low-<26='1<=25%' 26-<51='2<50%' 51-<75='3<=75%' 75-high='4>75%';
run;

/*data tmplib.homeless_recs;*/
/*set der.pwoud_comorb_records;*/
/*if cuml_nfa >=3 or home_icd=1;*/
/*run;*/


proc sql;
create table past_homeless as
select w.studyid_rev, w.episode, w.myepisode, w.time0,
	max(case when i.cuml_nfa >= 3 or i.home_icd=1 then i.dt else . end) as l_homeless format=date9.
from tmplib1.TT_cohort_m as w left join tmplib.homeless_recs as i
on w.studyid_rev=i.studyid_rev and w.time0-1826 <= i.dt < w.time0
group by w.studyid_rev, w.episode, w.myepisode, w.time0; *1826 days = 5 years;
quit;


data past_homeless;
set past_homeless;
length home_grp $15.;
if l_homeless=. then home_grp='0_neveror>5yr';
else if l_homeless <= time0-366 then home_grp='1_1to5years';
else home_grp='2_pastyear';
year=year(time0);
run;


/*data tmplib.lha_yr;*/
/*set der.lha_yr;*/
/*by studyid_rev year;*/
/*retain lha_pre;*/
/*if first.studyid_rev then lha_pre=lha;*/
/*if lha ne . and lha notin (888,999) then lha_pre=lha;*/
/*run;*/
/**/

*impute lha from previous nonmissing one if oop_clr=0;
proc sql;
create table lha as
select a.*, b.oop,
case when b.oop=0 and (b.lha=. or b.lha in (888,999)) then b.lha_pre else b.lha end as lha
from past_homeless as a left join tmplib.lha_yr as b on
a.studyid_rev=b.studyid_rev and a.year=b.year
order by a.studyid_rev, a.episode, a.myepisode, a.time0; 
quit; 



PROC IMPORT OUT=geo
            DATAFILE= "W:\SUD\Documents\2020\Documents\geo_rural.xlsx"
            DBMS=xlsx REPLACE ;
			GETNAMES=YES;
			sheet='sheet1';
RUN;


proc sql;
create table current_rural_status as
select dat.*, geo.rural, geo.ha_no
from lha as dat left join geo on (dat.lha=geo.lha_no)
order by studyid_rev, episode, myepisode, time0;
quit;

proc freq data=current_rural_status;
  tables rural / missing;
run;


/**virtual care;*/
/*data tmplib.msp_virtual;*/
/*set der.msp_virtual;*/
/*array diag diagcd diagcd2 diagcd3;*/
/*oud=0; */
/*if fitm in (39, 15039,13013,13014,36521) then oud=1;*/
/*do i=1 to 3;*/
/*if diag(i) in: ("304","3040","3047","3055","9650",*/
/*					"E8500","F11") then oud=1;*/
/*end;/*note that icd=304 was added*/*/
/*run;*/;

/*************************************************
	any virtual care in the past month
***************************************************/
proc sql;
create table past_virtual_care as
select e.studyid_rev, e.episode, e.myepisode, e.time0, max(p.servdt ne .) as virtual_care,
	max(p.oud=1) as virtual_oud
from tmplib1.TT_cohort_m as e left join tmplib.msp_virtual as p
on e.studyid_rev=p.studyid_rev and e.time0-30 <= p.servdt <= e.time0 
group by  e.studyid_rev, e.episode, e.myepisode, e.time0;
quit; /*note that wk_start=servdt was included as dispensation occurs after seeing a doctor*/


proc sql;
create table prac as
select e.studyid_rev, e.episode, e.myepisode, e.time0, p.ddose as f_ddose,
	p.prac_studyid1 as prac_studyid
from tmplib1.TT_cohort_m as e left join der.days_m as p
on e.studyid_rev=p.studyid_rev and e.time0=p.dt2 
order by e.studyid_rev, e.episode, e.myepisode;
quit;

*var related to prescriber;
proc sql;
create table prov2 as 
select e.studyid_rev, e.episode, e.myepisode, prov.prac_studyid_rev, e.time0, prov.chc, prov.chc_sa, prov.PRSCR_PRAC_PROF, clin_size, clin_size_sa
from prac as e left join der.prov_clinics_v2 as prov
on e.prac_studyid=prov.prac_studyid_rev and month(e.time0)=prov.month and
	year(e.time0)=prov.year
order by e.studyid_rev, e.episode, e.myepisode;
quit;

proc freq data=prov2; table clin_size/missing; run;


*pss at time0;
proc sql;
create table pss_time0 as
select e.studyid_rev, e.episode, e.myepisode, e.time0, p.dt2 ne . as pss
from tmplib1.TT_cohort_m as e left join pss.Pss_opioid_epi_days as p
on e.studyid_rev=p.studyid_rev and e.time0= p.dt2 
order by e.studyid_rev, e.episode, e.myepisode;
quit;



data past_comorb;
merge tmplib1.TT_cohort_m(in=a) der.cohort_comorb(keep=studyid_rev first_sud_date first_aud_date
	first_hiv_date first_hcv_date first_bpd_date first_mdd_date first_schiz_date
	first_canc_date first_palli_date first_asthma_date first_copd_date);
by studyid_rev;

sud=(first_sud_date ne . and first_sud_date < time0);
aud=(first_aud_date ne . and first_aud_date < time0);
hiv=(first_hiv_date ne . and first_hiv_date < time0);
hcv=(first_hcv_date ne . and first_hcv_date < time0);
asthma=(first_asthma_date ne . and first_asthma_date < time0);
copd=(first_copd_date ne . and first_copd_date < time0);
first_canc=(first_canc_date ne . and first_canc_date < time0);
first_palli=(first_palli_date ne . and first_palli_date < time0);
canc_palli= (first_canc or first_palli);
mdd=(first_mdd_date ne . and first_mdd_date < time0);
bpd=(first_bpd_date ne . and first_bpd_date < time0);
schiz=(first_schiz_date ne . and first_schiz_date < time0);
smd=max(mdd, bpd, schiz);
if a;
keep studyid_rev episode myepisode time0 sud--canc_palli smd;
run;

proc sql;
create table past_sed_psy as
select w.studyid_rev, w.episode, w.myepisode, w.time0,
	max(sed=1) as sed_med, max(psy=1) as psy_med
from tmplib1.TT_cohort_m as w 
left join der.Sed_psy_days(where=(dt >= '01DEC2009'd)) as b
on w.studyid_rev=b.studyid_rev and (
	w.time0-30 <= b.dt < w.time0 or
	w.time0-30 <= b.dt+b.days-1 < w.time0 or
	(b.dt<=w.time0-30 and w.time0<b.dt+b.days-1))
group by w.studyid_rev, w.episode, w.myepisode, w.time0;
quit;



proc sql;
create table past_opioid as
select w.studyid_rev, w.episode, w.myepisode, w.time0, 
	max(p.dt2 ne .) as po, count(distinct p.dt2) as d_po
from tmplib1.TT_cohort_m as w left join der.opioid_rx_days(where=(dt2 >= '01DEC2009'd)) as p
on w.studyid_rev=p.studyid_rev and w.time0-30 <= p.dt2 < w.time0
group by w.studyid_rev, w.episode, w.myepisode, w.time0;
quit;



proc sql;
create table past_hosp as
select w.studyid_rev, w.episode, w.myepisode, w.time0, 
	max(z.main_sud=1) as drug_acute, max(z.main_psych=1) as psy_hosp
from tmplib1.TT_cohort_m as w 
left join der.Sud_psych_acute(where=(dad=1 or nacrs=1)) as z
on w.studyid_rev=z.studyid_rev and w.time0-30 <= z.dt < w.time0
group by w.studyid_rev, w.episode, w.myepisode, w.time0;
quit;

proc sql;
create table past_comorb2 as
select w.studyid_rev, w.episode, w.myepisode, w.time0, 
	max(z.mh=1) as mh, max(z.cp=1) as cp, max(z.tob=1) as tob, 
	max(z.canc=1) as canc, max(z.palli=1) as palli, max(z.social=1) as social_assist
from tmplib1.TT_cohort_m as w left join der.pwoud_comorb_records(where=(dt >= '01JAN2009'd)) as z
on w.studyid_rev=z.studyid_rev and w.time0-365 <= z.dt < w.time0
group by w.studyid_rev, w.episode, w.myepisode, w.time0;
quit;

/*proc sql;*/
/*create table phy_att_m as*/
/*select e.studyid_rev, e.episode, e.myepisode, e.time0, w.tot_visits, w.max_visits, w.tot_disp, w.max_disp*/
/*from tmplib1.TT_cohort_m as e inner join der.TT_weekly_m as w */
/*on e.studyid_rev = w.studyid_rev and (w.wk_start between (e.time0 - 14) and e.time0)*/
/*group by e.studyid_rev, e.episode, e.myepisode, e.time0*/
/*having max(wk_start)=wk_start;*/
/*quit;*/ 
 

/*proc freq data=phy_att_m; table phy_att phy_att_pnet/missing; run; */
/*18491 missing values are found in phy_att, and phy_att_pnet*/

/*Trying differently*/

/**/
/*proc sql;*/
/*    create table phy_att_dat as*/
/*    select e.studyid_rev, e.time0, max(w.tot_visits) as tot_visits, max(w.max_visits) as max_visits, */
/*	max(w.tot_disp) as tot_disp, max(w.max_disp) as max_disp */
/*    from tmplib1.TT_cohort_m as e inner join der.TT_weekly_m as w */
/*    on e.studyid_rev = w.studyid_rev and (w.wk_start between (e.time0 - 20) and e.time0)*/
/*    group by e.studyid_rev, e.time0;*/
/*quit;*/


/*%include 'W:\SUD\Documents\2020\SAS code\cci_macro.sas';*/
/** for macro "_CharlsonICD9CM";*/
/*%include 'W:\SUD\Documents\Apr2019\SAS code\MCHP_SAS_macro_charlson_icd9CM_index.sas';*/
/** for macro "_CharlsonICD10";*/
/*%include 'W:\SUD\Documents\Apr2019\SAS code\MCHP_SAS_macro_charlson_icd10_index.sas';*/
/**/
/**/
/*%cci_prior_time0(dad=der.dad_rev, dat=tmplib.TT_analytic_b, id=studyid_rev, dt=time0,*/
/*	mo=12, output=cci_m);*/
/**/


data TT_cohort_m;
merge tmplib1.TT_cohort_m(in=a) tmplib.cci_m;
by studyid_rev time0;
cci_gt0 = cci > 0;
if a;
run;

proc freq data=TT_cohort_m; table cci cci_gt0/missing; run;


data TT_analytic_m;
merge TT_cohort_m (in=a) past_month_udt past_homeless lha current_rural_status past_virtual_care prov2 pss_time0 
past_comorb past_sed_psy past_opioid past_hosp past_hosp past_comorb2;
by studyid_rev episode myepisode time0;
if udt=0 then do; poc_udt=0; lab_udt=0; end;
if rural=. then rural=0;
*if clin_size=. then clin_size='1_small_<=13';
rename lha=clnt_lha ha_no=clnt_ha;
run;


proc freq data=TT_analytic_m; table clin_size /missing; run;


data tmplib1.TT_analytic_m1;
set TT_analytic_m;
by studyid_rev episode myepisode;
age_st=(time0-dob)/365.25;
age_eli=(age_st >= 18);
preg_eli = 0;
clin_size_b = clin_size = '3_large_>=278';
if f_preg_dt ne . then preg_eli = 1;
run;

proc freq data=tmplib1.TT_analytic_m1; table clin_size_b clin_size cci/missing; run;


/*data tmplib1.TT_analytic_m1;*/
/*set tmplib1.TT_analytic_m1;*/
/*by studyid_rev;*/
/**end of episode+gap;*/
/*if drug_end in ('M','MS','S') then dt_end_gap30=min(end_fu_dt, dt_end+30);*/
/*else if drug_end='I' then dt_end_gap30=min(end_fu_dt, dt_end+30);*/
/*else dt_end_gap30=min(end_fu_dt, dt_end+30);*/
/*format dt_end_gap30 date9.;*/
/*if drug_end in ('M','MS','S') then dt_end_gap14=min(end_fu_dt, dt_end+14);*/
/*else if drug_end='I' then dt_end_gap14=min(end_fu_dt, dt_end+14);*/
/*else dt_end_gap14=min(end_fu_dt, dt_end+14);*/
/*format dt_end_gap14 date9.;*/
/*run;*/



/*proc sql;*/
/*create table gp_msp as*/
/*select dat.studyid_rev, dat.time0, dat.prac_studyid_rev, */
/*	count(distinct msp.servdt) as n_visits*/
/*from tmplib1.TT_analytic_m1 as dat left join lib.msp (where=(clmspec in (0,21,87) and servdt>='1jan2009'd)) as msp*/
/*on dat.studyid_rev=msp.studyid_rev and dat.time0-365<= msp.servdt <dat.time0*/
/*group by dat.studyid_rev, dat.time0, dat.prac_studyid_rev; quit;*/
/**/
/*data tmplib.gp_msp_m; set gp_msp; run;*/

*find total # GP visits;
proc sql;
create table gp_msp_total as
select studyid_rev, time0, sum(n_visits) as tot_visits, max(n_visits) as max_visits
from tmplib.gp_msp_m group by studyid_rev, time0;
quit;

data gp_msp_total2;
set gp_msp_total;
if tot_visits<3 then phy_att=0;
else if max_visits<=tot_visits/2 then phy_att=1;
else if max_visits>tot_visits/2 then phy_att=2;
*more than 50% of the visits were to a single GP;
run;


data TT_analytic_m;
merge tmplib1.TT_analytic_m1 (in=a) gp_msp_total2;
by studyid_rev time0;
if a;
run;


/*#################
    Exposure and outcomes category 
################*/

data TT_analytic_m;
    set TT_analytic_m;
    length trial_cat $10 trt_p1 trt_p2 trt_p3 trt_s456 $20;

/*     Trial Category */
    if 1 <= trial_gap <= 3 then trial_cat = "1-3";
    else if trial_gap = 4 then trial_cat = "4";
    else if 5 <= trial_gap <= 14 then trial_cat = "5-14";
    else trial_cat = "";

/*     Trial 1 Treatment */
    if ddose_rev = 0 then trt_p1 = "no reversion";
    else if ddose_rev > 0 then trt_p1 = "reversion";
    else if ddose_rev < 0 then trt_p1 = "increase";
    else trt_p1 = "";

/*     Trial 2 Treatment */
    if ddose_rev = 0 then trt_p2 = "no reversion";
    else if 0 < ddose_rev < 50 then trt_p2 = "1-<50%";
    else if ddose_rev >= 50 or (ddose ne . and ddose < ddose_prev and ddose <= 40) then trt_p2 = ">=50%/<=40mg";
    else trt_p2 = "";

/*     Trial 3 Treatment */
    if ddose_rev = 0 then trt_p3 = "no reversion";
    else if 0 < ddose_rev < 50 then trt_p3 = "1-<50%";
    else if ddose_rev >= 50 or (ddose ne . and ddose < ddose_prev and ddose <= 40) then trt_p3 = ">=50%/<=40mg";
    else trt_p3 = "";

/*     Sensitivity for Trials 4, 5, 6 */
    if ddose_rev = 0 then trt_s456 = "no reversion";
    else if 0 < ddose_rev < 25 then trt_s456 = "1-<25%";
    else if 25 <= ddose_rev < 50 then trt_s456 = "25-<50%";
    else if ddose_rev >= 50 or (ddose ne . and ddose < ddose_prev and ddose <= 40) then trt_s456 = ">=50%/<=40mg";
    else if ddose_rev < 0 then trt_s456 = "increase";
    else trt_s456 = "";
run;


data TT_analytic_m;
set TT_analytic_m;
by studyid_rev episode myepisode;
*JM2: revised one from your code;
s_discont = 1;
if rev_dod ne . and rev_dod <= dt_end_gap then s_discont = 0; 
if f_dt_oop ne . and f_dt_oop >= time0 and f_dt_oop < dt_end_gap then s_discont = 0;
if dt_end_gap='31DEC2022'd  then s_discont = 0; 

/*if s_discont = 0 then t_discont = min(end_fu_dt, f_dt_oop)-time0;*/
if s_discont = 0 then do;
   if f_dt_oop >= time0 then t_discont = min(end_fu_dt, f_dt_oop) - time0;
   else t_discont = end_fu_dt - time0;
end;
else if s_discont = 1 then t_discont = dt_end_gap-time0;  


/*14 days gap*/
s_discont_14 = 1;
if rev_dod ne . and rev_dod <= dt_end_gap14 then s_discont_14 = 0; 
if f_dt_oop ne . and f_dt_oop >= time0 and f_dt_oop < dt_end_gap14 then s_discont_14 = 0;
if dt_end_gap14='31DEC2022'd  then s_discont_14 = 0; 

/*if s_discont = 0 then t_discont = min(end_fu_dt, f_dt_oop)-time0;*/
if s_discont_14 = 0 then do;
   if f_dt_oop >= time0 then t_discont_14 = min(end_fu_dt, f_dt_oop) - time0;
   else t_discont_14 = end_fu_dt - time0;
end;
else if s_discont_14 = 1 then t_discont_14 = dt_end_gap14-time0;  

/*30 days gap*/
s_discont_30 = 1;
if rev_dod ne . and rev_dod <= dt_end_gap30 then s_discont_30 = 0; 
if f_dt_oop ne . and f_dt_oop >= time0 and f_dt_oop < dt_end_gap30 then s_discont_30 = 0;
if dt_end_gap30='31DEC2022'd  then s_discont_30 = 0; 

/*if s_discont = 0 then t_discont = min(end_fu_dt, f_dt_oop)-time0;*/
if s_discont_30 = 0 then do;
   if f_dt_oop >= time0 then t_discont_30 = min(end_fu_dt, f_dt_oop) - time0;
   else t_discont_30 = end_fu_dt - time0;
end;
else if s_discont_30 = 1 then t_discont_30 = dt_end_gap30-time0;  



s_od_acute=(f_od_acute_dt ne .);
if f_od_acute_dt>end_fu_dt then s_od_acute = 0;
if rev_dod ne . and  rev_dod < f_od_acute_dt then s_od_acute = 0;
if f_dt_oop ne . and f_dt_oop >= time0 and f_dt_oop < f_od_acute_dt then s_od_acute = 0;

/*if s_od_acute = 0 then t_od_acute = min(end_fu_dt, f_dt_oop)- time0;*/
if s_od_acute = 0 then do;
   if f_dt_oop >= time0 then t_od_acute = min(end_fu_dt, f_dt_oop) - time0;
   else t_od_acute = end_fu_dt - time0;
end;
else if s_od_acute = 1 then t_od_acute = f_od_acute_dt-time0; 

s_all_acute=(f_acute_dt ne .);
if f_acute_dt>end_fu_dt then s_all_acute = 0;
if rev_dod ne . and rev_dod < f_acute_dt then s_all_acute = 0;
if f_dt_oop ne . and f_dt_oop >= time0 and f_dt_oop < f_acute_dt then s_all_acute = 0;

/*if s_all_acute = 0 then t_all_acute = min(end_fu_dt, f_dt_oop)- time0;*/
if s_all_acute = 0 then do;
   if f_dt_oop >= time0 then t_all_acute = min(end_fu_dt, f_dt_oop) - time0;
   else t_all_acute = end_fu_dt - time0;
end;
else if s_all_acute = 1 then t_all_acute = f_acute_dt-time0;  

*JM2: note that time variables can be negative if oop_bl=1;
if oop_bl=1 then do; s_discont=0;s_od_acute=0;s_all_acute=0; end;

*JM2: my version;
enddt_discont=min(dt_end_gap, end_fu_dt, f_dt_oop);
s_discont2=(enddt_discont=dt_end_gap and dt_end_gap ne '31dec2022'd and dt_end_gap ne rev_dod);
t_discont2=enddt_discont-time0;

enddt_od_acute=min(f_od_acute_dt, end_fu_dt, f_dt_oop);
s_od_acute2=(enddt_od_acute=f_od_acute_dt);
t_od_acute2=enddt_od_acute-time0;

enddt_all_acute=min(f_acute_dt,end_fu_dt, f_dt_oop);
s_all_acute2=(enddt_all_acute=f_acute_dt);
t_all_acute2=enddt_all_acute-time0;


format enddt_discont enddt_od_acute enddt_all_acute date9.;
run;

*JM2: compare two versions;
proc freq data=TT_analytic_m;
table s_discont*s_discont2  s_od_acute*s_od_acute2 s_all_acute*s_all_acute2/missing;
run;

proc freq data=TT_analytic_m;
table s_discont*s_discont_14  s_discont*s_discont_30/missing;
run;

proc means data=TT_analytic_m;
 var t_od_acute2 t_discont2 t_od_acute t_discont;
 run;


data tmplib1.TT_analytic_m1;
	set TT_analytic_m;
	by studyid_rev episode myepisode;

/*	Primary outcomes*/
	sens_discont = s_discont_14;
	time_discont = min(t_discont_14, 365); 
    if t_discont_14 >= 365 then sens_discont = 0;

	sens_od_acute = s_od_acute;
	time_od_acute = min(t_od_acute, 365); 
    if t_od_acute >= 365 then sens_od_acute = 0;
    
	sens_all_acute = s_all_acute;
	time_all_acute = min(t_all_acute, 365); 
    if t_all_acute >= 365 then sens_all_acute = 0;
 
/*	Sensitivity outcomes*/
	sens_discont30 = s_discont_30;
	time_discont30 = min(t_discont_30, 365); 
    if t_discont_30 >= 365 then sens_discont30 = 0;

run;


proc means data=tmplib1.TT_analytic_m1;
 var t_od_acute t_discont time_discont t_discont_30;
 run;



/*
Saving the whole data set
*/

data tmplib1.TT_analytic_all_m; set tmplib1.TT_analytic_m1; run;

proc freq data=tmplib1.TT_analytic_all_m;
where ddose ne . or predrug in ("M", "") and trial_cat="5-14";
table trt_p3*trt_s456/missing;*missing in trt_p3 should be only due to "increase";
run;


proc freq data=tmplib1.TT_analytic_all_m;
where ddose ne . or predrug in ("M", "") and trial_cat="4";
table trt_p2*trt_s456/missing;*missing in trt_p2 should be only due to "increase";
run;


data tmplib1.TT_analytic_all_m; set tmplib1.TT_analytic_all_m;
	if age_st < 30 then age_cat = '[-Inf,30)';
	else if age_st < 40 then age_cat = '[30,40)';
	else if age_st < 50 then age_cat = '[40,50)';
	else age_cat = '[50,Inf)';

	sexM = (gender = 'M');
	age_lt30=(age_st<30);
	age_30_39=(30<=age_st<40);
	age_40_49=(40<=age_st<50);
	age_ge50=(age_st>=50);

	rural=(rural=1);
	po_preOUD=(po_preOUD=1);
	phy_att_gr2=(phy_att=2);
	social_assist=(social_assist=1);
	 *homeless;
	home_grp_1=(home_grp='1_1to5years');
	home_grp_2=(home_grp='2_pastyear');
	home_grp_b=(home_grp_1 or home_grp_2);
	past_incar_grp=(incar_grp_bl='1_pastyear');

	clin_size_b=(clin_size='3_large_>=278');
	sed_psy_med=(sed_med = 1 or psy_med = 1);
	copd_asthma=(copd = 1 or asthma = 1);
	ret_length=time0-dt_start;

run;


data maint_2wkdts_ed;
    set der.maint_2wkdts;
	rename ddose = main_ddose;
    by studyid_rev episode f_maint_dt_end;
run;


proc sql;
create table f_induction as
select e.studyid_rev, e.episode, e.myepisode, e.time0,
	min(p.f_maint_dt_end) format=date9. as f_induction_comp_dt, min(p.main_ddose) as maint_ddose,
	min(p.maint_dur) as maint_dur, min(p.drug_maint) as drug_maint
from tmplib1.TT_analytic_all_m as e left join maint_2wkdts_ed as p
on e.studyid_rev=p.studyid_rev and (e.dt_start <= p.f_maint_dt_end <= e.time0) and p.drug_maint = "M"
group by e.studyid_rev, e.episode, e.myepisode, e.time0;
quit;



data tmplib1.TT_analytic_alls_m;
merge tmplib1.TT_analytic_all_m(in=a) f_induction;
by studyid_rev episode myepisode;
if a; 
induction = f_induction_comp_dt ne .;
run;



data data_ms;
set tmplib1.TT_analytic_alls_m;
by studyid_rev episode myepisode;

s_od_death=(od_death_dt ne .);
if od_death_dt > end_fu_dt then s_od_death = 0;
if f_dt_oop ne . and f_dt_oop >= time0 and f_dt_oop < od_death_dt then s_od_death = 0;

if s_od_death = 0 then do;
   if f_dt_oop >= time0 then t_od_death = min(end_fu_dt, f_dt_oop) - time0;
   else t_od_death = end_fu_dt - time0;
end;
else if s_od_death = 1 then t_od_death = od_death_dt-time0; 

sens_od_death = s_od_death;
time_od_death = min(t_od_death, 365); 
if t_od_death >= 365 then sens_od_death = 0;
run;

run;



data data_ms;
set data_ms;
by studyid_rev episode myepisode;
where dt_start >= '01JAN2010'd and dt_start < '31DEC2022'd and drug_start = "M" and predrug in ("M", "") and ddose_rev ne .;
if trial_gap <= 14;
run;


data tmplib1.TT_analytic_rprev_sens_nw_m;
set data_ms;
 if age_eli=1 and palli_excl=0 and canc_excl=0 and preg_eli=0 and incar_grp_bl ne '2_current';
run;



data data_m1s;
set tmplib1.TT_analytic_alls_m;
by studyid_rev episode;
where dt_start >= '01JAN2010'd and dt_start < '31DEC2022'd and drug_start = "M" and predrug in ("M", "") and ddose_rev ne .;
if first.episode and trial_gap <= 14;
run;


data data_m1s;
set data_m1s;
by studyid_rev episode myepisode;

s_od_death=(od_death_dt ne .);
if od_death_dt > end_fu_dt then s_od_death = 0;
if f_dt_oop ne . and f_dt_oop >= time0 and f_dt_oop < od_death_dt then s_od_death = 0;

if s_od_death = 0 then do;
   if f_dt_oop >= time0 then t_od_death = min(end_fu_dt, f_dt_oop) - time0;
   else t_od_death = end_fu_dt - time0;
end;
else if s_od_death = 1 then t_od_death = od_death_dt-time0; 

sens_od_death = s_od_death;
time_od_death = min(t_od_death, 365); 
if t_od_death >= 365 then sens_od_death = 0;
run;


data tmplib1.TT_analytic_prev_sens_nw_m;
set data_m1s;
 if age_eli=1 and palli_excl=0 and canc_excl=0 and preg_eli=0 and incar_grp_bl ne '2_current';
run;


data tmplib1.TT_analytic_incid_sens_nw_m;
set tmplib1.TT_analytic_prev_sens_m;
if episode = 1;
run;

 
