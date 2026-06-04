proc datasets kill;
quit;

options compress=y;

libname der "W:\SUD\Data\2022\derived";
libname lib "W:\SUD\Data\2022\SAS data";
libname pss "W:\SUD\Data\2022\derived\PSS";
libname tmplib "E:\SAS Data\Missed dose protocol_v2";



/*Creating IV variables */
proc sql;
create table iv_cohort_alls_m as
select a.*, b.com, b.comsize, b.n_pracs, b.n_pats 
from tmplib.TT_analytic_alls_m a left join der.prov_clinics_v2 b
on year(a.time0)=b.year and month(a.time0)=b.month and a.prac_studyid_rev=b.prac_studyid_rev
order by a.studyid_rev, a.time0;
quit; *395201;

data iv_cohort_alls_m;
set iv_cohort_alls_m;
trt_p23=trt_s456;
if trt_s456 in ("1-<25%", "25-<50%") then trt_p23 = "1-<50%";
incar_grp_bls=(incar_grp_bl="1_pastyear");
grp1 = (trt_s456="no reversion");
grp2 = (trt_s456="1-<50%" or trt_s456=">=50%/<=40mg");
grp3 = (trt_s456="1-<50%");
grp4 = (trt_s456=">=50%/<=40mg");
grp5 = (trt_s456="increase");
run;


data cohort_alls_m;
set iv_cohort_alls_m;
by studyid_rev episode; 
if first.episode;
if trial_gap > 2;
run;

proc freq data=cohort_alls_m; table trt_s456; run;

data dat_prev_m;
set tmplib.TT_analytic_prev_sens_nw_m;
time0_12m_prior=INTNX('MONTH', time0, -12, "same");
format time0_12m_prior date9.;
trt_p23=trt_s456;
if trt_s456 in ("1-<25%", "25-<50%") then trt_p23 = "1-<50%";
incar_grp_bls=(incar_grp_bl="1_pastyear");
grp1 = (trt_s456="no reversion");
grp2 = (trt_s456="1-<50%" or trt_s456=">=50%/<=40mg");
grp3 = (trt_s456="1-<50%");
grp4 = (trt_s456=">=50%/<=40mg");
grp5 = (trt_s456="increase");
if trial_gap > 2;
run; *10871;



/*Prescriber level preference*/
%macro praciv_prep_met(cond, trt_name, k, trial, cohort);

proc sql;
create table prep_&cohort._&k. as
select epi.prac_studyid_rev, epi.time0, 
	count(distinct catx(",", oat.studyid_rev, oat.time0)) as n_eps,

	count(distinct case when oat.&trt_name.="no reversion" then catx("", oat.studyid_rev, oat.time0) else "" end) as grp1_eps,
	(calculated grp1_eps/calculated n_eps*100) as grp1_eps_pct format 8.1,

	count(distinct case when oat.&trt_name.="1-<50%" or oat.&trt_name.=">=50%/<=40mg" then catx("", oat.studyid_rev, oat.time0) else "" end) as grp2_eps,
	(calculated grp2_eps/calculated n_eps*100) as grp2_eps_pct format 8.1,

	count(distinct case when oat.&trt_name.="1-<50%"  then catx("", oat.studyid_rev, oat.time0) else "" end) as grp3_eps,
	(calculated grp3_eps/calculated n_eps*100) as grp3_eps_pct format 8.1,

	count(distinct case when oat.&trt_name.=">=50%/<=40mg"  then catx("", oat.studyid_rev, oat.time0) else "" end) as grp4_eps,
	(calculated grp4_eps/calculated n_eps*100) as grp4_eps_pct format 8.1,

	count(distinct case when oat.&trt_name.="increase"  then catx("", oat.studyid_rev, oat.time0) else "" end) as grp5_eps,
	(calculated grp5_eps/calculated n_eps*100) as grp5_eps_pct format 8.1


from dat_&cohort._&k. (where=(&cond.)) as epi left join iv_cohort_alls_&k.(where=(&cond.)) as oat
on epi.prac_studyid_rev=oat.prac_studyid_rev and epi.time0_12m_prior<=oat.time0<epi.time0 
group by epi.prac_studyid_rev, epi.time0
order by epi.prac_studyid_rev, epi.time0;
quit; 

proc sort data=dat_&cohort._&k.(where=(&cond.)) out=dat;
by  prac_studyid_rev time0;
run;

data praciv_&cohort._&k._tr&trial.;
merge dat prep_&cohort._&k.;
by prac_studyid_rev time0;
run;
%mend;





%macro praciv_prep_met_fixed(cond, trt_name, k, trial, cohort);

proc sql;
create table praciv_&cohort._&k._tr&trial. as
select 
    d.*,

    /* Total episodes in past 12 months by this prescriber */
    (select count(*) 
     from iv_cohort_alls_&k. oat
     where oat.prac_studyid_rev = d.prac_studyid_rev
       and oat.time0 >= d.time0_12m_prior
       and oat.time0 < d.time0
       and (&cond.)) as n_eps,

    /* Count of each treatment group in past 12 months */
    (select count(*) 
     from iv_cohort_alls_&k. oat
     where oat.prac_studyid_rev = d.prac_studyid_rev
       and oat.time0 >= d.time0_12m_prior
       and oat.time0 < d.time0
       and oat.&trt_name = "no reversion"
       and (&cond.)) as grp1_eps,

    (calculated grp1_eps / calculated n_eps * 100) as grp1_eps_pct format 8.1,

    (select count(*) 
     from iv_cohort_alls_&k. oat
     where oat.prac_studyid_rev = d.prac_studyid_rev
       and oat.time0 >= d.time0_12m_prior
       and oat.time0 < d.time0
       and oat.&trt_name in ("1-<50%", ">=50%/<=40mg")
       and (&cond.)) as grp2_eps,

    (calculated grp2_eps / calculated n_eps * 100) as grp2_eps_pct format 8.1,

    (select count(*) 
     from iv_cohort_alls_&k. oat
     where oat.prac_studyid_rev = d.prac_studyid_rev
       and oat.time0 >= d.time0_12m_prior
       and oat.time0 < d.time0
       and oat.&trt_name = "1-<50%"
       and (&cond.)) as grp3_eps,

    (calculated grp3_eps / calculated n_eps * 100) as grp3_eps_pct format 8.1,

    (select count(*) 
     from iv_cohort_alls_&k. oat
     where oat.prac_studyid_rev = d.prac_studyid_rev
       and oat.time0 >= d.time0_12m_prior
       and oat.time0 < d.time0
       and oat.&trt_name = ">=50%/<=40mg"
       and (&cond.)) as grp4_eps,

    (calculated grp4_eps / calculated n_eps * 100) as grp4_eps_pct format 8.1,

    (select count(*) 
     from iv_cohort_alls_&k. oat
     where oat.prac_studyid_rev = d.prac_studyid_rev
       and oat.time0 >= d.time0_12m_prior
       and oat.time0 < d.time0
       and oat.&trt_name = "increase"
       and (&cond.)) as grp5_eps,

    (calculated grp5_eps / calculated n_eps * 100) as grp5_eps_pct format 8.1

from dat_&cohort._&k. d
where &cond.
order by prac_studyid_rev, time0;
quit;

%mend;



%praciv_prep_met(cond=trial_cat="1-3", trt_name=trt_p23, trial=1, k=m, cohort=incid)
%praciv_prep_met(cond=trial_cat="1-3", trt_name=trt_p23, trial=1, k=m, cohort=prev)
%praciv_prep_met(cond=trial_cat="1-3", trt_name=trt_p23, trial=1, k=m, cohort=rprev) 


%praciv_prep_met(cond=trial_cat="4", trt_name=trt_p23, trial=2, k=m, cohort=incid)
%praciv_prep_met(cond=trial_cat="4", trt_name=trt_p23, trial=2, k=m, cohort=prev)
%praciv_prep_met(cond=trial_cat="4", trt_name=trt_p23, trial=2, k=m, cohort=rprev)


%praciv_prep_met(cond=trial_cat="5-14", trt_name=trt_p23, trial=3, k=m, cohort=incid) 
%praciv_prep_met(cond=trial_cat="5-14", trt_name=trt_p23, trial=3, k=m, cohort=prev)
%praciv_prep_met(cond=trial_cat="5-14", trt_name=trt_p23, trial=3, k=m, cohort=rprev)


/*facility level preference*/

proc sql;
create table dat_incid_b_v2 as
select a.*, b.com, b.comsize, b.n_pracs, b.n_pats 
from dat_incid_b as a left join der.prov_clinics_v2 as b
on year(a.time0)=b.year and month(a.time0)=b.month and a.prac_studyid_rev=b.prac_studyid_rev
order by a.studyid_rev;
quit; *10871;


proc sql;
create table dat_prev_b_v2 as
select a.*, b.com, b.comsize, b.n_pracs, b.n_pats 
from dat_prev_b as a left join der.prov_clinics_v2 as b
on year(a.time0)=b.year and month(a.time0)=b.month and a.prac_studyid_rev=b.prac_studyid_rev
order by a.studyid_rev, a.time0;
quit; *119886;


proc sql;
create table dat_rprev_b_v2 as
select a.*, b.com, b.comsize, b.n_pracs, b.n_pats 
from dat_rprev_b as a left join der.prov_clinics_v2 as b
on year(a.time0)=b.year and month(a.time0)=b.month and a.prac_studyid_rev=b.prac_studyid_rev
order by a.studyid_rev, a.time0;
quit; *119886;



%macro faciiv_prep_bnx(cond, trt_name, k, trial, cohort);
proc sql;
create table prep_&cohort._&k._v2 as
select epi.com, epi.time0, 
	count(distinct catx(",", oat.studyid_rev, oat.time0)) as n_eps,

	count(distinct case when oat.&trt_name.="no reversion" then catx("", oat.studyid_rev, oat.time0) else "" end) as grp1_eps,
	(calculated grp1_eps/calculated n_eps*100) as grp1_eps_pct format 8.1,

	count(distinct case when oat.&trt_name.="1-<50%" or oat.&trt_name.=">=50%/<=40mg" then catx("", oat.studyid_rev, oat.time0) else "" end) as grp2_eps,
	(calculated grp2_eps/calculated n_eps*100) as grp2_eps_pct format 8.1,

	count(distinct case when oat.&trt_name.="1-<50%"  then catx("", oat.studyid_rev, oat.time0) else "" end) as grp3_eps,
	(calculated grp3_eps/calculated n_eps*100) as grp3_eps_pct format 8.1,

	count(distinct case when oat.&trt_name.=">=50%/<=40mg"  then catx("", oat.studyid_rev, oat.time0) else "" end) as grp4_eps,
	(calculated grp4_eps/calculated n_eps*100) as grp4_eps_pct format 8.1,

	count(distinct case when oat.&trt_name.="increase"  then catx("", oat.studyid_rev, oat.time0) else "" end) as grp5_eps,
	(calculated grp5_eps/calculated n_eps*100) as grp5_eps_pct format 8.1

from dat_&cohort._&k._v2 (where=(&cond.)) as epi left join iv_cohort_alls_&k.(where=(&cond.)) as oat
on epi.com=oat.com and epi.time0_12m_prior<=oat.time0<epi.time0 
group by epi.com, epi.time0
order by epi.com, epi.time0;
quit; 

proc sort data=dat_&cohort._&k._v2(where=(&cond.)) out=dat;
by com time0;
run;

data faciiv_&cohort._&k._tr&trial.;
merge dat prep_&cohort._&k._v2;
by com time0;
run;
%mend;


%faciiv_prep_bnx(cond=trial_cat="1-5", trt_name=trt_p45, trial=4, k=b, cohort=incid)
%faciiv_prep_bnx(cond=trial_cat="1-5", trt_name=trt_p45, trial=4, k=b, cohort=prev)
%faciiv_prep_bnx(cond=trial_cat="1-5", trt_name=trt_p45, trial=4, k=b, cohort=rprev)
 
%faciiv_prep_bnx(cond=trial_cat="6-14", trt_name=trt_p45, trial=5, k=b, cohort=incid)
%faciiv_prep_bnx(cond=trial_cat="6-14", trt_name=trt_p45, trial=5, k=b, cohort=prev)
%faciiv_prep_bnx(cond=trial_cat="6-14", trt_name=trt_p45, trial=5, k=b, cohort=rprev)




