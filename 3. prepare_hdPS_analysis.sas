
ods html close; /* close previous */
ods html; /* open new */

option compress=yes;

libname der "W:\SUD\Data\2022\derived";
libname lib "W:\SUD\Data\2022\SAS data";
libname tmplib1 "E:\SAS Data\Missed dose protocol_v2";

*******************************************************;
*Incident users*;
*******************************************************;
data TT_analytic_incid_sens_m;
 set tmplib1.TT_analytic_incid_sens_nw_m;
run;



%macro domain_incid_m(input, dt, code, dur_d, maxdigit, output);
  *One year prior to current week;
  proc sql;
    create table domain as
    select dat.studyid_rev, dat.time0, domain_dat.*
    from &input(keep=studyid_rev &dt &code:) as domain_dat,
         TT_analytic_incid_sens_m as dat
    where domain_dat.studyid_rev = dat.studyid_rev and
          &dt >= dat.time0 - &dur_d and &dt < dat.time0
    order by dat.studyid_rev, dat.time0, domain_dat.&dt;
  quit;

  *Create a visit number within the same date;
  data domain;
    set domain;
    by studyid_rev time0 &dt;
    if first.&dt then num = 0;
    num + 1;
  run;

  *Transpose;
  proc transpose data=domain(keep=studyid_rev time0 &dt num &code:)
                 out=&output(rename=(col1=&code) keep=studyid_rev time0 &dt col1);
    by studyid_rev time0 &dt num;
    var &code:;
  run;

  *Truncate;
  data &output;
    set &output;
    &code._trunc = substr(&code, 1, &maxdigit);
    if missing(&code) = 0;
  run;
%mend;




*admin data within the last 6 months (180 days);
*DAD diag (ICD10 up to 4 digits);
%domain_incid_m(input=der.dad_rev, dt=addate_rev, code=diagx, dur_d=180, maxdigit=4, output=dad_diagx_incid_m);
*DAD procedure (CCI up to 5 digits from field 1 to 3);
%domain_incid_m(input=der.dad_rev, dt=addate_rev, code=icode, dur_d=180, maxdigit=5, output=dad_proc_incid_m);
*MSP diag (ICD9 up to 3 digits);
%domain_incid_m(input=der.Msp_virtual, dt=servdt, code=diagcd, dur_d=180, maxdigit=3, output=msp_diag_incid_m);
*MSP fee item (all up to 5 digits);
%domain_incid_m(input=der.Msp_virtual, dt=servdt, code=fitm, dur_d=180, maxdigit=5, output=msp_fee_incid_m);
*PNET drug dispensation (TC5_CD up to 8 digits) in the past 6 months;
%domain_incid_m(input=der.Pnet_oat, dt=srv_date, code=TC5_CD, dur_d=180, maxdigit=8, output=pnet_drug_incid_m);
*get distinct group per day;
proc sql;
create table tmplib1.pnet_drug_rev_incid_m as
select distinct studyid_rev, time0, srv_date, TC5_CD_trunc
from pnet_drug_incid_m;
quit;



*******************************************************;
*  Prevalent new users*;
*******************************************************;
data TT_analytic_prev_sens_m;
 set tmplib1.TT_analytic_prev_sens_nw_m;
run;



%macro domain_prev_m(input, dt, code, dur_d, maxdigit, output);
  *One year prior to current week;
  proc sql;
    create table domain as
    select dat.studyid_rev, dat.time0, domain_dat.*
    from &input(keep=studyid_rev &dt &code:) as domain_dat,
         TT_analytic_prev_sens_m as dat
    where domain_dat.studyid_rev = dat.studyid_rev and
          &dt >= dat.time0 - &dur_d and &dt < dat.time0
    order by dat.studyid_rev, dat.time0, domain_dat.&dt;
  quit;

  *Create a visit number within the same date;
  data domain;
    set domain;
    by studyid_rev time0 &dt;
    if first.&dt then num = 0;
    num + 1;
  run;

  *Transpose;
  proc transpose data=domain(keep=studyid_rev time0 &dt num &code:)
                 out=&output(rename=(col1=&code) keep=studyid_rev time0 &dt col1);
    by studyid_rev time0 &dt num;
    var &code:;
  run;

  *Truncate;
  data &output;
    set &output;
    &code._trunc = substr(&code, 1, &maxdigit);
    if missing(&code) = 0;
  run;
%mend;




*admin data within the last 6 months (180 days);
*DAD diag (ICD10 up to 4 digits);
%domain_prev_m(input=der.dad_rev, dt=addate_rev, code=diagx, dur_d=180, maxdigit=4, output=dad_diagx_prev_m);
*DAD procedure (CCI up to 5 digits from field 1 to 3);
%domain_prev_m(input=der.dad_rev, dt=addate_rev, code=icode, dur_d=180, maxdigit=5, output=dad_proc_prev_m);
*MSP diag (ICD9 up to 3 digits);
%domain_prev_m(input=der.Msp_virtual, dt=servdt, code=diagcd, dur_d=180, maxdigit=3, output=msp_diag_prev_m);
*MSP fee item (all up to 5 digits);
%domain_prev_m(input=der.Msp_virtual, dt=servdt, code=fitm, dur_d=180, maxdigit=5, output=msp_fee_prev_m);
*PNET drug dispensation (TC5_CD up to 8 digits) in the past 6 months;
%domain_prev_m(input=der.Pnet_oat, dt=srv_date, code=TC5_CD, dur_d=180, maxdigit=8, output=pnet_drug_prev_m);
*get distinct group per day;
proc sql;
create table tmplib1.pnet_drug_rev_prev_m as
select distinct studyid_rev, time0, srv_date, TC5_CD_trunc
from pnet_drug_prev_m;
quit;




*******************************************************;
*  All prevalent users*;
*******************************************************;

data TT_analytic_rprev_sens_m;
 set tmplib1.TT_analytic_rprev_sens_nw_m;
run;



%macro domain_rprev_m(input, dt, code, dur_d, maxdigit, output);
  *One year prior to current week;
  proc sql;
    create table domain as
    select dat.studyid_rev, dat.time0, domain_dat.*
    from &input(keep=studyid_rev &dt &code:) as domain_dat,
         TT_analytic_rprev_sens_m as dat
    where domain_dat.studyid_rev = dat.studyid_rev and
          &dt >= dat.time0 - &dur_d and &dt < dat.time0
    order by dat.studyid_rev, dat.time0, domain_dat.&dt;
  quit;

  *Create a visit number within the same date;
  data domain;
    set domain;
    by studyid_rev time0 &dt;
    if first.&dt then num = 0;
    num + 1;
  run;

  *Transpose;
  proc transpose data=domain(keep=studyid_rev time0 &dt num &code:)
                 out=&output(rename=(col1=&code) keep=studyid_rev time0 &dt col1);
    by studyid_rev time0 &dt num;
    var &code:;
  run;

  *Truncate;
  data &output;
    set &output;
    &code._trunc = substr(&code, 1, &maxdigit);
    if missing(&code) = 0;
  run;
%mend;




*admin data within the last 6 months (180 days);
*DAD diag (ICD10 up to 4 digits);
%domain_rprev_m(input=der.dad_rev, dt=addate_rev, code=diagx, dur_d=180, maxdigit=4, output=dad_diagx_rprev_m);
*DAD procedure (CCI up to 5 digits from field 1 to 3);
%domain_rprev_m(input=der.dad_rev, dt=addate_rev, code=icode, dur_d=180, maxdigit=5, output=dad_proc_rprev_m);
*MSP diag (ICD9 up to 3 digits);
%domain_rprev_m(input=der.Msp_virtual, dt=servdt, code=diagcd, dur_d=180, maxdigit=3, output=msp_diag_rprev_m);
*MSP fee item (all up to 5 digits);
%domain_rprev_m(input=der.Msp_virtual, dt=servdt, code=fitm, dur_d=180, maxdigit=5, output=msp_fee_rprev_m);
*PNET drug dispensation (TC5_CD up to 8 digits) in the past 6 months;
%domain_rprev_m(input=der.Pnet_oat, dt=srv_date, code=TC5_CD, dur_d=180, maxdigit=8, output=pnet_drug);
*get distinct group per day;
proc sql;
create table tmplib1.pnet_drug_rev_rprev_m as
select distinct studyid_rev, time0, srv_date, TC5_CD_trunc
from pnet_drug;
quit;



/*
For bnx
*/

data TT_analytic_rprev_sens_b;
 set tmplib1.TT_analytic_rprev_sens_nw_b;
run;


*obtain admin data before the each week (for &dur_d) &code up to &maxdigit;

%macro domain_bnx(input, dt, code, dur_d, maxdigit, output);
  *One year prior to current week;
  proc sql;
    create table domain as
    select dat.studyid_rev, dat.time0, domain_dat.*
    from &input(keep=studyid_rev &dt &code:) as domain_dat,
         TT_analytic_rprev_sens_b as dat
    where domain_dat.studyid_rev = dat.studyid_rev and
          &dt >= dat.time0 - &dur_d and &dt < dat.time0
    order by dat.studyid_rev, dat.time0, domain_dat.&dt;
  quit;

  *Create a visit number within the same date;
  data domain;
    set domain;
    by studyid_rev time0 &dt;
    if first.&dt then num = 0;
    num + 1;
  run;

  *Transpose;
  proc transpose data=domain(keep=studyid_rev time0 &dt num &code:)
                 out=&output(rename=(col1=&code) keep=studyid_rev time0 &dt col1);
    by studyid_rev time0 &dt num;
    var &code:;
  run;

  *Truncate;
  data &output;
    set &output;
    &code._trunc = substr(&code, 1, &maxdigit);
    if missing(&code) = 0;
  run;
%mend;




*admin data within the last 6 months (180 days);
*DAD diag (ICD10 up to 4 digits);
%domain_bnx(input=der.dad_rev, dt=addate_rev, code=diagx, dur_d=180, maxdigit=4, output=dad_diagx_b);
*DAD procedure (CCI up to 5 digits from field 1 to 3);
%domain_bnx(input=der.dad_rev, dt=addate_rev, code=icode, dur_d=180, maxdigit=5, output=dad_proc_b);
*MSP diag (ICD9 up to 3 digits);
%domain_bnx(input=der.Msp_virtual, dt=servdt, code=diagcd, dur_d=180, maxdigit=3, output=msp_diag_b);
*MSP fee item (all up to 5 digits);
%domain_bnx(input=der.Msp_virtual, dt=servdt, code=fitm, dur_d=180, maxdigit=5, output=msp_fee_b);
*PNET drug dispensation (TC5_CD up to 8 digits) in the past 6 months;
%domain_bnx(input=der.Pnet_oat, dt=srv_date, code=TC5_CD, dur_d=180, maxdigit=8, output=pnet_drug_b);
*get distinct group per day;
proc sql;
create table tmplib1.pnet_drug_rev_b as
select distinct studyid_rev, time0, srv_date, TC5_CD_trunc
from pnet_drug_b;
quit;





*******************************************************;
*  All prevalent users*;
*******************************************************;

data TT_analytic_rprev_sens_b;
 set tmplib1.TT_analytic_rprev_sens_nw_b;
run;



%macro domain_rprev_b(input, dt, code, dur_d, maxdigit, output);
  *One year prior to current week;
  proc sql;
    create table domain as
    select dat.studyid_rev, dat.time0, domain_dat.*
    from &input(keep=studyid_rev &dt &code:) as domain_dat,
         TT_analytic_rprev_sens_b as dat
    where domain_dat.studyid_rev = dat.studyid_rev and
          &dt >= dat.time0 - &dur_d and &dt < dat.time0
    order by dat.studyid_rev, dat.time0, domain_dat.&dt;
  quit;

  *Create a visit number within the same date;
  data domain;
    set domain;
    by studyid_rev time0 &dt;
    if first.&dt then num = 0;
    num + 1;
  run;

  *Transpose;
  proc transpose data=domain(keep=studyid_rev time0 &dt num &code:)
                 out=&output(rename=(col1=&code) keep=studyid_rev time0 &dt col1);
    by studyid_rev time0 &dt num;
    var &code:;
  run;

  *Truncate;
  data &output;
    set &output;
    &code._trunc = substr(&code, 1, &maxdigit);
    if missing(&code) = 0;
  run;
%mend;




*admin data within the last 6 months (180 days);
*DAD diag (ICD10 up to 4 digits);
%domain_rprev_b(input=der.dad_rev, dt=addate_rev, code=diagx, dur_d=180, maxdigit=4, output=dad_diagx_rprev_b);
*DAD procedure (CCI up to 5 digits from field 1 to 3);
%domain_rprev_b(input=der.dad_rev, dt=addate_rev, code=icode, dur_d=180, maxdigit=5, output=dad_proc_rprev_b);
*MSP diag (ICD9 up to 3 digits);
%domain_rprev_b(input=der.Msp_virtual, dt=servdt, code=diagcd, dur_d=180, maxdigit=3, output=msp_diag_rprev_b);
*MSP fee item (all up to 5 digits);
%domain_rprev_b(input=der.Msp_virtual, dt=servdt, code=fitm, dur_d=180, maxdigit=5, output=msp_fee_rprev_b);
*PNET drug dispensation (TC5_CD up to 8 digits) in the past 6 months;
%domain_rprev_b(input=der.Pnet_oat, dt=srv_date, code=TC5_CD, dur_d=180, maxdigit=8, output=pnet_drug);
*get distinct group per day;
proc sql;
create table tmplib1.pnet_drug_rev_rprev_b as
select distinct studyid_rev, time0, srv_date, TC5_CD_trunc
from pnet_drug;
quit;
