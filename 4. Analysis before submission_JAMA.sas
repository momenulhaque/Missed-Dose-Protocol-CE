options compress=y;

libname der "W:\SUD\Data\2022\derived";
libname lib "W:\SUD\Data\2022\SAS data";
libname pss "W:\SUD\Data\2022\derived\PSS";
libname tmplib "E:\SAS Data\Missed dose protocol";
libname tmplib1 "E:\SAS Data\Missed dose protocol_v1";
libname tmplib2 "E:\SAS Data\Missed dose protocol_v2";

data dat1_rprev;
    set tmplib2.TT_analytic_rprev_sens_nw_m;

    length ret_length_cat $5;
    if ret_length <= 30 then ret_length_cat = "0-30";
    else if 31 <= ret_length <= 90 then ret_length_cat = "31-90";
    else if ret_length > 90 then ret_length_cat = ">90";

    /* Step 2: Create high_ddose */
    length high_ddose $3;
    if ddose_prev >= 80 then high_ddose = ">=80";
    else high_ddose = "<80";
run;

/* Step 3: Filter rows */
data dat1_rprev1;
    set dat1_rprev;
    if trial_gap in (1,2);
run;

data dat1_rprev1;
    set dat1_rprev1;
    trt = trt_p1; 
    if trt in ('no reversion','reversion','increase'); 
run;

proc freq data=dat1_rprev1; table trt; run;

proc means data=dat1_rprev1; vars prior_gap; run;



proc sort data=dat1_rprev1;
    by studyid_rev episode;
run;




data dat1_rprev1;
    set dat1_rprev1;
    by studyid_rev episode;

    retain prior_gap 0 prior_gap_num 0;

    if first.episode then do;
        prior_gap = 0;
        prior_gap_num = 0;
    end;

    prior_gap + gap_days;       /* cumulative sum */
    prior_gap_num + 1;          /* row number within group */
run;

proc logistic data=dat1_rprev1 noprint;
  class trt (ref="no reversion" order=data) sexM age_lt30 age_30_39 age_40_49 age_ge50 rural 
                po_preOUD phy_att_gr2 social_assist home_grp_b
                incar_grp_bl drug_acute cci_gt0 copd_asthma sed_psy_med
                po sud aud smd hcv tob udt ddose_prev clin_size_b virtual_care pss/ param=ref;
  model trt = sexM age_st age_lt30 age_30_39 age_40_49 age_ge50 rural 
                year po_preOUD phy_att_gr2 social_assist ret_length home_grp_b 
                incar_grp_bl drug_acute cci_gt0 copd_asthma sed_psy_med
                po sud aud smd hcv tob udt  clin_size_b virtual_care pss /* ddose_prev prior_gap *// link=glogit;
  output out=tx_d predprobs=i;
run;


proc logistic data=dat1_rprev1 noprint;
  class trt (ref="no reversion" order=data)/ param=ref;
  model trt = / link=glogit;
  output out=tx_n predprobs=i;
run;



  proc sort data=tx_d; by studyid_rev episode time0; run;
  proc sort data=tx_n; by studyid_rev episode time0; run;

 data _tx0;
  merge 
      tx_d(in=a 
           keep=studyid_rev episode time0 trt IP_: 
           rename=(IP_no_reversion=pT_IP_no_reversion_d 
                   IP_reversion=pT_reversion_d
                   IP_increase=pT_increase_d))
      tx_n(in=b 
           keep=studyid_rev episode time0 IP_: 
           rename=(IP_no_reversion=pT_IP_no_reversion_n 
                   IP_reversion=pT_reversion_n
                   IP_increase=pT_increase_n))
  ;   /* <-- THIS SEMICOLON WAS MISSING */

  by studyid_rev episode time0;

  if a;
  keep studyid_rev episode time0 pT_: trt;
run;


  data _tx_all;
    merge dat1_rprev1(in=a) _tx0(in=beg);
    by studyid_rev episode time0;
    if a;

    length ps_d ps_n wt_tx 8.;
    if trt="no reversion" then do; ps_d=pT_IP_no_reversion_d; ps_n=pT_IP_no_reversion_n; end;
    else if trt="reversion" then do; ps_d=pT_reversion_d; ps_n=pT_reversion_n; end;
    else if trt="increase" then do; ps_d=pT_increase_d; ps_n=pT_increase_n; end;
    else do; ps_d=.; ps_n=.; end;

    wt_tx = ps_n / ps_d;
  run;

  proc univariate data=_tx_all noprint;
    var wt_tx;
    output out=_tx_stat pctlpts=0.01 0.99 pctlpre=p pctlname=low high;
  run;

  data _null_;
    set _tx_stat;
    call symputx('tx_p1', plow);
    call symputx('tx_p99', phigh);
  run;

  data _tx_all;
    set _tx_all;
    wt_tx_trunc = min(max(wt_tx, &tx_p1), &tx_p99);
  run;

proc sort data=_tx_all;
    by trt;
run;

proc means data=_tx_all;
    by trt;
    var wt_tx wt_tx_trunc;
run;


data long_dat_od;
    set _tx_all;
    interval_count = ceil(time_od_acute / 14);
    do interval = 1 to interval_count;
        start = (interval-1)*14;
        stop  = interval*14;
        event = (interval = interval_count) * sens_od_acute;
        output;
    end;
run;


proc genmod data=long_dat_od descending;
    class trt (ref="no reversion") studyid_rev;
    model event = trt year ret_length ddose_prev interval
        / dist=bin link=logit;
    weight wt_tx;
    repeated subject=studyid_rev / type=ind;

    ods output ParameterEstimates=beta
               CovB=covb;
run;











%macro predict_cuminc_ow(
    Dat_expanded=,
    trt_var=,
    covs=,
    horizons=,
    cluster_var=,
    samples=1000,
    conf_int=1,
    ref_group=no reversion,
    wt_type=wt_trunc
);

%local i h n_horizons horizon_list;

/* Split horizons into macro variables */
%let n_horizons = %sysfunc(countw(&horizons, %str( )));
%let horizon_list = &horizons;

/* 1. Fit weighted pooled logistic regression with cluster-robust SE */
proc genmod data=&Dat_expanded descending;
    class &trt_var / param=ref ref=first;
    model event = &trt_var &covs ns(interval, df=3) / dist=bin link=logit;
    weight &wt_type;
    repeated subject=&cluster_var / type=ind; /* robust SE */
    ods output ParameterEstimates=coefs;
run;

/* 2. Extract covariance matrix for simulation */
proc iml;
    use coefs;
    read all var {Estimate} into beta;
    close coefs;

    /* For simplicity, use identity for covariance if not available */
    cov = diag( repeat(0.01, nrow(beta)) );

    /* Simulate coefficients if CI requested */
    %if &conf_int %then %do;
        call randseed(12345);
        call randnormal(samples, sim_betas, beta, cov);
    %end;
    %else %do;
        sim_betas = beta`;
    %end;

    /* Prepare for horizon computation */
    use &Dat_expanded;
    read all var {&trt_var &cluster_var interval stop event &wt_type &covs} into X[colname=names];
    close &Dat_expanded;

    n_horizons = &n_horizons;
    horizons = {&horizon_list};

    results = {};
    contrasts = {};

    do i = 1 to nrow(sim_betas);
        do h = 1 to n_horizons;
            horizon = horizons[h];

            /* select rows with stop <= horizon */
            idx = loc(X[, "stop"] <= horizon);
            Xh = X[idx,];

            trt_levels = unique(Xh[, "&trt_var"]);
            do t = 1 to ncol(trt_levels);
                trt = trt_levels[t];

                Xh_t = Xh;
                Xh_t[, "&trt_var"] = trt;

                /* Design matrix (intercept + covariates) */
                mm = j(nrow(Xh_t), ncol(sim_betas), 1);
                /* Here you would fill mm with proper covariates */
                /* Simplified: assume Xh_t columns match beta order */

                lp = mm * sim_betas[i,]`;
                p = 1 / (1 + exp(-lp));

                /* Compute cluster survival and risk */
                clusters = unique(Xh_t[, "&cluster_var"]);
                surv = j(ncol(clusters), 1, 1);

                do c = 1 to ncol(clusters);
                    cid = clusters[c];
                    idx_c = loc(Xh_t[, "&cluster_var"] = cid);
                    surv[c] = prod(1 - p[idx_c]);
                end;

                risk = 1 - mean(surv);
                results = results // (trt || horizon || risk);
            end;
        end;
    end;

    /* Convert results to dataset */
    create risks from results[colname={"trt" "horizon" "risk"}];
    append from results;
    close risks;
quit;

/* 3. Optional: Compute RD and RR relative to reference group */
%if &conf_int %then %do;
    proc sql;
        create table contrasts as
        select a.trt, a.horizon,
               a.risk - b.risk as RD,
               a.risk / b.risk as RR
        from risks a
        left join risks b
        on a.horizon = b.horizon and b.trt = "&ref_group"
        where a.trt ne "&ref_group";
    quit;
%end;

%mend predict_cuminc_ow;






/* Example usage */
%predict_cuminc_ow(
    Dat_expanded=Dat_expanded,
    trt_var=trt,
    covs=x1 x2 x3,
    horizons=30 60 90,
    cluster_var=cluster_id,
    samples=1000,
    conf_int=1,
    ref_group=no reversion,
    wt_type=wt_trunc
);
