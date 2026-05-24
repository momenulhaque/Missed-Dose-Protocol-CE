
# Figure 2: Panel A ######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
res_df <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs_iptw_ow.rds")

res_df2 <- res_df %>% filter(Design != "Prevalent new users") %>% 
  mutate(
    Comparison = recode(Comparison, 
                        "Dose reduction\nvs\nNo dose change" = "Dose reduction vs no change",
                        "Dose increase\nvs\nNo dose change" = "Dose increase vs no change",
                        "1-<50% reduction\nvs\nNo dose change" = "Dose reduction of <50% vs no change",
                        ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reduction of >=50%/initial dose* vs no change",
                        ">=50% reduction/<=4mg\nvs\nNo dose change" = "Dose reduction of >=50%/initial dose* vs no change"),
    Design = recode(Design,
                    "All prevalent users" = "Prevalent new users")
  ) %>% arrange(res_type, Model, Outcome, Trial, Comparison, Design)



res_df2 %>% count(res_type)

#### HRs- met forestplot ##########################

#### Trial 1: od: iptw #######
header1 <- data.frame(Comparison = "Trial 1: After MET interruption of 1-2 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr1_full_iptw_od_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                           Model == "IPTW",
                                           Outcome == "Overdose-related acute care visits",
                                           Trial == "Trial 1: MET interruption of 1-2 days",
                                           Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr1_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 1: MET interruption of 1-2 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr1_full_iptw_od <- bind_rows(header1, tr1_full_iptw_od_red, tr1_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: od: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 2: MET interruption of 3-4 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))



tr2_full_iptw_od <- bind_rows(header2, tr2_full_iptw_od_redlt50, tr2_full_iptw_od_redgt50, tr2_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)
  ) %>%
  ungroup()

#### Trial 3: od: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 3: MET interruption of 5-14 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr3_full_iptw_od <- bind_rows(header3, tr3_full_iptw_od_redlt50, tr3_full_iptw_od_redgt50, tr3_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()



format_num <- function(x) {
  ifelse(
    abs(x) >= 10,
    sprintf("%.1f", x),
    sprintf("%.2f", x)
  )
}

met_full_iptw_od <- bind_rows(tr1_full_iptw_od,
                              tr2_full_iptw_od,
                              tr3_full_iptw_od) %>%
  mutate(
    HR_display = dplyr::case_when(
      is.na(HR) ~ " ",
      TRUE ~ sprintf(
        "%s (%s, %s)",
        format_num(HR),
        format_num(lower),
        format_num(upper)
      )
    )
  )




met_full_iptw_od$`overdose-related acute care visits` <- "" 
met_full_iptw_od$`Adjusted hazard\nratio (95% CI)` <- met_full_iptw_od$HR_display
met_full_iptw_od$`Dose adjustment after interruption` <- met_full_iptw_od$Comparison
met_full_iptw_od$`User design` <- met_full_iptw_od$Design


#### Risk- met forestplot ###############
extract_365 <- function(res_obj, type = c("od", "dis")) {
  type <- match.arg(type)
  
  pops <- names(res_obj$points_ci[[type]])
  
  out <- lapply(pops, function(p) {
    pop_obj <- res_obj$points_ci[[type]][[p]]
    
    # skip if NULL, NA, or not a list
    if (is.null(pop_obj) || all(is.na(pop_obj)) || !is.list(pop_obj))
      return(NULL)
    
    # skip if contrasts is missing or NULL
    if (is.null(pop_obj$contrasts) || all(is.na(pop_obj$contrasts)))
      return(NULL)
    
    df <- pop_obj$contrasts
    
    # skip if empty or horizon column missing
    if (nrow(df) == 0 || !("horizon" %in% names(df)))
      return(NULL)
    
    df %>%
      dplyr::filter(horizon == 365) %>%
      dplyr::mutate(pop = p, type = type)
  })
  
  # bind only non-null results
  dplyr::bind_rows(out)
}


res_1 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk1_12d_fixed_iptw.rds")
res_2 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk2_34d_iptw.rds")
res_3 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk3_iptw.rds")
res_4 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk4_iptw.rds")
res_5 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk5_iptw.rds")

#### Trial 1: 1-2 days ##############
header1 <- data.frame(Comparison = "Trial 1: After MET interruption of 1-2 days",
                      Design = "",
                      RD = NA, RD_lower = NA, RD_upper = NA)

tr1_risk_iptw_od_red <- extract_365(res_1, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "reversion" = "Dose reduction vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  filter(Design != "Prevalent new users") %>% 
  mutate(Design = recode(Design, "Prevalent new users" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose reduction vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr1_risk_iptw_od_incr <- extract_365(res_1, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "reversion" = "Dose reduction vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  filter(Design != "Prevalent new users") %>% 
  mutate(Design = recode(Design, "Prevalent new users" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose increase vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr1_risk_iptw_od <- bind_rows(header1, tr1_risk_iptw_od_red, tr1_risk_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()




#### Trial 2: 3-4 days #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      RD = NA, RD_lower = NA, RD_upper = NA)
lt50_risk <- extract_365(res_2, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction of <50% vs no change",
                             ">=50%/<=40mg" = "Dose reduction of >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper)

tr2_risk_iptw_od_redlt50 <- bind_rows(lt50_risk, data.frame(
  Comparison = c("Dose reduction of <50% vs no change",
                 "Dose reduction of >=50%/initial dose* vs no change",
                 "Dose increase vs no change"),
  Design = "Incident users**",
  RD = NA, RD_lower = NA, RD_upper = NA
  
)) %>% mutate(Design = factor(Design, levels = c("Incident users**", "Prevalent new users",
                                                 "All prevalent users"))) %>% 
  arrange(Design) %>% 
  filter(Comparison == "Dose reduction of <50% vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr2_risk_iptw_od_redgt50 <- bind_rows(lt50_risk, data.frame(
  Comparison = c("Dose reduction of <50% vs no change",
                 "Dose reduction of >=50%/initial dose* vs no change",
                 "Dose increase vs no change"),
  Design = "Incident users**",
  RD = NA, RD_lower = NA, RD_upper = NA
  
)) %>% mutate(Design = factor(Design, levels = c("Incident users**", "Prevalent new users",
                                                 "All prevalent users"))) %>% 
  arrange(Design) %>% 
  filter(Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr2_risk_iptw_od_incr <- bind_rows(lt50_risk, data.frame(
  Comparison = c("Dose reduction of <50% vs no change",
                 "Dose reduction of >=50%/initial dose* vs no change",
                 "Dose increase vs no change"),
  Design = "Incident users**",
  RD = NA, RD_lower = NA, RD_upper = NA
  
)) %>% mutate(Design = factor(Design, levels = c("Incident users**", "Prevalent new users",
                                                 "All prevalent users"))) %>% 
  arrange(Design) %>% 
  filter(Comparison == "Dose increase vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr2_risk_iptw_od <- bind_rows(header2, 
                              tr2_risk_iptw_od_redlt50, 
                              tr2_risk_iptw_od_redgt50, 
                              tr2_risk_iptw_od_incr) %>%
  filter(Design != "Prevalent new users") %>% 
  mutate(Design = recode(Design, "All prevalent users"="Prevalent new users")) %>% 
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")
  ) %>%
  ungroup()




#### Trial 3  #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      RD = NA, RD_lower = NA, RD_upper = NA)
tr3_risk_iptw_od_redlt50 <- extract_365(res_3, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction of <50% vs no change",
                             ">=50%/<=40mg" = "Dose reduction of >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose reduction of <50% vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr3_risk_iptw_od_redgt50 <- extract_365(res_3, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction of <50% vs no change",
                             ">=50%/<=40mg" = "Dose reduction of >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr3_risk_iptw_od_incr <- extract_365(res_3, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction of <50% vs no change",
                             ">=50%/<=40mg" = "Dose reduction of >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose increase vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr3_risk_iptw_od <- bind_rows(header3, 
                              tr3_risk_iptw_od_redlt50, 
                              tr3_risk_iptw_od_redgt50, 
                              tr3_risk_iptw_od_incr) %>%
  filter(Design != "Prevalent new users") %>% 
  mutate(Design = recode(Design, "All prevalent users"="Prevalent new users")) %>% 
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")
  ) %>%
  ungroup()


format_num <- function(x) {
  ifelse(
    abs(x) >= 10,
    sprintf("%.1f", x),
    sprintf("%.2f", x)
  )
}
met_risk_iptw_od <- bind_rows(
  tr1_risk_iptw_od,
  tr2_risk_iptw_od,
  tr3_risk_iptw_od
) %>%
  mutate(
    RD_display = dplyr::case_when(
      is.na(RD) ~ " ",
      TRUE ~ sprintf(
        "%s (%s, %s)",
        format_num(RD * 100),
        format_num(RD_lower * 100),
        format_num(RD_upper * 100)
      )
    )
  )



met_risk_iptw_od$`overdose-related acute care visits` <- "" 
met_risk_iptw_od$`Cumulative incidence difference (95% CI)%` <- met_risk_iptw_od$RD_display
met_risk_iptw_od$`Dose adjustment after interruption` <- met_risk_iptw_od$Comparison
met_risk_iptw_od$`User design` <- met_risk_iptw_od$Design

#### creating the graph ######

met_risk_iptw_od_v1 <- met_risk_iptw_od %>% select(`Cumulative incidence difference (95% CI)%`, 
                                                   `overdose-related acute care visits`) %>% 
  rename(`% difference in cumulative\nincidence at 12 months (95% CI)`=`Cumulative incidence difference (95% CI)%`)
met_full_risk_iptw <- bind_cols(met_full_iptw_od, met_risk_iptw_od_v1)

names(met_full_risk_iptw)[c(7, 12)] <- c("                             ", 
                                         "                         ")

met_full_risk_iptw <- met_full_risk_iptw %>% 
  mutate(`Dose adjustment after interruption`= 
           recode(`Dose adjustment after interruption`,
                  "   Dose reduction of >=50%/initial dose* vs no change"=
                    "   Dose reduction of >=50%/initial dose* vs no change "))


library(forestploter)
library(gtable)


p_od_dis_met <- forest(
  met_full_risk_iptw[, c(9,10,8,7,11,12)],
  est   = list(met_full_risk_iptw$HR, met_risk_iptw_od$RD*100),
  lower = list(met_full_iptw_od$lower, met_risk_iptw_od$RD_lower*100),
  upper = list(met_full_iptw_od$upper, met_risk_iptw_od$RD_upper*100),
  ci_column = c(4, 6),
  ref_line = c(1,0),
  sizes = .6,
  x_trans = c("log", "none"),
  xlim = list(c(0.35, 4), c(-5, 5)),
  ticks_at = list(c(0.5, 1, 1.5, 3), c(-2.5, 0, 2.5)),
  xlab = NULL,
  arrow_lab = c("Beneficial", "Harmful"),
  title = "A. Methadone (MET) trials",
  footnote = NULL,
  nudge_y = 0.02,
  theme =
    forest_theme(
      base_size = 14,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8,
                         fontface="plain", col = "black"),
      arrow_length = 0.05,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))


txt <- "<sup>*</sup> <=40mg MET for Trial 2 and 3 <br><sup>**</sup> NA for few events<br>CI. compatibility intervals"


p_od_dis_met1 <- add_grob(p_od_dis_met, 
                          row = 20, 
                          col = 1:2,
                          order = "background",
                          gb_fn = gridtext::richtext_grob,
                          text = txt,
                          gp = gpar(fontsize = 11),
                          hjust = 0, vjust = 1, halign = 0, valign = 1,
                          x = unit(0, "npc"), y = unit(1, "npc"))


name = "Fig2_od_risk_A_fixed_Jan27_MK"
figure=paste0(folder, name, ".jpg")
ppi = 300
jpeg(figure, width=4450, height=2200, res=300)
p_od_dis_met1
dev.off()




### Figure 2: Panel B ######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
res_df <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs_iptw_ow.rds")

res_df2 <- res_df %>%
  mutate(
    Comparison = recode(Comparison, 
                        "Dose reduction\nvs\nNo dose change" = "Dose reduction vs no change",
                        "Dose increase\nvs\nNo dose change" = "Dose increase vs no change",
                        "1-<50% reduction\nvs\nNo dose change" = "Dose reduction of <50% vs no change",
                        ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reduction of >=50%/initial dose* vs no change",
                        ">=50% reduction/<=4mg\nvs\nNo dose change" = "Dose reduction of >=50%/initial dose* vs no change")
  ) %>% arrange(res_type, Model, Outcome, Trial, Comparison, Design) %>% 
filter(Design != "Prevalent new users") %>%
  mutate(Design = recode(Design,
                         "All prevalent users" = "Prevalent new users")) 


res_df2 %>% count(res_type)

#### HRs- bnx forestplot ##########################

#### Trial 4: od: iptw #######
header4 <- data.frame(Comparison = "Trial 4: After BNX interruption of 1-5 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr4_full_iptw_od_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                           Model == "IPTW",
                                           Outcome == "Overdose-related acute care visits",
                                           Trial == "Trial 4: BNX interruption of 1-5 days",
                                           Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr4_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 4: BNX interruption of 1-5 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr4_full_iptw_od <- bind_rows(header4, tr4_full_iptw_od_red, tr4_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: od: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 5: BNX interruption of 6-14 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr5_full_iptw_od <- bind_rows(header5, tr5_full_iptw_od_redlt50, 
                              tr5_full_iptw_od_redgt50, 
                              tr5_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")
  ) %>%
  ungroup()



# bnx_full_iptw_od <- bind_rows(tr4_full_iptw_od,
#                               tr5_full_iptw_od) %>%
#   mutate(HR_display = ifelse(
#     is.na(HR),
#     " ",
#     sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
#   )

format_num <- function(x) {
  ifelse(
    abs(x) >= 10,
    sprintf("%.1f", x),
    sprintf("%.2f", x)
  )
}


bnx_full_iptw_od <- bind_rows(tr4_full_iptw_od,
                              tr5_full_iptw_od) %>%
  mutate(
    HR_display = dplyr::case_when(
      is.na(HR) ~ " ",
      TRUE ~ sprintf(
        "%s (%s, %s)",
        format_num(HR),
        format_num(lower),
        format_num(upper)
      )
    )
  )



bnx_full_iptw_od$`overdose-related acute care visits` <- "" 
bnx_full_iptw_od$`Adjusted hazard\nratio (95% CI)` <- bnx_full_iptw_od$HR_display
bnx_full_iptw_od$`Dose adjustment after interruption` <- bnx_full_iptw_od$Comparison
bnx_full_iptw_od$`User design` <- bnx_full_iptw_od$Design


#### Risk- bnx forestplot ###############
extract_365 <- function(res_obj, type = c("od", "dis")) {
  type <- match.arg(type)
  
  pops <- names(res_obj$points_ci[[type]])
  
  out <- lapply(pops, function(p) {
    pop_obj <- res_obj$points_ci[[type]][[p]]
    
    # skip if NULL, NA, or not a list
    if (is.null(pop_obj) || all(is.na(pop_obj)) || !is.list(pop_obj))
      return(NULL)
    
    # skip if contrasts is missing or NULL
    if (is.null(pop_obj$contrasts) || all(is.na(pop_obj$contrasts)))
      return(NULL)
    
    df <- pop_obj$contrasts
    
    # skip if empty or horizon column missing
    if (nrow(df) == 0 || !("horizon" %in% names(df)))
      return(NULL)
    
    df %>%
      dplyr::filter(horizon == 365) %>%
      dplyr::mutate(pop = p, type = type)
  })
  
  # bind only non-null results
  dplyr::bind_rows(out)
}

res_4 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk4_iptw.rds")
res_5 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk5_iptw.rds")

#### Trial 1: 1-2 days ##############
header4 <- data.frame(Comparison = "Trial 4: After BNX interruption of 1-5 days",
                      Design = "",
                      RD = NA, RD_lower = NA, RD_upper = NA)

tr4_risk_iptw_od_red <- extract_365(res_4, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "reversion" = "Dose reduction vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose reduction vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr4_risk_iptw_od_incr <- extract_365(res_4, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "reversion" = "Dose reduction vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose increase vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr4_risk_iptw_od <- bind_rows(header4, tr4_risk_iptw_od_red, tr4_risk_iptw_od_incr) %>%
  filter(Design != "Prevalent new users") %>%
  mutate(Design = recode(Design,
                         "All prevalent users" = "Prevalent new users")) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()




#### Trial 5: 6-14 days #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      RD = NA, RD_lower = NA, RD_upper = NA)
tr5_risk_od_redlt50 <- extract_365(res_5, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction of <50% vs no change",
                             ">=50%/<=4mg" = "Dose reduction of >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose reduction of <50% vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr5_risk_od_redgt50 <- extract_365(res_5, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction of <50% vs no change",
                             ">=50%/<=4mg" = "Dose reduction of >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>%  
  filter(Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_risk_od_incr <- extract_365(res_5, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction of <50% vs no change",
                             ">=50%/<=4mg" = "Dose reduction of >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>%  
  filter(Comparison == "Dose increase vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr5_risk_iptw_od <- bind_rows(header5, 
                              tr5_risk_od_redlt50, 
                              tr5_risk_od_redgt50, 
                              tr5_risk_od_incr) %>%
  filter(Design != "Prevalent new users") %>%
  mutate(Design = recode(Design,
                         "All prevalent users" = "Prevalent new users")) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")
  ) %>%
  ungroup()
# 
# names(bnx_risk_iptw_od)
# bnx_risk_iptw_od <- bind_rows(tr4_risk_iptw_od,
#                               tr5_risk_iptw_od) %>%
#   mutate(RD_display = ifelse(
#     is.na(RD),
#     " ",
#     sprintf("%.2f(%.2f,%.2f)", RD*100, RD_lower*100, RD_upper*100))
#   )


format_num <- function(x) {
  ifelse(
    abs(x) >= 10,
    sprintf("%.1f", x),
    sprintf("%.2f", x)
  )
}

bnx_risk_iptw_od <- bind_rows(tr4_risk_iptw_od,
                              tr5_risk_iptw_od) %>%
  mutate(
    RD_display = dplyr::case_when(
      is.na(RD) ~ " ",
      TRUE ~ sprintf(
        "%s (%s, %s)",
        format_num(RD * 100),
        format_num(RD_lower * 100),
        format_num(RD_upper * 100)
      )
    )
  )



bnx_risk_iptw_od$`overdose-related acute care visits` <- "" 
bnx_risk_iptw_od$`Cumulative incidence difference (95% CI)%` <- bnx_risk_iptw_od$RD_display
bnx_risk_iptw_od$`Dose adjustment after interruption` <- bnx_risk_iptw_od$Comparison
bnx_risk_iptw_od$`User design` <- bnx_risk_iptw_od$Design

#### creating the graph ######

bnx_risk_iptw_od_v1 <- bnx_risk_iptw_od %>% select(`Cumulative incidence difference (95% CI)%`, 
                                                   `overdose-related acute care visits`) %>% 
  rename(`% difference in cumulative\nincidence at 12 months (95% CI)`=`Cumulative incidence difference (95% CI)%`)
bnx_full_risk_iptw <- bind_cols(bnx_full_iptw_od, bnx_risk_iptw_od_v1)

names(bnx_full_risk_iptw)[c(7, 12)] <- c("                             ", 
                                         "                         ")

bnx_full_risk_iptw <- bnx_full_risk_iptw %>% 
  mutate(`Dose adjustment after interruption`= 
           recode(`Dose adjustment after interruption`,
                  "   Dose reduction of >=50%/initial dose* vs no change"=
                    "   Dose reduction of >=50%/initial dose* vs no change "))



p_od_bnx <- forest(
  bnx_full_risk_iptw[, c(9,10,8,7,11,12)],
  est   = list(bnx_full_risk_iptw$HR, bnx_risk_iptw_od$RD*100),
  lower = list(bnx_full_iptw_od$lower, bnx_risk_iptw_od$RD_lower*100),
  upper = list(bnx_full_iptw_od$upper, bnx_risk_iptw_od$RD_upper*100),
  ci_column = c(4, 6),
  ref_line = c(1,0),
  sizes = .6,
  x_trans = c("log", "none"),
  xlim = list(c(0.35, 4), c(-5, 5)),
  ticks_at = list(c(0.5, 1, 1.5, 3), c(-2.5, 0, 2.5)),
  xlab = NULL,
  arrow_lab = c("Beneficial", "Harmful"),
  title = "B. Buprenorphine/naloxone (BNX) trials",
  footnote = NULL,
  nudge_y = 0.02,
  theme =
    forest_theme(
      base_size = 14,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8,
                         fontface="plain", col = "black"),
      arrow_length = 0.05,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))


txt <- "<sup>*</sup> <=40mg BNX for Trial 5<br>CI. compatibility intervals"

p_od_bnx1 <- add_grob(p_od_bnx, 
                      row = 13, 
                      col = 1:2,
                      order = "background",
                      gb_fn = gridtext::richtext_grob,
                      text = txt,
                      gp = gpar(fontsize = 11),
                      hjust = 0, vjust = 1, halign = 0, valign = 1,
                      x = unit(0, "npc"), y = unit(1, "npc"))


name = "Fig2_od_risk_B_fixed_Jan27_MK"
figure=paste0(folder, name, ".jpg")
ppi = 300
jpeg(figure, width=4450, height=1600, res=300)
p_od_bnx1
dev.off()






