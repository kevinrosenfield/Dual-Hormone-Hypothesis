library(pacman); pacman::p_load(lmerTest, psych, sjPlot, sjmisc, sjlabelled,
                                data.table, reshape2, ggcorrplot, stargazer, gridExtra, htmlTable, tidyverse)

###############################################################################################################################
################################################                               ################################################
################################################        DATA PREPERATION       ################################################
################################################                               ################################################
###############################################################################################################################

lapply(c("tidyverse", "dplyr", "lmerTest", "psych", "sjPlot", "sjmisc", "sjlabelled",
         "data.table", "reshape2", "ggcorrplot", "stargazer", "gridExtra", "htmlTable"), FUN = citation)

# read in data

df <- as.data.frame(read.csv("df.csv"))

# change variable types

df <- df %>%
  filter(study == "ihh") %>%
  select(-c(study, ihh_race, msu_ethnicity, ihh_ethnicity, AM_T, PM_T, Kinsey_Attr:cgnq24b4_score)) %>%
  mutate_at(vars(ID:IUD), as.factor) %>%
  mutate_at(vars(age:sdi_mean), as.character) %>%
  mutate_at(vars(age:sdi_mean), as.numeric)

# calculate mean soi and soi psychology scores

df <- df %>%
  mutate(soi_mean = rowMeans(select(., c(soit1:soit9)), na.rm = TRUE),
         soi_psych = rowMeans(select(., c(soi_des, soi_attit)), na.rm = TRUE))

# calculate subject means and subject mean-center hormone variables; center subject-means; scale subject means and subject-centered data within each sex

df <- df %>%
  group_by(subID) %>%
  mutate(ihh_C.cm = mean(ihh_C, na.rm =TRUE),
         ihh_T.cm = mean(ihh_T, na.rm =TRUE),
         ihh_P.cm = mean(ihh_P, na.rm =TRUE),
         ihh_E.cm = mean(ihh_E, na.rm =TRUE),
         ihh_C.cwc = ihh_C - mean(ihh_C, na.rm =TRUE),
         ihh_T.cwc = ihh_T - mean(ihh_T, na.rm =TRUE),
         ihh_P.cwc = ihh_P - mean(ihh_P, na.rm =TRUE),
         ihh_E.cwc = ihh_E - mean(ihh_E, na.rm =TRUE),
         binary_T = ifelse(ihh_T == min(ihh_T), "low", "high"),
         binary_C = ifelse(ihh_C == min(ihh_C), "low", "high")) %>%
  ungroup() %>%
  group_by(sex) %>%
  mutate(ihh_C.cwc = ihh_C.cwc / sd(ihh_C.cwc, na.rm = TRUE),
         ihh_T.cwc = ihh_T.cwc / sd(ihh_T.cwc, na.rm = TRUE),
         ihh_P.cwc = ihh_P.cwc / sd(ihh_P.cwc, na.rm = TRUE),
         ihh_E.cwc = ihh_E.cwc / sd(ihh_E.cwc, na.rm = TRUE),
         ihh_C.cm = (ihh_C.cm - mean(ihh_C.cm, na.rm = TRUE)) / sd(ihh_C.cm, na.rm = TRUE),
         ihh_T.cm = (ihh_T.cm - mean(ihh_T.cm, na.rm = TRUE)) / sd(ihh_T.cm, na.rm = TRUE),
         ihh_P.cm = (ihh_P.cm - mean(ihh_P.cm, na.rm = TRUE)) / sd(ihh_P.cm, na.rm = TRUE),
         ihh_E.cm = (ihh_E.cm - mean(ihh_E.cm, na.rm = TRUE)) / sd(ihh_E.cm, na.rm = TRUE)) %>%
  ungroup()

# log-transformed

df.log <- df %>%
  group_by(subID) %>%
  mutate(ihh_C.log = log(ihh_C),
         ihh_T.log = log(ihh_T),
         ihh_P.log = log(ihh_P),
         ihh_E.log = log(ihh_E)) %>%
  mutate(ihh_C.log.cm = mean(ihh_C.log, na.rm =TRUE),
         ihh_T.log.cm = mean(ihh_T.log, na.rm =TRUE),
         ihh_P.log.cm = mean(ihh_P.log, na.rm =TRUE),
         ihh_E.log.cm = mean(ihh_E.log, na.rm =TRUE),
         ihh_C.log.cwc = ihh_C.log - mean(ihh_C.log, na.rm =TRUE),
         ihh_T.log.cwc = ihh_T.log - mean(ihh_T.log, na.rm =TRUE),
         ihh_P.log.cwc = ihh_P.log - mean(ihh_P.log, na.rm =TRUE),
         ihh_E.log.cwc = ihh_E.log - mean(ihh_E.log, na.rm =TRUE)) %>%
  ungroup() %>%
  group_by(sex) %>%
  mutate(ihh_C.log.cwc = ihh_C.log.cwc / sd(ihh_C.log.cwc, na.rm = TRUE),
         ihh_T.log.cwc = ihh_T.log.cwc / sd(ihh_T.log.cwc, na.rm = TRUE),
         ihh_P.log.cwc = ihh_P.log.cwc / sd(ihh_P.log.cwc, na.rm = TRUE),
         ihh_E.log.cwc = ihh_E.log.cwc / sd(ihh_E.log.cwc, na.rm = TRUE),
         ihh_C.log.cm = (ihh_C.log.cm - mean(ihh_C.log.cm, na.rm = TRUE)) / sd(ihh_C.log.cm, na.rm = TRUE),
         ihh_T.log.cm = (ihh_T.log.cm - mean(ihh_T.log.cm, na.rm = TRUE)) / sd(ihh_T.log.cm, na.rm = TRUE),
         ihh_P.log.cm = (ihh_P.log.cm - mean(ihh_P.log.cm, na.rm = TRUE)) / sd(ihh_P.log.cm, na.rm = TRUE),
         ihh_E.log.cm = (ihh_E.log.cm - mean(ihh_E.log.cm, na.rm = TRUE)) / sd(ihh_E.log.cm, na.rm = TRUE))


# break out datasets by sex/study/hormonal contraception use

df_ihh.m <- df %>% as.data.frame(filter(sex == "m"))
df_ihh.f <- df %>% filter(sex == "f")
df_ihh.hc <- df_ihh.f[df_ihh.f$HC == 1,]
df_ihh.nc <- df_ihh.f[df_ihh.f$HC == 0,]
df.log_ihh.m <- df.log %>% filter(sex == "m")
df.log_ihh.f <- df.log %>% filter(sex == "f")

# create function for easier viewing of model output

table.model <- function(model) {
  model$coefficients %>%
    as.data.frame() %>%
    rownames_to_column(var = "Effect") %>%
    mutate_if(is.numeric, round, 3) %>%
    rename(p = `Pr(>|t|)`)
}

###############################################################################################################################
################################################                               ################################################
################################################     DESCRIPTIVE STATISTICS    ################################################
################################################                               ################################################
###############################################################################################################################

# inter-item reliability tests (cronbach's alpha for SOI and SDI questions)

df_sdi <- df %>%
  dplyr::select(soit1:soit9)

df_sdi <- df %>%
  dplyr::select(sdi1:sdi14)

library(ltm)
cronbach.alpha(df_soi, na.rm = T)
cronbach.alpha(df_sdi, na.rm = T)
detach(package:ltm)
remove(df_sdi)
remove(df_soi)

# Table 1 - age, total n, and n of subs with at least one set of T, C, and either SDI or SOI samples

options(table_counter = TRUE)

df %>%
  filter(session == 1) %>%
  group_by(sex, HC) %>%
  summarise("Age (mean)" = round(mean(age), 2), "SD" = round(sd(age), 2), n = n()) %>%
  filter(!is.na(HC) | sex == "m") %>%
  ungroup() %>%
  dplyr::select(-c(HC, sex)) %>%
  htmlTable(align = "c",
            caption = "Age statistics broken out by subgroup",
            rnames = c("Women (NC)", "Women (HC)", "Men"))

# Table 2 - n of subs with one or twos set of T, C, and either SDI or SOI samples

df %>%
  filter(session == 1) %>%
  dplyr::select(n_sessions, sex, HC, ihh_T, ihh_C, sdi_mean, soi_mean) %>%
  group_by(n_sessions, sex, HC) %>%
  summarise(n = n(),
            "n (SOI)" = sum(!is.na(ihh_T) & !is.na(ihh_C) & !is.na(soi_mean)),
            "n (SDI)" = sum(!is.na(ihh_T) & !is.na(ihh_C) & !is.na(sdi_mean))) %>%
  filter(!is.na(HC) | sex == "m") %>%
  ungroup() %>%
  dplyr::select(-c(HC, sex, n_sessions)) %>%
  htmlTable(align = "c",
            col.rgroup = c("none", "#F7F7F7"),
            tspanner = c("Single session","2 Sessions"), n.tspanner = c(3,3),
            caption = "Sample sizes broken out by number of sessions",
            rnames = c("Women (NC)", "Women (HC)", "Men",
                       "Women (NC)", "Women (HC)", "Men"))


# Table 3 - SOI-R and SDI-2 Statistics

table3Rows = c("SOI-R ", "SOI: Attitudes", "SOI: Behavior", "SOI: Desire", "SDI-2",
               "SDI: Solitary", "SDI: Dyadic")
table3Cols = c("n", "Mean", "SD")

table3 <- df_ihh.m %>%
  group_by(subID) %>%
  dplyr::select(soi_mean, soi_attit, soi_behav, soi_des, sdi_mean, sdi_sol, sdi_dyad, ihh_T, ihh_C) %>%
  filter(!is.na(ihh_T), !is.na(ihh_C)) %>%
  summarise_all(funs(mean)) %>%
  describe %>%
  bind_cols(as.numeric(rep("",10)), df_ihh.nc %>%
              group_by(subID) %>%
              dplyr::select(soi_mean, soi_attit, soi_behav, soi_des, sdi_mean, sdi_sol, sdi_dyad, ihh_T, ihh_C) %>%
              filter(!is.na(ihh_T), !is.na(ihh_C)) %>%
              summarise_all(funs(mean)) %>%
              describe %>%
              bind_cols(as.numeric(rep("",10)), df_ihh.hc %>%
                          group_by(subID) %>%
                          dplyr::select(soi_mean, soi_attit, soi_behav, soi_des, sdi_mean, sdi_sol, sdi_dyad, ihh_T, ihh_C) %>%
                          filter(!is.na(ihh_T), !is.na(ihh_C)) %>%
                          summarise_all(funs(mean)) %>%
                          describe))

htmlTable(round(table3[2:8,c(2:4, 14, 16:18, 28, 30:32)], 2),
          align = "c",
          caption = "Responses to SOI-R, SDI-2, and subscales",
          rnames = table3Rows, 
          header = c(table3Cols, "", table3Cols, "", table3Cols),
          col.rgroup = c("none", "#F7F7F7"),
          cgroup = c("Men", "", "Women (NC)", "", "Women (HC)"), n.cgroup = c(3, 1, 3, 1, 3))

# Table 4 - Raw hormone average levels

table4Rows = c("T", "C", "P", "E")
table4Cols = c("n", "Mean", "SD")

table4 <- df_ihh.m %>%
  group_by(subID) %>%
  dplyr::select(ihh_T, ihh_C, ihh_P, ihh_E) %>%
  summarise_all(funs(mean)) %>%
  describe %>%
  bind_cols(as.numeric(rep("",5)), df_ihh.hc %>%
              group_by(subID) %>%
              dplyr::select(ihh_T, ihh_C, ihh_P, ihh_E) %>%
              summarise_all(funs(mean)) %>%
              describe, as.numeric(rep("",5)), df_ihh.nc %>%
              group_by(subID) %>%
              dplyr::select(ihh_T, ihh_C, ihh_P, ihh_E) %>%
              summarise_all(funs(mean)) %>%
              describe)
table4 <- round(table4, 2)
table4[4:5, 2:4] <- "NA"

htmlTable(table4[2:5,c(14, 2:4, 14, 16:18, 28, 30:32)],
          align = "c",
          caption = "Raw hormone descriptive statistics",
          rnames = table4Rows, 
          header = c("", table4Cols, "", table4Cols, "", table4Cols),
          col.rgroup = c("none", "#F7F7F7"),
          cgroup = c("", "Men", "", "Women (HC)", "", "Women (NC)"), n.cgroup = c(1, 3, 1, 3, 1, 3))

# sex differences on SOI-R, SDI-2, and subject-mean hormones

ihh.sex <- df %>%
  filter(!is.na(soi_mean), !is.nan(soi_mean), !is.na(sex)) %>%
  dplyr::select(subID, soi_mean, sex, sdi_mean, ihh_C, ihh_T) %>%
  group_by(subID) %>%
  mutate(soi_mean = mean(soi_mean),
         sdi_mean = mean(sdi_mean),
         C_mean = mean(soi_mean),
         T_mean = mean(sdi_mean))

wilcox.test(ihh.sex$soi_mean ~ ihh.sex$sex)
t.test(ihh.sex$sdi_mean ~ ihh.sex$sex)
t.test(ihh.sex$T_mean ~ ihh.sex$sex)
t.test(ihh.sex$C_mean ~ ihh.sex$sex)
remove(ihh_sex)

###############################################################################################################################
##############################################                                     ############################################
##############################################          Data visualization         ############################################
##############################################                                     ############################################
###############################################################################################################################


# Raw hormone data plots; data for both sexes

par(mfrow=c(3,2))

plot(df$ihh_C[!is.na(df$ihh_C)], ylab = "Cortisol")
plot(df$ihh_T[!is.na(df$ihh_T)], ylab = "Testosterone")
plot(df$ihh_P[!is.na(df$ihh_P)], ylab = "Progesterone")
plot(df$ihh_E[!is.na(df$ihh_E)], ylab = "Estradiol")

plot(df$ihh_C[!is.na(df$ihh_C) & df$ihh_C < 40], ylab = "Cortisol - 1 Outlier") # cortisol w/out 20 sd outlier
plot(df$ihh_T[!is.na(df$ihh_T) & df$ihh_T < 145], ylab = "Testosterone - 3 Outliers") # testosterone w/out 20 sd outlier

# Hormone data: Each point is a subject mean, global-mean centered on zero, scaled (to sd = 1) within sex; data for both sexes

plot(df$ihh_C.cm[!is.na(df$ihh_C.cm)], ylab = "Cortisol")
plot(df$ihh_T.cm[!is.na(df$ihh_T.cm)], ylab = "Testosterone")
plot(df$ihh_P.cm[!is.na(df$ihh_P.cm)], ylab = "Progesterone")
plot(df$ihh_E.cm[!is.na(df$ihh_E.cm)], ylab = "Estradiol")

plot(df$ihh_C.cm[!is.na(df$ihh_C.cm) & df$ihh_C.cm < 19.5], ylab = "Cortisol") # cortisol without ~ 20 sd outlier
plot(df$ihh_T.cm[!is.na(df$ihh_T.cm) & df$ihh_T.cm < 9.5], ylab = "Testosterone") # testosterone without ~ 20 sd outlier

# Hormone data: Each point is a session, subject-mean centered on zero, scaled (to sd = 1) within sex; data for both sexes

par(mfrow=c(2,2))

plot(df$ihh_C.cwc[!is.na(df$ihh_C.cwc)])
plot(df$ihh_T.cwc[!is.na(df$ihh_T.cwc)])
plot(df$ihh_P.cwc[!is.na(df$ihh_P.cwc)])
plot(df$ihh_E.cwc[!is.na(df$ihh_E.cwc)])

# Supplementary Figures 1-4 : Correlations of scales and their subscales; each scale/sex combination presented separately

pairs.panels(
  df_ihh.m %>% select(sdi_mean, sdi_sol, sdi_dyad), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, cex.cor = 2, ellipses = FALSE)

pairs.panels(
  df_ihh.f %>% select(sdi_mean, sdi_sol, sdi_dyad), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, cex.cor = 2, ellipses = FALSE)

pairs.panels(
  df_ihh.m %>% select(soi_mean, soi_attit, soi_behav, soi_des), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, cex.cor = 2, ellipses = FALSE)

pairs.panels(
  df_ihh.f %>% select(soi_mean, soi_attit, soi_behav, soi_des), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, cex.cor = 2, ellipses = FALSE)

# Supplementary Figures 5-6: Correlations of ALL scales and subscales, seperately by sex

pairs.panels(
  df_ihh.m %>% select(soi_mean, soi_attit, soi_behav, soi_des, sdi_mean, sdi_sol, sdi_dyad), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, cex.cor = 2, ellipses = FALSE)

pairs.panels(
  df_ihh.f %>% select(soi_mean, soi_attit, soi_behav, soi_des, sdi_mean, sdi_sol, sdi_dyad), stars= TRUE, method = "pearson", hist.col = "gray", cex.labels=1.9, cex.cor = 2, ellipses = FALSE)

# Supplementary Figure 7: Log-transformed raw hormone data plots; data for both sexes

par(mfrow=c(2,2))

plot(df.log$ihh_C.log[!is.na(df.log$ihh_C.log)], ylab = "Log Cortisol")
plot(df.log$ihh_T.log[!is.na(df.log$ihh_T.log)], ylab = "Log Testosterone")
plot(df.log$ihh_P.log[!is.na(df.log$ihh_P.log)], ylab = "Log Progesterone")
plot(df.log$ihh_E.log[!is.na(df.log$ihh_E.log)], ylab = "Log Estradiol")

# Supplementary Figure 8: Log-transformed hormone data: Each point is a subject mean, global-mean centered on zero, scaled (to sd = 1) within sex; data for both sexes

plot(df.log$ihh_C.log.cm[!is.na(df.log$ihh_C.log.cm)], ylab = "Log Cortisol")
plot(df.log$ihh_T.log.cm[!is.na(df.log$ihh_T.log.cm)], ylab = "Log Testosterone")
plot(df.log$ihh_P.log.cm[!is.na(df.log$ihh_P.log.cm)], ylab = "Log Progesterone")
plot(df.log$ihh_E.log.cm[!is.na(df.log$ihh_E.log.cm)], ylab = "Log Estradiol")

# Supplementary Figure 9: Log-transformed hormone data: Each point is a session, subject-mean centered on zero, scaled (to sd = 1) w/in sex; data for both sexes

plot(df.log$ihh_C.log.cwc[!is.na(df.log$ihh_C.log.cwc)], ylab = "Log Cortisol")
plot(df.log$ihh_T.log.cwc[!is.na(df.log$ihh_T.log.cwc)], ylab = "Log Testosterone")
plot(df.log$ihh_P.log.cwc[!is.na(df.log$ihh_P.log.cwc)], ylab = "Log Progesterone")
plot(df.log$ihh_E.log.cwc[!is.na(df.log$ihh_E.log.cwc)], ylab = "Log Estradiol")

# wrangle data for surface plots (Figures 1-2; produced in Excel: Insert->Charts->Surface while highlighting all data in .csv file(s) produced)

surface.csv <- function(df, yVarLim, zVarLim) {
  
  df.surface <- df %>%
    select(subID, session, ihh_C, ihh_T, sdi_sol, HC, ihh_T.cwc, ihh_C.cwc, ihh_T.cm, ihh_C.cm)
  
  surface_model <- lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1|subID), data  = df.surface)
  grid <- data.frame("ihh_T.cwc" = seq(0, yVarLim, yVarLim / 100), "ihh_C.cwc" = round(seq(0, zVarLim, zVarLim / 100), 2))
  #grid <- data.frame("ihh_T.cwc" = seq(0 - yVarLim, yVarLim, yVarLim / 50), "ihh_C.cwc" = round(seq(0 - zVarLim, zVarLim, zVarLim / 50), 2))
  grid <- data.frame(expand.grid(grid)) %>%
    bind_cols('subID' = NA, 'ihh_T.cm' = rep(0,101*101), 'ihh_C.cm' = rep(0,101*101)) %>%
    mutate(predicted = predict(surface_model, ., re.form = NA)) %>%
    select(-c(ihh_T.cm, ihh_C.cm, subID)) %>%
    pivot_wider(values_from = predicted, names_from = ihh_T.cwc)

  return(grid)
}

grid.m <- surface.csv(df_ihh.m, 5, 5)
write_csv(grid.m, path = '/Users/kevinrosenfield/Box/PSU/Projects/CGN/Repo/CGN/surface_men.csv')

grid.f.nc <- surface.csv(df_ihh.nc, 5, 5)
write_csv(grid.f.nc, path = '/Users/kevinrosenfield/Box/PSU/Projects/CGN/Repo/CGN/surface_women.nc.csv')

grid.f.hc <- surface.csv(df_ihh.hc, 5, 5)
write_csv(grid.f.hc, path = '/Users/kevinrosenfield/Box/PSU/Projects/CGN/Repo/CGN/surface_women.hc.csv')

###############################################################################################################################
##############################################                                     ############################################
##############################################        Preliminary Analyses         ############################################
##############################################                                     ############################################
###############################################################################################################################

# are SOI and/or SDI higher in higher T sessions in women?

soi_binary_T_f <- df_ihh.f %>%
  filter(n_sessions == 2, !is.na(soi_mean) & !is.na(binary_T)) %>%
  select(subID, binary_T, binary_C, soi_mean) %>%
  pivot_wider(id_cols = subID, names_from = c(binary_T), values_from = c(soi_mean, binary_C))

sdi_binary_T_f <- df_ihh.f %>%
  filter(n_sessions == 2, !is.na(sdi_mean) & !is.na(binary_T)) %>%
  select(subID, binary_T, binary_C, sdi_mean) %>%
  pivot_wider(id_cols = subID, names_from = c(binary_T), values_from = c(sdi_mean, binary_C))

t.test(soi_binary_T_f$soi_mean_low, soi_binary_T_f$soi_mean_high, paired = TRUE)
t.test(sdi_binary_T_f$sdi_mean_low, sdi_binary_T_f$sdi_mean_high, paired = TRUE)

# are SOI and/or SDI higher in higher T sessions in men?

soi_binary_T_m <- df_ihh.m %>%
  filter(n_sessions == 2, !is.na(soi_mean) & !is.na(binary_T)) %>%
  select(subID, binary_T, binary_C, soi_mean) %>%
  pivot_wider(id_cols = subID, names_from = c(binary_T), values_from = c(soi_mean, binary_C))

sdi_binary_T_m <- df_ihh.m %>%
  filter(n_sessions == 2, !is.na(sdi_mean) & !is.na(binary_T)) %>%
  select(subID, binary_T, binary_C, sdi_mean) %>%
  pivot_wider(id_cols = subID, names_from = c(binary_T), values_from = c(sdi_mean, binary_C))

t.test(soi_binary_T_m$soi_mean_low, soi_binary_T_m$soi_mean_high, paired = TRUE)
t.test(sdi_binary_T_m$sdi_mean_low, sdi_binary_T_m$sdi_mean_high, paired = TRUE)

# are SOI and/or SDI higher in higher T sessions that are also higher/lower C sessions (additive/antagonistic effects) in women?

soi_add_C_f <- soi_binary_T_f %>%
  filter(binary_C_low == "low")

soi_ant_C_f <- soi_binary_T_f %>%
  filter(binary_C_low == "high")

sdi_add_C_f <- sdi_binary_T_f %>%
  filter(binary_C_low == "low")

sdi_ant_C_f <- sdi_binary_T_f %>%
  filter(binary_C_low == "high")

t.test(soi_add_C_f$soi_mean_low, soi_add_C_f$soi_mean_high, paired = TRUE)
t.test(soi_ant_C_f$soi_mean_low, soi_ant_C_f$soi_mean_high, paired = TRUE)

t.test(sdi_add_C_f$sdi_mean_low, sdi_add_C_f$sdi_mean_high, paired = TRUE)
t.test(sdi_ant_C_f$sdi_mean_low, sdi_ant_C_f$sdi_mean_high, paired = TRUE)

# are SOI and/or SDI higher in higher T sessions that are also higher/lower C sessions (additive/antagonistic effects) in men?

soi_add_C_m <- soi_binary_T_m %>%
  filter(binary_C_low == "low")

soi_ant_C_m <- soi_binary_T_m %>%
  filter(binary_C_low == "high")

sdi_add_C_m <- sdi_binary_T_m %>%
  filter(binary_C_low == "low")

sdi_ant_C_m <- sdi_binary_T_m %>%
  filter(binary_C_low == "high")

t.test(soi_add_C_m$soi_mean_low, soi_add_C_m$soi_mean_high, paired = TRUE)
t.test(soi_ant_C_m$soi_mean_low, soi_ant_C_m$soi_mean_high, paired = TRUE)

t.test(sdi_add_C_m$sdi_mean_low, sdi_add_C_m$sdi_mean_high, paired = TRUE)
t.test(sdi_ant_C_m$sdi_mean_low, sdi_ant_C_m$sdi_mean_high, paired = TRUE)

### Preliminarry analysis: HORMONES UNCENTERED: Within- and between-individual effects ARE NOT isolated

# men

m.d <- table.model(summary(lmer(sdi_mean ~ ihh_T*ihh_C + (1 | subID), data = df_ihh.m))) # sdi-2

m.sol <- table.model(summary(lmer(sdi_sol ~ ihh_T*ihh_C + (1 | subID), data = df_ihh.m))) # sdi-2[solitary]

m.dyad <- table.model(summary(lmer(sdi_dyad ~ ihh_T*ihh_C + (1 | subID), data = df_ihh.m))) # sdi-2[dyadic]

m.s <- table.model(summary(lmer(soi_mean ~ ihh_T*ihh_C + (1 | subID), data = df_ihh.m))) # soiR

m.att <- table.model(summary(lmer(soi_attit ~ ihh_T*ihh_C + (1 | subID), data = df_ihh.m))) # soi-R[attitude]

m.beh <- table.model(summary(lmer(soi_behav ~ ihh_T*ihh_C + (1 | subID), data = df_ihh.m))) # soi-R[behavior]

m.des <- table.model(summary(lmer(soi_des ~ ihh_T*ihh_C + (1 | subID), data = df_ihh.m))) # soi-R[desire]


# women

w.d <- table.model(summary(lmer(sdi_mean ~ ihh_T*ihh_C*HC + (1 | subID), data = df_ihh.f))) # sdi-2

w.sol <- table.model(summary(lmer(sdi_sol ~ ihh_T*ihh_C*HC + (1 | subID), data = df_ihh.f))) # sdi-2[solitary]

w.dyad <- table.model(summary(lmer(sdi_dyad ~ ihh_T*ihh_C*HC + (1 | subID), data = df_ihh.f))) # sdi-2[dyadic]

w.s <- table.model(summary(lmer(soi_mean ~ ihh_T*ihh_C*HC + (1 | subID), data = df_ihh.f))) # soiR

w.att <- table.model(summary(lmer(soi_attit ~ ihh_T*ihh_C*HC + (1 | subID), data = df_ihh.f))) # soi-R[attitude]

w.beh <- table.model(summary(lmer(soi_behav ~ ihh_T*ihh_C*HC + (1 | subID), data = df_ihh.f))) # soi-R[behavior]

w.des <- table.model(summary(lmer(soi_des ~ ihh_T*ihh_C*HC + (1 | subID), data = df_ihh.f))) # soi-R[desire]


###############################################################################################################################
##############################################                                     ############################################
##############################################    DUAL-HORMONE HYPOTHESIS TESTS    ############################################
##############################################                                     ############################################
###############################################################################################################################

### Main Analyses: hormones centered; Within(cwc)- and between(cm)-individual effects ARE isolated

# men

m.d_w <- table.model(summary(lmer(sdi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.m))) # sdi-2

m.sol_w <- table.model(summary(lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.m))) # sdi-2[solitary]

m.dyad_w <- table.model(summary(lmer(sdi_dyad ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.m))) # sdi-2[dyadic]

m.s_w <- table.model(summary(lmer(soi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.m))) # soiR

m.att_w <- table.model(summary(lmer(soi_attit ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.m))) # soi-R[attitude]

m.beh_w <- table.model(summary(lmer(soi_behav ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.m))) # soi-R[behavior]

m.des_w <- table.model(summary(lmer(soi_des ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.m))) # soi-R[desire]

# women

w.d_w <- table.model(summary(lmer(sdi_mean ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + (1 | subID), data = df_ihh.f))) # sdi-2

w.sol_w <- table.model(summary(lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + (1 | subID), data = df_ihh.f))) # sdi-2[solitary]

w.dyad_w <- table.model(summary(lmer(sdi_dyad ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + (1 | subID), data = df_ihh.f))) # sdi-2[dyadic]

w.s_w <- table.model(summary(lmer(soi_mean ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + (1 | subID), data = df_ihh.f))) # soiR

w.att_w <- table.model(summary(lmer(soi_attit ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + (1 | subID), data = df_ihh.f))) # soi-R[attitude]

w.beh_w <- table.model(summary(lmer(soi_behav ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + (1 | subID), data = df_ihh.f))) # soi-R[behavior]

w.des_w <- table.model(summary(lmer(soi_des ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + (1 | subID), data = df_ihh.f))) # soi-R[desire]

###############################################################################################################################
##############################################                                     ############################################
##############################################         Robustness Analyses         ############################################
##############################################                                     ############################################
###############################################################################################################################

# women taking hormonal contraception

hc.d_w <- table.model(summary(lmer(sdi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.hc))) # sdi-2

hc.sol_w <- table.model(summary(lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.hc))) # sdi-2[solitary]

hc.dyad_w <- table.model(summary(lmer(sdi_dyad ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.hc))) # sdi-2[dyadic]

hc.s_w <- table.model(summary(lmer(soi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.hc))) # soiR

hc.att_w <- table.model(summary(lmer(soi_attit ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.hc))) # soi-R[attitude]

hc.beh_w <- table.model(summary(lmer(soi_behav ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.hc))) # soi-R[behavior]

hc.des_w <- table.model(summary(lmer(soi_des ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.hc))) # soi-R[desire]

# women NOT taking hormonal contraception

nc.d_w <- table.model(summary(lmer(sdi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.nc))) # sdi-2

nc.sol_w <- table.model(summary(lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.nc))) # sdi-2[solitary]

nc.dyad_w <- table.model(summary(lmer(sdi_dyad ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.nc))) # sdi-2[dyadic]

nc.s_w <- table.model(summary(lmer(soi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.nc))) # soiR

nc.att_w <- table.model(summary(lmer(soi_attit ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.nc))) # soi-R[attitude]

nc.beh_w <- table.model(summary(lmer(soi_behav ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.nc))) # soi-R[behavior]

nc.des_w <- table.model(summary(lmer(soi_des ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + (1 | subID), data = df_ihh.nc))) # soi-R[desire]

##### ESM Table 1 - Sociosexual orientation (women broken out by hormonal contraception use)

table3cols <- c("Estimate", "df", "t", "p")
table3rows <- c("Intercept", "T (within)", "C (within)", "T (between)",
                "C (between)", "T x C (within)", "T x C (between)")
table3 <- round(bind_rows(bind_cols(m.s_w[-c(1,3)], m.att_w[-c(1,3)], m.beh_w[-c(1,3)], m.des_w[-c(1,3)]),
                          bind_cols(nc.s_w[-c(1,3)], nc.att_w[-c(1,3)], nc.beh_w[-c(1,3)], nc.des_w[-c(1,3)]),
                          bind_cols(hc.s_w[-c(1,3)], hc.att_w[-c(1,3)], hc.beh_w[-c(1,3)], hc.des_w[-c(1,3)])),2)

htmlTable(table3,
          align = "c",
          caption = "SOI-R Model results",
          rnames = c(table3rows, table3rows, table3rows),
          header = c(table3cols, table3cols, table3cols, table3cols),
          col.rgroup = c("none", "#F7F7F7"),
          cgroup = list(c("SOI-R"), c("General", "Attitudes", "Behavior", "Desire")),
          n.cgroup = list(c(4), c(4, 4, 4, 4)),
          tspanner = c("Men","Women not usng hormonal contraception", "Women using hormonal contraception"), n.tspanner = c(7,7),
          css.tspanner = "padding-top: 0.5em; font-weight:bold",
          css.cell = "padding-right: 0.5em")

##### ESM Table 2 - Sexual desire (women broken out by hormonal contraception use)

table4cols <- c("Estimate", "df", "t", "p")
table4rows <- c("Intercept", "T (within)", "C (within)", "T (between)",
                "C (between)", "T x C (within)", "T x C (between)")
table4 <- round(bind_rows(bind_cols(m.d_w[-c(1,3)], m.sol_w[-c(1,3)], m.dyad_w[-c(1,3)]),
                          bind_cols(nc.d_w[-c(1,3)], nc.sol_w[-c(1,3)], nc.dyad_w[-c(1,3)]),
                          bind_cols(hc.d_w[-c(1,3)], hc.sol_w[-c(1,3)], hc.dyad_w[-c(1,3)])),2)

htmlTable(table4,
          align = "c",
          caption = "SDI-2 Model results",
          rnames = c(table4rows, table4rows, table4rows),
          header = c(table4cols, table4cols, table4cols),
          col.rgroup = c("none", "#F7F7F7"),
          cgroup = list(c("SDI-2"), c("General", "Solitary", "Dyadic")),
          n.cgroup = list(c(3), c(4, 4, 4)),
          tspanner = c("Men","Women not usng hormonal contraception", "Women using hormonal contraception"), n.tspanner = c(7,7),
          css.tspanner = "padding-top: 0.5em; font-weight:bold",
          css.cell = "padding-right: 0.5em")

#####

# do effects of C and T remain after accounting for effects of ovarian hormones on SDI and/or SOI in non-contracepting women?

ov.d_w <- table.model(summary(lmer(sdi_mean ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                 ihh_C.cm*ihh_T.cm + (1 | subID), data = df_ihh.nc))) # sdi-2

ov.sol_w <- table.model(summary(lmer(sdi_sol ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                   ihh_C.cm*ihh_T.cm + (1 | subID), data = df_ihh.nc))) # sdi-2[solitary]

ov.dyad_w <- table.model(summary(lmer(sdi_dyad ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                    ihh_C.cm*ihh_T.cm + (1 | subID), data = df_ihh.nc))) # sdi-2[dyadic]

ov.s_w <- table.model(summary(lmer(soi_mean ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                 ihh_C.cm*ihh_T.cm + (1 | subID), data = df_ihh.nc))) # soiR

ov.att_w <- table.model(summary(lmer(soi_attit ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                   ihh_C.cm*ihh_T.cm + (1 | subID), data = df_ihh.nc))) # soi-R[attitude]

ov.beh_w <- table.model(summary(lmer(soi_behav ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                   ihh_C.cm*ihh_T.cm + (1 | subID), data = df_ihh.nc))) # soi-R[behavior]

ov.des_w <- table.model(summary(lmer(soi_des ~ ihh_P.cwc*ihh_E.cwc + ihh_P.cm*ihh_E.cm + ihh_C.cwc*ihh_T.cwc +
                   ihh_C.cm*ihh_T.cm + (1 | subID), data = df_ihh.nc))) # soi-R[desire]

# With predictor hormones log-transformed before centering and scaling

### HORMONES ARE CENTERED: Within(cwc)- and between(cm)-individual effects ARE isolated; hormones are log-transformed

# men

m.d.log <- table.model(summary(lmer(sdi_mean ~ ihh_T.log.cwc*ihh_C.log.cwc +
                                      ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df.log_ihh.m))) # sdi-2

m.sol.log <- table.model(summary(lmer(sdi_sol ~ ihh_T.log.cwc*ihh_C.log.cwc +
                                        ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df.log_ihh.m))) # sdi-2[solitary]

m.dyad.log <- table.model(summary(lmer(sdi_dyad ~ ihh_T.log.cwc*ihh_C.log.cwc +
                                         ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df.log_ihh.m))) # sdi-2[dyadic]

m.s.log <- table.model(summary(lmer(soi_mean ~ ihh_T.log.cwc*ihh_C.log.cwc +
                                      ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df.log_ihh.m))) # soiR

m.att.log <- table.model(summary(lmer(soi_attit ~ ihh_T.log.cwc*ihh_C.log.cwc +
                                        ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df.log_ihh.m))) # soi-R[attitude]

m.beh.log <- table.model(summary(lmer(soi_behav ~ ihh_T.log.cwc*ihh_C.log.cwc +
                                        ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df.log_ihh.m))) # soi-R[behavior]

m.des.log <- table.model(summary(lmer(soi_des ~ ihh_T.log.cwc*ihh_C.log.cwc +
                                        ihh_T.log.cm*ihh_C.log.cm + (1 | subID), data = df.log_ihh.m))) # soi-R[desire]

# women

w.d.log <- table.model(summary(lmer(sdi_mean ~ ihh_T.log.cwc*ihh_C.log.cwc*HC +
                                      ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df.log_ihh.f))) # sdi-2

w.sol.log <- table.model(summary(lmer(sdi_sol ~ ihh_T.log.cwc*ihh_C.log.cwc*HC +
                                        ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df.log_ihh.f))) # sdi solitary

w.dyad.log <- table.model(summary(lmer(sdi_dyad ~ ihh_T.log.cwc*ihh_C.log.cwc*HC + 
                                         ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df.log_ihh.f))) # sdi dyadic

w.s.log <- table.model(summary(lmer(soi_mean ~ ihh_T.log.cwc*ihh_C.log.cwc*HC +
                                      ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df.log_ihh.f))) # soiR

w.att.log <- table.model(summary(lmer(soi_attit ~ ihh_T.log.cwc*ihh_C.log.cwc*HC +
                                        ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df.log_ihh.f))) # attitude

w.beh.log <- table.model(summary(lmer(soi_behav ~ ihh_T.log.cwc*ihh_C.log.cwc*HC +
                                        ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df.log_ihh.f))) # behavior

w.des.log <- table.model(summary(lmer(soi_des ~ ihh_T.log.cwc*ihh_C.log.cwc*HC +
                                        ihh_T.log.cm*ihh_C.log.cm*HC + (1 | subID), data = df.log_ihh.f))) # soi desire

# men

m.d_w_age <- table.model(summary(lmer(sdi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + age + (1 | subID), data = df_ihh.m))) # sdi-2

m.sol_w_age <- table.model(summary(lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + age + (1 | subID), data = df_ihh.m))) # sdi-2[solitary]

m.dyad_w_age <- table.model(summary(lmer(sdi_dyad ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + age + (1 | subID), data = df_ihh.m))) # sdi-2[dyadic]

m.s_w_age <- table.model(summary(lmer(soi_mean ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + age + (1 | subID), data = df_ihh.m))) # soiR

m.att_w_age <- table.model(summary(lmer(soi_attit ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + age + (1 | subID), data = df_ihh.m))) # soi-R[attitude]

m.beh_w_age <- table.model(summary(lmer(soi_behav ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + age + (1 | subID), data = df_ihh.m))) # soi-R[behavior]

m.des_w_age <- table.model(summary(lmer(soi_des ~ ihh_T.cwc*ihh_C.cwc + ihh_T.cm*ihh_C.cm + age + (1 | subID), data = df_ihh.m))) # soi-R[desire]

# women

w.d_w_age <- table.model(summary(lmer(sdi_mean ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + age + (1 | subID), data = df_ihh.f))) # sdi-2

w.sol_w_age <- table.model(summary(lmer(sdi_sol ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + age + (1 | subID), data = df_ihh.f))) # sdi-2[solitary]

w.dyad_w_age <- table.model(summary(lmer(sdi_dyad ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + age + (1 | subID), data = df_ihh.f))) # sdi-2[dyadic]

w.s_w_age <- table.model(summary(lmer(soi_mean ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + age + (1 | subID), data = df_ihh.f))) # soiR

w.att_w_age <- table.model(summary(lmer(soi_attit ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + age + (1 | subID), data = df_ihh.f))) # soi-R[attitude]

w.beh_w_age <- table.model(summary(lmer(soi_behav ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + age + (1 | subID), data = df_ihh.f))) # soi-R[behavior]

w.des_w_age <- table.model(summary(lmer(soi_des ~ ihh_T.cwc*ihh_C.cwc*HC + ihh_T.cm*ihh_C.cm*HC + age + (1 | subID), data = df_ihh.f))) # soi-R[desire]

