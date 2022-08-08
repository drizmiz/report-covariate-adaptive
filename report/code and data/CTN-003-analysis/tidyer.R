
library(tidyverse)

# treatment arm, sex

dm = read_csv("./dm.csv")

cc = dm %>% select(USUBJID, SEX, ARMCD)

# dose (stratification)

sc = read_csv("./sc.csv")

sc_shorter = sc %>%
  filter(SCTEST == "BUP DOSE AT RANDOMIZATION") %>%
  select(USUBJID, SCORRES)

cc = cc %>% full_join(sc_shorter, by = "USUBJID")

# ARSW, COWS and VAS

qs = read_csv("./qs.csv") %>% filter(VISITNUM == 0)

# ARSW

arsw = qs %>%
  filter(QSCAT == "ADJECTIVE RATING SCALE FOR WITHDRAWAL") %>%
  select(USUBJID, QSTESTCD, QSORRES)

translate = function(scl) {
  str_match(scl, "[0-9]")[[1]] %>% utf8ToInt() - utf8ToInt("0")
}

arsw$QSORRES = arsw$QSORRES %>% map(translate)

arsw = arsw %>%
  pivot_wider(names_from = QSTESTCD, values_from = QSORRES, values_fill=list(as.integer(NA)))

calc_arsw = function(...) {
  x = list(...)
  
  x$ARSWTOTAL = sum(x[2:17] %>% as_vector())
  
  return(x)
}

arsw = arsw %>% pmap_dfr(calc_arsw) %>% select(USUBJID, ARSWTOTAL)

cc = cc %>% full_join(arsw, by = "USUBJID")

# COWS

cows = qs %>%
  filter(QSTESTCD == "COWS012") %>%
  select(USUBJID, QSORRES) %>%
  rename(COWSTOTAL = QSORRES)

cc = cc %>% full_join(cows, by = "USUBJID")

# VAS

vas =  qs %>%
  filter(QSCAT == "VISUAL ANALOG SCALE") %>%
  select(USUBJID, QSTESTCD, QSORRES)

vas = vas %>%
  pivot_wider(names_from = QSTESTCD, values_from = QSORRES)

cc = cc %>% full_join(vas, by = "USUBJID") 

# toxicology

urine_res = read_csv("su.csv") %>% filter(SUCAT == "CLINIC VISIT DRUG USE")

tox = urine_res %>%
  filter(SUSCAT == "TOXICOLOGY RESULTS", VISITNUM == 0) %>% 
  select(USUBJID, SUTRT, SUOCCUR) %>% 
  pivot_wider(names_from = SUTRT, values_from = SUOCCUR)

cc = cc %>% full_join(tox, by = "USUBJID") 

# outcome

out = urine_res %>%
  filter(SUSCAT == "TOXICOLOGY RESULTS", EPOCH == "END OF TAPER") %>% 
  select(USUBJID, SUTRT, SUOCCUR) %>% 
  pivot_wider(names_from = SUTRT, values_from = SUOCCUR)

opi = c("MTD", "MOR", "OXY")

calc_out = function(...) {
  x = list(...)
  
  x[[x$name]] = (sum((x[opi] %>% as_vector()) == "Y") != 0) 
  
  return(x)
}

out = out %>% pmap_dfr(calc_out, name = "outcome") %>% select(USUBJID, outcome)

cc = cc %>% full_join(out, by = "USUBJID") 

cc = cc %>% type_convert()

# cc %>% write_excel_csv("dataset.csv")

calc_vas = function(...) {
  x = list(...)
  
  x$VAS = sum(x[c("VAS001", "VAS002", "VAS003")] %>% as_vector())
  
  return(x)
}

cc = cc %>% 
  pmap_dfr(calc_out, name = "baseline.OPI") %>%
  pmap_dfr(calc_vas) %>% 
  select(USUBJID, SEX, ARMCD, SCORRES, ARSWTOTAL, COWSTOTAL, VAS, baseline.OPI, outcome) %>% 
  rename(COWS = COWSTOTAL, ARSW = ARSWTOTAL, sex = SEX, 
         outcome.OPI = outcome, strata = SCORRES, arm = ARMCD)

cc = cc %>% select(-USUBJID) %>% filter(arm == "TAPER7" | arm == "TAPER28")

  cc %>% write_excel_csv("CTN03.csv")
