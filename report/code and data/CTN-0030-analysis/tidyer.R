
library(tidyverse)

dem = read_csv("T_FRDEM.csv")

dem = dem %>% select(patdeid, DEM001, DEM002)

dem$DEM001 = -dem$DEM001 / 365

dem = dem %>% rename(age = DEM001, sex = DEM002)

#dem$age = dem$age %>% replace_na(mean(dem$age, na.rm = TRUE))

cc = dem

rnd = read_csv("T_FRRND.csv")

rnd = rnd %>%
  filter(!is.na(RND005)) %>% 
  select(patdeid, RND005, RND007, RND008) %>% 
  rename(arm1 = RND005, heroin_history = RND007, chronic_pain = RND008)

cc = cc %>% right_join(rnd, by = "patdeid")

ur = read_csv("T_FRUDS.csv")

# no UDS014
# tests = c("UDS005", "UDS006", "UDS007", "UDS008", "UDS009", "UDS010", "UDS011", "UDS012", "UDS013")
tests = c("UDS011")

ur = ur %>% 
  filter(VISIT %in% c("BL", "P1Wk1A", "P1Wk1B", "P1Wk2", "P1Wk3", "P1Wk4"))

calc_out = function(...) {
  x = list(...)
  
  t_vec = x[tests] %>% as_vector()
  
  x$UDSPOS = any(t_vec > 0) 
  
  return(x)
}

ur = ur %>% pmap_dfr(calc_out)

ur = ur %>% 
  select(patdeid, VISIT, UDSPOS) %>% 
  pivot_wider(names_from = VISIT, values_from = UDSPOS)
  
cc = cc %>% left_join(ur, by = "patdeid")

cc %>% write_excel_csv("CTN30.csv")
