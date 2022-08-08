
library(tidyverse)

bfs = read_csv("BFS.csv")

bfs = bfs %>%
  select(patdeid, BFSEX, BFAGE)

enrb = read_csv("ENRB.csv")

enrb = enrb %>%
  select(STRATA, patdeid)

# enrb = enrb %>%
#   select(RAPRIMSU, RABLURBR, patdeid) %>%
#   mutate(STRATA = interaction(RAPRIMSU, RABLURBR)) %>% 
#   select(!c(RAPRIMSU, RABLURBR))

trt = read_csv("TRT.csv") %>% 
  select(trt, patdeid)

uds = read_csv("UDS.csv")

uds = uds %>%
  filter(UD1ADULT == 0 | UD2ADULT == 0)

calc_out = function(...) {
  x = list(...)

  x$UD1RES = x[15:24] %>% as.logical() %>% any()
  x$UD2RES = x[32:41] %>% as.logical() %>% any()
  if(x$UD1ADULT == 0) {
    x$UDRES = x$UD1RES
  } else {
    x$UDRES = x$UD2RES
  }

  return(x)
}

uds = uds %>%
  pmap_dfr(calc_out, name = "outcome") %>%
  filter(VISNO != "3MFU" & VISNO != "6MFU") %>%
  select(patdeid, VISNO, UDRES)

uds = uds %>%
  pivot_wider(names_from = VISNO, values_from = UDRES)

cc = enrb %>%
  left_join(bfs) %>% 
  left_join(trt) %>% 
  left_join(uds) %>% 
  rename(strata = STRATA, age = BFAGE, `0` = `00`,
         gender = BFSEX, arm = trt)

cc %>% write_excel_csv("CTN44.csv")
