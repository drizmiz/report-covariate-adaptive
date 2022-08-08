
library(tidyverse)

source("./CTN-003-analysis/analysis.R")

proc = function(x) {
J = x %>% as_tibble()
J$rest = ")"
J = J %>% 
  unite(CI.upper, rest, col = "CI.upper", sep = "") %>% 
  unite(CI.lower, CI.upper, col = "CI", sep = ",") %>% 
  unite(est, CI, col = "estCI", sep = " (")
}

res = proc(I)

source("./CTN-0030-analysis/analysis.R")
res = res %>% bind_cols(proc(I))
source("./CTN-0044-analysis/analysis.R")
res = res %>% bind_cols(proc(I))

res = res %>% 
  rename(`NIDA-CTN-0003` = `estCI...1`, 
         `NIDA-CTN-0030` = `estCI...3`, 
         `NIDA-CTN-0044` = `estCI...5`
         ) %>% 
  select(`NIDA-CTN-0003`, `NIDA-CTN-0030`, `NIDA-CTN-0044`)
