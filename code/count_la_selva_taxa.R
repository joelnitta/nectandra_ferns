# From the raw flora list of La Selva (Feb. 2017),
# subset to only pteridophytes and count the number of 
# unique taxa.

library(tidyverse)

la_selva_pteridos <- readxl::read_excel("data_raw/Lista_especies_LS_feb2017.xlsx", skip = 12) %>%
  janitor::clean_names() %>%
  rename(family = familia, genus = genero, specific_epithet = especie, 
         author = autor, notes = historia_taxonomica,
         habit = habito, habit_atr = habito_atributo) %>%
  left_join(taxize::apgFamilies()) %>%
  filter(order %in% c(
    "Cyatheales",
    "eupolypod II",
    "Gleicheniales" ,
    "Hymenophyllales",
    "Lycopodiales",
    "Marattiales",
    "Ophioglossales",
    "Polypodiales",
    "Polypodiales-eupolypod I",
    "Polypodiales-eupolypod II"
  )) %>%
  mutate(taxon = jntools::paste3(genus, specific_epithet))

# 197
n_distinct(la_selva_pteridos$taxon)
