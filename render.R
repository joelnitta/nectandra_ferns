rmarkdown::render(
  here::here("ms/nectandra_pteridos.Rmd"),
  output_file = here::here("ms/nectandra_pteridos.pdf"),
  quiet = TRUE)