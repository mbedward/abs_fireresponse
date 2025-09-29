# Generate PDF versions of the package vignettes

library(here)
library(fs)

vigns <- tools::pkgVignettes(dir = ".")$docs

for (v in vigns) {
  pdf_fname <- fs::path_file(v) |>
    fs::path_ext_set("pdf")

  pdf_path <- here(pdf_fname)

  rmarkdown::render(v, output_format = "pdf_document", output_file = pdf_path)
}
