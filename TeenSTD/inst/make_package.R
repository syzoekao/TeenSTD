# Set working directory in `COVID19pack/` folder
# setwd("..")
tmp_path <- "inst/raw_data/"
MDH_data <- readRDS(paste0(tmp_path, "MDH.yr.rds"))
preg_data <- readRDS(paste0(tmp_path, "Preg.yr.rds"))
sex_behave_data <- get(load(paste0(tmp_path, "sex.behave.RData")))
posterior_set <- readRDS(paste0("test/posterior_set.RDS"))

usethis::use_data(MDH_data, overwrite = T)
usethis::use_data(preg_data, overwrite = T)
usethis::use_data(sex_behave_data, overwrite = T)
usethis::use_data(posterior_set, overwrite = T)

devtools::document()
package_loc <- devtools::build()
install.packages(package_loc, repos = NULL)

#### Removing tar.gz file
# setwd("..")
# file.remove("TeenSTD_0.1.0.tar.gz")
