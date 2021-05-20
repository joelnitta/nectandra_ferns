# Rename 2013 photos exported from Photos app
#
# Photos from 2013 are currently named by default number when taken by camera ("IMGP4758.jpeg", etc)
# These need to be renamed to include species and specimen collection number before uploading to
# Ferns of the World (FOW)
#
# Photos from 2008 and 2011 were already renamed to include species and specimen collection number
# manually, so they don't need to be renamed, just exported from Photos.
# 
# Use the following settings for exporting photos from Photos:
# export settings:
#   - quality: medium
#   - size: large

# Setup ----
# Set working directory
setwd(here::here())

# Load packages
source("code/packages.R")

# Load data ----
# - Nectandra specimen data used for data analysis
nectandra_specimens <- read_csv("data/nectandra_specimens.csv")
# - Nectandra specimen data uploaded to Ferns of the World (FOW)
fow_coll_data <- read_csv("results/other/fow_nectandra_collection_data.csv") %>%
  mutate(fow_collection_number = as.character(fow_collection_number))
# - Photo metadata (collector name and number, file name, specimen ID)
photos <- read_csv("data_raw/photos.csv") %>%
  mutate(file = paste0(file, ".jpeg")) %>%
  select(specimen_id = specimenID, file) %>%
  # Exclude duplicate entries (or will create duplicate renamed files)
  unique()

# Make new file names ----

# Make a tibble linking collection name and number to specimen ID
ids <-
select(nectandra_specimens, specimen_id, specimen) %>%
  separate(specimen, c("primarycollectors", "fow_collection_number")) %>% 
  mutate(primarycollectors = str_replace_all(primarycollectors, "Nitta", "J.H. Nitta"))

# Make a tibble with the original file name and new file name (with specimen info)
photo_data_2013 <-
fow_coll_data %>%
  left_join(ids, by = c("primarycollectors", "fow_collection_number")) %>%
  left_join(photos, by = "specimen_id") %>%
  filter(str_detect(fow_date, "2013")) %>% 
  # filter(is.na(file)) %>% # 3 specimens with no field photo
  mutate(
    new_file = paste(title, fow_collection_number) %>% str_replace_all(" ", "_") %>% str_remove_all("\\."),
    new_file = make.unique(new_file, sep = "-") %>% paste0(".jpg")) %>%
  select(title, collection_number = fow_collection_number, file, new_file)

# Make tibble for moving files to folder with new names
files_to_move <- 
  tibble(file = list.files("fow_images/2013_photos/")) %>%
  inner_join(photo_data_2013, by = "file") %>%
  mutate(
    file = fs::path_abs(glue::glue("fow_images/2013_photos/{file}")),
    new_file = fs::path_abs(glue::glue("fow_images/2013_photos_renamed/{new_file}"))
    ) %>%
  select(path = file, new_path = new_file)

# Copy (rename) the files ----
pwalk(files_to_move, fs::file_copy)
