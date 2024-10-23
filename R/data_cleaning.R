library(here)
library(tidyverse)
library(readxl)
library(plater)
library(janitor)

source(here("R", "functions.R"))

############################
# download DRC data from dropbox
############################
db_dest <- here("data", "dropbox")
dropbox_downloader("https://www.dropbox.com/scl/fo/5euixwis9o4bz95r0grd6/AP989kq5uM_tntZDXOCBdf4?rlkey=ibgdxnxql0onmovsh4fq9wpkj&dl=1", db_dest)

############################
# import ctg data
############################
#list all files you copied over from dropbox
all_paths <- list.files(db_dest, recursive = TRUE, full.names = TRUE)
#get the paths to excel files of ctg data
ctg_paths <- all_paths[grepl("ctg", all_paths)]
#get the paths to excel files of nn (neural network) data, 
    #to be used in neural network section
nn_paths <- all_paths[grepl("NN", all_paths)]
#remove all_paths to save space
rm(all_paths)
#run ctg2plater on the ctg data
data_path <- here("data", "experiment_output_data")
if (!dir.exists(data_path)) {
  dir.create(data_path)
}
ctg_df <- ctg_paths %>% 
  map(ctg2plater, collect_metadata = TRUE, data_path = data_path) %>%
  bind_rows()

############################
# import nn data
############################
nn_df <- nn_paths %>%
  map(microsNNreader, write_meta = TRUE) %>%
  bind_rows() %>%
  clean_names()

#deal with weird well names for NewNN data
  #for the old NN data, Des was able to identify wells 
nn_df2 <- nn_df %>%
  #isolate the weird well numbers from the NewNN data
  mutate(new_well = case_when(is.na(well) & nn == "NewNN" ~ str_remove_all(str_extract(document, "_[[:digit:]]+.vsi"), "_|.vsi"),
                              TRUE ~ NA)) %>%
  mutate(new_well = as.numeric(new_well))
#next step is to match up the real well numbers with the weird ones
normal_wells <- sort(unique(nn_df2$well[!is.na(nn_df2$well)]))
weird_wells <- sort(unique(nn_df2$new_well[!is.na(nn_df2$new_well)]))
#these lists need to be the exact same length
if (length(normal_wells) == length(weird_wells)) {
  well_lookup <- data.frame(normal_wells, weird_wells)
  
  #left merge nn_df2 on well_lookup
  nn_df3 <- merge(x = nn_df2, y = well_lookup, by.x = "new_well", by.y = "weird_wells", all.x = TRUE)
  
  #now set the values for those previously NA wells!
  nn_df <- nn_df3 %>%
    mutate(well = case_when(is.na(well) ~ normal_wells,
                            TRUE ~ well))
  
} else {
  stop("Wells don't match up between NewNN and old NN data!! Look at your raw data!")
}

############################
#merge on well and merger id, calculate stats and normalize
############################
#merge in metadata from ctg_df (no signal yet)
df0 <- merge(ctg_df[!(names(ctg_df) %in% c("signal", "treatment_duration"))], nn_df, 
             by.x = c("sheet_name", "well"), 
             by.y = c("ctg_merge_id", "well"), all = TRUE,
             suffixes = c("_ctg", "_nn")) %>%
  #filter out empty outer wells
  filter(!(grepl("A|H|01|12", well)))
#now merge in the signal column
df <- merge(df0, ctg_df[,c("well", "sheet_name", "treatment_duration", "signal")], 
            by.x = c("sheet_name", "well", "treatment_duration"), 
            by.y = c("sheet_name", "well", "treatment_duration"), all.x = TRUE,
            suffixes = c("", "_ctg")) %>%
  #group_by plate
    #and the nn u used to analyze that plate!
  group_by(sheet_name, operator, cell_line, treatment_duration, nn) %>%
  mutate(min_hctrl_conc = min(drug_concentration[condition == "DMSO"], na.rm = TRUE)) %>%
  ungroup() %>%
  #make total_live_area column
  mutate(total_live_area_nn = mean_area_mm2_live*object_count_number_live,
         #make control column
         control = case_when(condition == "DMSO" & drug_concentration == min_hctrl_conc ~ "high control",
                             condition == "Staurosporin" ~ "low control",
                             TRUE ~ "sample"),
         control = factor(control, levels = c("sample", "high control", "low control"))) %>%
  #group_by plate
  #and the nn u used to analyze that plate!
  group_by(sheet_name, operator, cell_line, treatment_duration, nn) %>%
  #ctg z'factor calculations
  mutate(hctrl_avg_ctg = mean(signal[control == "high control"], na.rm = TRUE), 
         lctrl_avg_ctg = mean(signal[control == "low control"], na.rm = TRUE),
         hctrl_sd_ctg = sd(signal[control == "high control"], na.rm = TRUE), 
         lctrl_sd_ctg = sd(signal[control == "low control"], na.rm = TRUE),
         #nn z'factor calculations
         hctrl_avg_nn_tla = mean(total_live_area_nn[control == "high control"], na.rm = TRUE), 
         lctrl_avg_nn_tla = mean(total_live_area_nn[control == "low control"], na.rm = TRUE),
         hctrl_sd_nn_tla = sd(total_live_area_nn[control == "high control"], na.rm = TRUE), 
         lctrl_sd_nn_tla = sd(total_live_area_nn[control == "low control"], na.rm = TRUE)) %>%
  ungroup %>%
  #calculate hctrl_avg_ctg and hctrl_avg_nn_tla at day 0
    #so group by plate (and nn), but not day
  group_by(plate_id, operator, cell_line, nn) %>%
  mutate(d0_avg_nn_tla = mean(total_live_area_nn[treatment_duration == 0], na.rm = TRUE)) %>%
  #group by plate AND condition AND the nn u used to analyze the data
  group_by(plate_id, operator, treatment_duration, condition, drug_concentration, cell_line, nn) %>%
  #calculate percent cv
          #ctg %cv calculations
  mutate(mean_signal_ctg = mean(signal, na.rm = TRUE),
         sd_signal_ctg = sd(signal, na.rm = TRUE),
         percent_cv_ctg = (sd_signal_ctg/mean_signal_ctg)*100,
          #nn %cv calculations
         mean_total_live_area_nn = mean(total_live_area_nn, na.rm = TRUE),
         sd_total_live_area_nn = sd(total_live_area_nn, na.rm = TRUE),
         percent_cv_nn_tla = (sd_total_live_area_nn/mean_total_live_area_nn)*100) %>%
  #group by plate and nn
  group_by(sheet_name, operator, cell_line, treatment_duration, nn) %>%
  #now get plate_percent_cv for both ctg and nn
  mutate(plate_percent_cv_ctg = mean(percent_cv_ctg, na.rm = TRUE),
         plate_percent_cv_nn_tla = mean(percent_cv_nn_tla, na.rm = TRUE),
         gr_nn_tla = 2^(log2(total_live_area_nn/d0_avg_nn_tla)/log2(hctrl_avg_nn_tla/d0_avg_nn_tla)) - 1) %>%
  ungroup() %>%
          #ctg percent control, plate zprime
  mutate(percent_control_ctg = (signal - lctrl_avg_ctg)/(hctrl_avg_ctg - lctrl_avg_ctg)*100,
         plate_zprime_ctg = 1 - (3*hctrl_sd_ctg + 3*lctrl_sd_ctg)/abs(hctrl_avg_ctg - lctrl_avg_ctg),
          #nn percent control, plate zprime
         percent_control_nn_tla = (total_live_area_nn - lctrl_avg_nn_tla)/(hctrl_avg_nn_tla - lctrl_avg_nn_tla)*100,
         plate_zprime_nn_tla = 1 - (3*lctrl_sd_nn_tla + 3*hctrl_sd_nn_tla)/abs(lctrl_avg_nn_tla - hctrl_avg_nn_tla)) %>%
  #group by plate and condition and nn
  group_by(sheet_name, condition, drug_concentration, treatment_duration, nn) %>%
  mutate(zfactor_ctg = 1 - (3*sd_signal_ctg + 3*lctrl_sd_ctg)/abs(mean_signal_ctg - lctrl_avg_ctg),
         zfactor_nn_tla = 1 - (3*sd_total_live_area_nn + 3*hctrl_sd_nn_tla)/abs(mean_total_live_area_nn - hctrl_avg_nn_tla)) %>%
  ungroup() %>%
  arrange()

#i want to add suffixes to the columns that came from the nn data vs the ctg data
#because otherwise it's kind of a headache
#i defined a function so I could use lapply and apply it to the list of names and get out a new list 
  #to set the names to
add_suffixes <- function(colname) {
  #add the _ctg suffix if:
    #that name comes from the ctg_df,
  if (colname %in% names(ctg_df) & 
      #it is not in the nn_df, 
      !(colname %in% names(nn_df)) &
      #and it doesn't already have "_ctg" in it
      !grepl("ctg", colname)) {
    output_colname <- paste0(colname, "_ctg")
  #add the _nn suffix if:
              #that name comes from the nn_df,
  } else if (colname %in% names(nn_df) & 
             #it is not in the ctg_df,
             !(colname %in% names(ctg_df)) &
             #and it doesn't already have "_nn" in it
             !grepl("nn", colname)) {
    output_colname <- paste0(colname, "_nn")
  } else {
    #if it's not in the column titles of either of these dataframes,
      #i must have made it after merging, and then I don't want to add a suffix
    output_colname <- colname
  }
  return(output_colname)
}
new_names <- unlist(lapply(names(df), add_suffixes))
names(df) <- new_names
############################
#check that df is merged correctly
############################
checker_df <- df %>%
  dplyr::select(well, sheet_name_ctg, treatment_duration, nn, signal_ctg, object_count_number_dead_nn) %>%
  distinct() %>%
  arrange(sheet_name_ctg, well, treatment_duration, nn)
checker_nn_tla_df <- nn_df %>%
  dplyr::select(well, ctg_merge_id, treatment_duration, nn) %>%
  distinct() %>%
  arrange(ctg_merge_id, well, treatment_duration, nn)
checker_zprime <- df %>%
  dplyr::select(sheet_name_ctg, treatment_duration, seeding_date_ctg, seeding_date_nn, nn, plate_zprime_ctg, plate_zprime_nn_tla)
checker_cv <- df %>%
  dplyr::select(sheet_name_ctg, treatment_duration, , nn, condition_ctg, percent_cv_ctg, plate_percent_cv_ctg, percent_cv_nn_tla, plate_percent_cv_nn_tla)
checker_gr <- df %>%
  dplyr::select(sheet_name_ctg, treatment_duration, nn, d0_avg_nn_tla, gr_nn_tla)
#remove seeding_date_ctg
#because i checked that it was the same as seeding_date_nn
#both columns have the same data, but seeding_date_nn is formatted nice and seeding_date_ctg is not
df <- df %>%
  dplyr::select(-seeding_date_ctg) %>%
  #get nM drug concentration
  #(in the raw data, they are entered as uM)
  mutate(drug_concentration_nm = 1000*drug_concentration_ctg)
#save as RDS
if (!dir.exists(here("data", "cleaned"))) {
  dir.create(here("data", "cleaned"))
}
saveRDS(df, file = here("data", "cleaned", paste0("data_cleaning_output_", Sys.Date(), ".RDS")))

