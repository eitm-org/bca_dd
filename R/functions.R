library(curl)
library(zip)
library(here)
library(maditr)

#you can put dropbox_downloader.R in the scripts folder of your project
#and use it to copy over dropbox files for your project
#if you do it like this, you won't have to open up old files to read them in!
#it's similar to rdrop2, but rdrop2 needs a manager and might die soon :/
#also imo this is marginally easier to do


#first, get the link to the dropbox folder you want to read in 
#it will look like this: https://www.dropbox.com/sh/v8526atx7vfj7vb/AABFVp2nSFV1cCzNUppo34pHa?dl=0
#change the 0 at the end of your link to a 1
#like this: https://www.dropbox.com/sh/v8526atx7vfj7vb/AABFVp2nSFV1cCzNUppo34pHa?dl=1
#this means that instead of taking you to the folder, this link will automatically download the folder
#this example is the dropbox link to  'Dropbox (EITM)', 'EITM AR SPRC 2022 Docs', 'Data' folder
# dropbox_link <- "https://www.dropbox.com/sh/v8526atx7vfj7vb/AABFVp2nSFV1cCzNUppo34pHa?dl=1"

#local dest must be a character string that refers to a directory within your project
dropbox_downloader <- function(dropbox_link = dropbox_link, local_dest = local_dest) {
  
  #check that the dropbox_link ends with a 1
  if (substr(dropbox_link, nchar(dropbox_link), nchar(dropbox_link) + 1) != "1") {
    stop("ERROR: Your dropbox link does not end in a 1!\nRemember to change the 0 at the end of your dropbox link to a 1!")
    
  }
  
  #check if the local_dest directory exists
  if(!dir.exists(local_dest)) {
    #check and see if it exists again
    if (!dir.exists(local_dest)) {
      #and if it doesn't, create it
      dir.create(local_dest)
    }
  }
  
  zip_file_path <- file.path(local_dest, "db_download.zip")
  #this sets the folder where you want to download your dropbox files to
  #I usually put my dropbox input files in a subdirectory of the Data_Input folder
  #this line downloads the dropbox folder you designated in the link into the folder you designated in the destination_dropbox statement
  message(curl::multi_download(url = dropbox_link, destfile = zip_file_path))
  
  #unzip the file
  message(zip::unzip(zipfile = zip_file_path, exdir = here(local_dest, "unzipped")))
  
}

# read in CTG xls files, create a .csv file. read in using plater
ctg2plater <- function(xls_path, collect_metadata = TRUE, range = "A1:M50", data_path = here()) {
  require(readxl)
  require(dplyr)
  require(openxlsx)
  
  #collect the sheets containing plate data
  sheets_names <- excel_sheets(path = xls_path)
  all_sheets <- lapply(sheets_names, function(x) read_excel(xls_path, sheet = x, range = range))
  names(all_sheets) <- sheets_names
  plate_names <- names(all_sheets[names(all_sheets) %in% c("Formatting Parameters", "Sheet1", "Sheet2", 
                                                           "metadata", "plate template", "template dictionary", 
                                                           # "^(Raw Data)", "^(Working Data)", 
                                                           "Layout", 
                                                           "Associated Files", "empty_plate template", "validation") == FALSE])
  
  #collect the sheets containing metadata
  meta_sheet_names <- names(all_sheets[names(all_sheets) %in% c("metadata", "Metadata")])
  #read in metadata & write to living csv
  meta_df <- read_xlsx(path = xls_path, sheet = meta_sheet_names)
  #in the first spreadsheet, 
  meta_df <- meta_df %>% 
    distinct(plate_id, .keep_all = TRUE)
  
  if(collect_metadata == TRUE) {
    if (!file.exists(file.path(data_path, "ctg_metadata_tracker.csv"))) {
      file.create(file.path(data_path, "ctg_metadata_tracker.csv"))
    }
    write_csv(meta_df, file = file.path(data_path, "ctg_metadata_tracker.csv"), append = TRUE)
  }
  #add csvs folder if it's not there
  if (!dir.exists(file.path(data_path, "csvs"))) {
    dir.create(file.path(data_path, "csvs"))
  }
  #write plater data to csv
  sapply(
    plate_names,
    function(x) write_csv(all_sheets[[x]], file = paste0(file.path(data_path, "csvs/"), x, ".csv"))
  )
  # read in all the csv that were just written using plater
  sapply(plate_names, function(x) check_plater_format(file = paste0(file.path(data_path, "csvs/"), x, ".csv")))
  plate_df <- lapply(plate_names, function(x) read_plate(
    file = paste0(file.path(data_path, "csvs/"), x, ".csv"), # full path to the .csv file
    well_ids_column = "Well" # name to give column of well IDs
  ))
  #extract plate_id from the sheet names so we can merge the plate data w metadata
  #force ImageQC_Corrected to be character type so that we can combine plates that include 'Omit's
  for (i in 1:length(plate_names)) {
    plate_df[[i]]$plate_id <- meta_df$plate_id[i]
    plate_df[[i]]$ImageQC_Corrected <- as.character(plate_df[[i]]$ImageQC)
    plate_df[[i]]$sheet_name <- plate_names[[i]]
  }
  #bind plates into one dataframe
  bound_data <- bind_rows(plate_df)
  #merge plate data with metadata
  experiment_df <- full_join(bound_data, meta_df)
  experiment_df <- experiment_df %>% janitor::clean_names()
  
  return(experiment_df)
}

#read in microscopy neural network data, extract info from document name, write metadata to living csv
microsNNreader <- function(xls_path, write_meta = TRUE) {
  #collect raw data from NN files
  sheets_names <- excel_sheets(path = xls_path)
  all_sheets <- lapply(sheets_names, function(x) read_excel(xls_path, sheet = x))
  names(all_sheets) <- sheets_names
  nn_names <- names(all_sheets[names(all_sheets) %in% c("Class Measurements D0", "Class Measurement D0","Class Measurements D3", "Class Measurements D5", "Class Measurements",
                                                        "Class Measurements D1")])
  
  nn_df <- lapply(
    nn_names,
    function(x) read_excel(path = xls_path, sheet = x)
  )
  
  #bind into large dataframe
  bound_data <- bind_rows(nn_df)
  
  #get filename
  filename <- basename(xls_path)
  #get seeding date from filename
  s_date <- stringr::str_extract(filename, "[[:digit:]]{8}")
  
  #clean large df. extract well, treatment treatment_duration, and experiment date info from the document name
  experiment_data <- bound_data %>%
    janitor::clean_names() %>%
    filter(!is.na(document)) %>%
    #plug in seeding date from filename
    mutate(seeding_date = s_date) %>%
    #extract well number from the "document" column
    extract(col = document, into = "well", regex = "([A-H]\\d{2})", remove = FALSE) %>%
    #extract metadata from the "document" column
    extract(col = document, into = "treatment_duration", regex = "Image(\\d)", remove = FALSE) %>%
    extract(col = document, into = "ctg_merge_id", regex = "_(\\d{6})_", remove = FALSE) %>%
    #clean up extracted columns
    mutate(treatment_duration = str_remove_all(treatment_duration, "Image"),
           ctg_merge_id = str_remove_all(ctg_merge_id, "_"),
           seeding_date = as.Date(seeding_date, format = "%Y%m%d"))
  #unmelt live/dead
  experiment_data2 <- dcast(experiment_data, ...~ object_class, 
                           #value.vars are all numeric variables
                           value.var = names(experiment_data[unlist(lapply(experiment_data, is.numeric), use.names = FALSE)]),
                           id.vars = names(experiment_data[!unlist(lapply(experiment_data, is.numeric), use.names = FALSE)]))
  #dcast used sum to aggregate the values, but there was only one value to aggregate
  #so we can remove the word sum from all the names in the df
  names(experiment_data2) <- str_remove_all(names(experiment_data2), "_sum")
  
  return(experiment_data2)
}

does_it_model <- function(data) {
  tryCatch({
    drm(yvar ~ drug_concentration_ctg, data = data, fct = LL.4(names = c("hill", "min_value", "max_value", "ec_50")))
    return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}

lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

#to plot nn and ctg data from multiple cell lines on the same plot
plot_drcs <- function(condition,
                      experiment_df,
                      colname = "percent_control_ctg",
                      y_limit = 200,
                      subtitle = "",
                      col_list = c()) {
  require(drc)
  #make a blank plot to fill in with dose response curves
  xmin <- min(experiment_df$drug_concentration_ctg, na.rm = TRUE)
  if (xmin == 0) {
    xmin <- min(experiment_df$drug_concentration_ctg[experiment_df$drug_concentration_ctg != 0], na.rm = TRUE)
    xmin <- xmin/2
  }
  # breaks <- sort(unique(experiment_df$drug_concentration_ctg))
  # breaks <- c(breaks[1], breaks[round(length(breaks)/3)], breaks[2*round(length(breaks)/3)], breaks[length(breaks)])
  the_plot <- ggplot() +
    theme_bw() +
    scale_x_continuous(trans='log2', limit = c(xmin, max(experiment_df$drug_concentration_ctg, na.rm = TRUE))) +
    scale_y_continuous(limit = c(0, y_limit))
  #rename colname in df
  experiment_df["yvar"] <- experiment_df[colname]
  df_nested <- experiment_df %>%
    filter(!is.na(yvar)) %>%
    #group by plate and cell line
    group_by(sheet_name_ctg, seeding_date_nn, cell_line_ctg) %>%
    nest() %>%
    mutate(does_it_model_col = map_lgl(data, function(data) does_it_model(data)))
  
  df_no_model <- df_nested %>% 
    filter(does_it_model_col == FALSE)  %>%
    mutate(model = NULL, 
           rse = NULL, 
           # emax = vapply(data, function(df) min(df$sig_mean, na.rm = TRUE), numeric(1)),
           plate_cv = vapply(data, function(df) mean(df$percent_cv, na.rm=TRUE), numeric(1)))
  
  df_nested <- df_nested %>%
    filter(does_it_model_col == TRUE) %>%
    mutate(model = lapply(data, function(df) drm(yvar ~ drug_concentration_ctg, fct = LL.4(names = c("hill", "min_value", 'max_value', "ec_50")), data = df, control = drmc(errorm = FALSE))),
           rse = vapply(model, function(model) summary(model)$rseMat[[1]], numeric(1)))
           # emax = vapply(data, function(df) min(df$sig_mean, na.rm = TRUE), numeric(1)),
           # plate_cv = vapply(data, function(df) mean(df$percent_cv, na.rm=TRUE), numeric(1)))
  
  abs_EDs <- lapply(df_nested$model, function(model) ED(object = model, respLev = 50, type = 'absolute', interval = 'delta', display = FALSE))
  df_nested[c("absolute_ec50", "absolute_ec50_se", "absolute_lower_ci", "absolute_upper_ci")] <- do.call(rbind, abs_EDs)
  
  rel_EDs <- lapply(df_nested$model, function(model) ED(object = model, respLev = 50, type = 'relative', interval = 'delta', display = FALSE))
  df_nested[c("rel_ec50", "rel_ec50_se", "relative_lower_ci", "relative_upper_ci")] <- do.call(rbind, rel_EDs)
  
  coefficients <- lapply(df_nested$model, function(model) summary(model)$coefficients[,1])
  df_nested[c("hill", "einf", "max", "rel_ec50_2")] <- do.call(rbind, coefficients)
  
  df_nested <- bind_rows(df_nested, df_no_model)
  
  yheight <- 10

  for (i in 1:nrow(df_nested)) {
    this_date <- unlist(df_nested[i,"seeding_date_nn"])
    this_cell <- df_nested[i, "cell_line_ctg"]
    this_exp <- df_nested[i, "sheet_name_ctg"]
    if (length(col_list) == 0) {
      this_color <- "black"
      print("do you want the curves to be different colors? then pass through col_list!")
    } else {
      this_color <- col_list[[toString(this_exp)]]
    }
    
    this_model <- df_nested[i, "model"][[1]][[1]]
    this_df <- df_nested[i, "data"][[1]][[1]]
    #if there's only one concentration of this drug, we don't have to plot it
    if (length(unique(this_df$drug_concentration_ctg)) <= 1) {
      next
    }
    this_ic50 <- unlist(df_nested[i, "absolute_ec50"])
    this_rse <- df_nested[i, "rse"]

    does_it_model <- as.logical(df_nested[i, "does_it_model_col"])
    
    point.data <- data.frame(
      y = this_df$yvar,
      x = this_df$drug_concentration_ctg
    )
    
    if (does_it_model == TRUE) {
      curve.data.drc <- PR(
        this_model,
        xVec = seq(to = min(this_df[this_df$drug_concentration_ctg > 0, "drug_concentration_ctg"], na.rm = TRUE), 
                   from = max(this_df$drug_concentration_ctg, na.rm = TRUE),
                   length.out = 1000)
      )
      curve.data <- data.frame(
        y = as.vector(unlist(curve.data.drc)),
        x = as.numeric(names(curve.data.drc))
      )
      
      the_plot <- the_plot +
        geom_line(data = curve.data, aes(x = x, y = y), color = this_color)
        # geom_vline(xintercept = this_ic50[[1]], color = this_color) +
        # annotate("text", x = this_ic50[[1]], y = -5, label = round(this_ic50[[1]], digits = 3), color = this_color, size = 7.5)       
    } else {
      the_plot <- the_plot +
        geom_smooth(data = this_df, linewidth = 1, mapping = aes(x = drug_concentration_ctg, y = yvar), color = this_color, se = FALSE)
    }
    
    the_plot <- the_plot +
      annotate("text", x = min(this_df$drug_concentration_ctg), y = yheight, label = this_exp, hjust = 0, color = this_color, size = 3) +
      geom_point(data = this_df, size = 2, mapping = aes(x = drug_concentration_ctg, y = yvar), color = this_color, alpha = .5)
    
    yheight <- yheight + 15
  }
  #make the column name nice to add it as the yaxis label
  nice_colname <- str_replace_all(str_to_title(str_replace_all(colname, "_", " ")), "Ctg", "CTG")
  #replace sstrings in nn title
  nice_colname <- str_replace_all(str_replace_all(nice_colname, "Nn", "NN"), "Tla", "Total Live Area")
  nice_colname <- paste(nice_colname, "(%)")
  the_plot <- the_plot +
    labs(x = parse(text = "Dose (μM)"),
         y = nice_colname,
         subtitle = paste0("Cell Line: ", unique(experiment_df$cell_line_ctg))) +
    ggtitle(paste0("Condition: ", condition)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    geom_hline(yintercept = 100, linetype = 'dashed', color = 'grey', size = 0.5) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey', size = 0.5)
  
  efficacy.table <- df_nested %>%
    dplyr::select(-c(data, model, does_it_model_col))
  names(efficacy.table) <- str_to_title(str_replace_all(names(efficacy.table), "_", " "))
  return_list <- list("condition" = condition,
                      "plot" = the_plot,
                      "table" = efficacy.table) 
  return(return_list)
}

#input = dataframe with only 4 variables
  #one variable for the x value (ctg)
  #one variable for the y value (neural network)
  #sheet_name
  #treatment duration
make_ctg_v_nn_plot <- function(df, x_var_nice = NA, y_var_nice = NA) {
  trt_dur <- max(df$treatment_duration)
  df <- df %>%
    filter(treatment_duration == trt_dur)
  #get nice variable names
  if (is.na(x_var_nice)) {
    x_var_nice <- str_to_title(str_replace_all(names(df)[[1]], "_", " "))
  }
  if (is.na(y_var_nice)) {
    y_var_nice <- str_to_title(str_replace_all(names(df)[[2]], "_", " "))
  }
  #rename variables
  df["xvar"] <- df[,1]
  df["yvar"] <- df[,2]

  #linear model
  linmod <- lm(yvar ~ xvar, data = df)
  linmod_summ <- summary(linmod)
  coef <-  round(linmod_summ$coefficients[2,1], 3)
  pval <-  linmod_summ$coefficients[2,4]
  rsq <- round(linmod_summ$r.squared, 3)
  if (pval < .001) {
    pval <- formatC(pval, format = "e", digits = 3)
  } else {
    pval <- round(pval, 3)
  }
  #plot
  ctgvnn_plot <- ggplot(data = df, aes(x = xvar, y = yvar)) +
    geom_smooth(color = "black", method = "lm") +
    geom_point(size = 3, alpha = .5, aes(color = sheet_name)) +
    theme_bw() +
    scale_color_viridis_d(end = .8) +
    ggtitle("Neural Network vs CTG") +
    labs(subtitle = paste("Treatment Duration =", trt_dur, "Days"),
         color = "Plate ID") +
    xlab(x_var_nice) +
    ylab(y_var_nice) +
    geom_label(aes(x = max(df$xvar, na.rm = TRUE), y = min(yvar, na.rm = TRUE)), label = paste("Pval:", pval), hjust = 1, vjust = .5) +
    geom_label(aes(x = max(df$xvar, na.rm = TRUE), y = min(yvar, na.rm = TRUE)), label = paste("R^2:", toString(rsq)), hjust = 1, vjust = -0.65, parse = TRUE) +
    geom_label(aes(x = max(df$xvar, na.rm = TRUE), y = min(yvar, na.rm = TRUE)), label = paste("Coef:", coef), hjust = 1, vjust = -2.3) +
    xlim(c(min(df$xvar, na.rm = TRUE), max(df$xvar, na.rm = TRUE))) +
    ylim(c(min(df$yvar, na.rm = TRUE), max(df$yvar, na.rm = TRUE)))
  
  return(ctgvnn_plot)
}
