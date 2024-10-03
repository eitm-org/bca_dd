library(curl)
library(zip)
library(here)

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