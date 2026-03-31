### Title: TOPX 3 column files
### Author: Becca Kazinka, edits by Jessica Arend
### Last updated: 2026.03.30

# Notes: 
# 37 = left button press; 38 = middle (or up) button press; 39 = right button press
# Therefore, correct button presses are:
# AX = cue:37, probe:38
# AY = cue:37, probe:39
# BX = cue:37, probe:39
# BY = cue:37, probe:39

# This script makes 3 column timing files for:
# A Cue
# B Cue
# AX
# AY
# BX
# BY
# Outcome: Correct
# Outcome: Incorrect
# Outcome: Missed
# Cue Errors
# Probe Errors
# Missed Cues
# Missed Probes

#### Prep workspace ####

# clear workspace
rm(list=ls())
# load packages
pacman::p_load(ggplot2, ggsignif, lme4, purrr, DescTools, dplyr, rstatix, psycho, tidyverse, jsonlite)

#### Create function to make 3 column files ####

makeTOPX3Columns <- function(filePath, savePath = "EVs"){
  
  #use jsonlite to read raw data file
  dat <- as.data.frame(stream_in(file(filePath)))
  
  # get the name of just the file
  strName <- strsplit(filePath,"/")
  strName <- unlist(strName)
  fileName <- strName[length(strName)]
  
  # from the file name, extract the subject and session names
  strings <- unlist(strsplit(fileName,"_"))
  sub <- strings[1]
  ses <- strings[3]
  
  # calculate duration
  # there will be NAs because the words "timeout" cannot be turned into numbers
  dat$duration <- dat$end_time - dat$start_time
  suppressWarnings({
    dat$response <- as.numeric(as.character(dat$response))
  })
  
  # calculate onset (in seconds) relative to system start
  dat$onset <- ifelse(dat$blocksCompleted == 0,
                      (dat$start_time - dat$System_Start_Time_1[1])/1000,
                      (dat$start_time - dat$System_Start_Time_2[1])/1000)
 
  # extract cue and probe names from image files
  dat$cue_type <- substring(dat$cue,6,6)
  dat$probe_type <- substring(dat$probe,6,6)
  
  # make savePath folder (location of output) if necessary
  if (!dir.exists(savePath)){
    dir.create(savePath)
    print("making new output folder")
  }   
  
  #### Loop through two blocks ####
  
  # looping skips having to split the blocks in half and process separately
  for (ii in c(0,1)) {
    run <- ifelse(ii==0, "run1","run2")
    task <- ifelse(ii==0,"task-topx1_rec-NORDIC_run-1_part-mag_bold",
                   "task-topx2_rec-NORDIC_run-2_part-mag_bold")

    # make new task folder if necessary
    if (!dir.exists(paste(savePath,task,sep="/"))){
      dir.create(paste(savePath,task,sep="/"))
      print("making new task folder")
    }
    
    # make new subject folder if necessary
    if (!dir.exists(paste(savePath,task, paste0("sub-",sub),sep="/"))){
      dir.create(paste(savePath,task,paste0("sub-",sub),sep="/"))
      print("making new subject folder")
    }
    
    # make new session folder if necessary
    if (!dir.exists(paste(savePath,task, paste0("sub-",sub),paste0("ses-",ses),sep="/"))){
      dir.create(paste(savePath,task, paste0("sub-",sub),paste0("ses-",ses),sep="/"))
      print("making new session folder")
    }
    
    #### Split into two halves ####
    
    # last outcome_focus might be in the wrong block, no biggie
    dat_half1 <- filter(dat, blocksCompleted == ii) 

    # remove the system start info and reshape to wider format
    dat_half <- dat_half1 %>%
      filter(!is.na(start_time))%>%
      select(c("response_time","response","response_state","duration","onset","outcome",
               "blocksCompleted","trialsCompleted","cue_type","probe_type")) %>%
      pivot_wider(names_from = response_state, values_from = c(response_time,response,duration, onset, outcome))
    
    # set column weights to 1
    dat_half$tcf_weight <- rep(c(1), times=nrow(dat_half))
    
    #### Correct Cue 3cols ####
    
    ##### Prep components #####
    # keep only trials with correct cue responses
    dat_correct <- filter(dat_half, response_cue == 37 | response_cue_focus == 37)
    # calculate duration (in seconds) for cue + cue focus
    dat_correct$tcf_dur <- (dat_correct$duration_cue + dat_correct$duration_cue_focus)/1000
   
    ###### A Cue ######
    dat_acue <- filter(dat_correct, cue_type == "A")
    # filter out probe errors
    dat_acue <- filter(dat_acue, (probe_type == "X" & (response_probe == 38 | 
                                                         response_probe_focus == 38)) | 
                                (probe_type == "Y" & (response_probe == 39 | 
                                                        response_probe_focus == 39)))
    # write 3col file for A cue
    dat_acue <- as.data.frame(cbind(dat_acue$onset_cue, dat_acue$tcf_dur, dat_acue$tcf_weight))
    write.table(dat_acue, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,
                                      "/",sub,"_",ses,"_Acue_",run,".txt"), 
                row.names = FALSE, col.names = FALSE)
    
    ######  B Cue  ###### 
    dat_bcue <- filter(dat_correct, cue_type == "B")
    # filter out probe errors
    dat_bcue <- filter(dat_bcue, response_probe == 39 | response_probe_focus == 39)
    # write 3col file for B cue
    dat_bcue <- as.data.frame(cbind(dat_bcue$onset_cue, dat_bcue$tcf_dur, dat_bcue$tcf_weight))
    write.table(dat_bcue, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,
                                      "/",sub,"_",ses,"_Bcue_",run,".txt"), 
                row.names = FALSE, col.names = FALSE)
    
    
    #### Correct Probe 3cols ####
    
    ##### Prep components #####
    # get the timing based on reaction time, not duration:
    # if responded during probe, keep probe RT
    # if responded during probe focus, add probe duration to probe focus RT 
    dat_half$response_time_probe_adj <- (ifelse(is.na(dat_half$response_time_probe), 
                                                    dat_half$duration_probe +
                                                      dat_half$response_time_probe_focus, 
                                                    dat_half$response_time_probe))/1000
    dat_ax <- filter(dat_half, cue_type == "A" & probe_type =="X")
    dat_ay <- filter(dat_half, cue_type == "A" & probe_type =="Y")
    dat_bx <- filter(dat_half, cue_type == "B" & probe_type =="X")
    dat_by <- filter(dat_half, cue_type == "B" & probe_type =="Y")
    
    # need correct trials only (calculate based on button presses, as "outcome" col is not always accurate)
    
    ###### AX ######
    # filter out errors
    dat_ax <- filter(dat_ax, (response_probe == 38 | response_probe_focus == 38) &
                       (response_cue == 37 | response_cue_focus == 37))
    # write 3col file for AX
    tcf_ax <- as.data.frame(cbind(dat_ax$onset_probe, dat_ax$response_time_probe_adj, dat_ax$tcf_weight))
    write.table(tcf_ax, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",
                                    sub,"_",ses,"_AX_",run,".txt"), row.names = FALSE, col.names = FALSE)
    
    ###### AY ###### 
    # filter out errors
    dat_ay<- filter(dat_ay, (response_probe == 39 | response_probe_focus == 39) &
                      (response_cue == 37 | response_cue_focus == 37))
    # write 3col file for AY
    tcf_ay <- as.data.frame(cbind(dat_ay$onset_probe, dat_ay$response_time_probe_adj, dat_ay$tcf_weight))
    write.table(tcf_ay, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",
                                    sub,"_",ses,"_AY_",run,".txt"), row.names = FALSE, col.names = FALSE)
    
    ###### BX ###### 
    # filter out errors
    dat_bx <- filter(dat_bx, (response_probe == 39 | response_probe_focus == 39) &
                       (response_cue == 37 | response_cue_focus == 37))
    # write 3col file for BX
    tcf_bx <- as.data.frame(cbind(dat_bx$onset_probe, dat_bx$response_time_probe_adj, dat_bx$tcf_weight))
    write.table(tcf_bx, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",
                                    sub,"_",ses,"_BX_",run,".txt"), row.names = FALSE, col.names = FALSE)
    
    ###### BY ###### 
    # filter out errors
    dat_by <- filter(dat_by, (response_probe == 39| response_probe_focus == 39) &
                       (response_cue == 37 | response_cue_focus == 37))
    # write 3col file for BY
    tcf_by <- as.data.frame(cbind(dat_by$onset_probe, dat_by$response_time_probe_adj, dat_by$tcf_weight))
    write.table(tcf_by, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",
                                    sub,"_",ses,"_BY_",run,".txt"), row.names = FALSE, col.names = FALSE)
    
    
    #### Outcome 3cols ####
    
    ##### Prep components ###### 
    dat_otc_corr <- filter(dat_half,  outcome_outcome == "Correct")
    dat_otc_incor <-  filter(dat_half,  outcome_outcome == "Incorrect")
    dat_otc_miss <-  filter(dat_half,  outcome_outcome == "Respond Faster")
    
    ###### Correct ######
    # write 3col file for correct outcomes ("+1 point")
    tcf_otc_corr <- as.data.frame(cbind(dat_otc_corr$onset_outcome, 
                                        dat_otc_corr$duration_outcome/1000, 
                                        dat_otc_corr$tcf_weight))
    write.table(tcf_otc_corr, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,
                                          "/",sub,"_",ses,"_Outcome_win_",run,".txt"), 
                row.names = FALSE, col.names = FALSE)
    
    ###### Incorrect ######
    # combine 3cols for incorrect outcomes ("+0 points")
    tcf_otc_incor <- as.data.frame(cbind(dat_otc_incor$onset_outcome, 
                                         dat_otc_incor$duration_outcome/1000, 
                                         dat_otc_incor$tcf_weight))
    # create dummy variable if no errors were made
    if(nrow(tcf_otc_incor) ==0) {tcf_otc_incor <- data.frame(0,0,0)}
    # write 3col file for incorrect outcomes ("+0 points")
    write.table(tcf_otc_incor, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,
                                           "/",sub,"_",ses,"_Outcome_loss_",run,".txt"),
                row.names = FALSE, col.names = FALSE)
    
    ###### Missed ###### 
    # combine 3cols for missed outcomes ("respond faster")
    tcf_otc_miss <- as.data.frame(cbind(dat_otc_miss$onset_outcome, 
                                        dat_otc_miss$duration_outcome/1000, 
                                        dat_otc_miss$tcf_weight))
    # create dummy variable if no trials were missed
    if(nrow(tcf_otc_miss) ==0) {tcf_otc_miss <- data.frame(0,0,0)}
    # write 3col file for missed outcomes ("respond faster")
    write.table(tcf_otc_miss, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",
                                          sub,"_",ses,"_Outcome_miss_",run,".txt"), 
                row.names = FALSE, col.names = FALSE)
    
    
    #### Error 3cols ####
    
    ###### Cue errors ######
    # Missed cues will be modeled separately from cue errors
    # Correct trials must have both cue and probe correct
    # So incorrect trials will have either cue or probe (or both) incorrect
    
    # cue trials with only cue errors
    # dat_cue_err <- filter(dat_half, response_cue == 38| response_cue == 39 |
                            #response_cue_focus == 38| response_cue_focus == 39)
    
    # cue trials with cue errors or probe errors;
    dat_cue_err <- filter(dat_half, (response_cue == 38 | response_cue == 39 |
                                       response_cue_focus == 38 | response_cue_focus == 39) | # cue error
                            (cue_type== "A" & probe_type=="X" & 
                               (response_probe == 37 | response_probe == 39 |
                                  response_probe_focus == 37 | response_probe_focus == 39)) | # probe error
                            (cue_type == "A" & probe_type=="Y" &
                               (response_probe == 37 | response_probe == 38 |
                                  response_probe_focus == 37 | response_probe_focus == 38)) | # probe error
                            (cue_type == "B" & probe_type=="X" &
                               (response_probe == 37 | response_probe == 38 |
                                  response_probe_focus == 37 | response_probe_focus == 38)) | # probe error
                            (cue_type == "B" & probe_type=="Y" &
                               (response_probe == 37 | response_probe == 38 |
                                  response_probe_focus == 37 | response_probe_focus == 38))) # probe error
    dat_cue_err$tcf_cue_err_dur <- ((dat_cue_err$duration_cue + dat_cue_err$duration_cue_focus)/1000)
    # combine 3cols for cue errors
    tcf_cue_err <- as.data.frame(cbind(dat_cue_err$onset_cue, 
                                       dat_cue_err$tcf_cue_err_dur, 
                                       dat_cue_err$tcf_weight))
    # order by onset time
    tcf_cue_err <- arrange(tcf_cue_err, V1)
    # create dummy variable if no errors were made
    if(nrow(tcf_cue_err) == 0) {tcf_cue_err <- data.frame(0,0,0)}
    # write 3col file for cue errors
    write.table(tcf_cue_err, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",
                                         sub,"_",ses,"_Error_cue_",run,".txt"), 
                row.names = FALSE, col.names = FALSE)
    
    ###### Probe errors ######
    # probe trials with probe errors (same errors as above, just during the probe time):
    dat_probe_err <- filter(dat_half, (cue_type== "A" & probe_type=="X" &
                                             (response_probe == 37 | response_probe == 39 |
                                                response_probe_focus == 37 | response_probe_focus == 39)) |
                                (cue_type == "A" & probe_type=="Y" &
                                   (response_probe == 37 | response_probe == 38 |
                                      response_probe_focus == 37 | response_probe_focus == 38)) |
                                (cue_type == "B" & probe_type=="X" &
                                   (response_probe == 37 | response_probe == 38 |
                                      response_probe_focus == 37 | response_probe_focus == 38)) |
                                (cue_type == "B" & probe_type=="Y" &
                                   (response_probe == 37 | response_probe == 38 |
                                      response_probe_focus == 37 | response_probe_focus == 38)))
    # combine 3cols for probe errors
    tcf_probe_err <- as.data.frame(cbind(dat_probe_err$onset_probe,
                                         dat_probe_err$response_time_probe_adj, 
                                         dat_probe_err$tcf_weight))
    # order by onset time
    tcf_probe_err <- arrange(tcf_probe_err, V1)
    # create dummy variable if no errors were made
    if(nrow(tcf_probe_err) ==0) {tcf_probe_err <-  data.frame(0,0,0)}
    # write 3col file for probe errors
    write.table(tcf_probe_err, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,
                                           "/",sub,"_",ses,"_Error_probe_",run,".txt"),
                row.names = FALSE, col.names = FALSE)
    
    #### Missed Trial 3cols ####
    
    ###### Missed Cues ######
    # select trials with missed cue or probe
    dat_cue_miss <- filter(dat_half, is.na(response_cue) & is.na(response_cue_focus) |
                             is.na(response_probe) & is.na(response_probe_focus))
    # catches trials where the two halves were not split properly
    dat_cue_miss <- filter(dat_cue_miss, !is.na(onset_cue))
    dat_cue_miss$tcf_cue_miss_dur <- ((dat_cue_miss$duration_cue + dat_cue_miss$duration_cue_focus)/1000)
    # combine 3cols for missed cues
    tcf_cue_miss <- as.data.frame(cbind(dat_cue_miss$onset_cue, 
                                        dat_cue_miss$tcf_cue_miss_dur, 
                                        dat_cue_miss$tcf_weight))
    # order by onset time
    tcf_cue_miss <- arrange(tcf_cue_miss, V1)
    # create dummy variable if no cues were missed
    if(nrow(tcf_cue_miss) ==0) {tcf_cue_miss <-  data.frame(0,0,0)}
    # write 3col file for missed cues
    write.table(tcf_cue_miss, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,
                                          "/",sub,"_",ses,"_Missed_cue_",run,".txt"), 
                row.names = FALSE, col.names = FALSE)
    
    ###### Missed Probes ######
    # select trials with missed cue or probe
    dat_probe_miss <- filter(dat_half, is.na(response_probe) & is.na(response_probe_focus) |
                             is.na(response_cue) & is.na(response_cue_focus))
    # catches trials where the two halves were not split properly
    dat_probe_miss <- filter(dat_probe_miss, !is.na(onset_probe))
    dat_probe_miss$tcf_probe_miss_dur  <- ((dat_probe_miss$duration_probe +
                                            dat_probe_miss$duration_probe_focus)/1000)
    # combine 3cols for missed probes
    tcf_probe_miss <- as.data.frame(cbind(dat_probe_miss$onset_probe,
                                            dat_probe_miss$tcf_probe_miss_dur,
                                            dat_probe_miss$tcf_weight))
    # order by onset time
    tcf_probe_miss <- arrange(tcf_probe_miss, V1)
    # create dummy variable if no probes were missed
    if(nrow(tcf_probe_miss) == 0) {tcf_probe_miss <-  data.frame(0,0,0)}
    # write 3col file for missed probes
    write.table(tcf_probe_miss, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,
                                              "/",sub,"_",ses,"_Missed_probe_",run,".txt"),
                row.names = FALSE, col.names = FALSE)
    
  } #end of run loop
} #end of function

#### Using the function - Example ####

# Port in variables of interest
# source("STEP_TOPX_make_3cols.R") # don't need this step unless running function in separate script
#filePath <- "TOPX behavioral RAW/SP1001_TOPX_6M_20221019.txt"
#makeTOPX3Columns(filePath)
# defaults output/"save path" to "EVs/task-topx1_rec-NORDIC_run-1_part-mag_bold" (run 1) & "EVs/task-topx2_rec-NORDIC_run-2_part-mag_bold", (run 2)
# but can enter an alternate savePath (replaces "EVs"), e.g., makeTOPX3Columns(filePath, savePath = newPath)

#input: ~/Library/CloudStorage/Box-Box/STEP-P50/Reporting & Data Submission/MRI Task Data/TOPX Scanner Behavioral Data + Analyses/TOPX behavioral RAW
# output: ~/Library/CloudStorage/Box-Box/STEP-P50/Reporting & Data Submission/MRI Task Data/TOPX Scanner Behavioral Data + Analyses/EVs

####  Using the function - Actual ######

# set working directory
setwd("/Users/arend103/Library/CloudStorage/Box-Box/STEP-P50/Reporting & Data Submission/MRI Task Data/TOPX Scanner Behavioral Data + Analyses")

# read in all the file names of TOPX raw data
txt_files <- list.files("/Users/arend103/Library/CloudStorage/Box-Box/STEP-P50/Reporting & Data Submission/MRI Task Data/TOPX Scanner Behavioral Data + Analyses/TOPX behavioral RAW",
  pattern = "\\.txt$",
  full.names = TRUE)

# filter to only keep 6-Month files
txt_files_6M <- grep("6M", txt_files, value = TRUE)

# Save to a text file, one per line
writeLines(subject_ids_sub, "subject_ids_sub.txt")

# run iteratively
for (filePath in txt_files_6M) {
  message("Processing: ", basename(filePath))
  makeTOPX3Columns(filePath, savePath = "EVs")
}

# save list of only subject IDs, for subject list in fMRI analyses 
# (e.g., subs_6m_YYYYMMDD_basic.txt)
fn <- basename(txt_files_6M)
(subject_ids <- sub("^(SP[0-9]{4}).*$", "\\1", fn))
subject_ids <- paste0("sub-", subject_ids)
writeLines(subject_ids,
           "/Users/arend103/Documents/umn_work/analyses/topx_analyses/topx_3col/topx_subject_ids_6m.txt")

# # as needed, remove single files that break the script
# txt_files_6M <- txt_files_6M[basename(txt_files_6M) != "SP1109_TOPX_6M_20240804 (1).txt"]
# txt_files_6M <- txt_files_6M[basename(txt_files_6M) != "SP2087_TOPX_6M_20231117.txt"]
# 
# # as needed, remove files already run by filtering by ID number
# fn <- basename(txt_files_6M)
# sp_num <- as.numeric(sub("SP([0-9]{4}).*", "\\1", fn))
# txt_files_6M <- txt_files_6M[sp_num >= 1109]