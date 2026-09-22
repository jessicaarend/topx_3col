pacman::p_load(ggplot2, ggsignif, lme4, purrr, DescTools, dplyr, rstatix, psycho, tidyverse, jsonlite)

#### load spreadsheet ####

#need to port in variables of interest
#example:   
#source("STEP_TOPX_make_3cols.R")
#filePath <- "TOPX behavioral RAW/SP1105_TOPX_BL2_20230913.txt"
#makeTOPX3Columns(filePath)
#defaults path to "EVs/task-topx1_rec-NORDIC_run-1_part-mag_bold" & "EVs/task-topx2_rec-NORDIC_run-2_part-mag_bold",
#but you can enter an alternate savePath (replaces "EVs"). 
#makeTOPX3Columns(filePath, savePath = newPath)

makeTOPX3Columns <- function(filePath, savePath = "EVs"){
  behdat.df <- as.data.frame(stream_in(file(filePath))) #use jsonlite to read raw data file
  
  #get the name of just the file
  strName <- strsplit(filePath,"/")
  strName <- unlist(strName)
  fileName <- strName[length(strName)]
  
  #from the file name, extract the subject and session names
  strings <- unlist(strsplit(fileName,"_"))
  sub <- strings[1]
  ses <- strings[3]
  
  behdat.df$duration <-behdat.df$end_time - behdat.df$start_time
  suppressWarnings({
    behdat.df$response <- as.numeric(as.character(behdat.df$response))
  }) #there will be NAs because the words "timeout" cannot be turned into numbers
  
  if(is.null(behdat.df$System_Start_Time_1)){
    if(!is.null(behdat.df$System_Start_Time)){
      behdat.df$System_Start_Time_1[1]<-behdat.df$System_Start_Time[1]
    }else{
      print("System start time 1 missing")
    behdat.df$System_Start_Time_1[1]<-behdat.df$start_time[1]-5000
    }
  }
  if (is.null(behdat.df$System_Start_Time_2)){
    print("System start time 2 missing")
    tmp<-subset(behdat.df,behdat.df$blocksCompleted ==1)
    if(tmp$response_state[1] =="cue"){
    behdat.df$System_Start_Time_2[1]<- tmp$start_time[1]-5000 
    }else if (tmp$response_state[1] =="outcome_focus"){
      behdat.df$System_Start_Time_2[1]<- tmp$start_time[2]-5000 
    }#outcome_focus always seems to hang on
  }
  start_datetime <- as.POSIXct(behdat.df$System_Start_Time_1[1]/1000, origin = "1970-01-01")
  print(start_datetime)
  #calculate onset relative to system start
  behdat.df$onset <- ifelse(behdat.df$blocksCompleted ==0,(behdat.df$start_time - behdat.df$System_Start_Time_1[1])/1000,
                            (behdat.df$start_time - behdat.df$System_Start_Time_2[1])/1000)
 
  behdat.df$cue_type <- substring(behdat.df$cue,6,6) #extract cue name from image file
  behdat.df$probe_type <- substring(behdat.df$probe,6,6) #extract probe name from image file
  
  #make savePath folder if necessary
  if (!dir.exists(savePath)){
    dir.create(savePath)
    print("making new output folder")
  }   
  
  #### loop through the two blocks ####
  for (ii in c(0,1)) {
    #skip the splitting in half
    run<-ifelse(ii==0, "run1","run2")
    task <- ifelse(ii==0,"task-topx1_rec-NORDIC_run-1_part-mag_bold","task-topx2_rec-NORDIC_run-2_part-mag_bold")

    #make new task folder if necessary
    if (!dir.exists(paste(savePath,task,sep="/"))){
      dir.create(paste(savePath,task,sep="/"))
      print("making new task folder")
    }
    #make new subject folder if necessary
    if (!dir.exists(paste(savePath,task, paste0("sub-",sub),sep="/"))){
      dir.create(paste(savePath,task,paste0("sub-",sub),sep="/"))
      print("making new subject folder")
    }
    #make new session folder if necessary
    if (!dir.exists(paste(savePath,task, paste0("sub-",sub),paste0("ses-",ses),sep="/"))){
      dir.create(paste(savePath,task, paste0("sub-",sub),paste0("ses-",ses),sep="/"))
      print("making new session folder")
    }
    
    #### split into first and second halves ####
    half1.dat <- filter(behdat.df, blocksCompleted ==ii) #last outcome_focus might be in the wrong block, no biggie
    
    #remove the system start info and reshape to wider format
    half.dat<-half1.dat %>%
      filter(!is.na(start_time))%>%
      select(c("response_time","response","response_state","duration","onset","outcome",
               "blocksCompleted","trialsCompleted","cue_type","probe_type"))%>%
      pivot_wider(names_from = response_state, values_from = c(response_time,response,duration, onset, outcome))
    
    half.dat$tcfweight <- rep(c(1), times=nrow(half.dat))
    
    #### ACue and BCue 3cols ####
    
    ###Prep components #####
  
    c.c.final.1 <- filter(half.dat, response_cue ==37 | response_cue_focus == 37)#correct responses
    c.c.final.1$tcfdur.1 <- (c.c.final.1$duration_cue + c.c.final.1$duration_cue_focus)/1000
    ###A Cue ####
    
    tcf_acue.prep.1 <- filter(c.c.final.1, cue_type == "A")
    tcf_acue.prep.1 <- filter(tcf_acue.prep.1, (probe_type == "X" & (response_probe == 38| response_probe_focus == 38)) | 
                                (probe_type == "Y" & (response_probe == 39| response_probe_focus == 39))) #filter out probe errors
    tcf_acue.1 <- as.data.frame(cbind(tcf_acue.prep.1$onset_cue, tcf_acue.prep.1$tcfdur.1, tcf_acue.prep.1$tcfweight))
    
    ##write.table(tcf_acue.1, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",sub,"_",ses,"_Acue_",run,".txt"), row.names = FALSE, col.names = FALSE)
    
    ###B Cue ####
    
    tcf_bcue.prep.1 <- filter(c.c.final.1, cue_type == "B")
    tcf_bcue.prep.1 <- filter(tcf_bcue.prep.1, response_probe == 39| response_probe_focus == 39) #filter out probe errors
    tcf_bcue.1 <- as.data.frame(cbind(tcf_bcue.prep.1$onset_cue, tcf_bcue.prep.1$tcfdur.1, tcf_bcue.prep.1$tcfweight))
    
    ##write.table(tcf_bcue.1, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",sub,"_",ses,"_Bcue_",run,".txt"), row.names = FALSE, col.names = FALSE)
    
    #### Probe 3Cols ####
  
    ###Prep#####
    #get the timing based on reaction time, not duration
    half.dat$response_time_probe_adjusted <- (ifelse(is.na(half.dat$response_time_probe), 
                                                    half.dat$duration_probe + half.dat$response_time_probe_focus, 
                                                    half.dat$response_time_probe))/1000
    
    ax.1 <- filter(half.dat, cue_type == "A" & probe_type =="X")
    ay.1 <- filter(half.dat, cue_type == "A" & probe_type =="Y")
    bx.1 <- filter(half.dat, cue_type == "B" & probe_type =="X")
    by.1 <- filter(half.dat, cue_type == "B" & probe_type =="Y")
    
    #need correct trials only (ignore feedback)
    
    ###AX####
    
    ax.final.1 <- filter(ax.1, (response_probe == 38| response_probe_focus == 38) & (response_cue == 37 | response_cue_focus == 37))
    tcf_ax.1 <- as.data.frame(cbind(ax.final.1$onset_probe, ax.final.1$response_time_probe_adjusted, ax.final.1$tcfweight))
    
    ##write.table(tcf_ax.1, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",sub,"_",ses,"_AX_",run,".txt"), row.names = FALSE, col.names = FALSE)
    
    ###AY####
    
    ay.final.1 <- filter(ay.1, (response_probe == 39| response_probe_focus == 39) & (response_cue == 37 | response_cue_focus == 37))
    tcf_ay.1 <- as.data.frame(cbind(ay.final.1$onset_probe, ay.final.1$response_time_probe_adjusted, ay.final.1$tcfweight))

    ##write.table(tcf_ay.1, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",sub,"_",ses,"_AY_",run,".txt"), row.names = FALSE, col.names = FALSE)
    
    ###BX####
    
    bx.final.1 <- filter(bx.1, (response_probe == 39| response_probe_focus == 39) & (response_cue == 37 | response_cue_focus == 37))
    tcf_bx.1 <- as.data.frame(cbind(bx.final.1$onset_probe, bx.final.1$response_time_probe_adjusted, bx.final.1$tcfweight))

    ##write.table(tcf_bx.1, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",sub,"_",ses,"_BX_",run,".txt"), row.names = FALSE, col.names = FALSE)
    
    ###BY####
    
    by.final.1 <- filter(by.1, (response_probe == 39| response_probe_focus == 39) & (response_cue == 37 | response_cue_focus == 37))
    tcf_by.1 <- as.data.frame(cbind(by.final.1$onset_probe, by.final.1$response_time_probe_adjusted, by.final.1$tcfweight))
    
    ##write.table(tcf_by.1, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",sub,"_",ses,"_BY_",run,".txt"), row.names = FALSE, col.names = FALSE)
    
    #### Outcome 3cols ####
    
    ###Prep#####
    
    correct.1 <- filter(half.dat,  outcome_outcome == "Correct")
    incorrect.1 <-  filter(half.dat,  outcome_outcome == "Incorrect")
    missed.1 <-  filter(half.dat,  outcome_outcome == "Respond Faster")
    
    correct.ax.1 <- filter(correct.1, cue_type == "A" & probe_type =="X")
    correct.ay.1 <- filter(correct.1, cue_type == "A" & probe_type =="Y")
    correct.bx.1 <- filter(correct.1, cue_type == "B" & probe_type =="X")
    correct.by.1 <- filter(correct.1, cue_type == "B" & probe_type =="Y")
    
    ###Correct####
    
    tcf_correct.1 <- as.data.frame(cbind(correct.1$onset_outcome, correct.1$duration_outcome/1000, correct.1$tcfweight))
    tcf_correct.ax.1 <- as.data.frame(cbind(correct.ax.1$onset_outcome, correct.ax.1$duration_outcome/1000, correct.ax.1$tcfweight))
    tcf_correct.ay.1 <- as.data.frame(cbind(correct.ay.1$onset_outcome, correct.ay.1$duration_outcome/1000, correct.ay.1$tcfweight))
    tcf_correct.bx.1 <- as.data.frame(cbind(correct.bx.1$onset_outcome, correct.bx.1$duration_outcome/1000, correct.bx.1$tcfweight))
    tcf_correct.by.1 <- as.data.frame(cbind(correct.by.1$onset_outcome, correct.by.1$duration_outcome/1000, correct.by.1$tcfweight))
    
    #write.table(tcf_correct.1, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",sub,"_",ses,"_Outcome_win_",run,".txt"), row.names = FALSE, col.names = FALSE)
    #write.table(tcf_correct.ax.1, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",sub,"_",ses,"_Outcome_win_AX_",run,".txt"), row.names = FALSE, col.names = FALSE)
    #write.table(tcf_correct.ay.1, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",sub,"_",ses,"_Outcome_win_AY_",run,".txt"), row.names = FALSE, col.names = FALSE)
    #write.table(tcf_correct.bx.1, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",sub,"_",ses,"_Outcome_win_BX_",run,".txt"), row.names = FALSE, col.names = FALSE)
    #write.table(tcf_correct.by.1, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",sub,"_",ses,"_Outcome_win_BY_",run,".txt"), row.names = FALSE, col.names = FALSE)
    
    
    ###Incorrect####
    
    tcf_incorrect.1 <- as.data.frame(cbind(incorrect.1$onset_outcome, incorrect.1$duration_outcome/1000, incorrect.1$tcfweight))
    #create dummy variable if no errors were made
    if(nrow(tcf_incorrect.1) ==0) {tcf_incorrect.1 <-  data.frame(0,0,0)}
    
    #write.table(tcf_incorrect.1, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",sub,"_",ses,"_Outcome_loss_",run,".txt"), row.names = FALSE, col.names = FALSE)
    
    ###Missed####
  
    tcf_missed.1 <- as.data.frame(cbind(missed.1$onset_outcome, missed.1$duration_outcome/1000, missed.1$tcfweight))
    #create dummy variable if no trials were missed
    if(nrow(tcf_missed.1) ==0) {tcf_missed.1 <-  data.frame(0,0,0)}
    
    #write.table(tcf_missed.1, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",sub,"_",ses,"_Outcome_miss_",run,".txt"), row.names = FALSE, col.names = FALSE)
    
    
    #### Error 3cols ####
    
    #Missed trials will be modeled separately
    
    ###Cue errors####
    #Note: cue errors = responded either with '38' or '39'.
    
    #cue.err.final.1 <- filter(half.dat, response_cue == 38| response_cue == 39 | response_cue_focus == 38| response_cue_focus == 39)# for only cue errors
    
    #also remove cue trials where there were probe errors
    cue.err.final.1 <- filter(half.dat, (cue_type== "A" & probe_type=="X" & (response_probe == 37| response_probe == 39|response_probe_focus == 37| response_probe_focus == 39)) |
                                   (cue_type == "A" & probe_type=="Y" & (response_probe == 37| response_probe == 38|response_probe_focus == 37| response_probe_focus == 38))  |
                                   (cue_type == "B" & probe_type=="X" & (response_probe == 37| response_probe == 38|response_probe_focus == 37| response_probe_focus == 38))  |
                                   (cue_type == "B" & probe_type=="Y" & (response_probe == 37| response_probe == 38|response_probe_focus == 37| response_probe_focus == 38)) )
    
    cue.err.final.1$tcfdur.cue.err <- ((cue.err.final.1$duration_cue + cue.err.final.1$duration_cue_focus)/1000)
    
    tcf_cue.err.1 <- as.data.frame(cbind(cue.err.final.1$onset_cue, cue.err.final.1$tcfdur.cue.err, cue.err.final.1$tcfweight))
    tcf_cue.err.1 <- arrange(tcf_cue.err.1, V1)
    
    #create dummy variable if no errors were made
    if(nrow(tcf_cue.err.1) == 0) {tcf_cue.err.1 <- data.frame(0,0,0)}
    
    #write.table(tcf_cue.err.1, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",sub,"_",ses,"_Error_cue_",run,".txt"), row.names = FALSE, col.names = FALSE)
    
    
    ##Probe errors####
    #same errors as above, just during the probe time
    probe.err.final.1 <- filter(half.dat, (cue_type== "A" & probe_type=="X" & (response_probe == 37| response_probe == 39|response_probe_focus == 37| response_probe_focus == 39)) |
                                (cue_type == "A" & probe_type=="Y" & (response_probe == 37| response_probe == 38|response_probe_focus == 37| response_probe_focus == 38))  |
                                (cue_type == "B" & probe_type=="X" & (response_probe == 37| response_probe == 38|response_probe_focus == 37| response_probe_focus == 38))  |
                                (cue_type == "B" & probe_type=="Y" & (response_probe == 37| response_probe == 38|response_probe_focus == 37| response_probe_focus == 38)) )
    
    tcf_probe.err.1 <- as.data.frame(cbind(probe.err.final.1$onset_probe, probe.err.final.1$response_time_probe_adjusted, probe.err.final.1$tcfweight))
    tcf_probe.err.1 <- arrange(tcf_probe.err.1, V1)
    
    #create dummy variable if no errors were made
    if(nrow(tcf_probe.err.1) ==0) {tcf_probe.err.1 <-  data.frame(0,0,0)}
    
    #write.table(tcf_probe.err.1, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",sub,"_",ses,"_Error_probe_",run,".txt"), row.names = FALSE, col.names = FALSE)
    
    #combine to model error trials
    probe.err.final.1$trialLength <- cue.err.final.1$tcfdur.cue.err + probe.err.final.1$response_time_probe_adjusted
    tcf_all.err.1 <- as.data.frame(cbind(probe.err.final.1$onset_cue, probe.err.final.1$trialLength, probe.err.final.1$tcfweight))
    tcf_all.err.1 <- arrange(tcf_all.err.1, V1)
    ##Missed Trial 3cols
    
    ###Missed Cues####
    
    #removes trials with missed cue or probe
    cue.miss.1 <- filter(half.dat, is.na(response_cue) & is.na(response_cue_focus) | is.na(response_probe) & is.na(response_probe_focus))
    cue.miss.1 <- filter(cue.miss.1, !is.na(onset_cue)) #catches trials where the two halves were not split properly
    
    cue.miss.1$tcfdur.cue.miss <- ((cue.miss.1$duration_cue + cue.miss.1$duration_cue_focus)/1000)

    tcf_cue.miss.1 <- as.data.frame(cbind(cue.miss.1$onset_cue, cue.miss.1$tcfdur.cue.miss, cue.miss.1$tcfweight))
    tcf_cue.miss.1 <- arrange(tcf_cue.miss.1, V1)
    
    if(nrow(tcf_cue.miss.1) ==0) {tcf_cue.miss.1 <-  data.frame(0,0,0)}
    
    #write.table(tcf_cue.miss.1, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",sub,"_",ses,"_Missed_cue_",run,".txt"), row.names = FALSE, col.names = FALSE)
    
    ##Missed Probes####
    
    #removes trials with missed cue or probe
    probe.miss.1 <- filter(half.dat, is.na(response_probe) & is.na(response_probe_focus) | is.na(response_cue) & is.na(response_cue_focus))
    probe.miss.1 <- filter(probe.miss.1, !is.na(onset_probe))#catches trials where the two halves were not split properly
    
    probe.miss.1$tcfdur.probe.miss  <- ((probe.miss.1$duration_probe + probe.miss.1$duration_probe_focus)/1000)
  
    tcf_probe.miss.1 <- as.data.frame(cbind(probe.miss.1$onset_probe, probe.miss.1$tcfdur.probe.miss, probe.miss.1$tcfweight))
    tcf_probe.miss.1 <- arrange(tcf_probe.miss.1, V1)
    
    if(nrow(tcf_probe.miss.1) == 0) {tcf_probe.miss.1 <-  data.frame(0,0,0)}
    
    #write.table(tcf_probe.miss.1, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",sub,"_",ses,"_Missed_probe_",run,".txt"), row.names = FALSE, col.names = FALSE)
    
    #combine to model missed trials
    probe.miss.1$trialLength <- cue.miss.1$tcfdur.cue.miss + probe.miss.1$tcfdur.probe.miss
    tcf_all.miss.1 <- as.data.frame(cbind(probe.miss.1$onset_cue, probe.miss.1$trialLength, probe.miss.1$tcfweight))
    tcf_all.miss.1 <- arrange(tcf_all.miss.1, V1)
    
    tcf_all.issue.1 <- rbind(tcf_all.err.1,tcf_all.miss.1) 
    tcf_all.issue.1 <- arrange(tcf_all.issue.1,V1) 
    write.table(tcf_all.issue.1, file=paste0(savePath,"/",task,"/sub-",sub, "/ses-", ses,"/",sub,"_",ses,"_issue_trial_",run,".txt"), row.names = FALSE, col.names = FALSE)
    
    
  } #end of run loop
} #end of function

####  Using the function - Actual ######

# set working directory
setwd("/Users/arend103/Library/CloudStorage/Box-Box/STEP-P50/Reporting & Data Submission/MRI Task Data/TOPX Scanner Behavioral Data + Analyses")

# read in all the file names of TOPX raw data
txt_files <- list.files("/Users/arend103/Library/CloudStorage/Box-Box/STEP-P50/Reporting & Data Submission/MRI Task Data/TOPX Scanner Behavioral Data + Analyses/TOPX behavioral RAW",
                        pattern = "\\.txt$",
                        full.names = TRUE)

# filter to only keep 6-Month files
txt_files_6M <- grep("6M", txt_files, value = TRUE)

# save list of only subject IDs, for subject list in fMRI analyses 
# (e.g., subs_6m_YYYYMMDD_basic.txt)
fn <- basename(txt_files_6M)
(subject_ids <- sub("^(SP[0-9]{4}).*$", "\\1", fn))
subject_ids <- paste0("sub-", subject_ids)
writeLines(subject_ids,
           "/Users/arend103/Documents/umn_work/analyses/topx_analyses/topx_3col/topx_subject_ids_6m.txt")

# run iteratively
for (filePath in txt_files_6M) {
  message("Processing: ", basename(filePath))
  makeTOPX3Columns(filePath, savePath = "EVs")
}

