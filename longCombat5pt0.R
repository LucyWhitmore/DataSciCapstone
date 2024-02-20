## Multi-Site effect harmonization using Longitudinal ComBat
# Beer JC, Tustison NJ, Cook PA, Davatzikos C, Sheline YI, Shinohara RT, Linn KA. (2020)
# Longitudinal ComBat: A method for harmonizing longitudinal multi-scanner imaging data. 
# NeuroImage. https://doi.org/10.1016/j.neuroimage.2020.117129.

#install.packages("devtools")
#devtools::install_github("jcbeer/longCombat")

library(dplyr)
library(longCombat)
# possible that longCombat needs to be run before loading tidyverse (unclear if that included dplyr)

?longCombat


# load data
# need to figure out which id will help match these across waves - I'll filter by that
# might need abcd_imgincl01_id
# alternately, could do the exclusion separately for waves and rejoin
# i.e. could rbind smri_baseline and smri_followup

datapath = "/Volumes/devbrainlab/ABCD_Data/ABCD5pt1/core/imaging/"
#directory for saved data is "/Volumes/devbrainlab/

# could filter the inclusion, then join by id & wave (drop other columns from inclusion)
struc_include5<-rio::import(paste0(datapath, "mri_y_qc_incl.csv")) %>%
  select(src_subject_id, eventname, imgincl_t1w_include) %>% 
  #mutate(imgincl_t1w_include = as.numeric(imgincl_t1w_include)) %>% 
  filter(imgincl_t1w_include == 1) 

#qc_info5 <- rio::import("/Volumes/devbrainlab/ABCD_Data/ABCD5pt0/core/imaging/mri_y_qc_man_post_t2w.csv") %>%
#  select(src_subject_id, eventname, imgincl_t1w_include)

#Load structural data
# Volume
smri_vol <- rio::import(paste0(datapath,"mri_y_smr_vol_dsk.csv")) # %>%
 # select(subjectkey, src_subject_id, interview_age, sex, eventname, matches("vol")) 
# note: might just select all, don't need to match for volume anymore  

# Area
smri_area <- rio::import(paste0(datapath,"mri_y_smr_area_dsk.csv"))# %>%
  #select(subjectkey, src_subject_id, interview_age, sex, eventname, matches("area")) 
  
# Thickness
smri_thickness <- rio::import(paste0(datapath,"mri_y_smr_thk_dsk.csv")) #%>%
  #select(subjectkey, src_subject_id, interview_age, sex, eventname, matches("thick"))
 
# Subcortical Volume 
smri_subvol <- rio::import(paste0(datapath,"mri_y_smr_vol_aseg.csv")) %>%
  select(-c(smri_vol_scs_lesionlh, smri_vol_scs_lesionrh))
 # select(subjectkey, src_subject_id, interview_age, sex, eventname, matches("vol"))

# Scanner Info
mri_info5 <- rio::import(paste0(datapath,"mri_y_adm_info.csv")) %>% 
  select(src_subject_id, eventname, mri_info_deviceserialnumber) #, mri_info_manufacturer, mri_info_manufacturersmn)

# Demographic Info (age)
demo <- rio::import("/Volumes/devbrainlab/ABCD_Data/ABCD5pt0/core/abcd-general/abcd_y_lt.csv") %>% 
  select(src_subject_id, eventname, interview_age) %>% 
  mutate(interview_age = interview_age/12)



# don't need to join all measures together because we run combat on them separately
#join filter for usable data

# Cortical Volume
# join cortical and subcortical smri
good_mri_vol <- left_join(struc_include5, smri_vol, by = c("src_subject_id", "eventname"))
good_mri_vol <- left_join(good_mri_vol, mri_info5, by = c("src_subject_id", "eventname")) %>% 
  mutate(src_subject_id = as.factor(src_subject_id),
         eventname = as.factor(eventname), 
         mri_info_deviceserialnumber = as.factor(mri_info_deviceserialnumber)) %>% 
  select(-imgincl_t1w_include)

volume <- left_join(good_mri_vol, demo, c("src_subject_id", "eventname"))
#need to decide whether to bind volume and subcortical volume before combat

# Subcortical Volume
good_mri_subvol <- left_join(struc_include5, smri_subvol, by = c("src_subject_id", "eventname"))
good_mri_subvol <- left_join(good_mri_subvol, mri_info5, by = c("src_subject_id", "eventname")) %>% 
  mutate(src_subject_id = as.factor(src_subject_id),
         eventname = as.factor(eventname), 
         mri_info_deviceserialnumber = as.factor(mri_info_deviceserialnumber)) %>% 
  select(-imgincl_t1w_include)

subvolume <- left_join(good_mri_subvol, demo, c("src_subject_id", "eventname")) %>% 
  select(-c(smri_vol_scs_wmhintlh, smri_vol_scs_wmhintrh))


# Cortical Area
good_mri_area <- left_join(struc_include5, smri_area, by = c("src_subject_id", "eventname"))
good_mri_area <- left_join(good_mri_area, mri_info5, by = c("src_subject_id", "eventname"))%>% 
  mutate(src_subject_id = as.factor(src_subject_id),
         eventname = as.factor(eventname), 
         mri_info_deviceserialnumber = as.factor(mri_info_deviceserialnumber)) %>% 
  select(-imgincl_t1w_include)

area <- left_join(good_mri_area, demo, c("src_subject_id", "eventname"))

# note: combat doesn't work well when there's large sample differences between sites


random ='(1|src_subject_id)'

#################################
# Model: area
#################################

#################################
# longCombat() -- apply longitudinal ComBat
#################################
brainname_area <- c(names(area)[3:73]) # names(area)[5:72]


area_combat <- longCombat(idvar='src_subject_id',
                        timevar='eventname',  # could also be smri_visitid
                        batchvar='mri_info_deviceserialnumber', # should be a factor
                        features=brainname_area,
                        formula='interview_age', # change formula
                        ranef=random,   # for random intercept or slope. Dani did a random intercept, (1|ID)
                        data=area)


# get the harmonized data
area_harmonized <- left_join(area_combat$data_combat, demo, by = c("src_subject_id", "eventname"))

save(area_harmonized, file="area_harmonized.Rda")

# save combat feature names
#featurenames.combat <- names(area_harmonized)[4:71] # change indices

# merge with original dataframe
#area_final <- merge(area, area_harmonized[,c(1,2,4:41)], by=c('ID', 'TP'))  # change indices



# check area boxplot before combat
area_df <- as.data.frame(smri_area)
batchBoxplot(idvar='src_subject_id', 
             batchvar='mri_info_deviceserialnumber', 
             feature='brainname_area', 
             formula='interview_age',
             ranef='(1|src_subject_id)',
             data=unlist(area_df),
             colors=1:31,
             title='Area before ComBat')

# check area boxplot after combat
batchBoxplot(idvar='src_subject_id', 
             batchvar='mri_info_deviceserialnumber', 
             feature='featurenames.combat', 
             formula='interview_age',
             ranef='(1|ID)',
             data=area_harmonized,
             colors=1:31,
             title='Area after combat')


#################################
# Model: volume
#################################
brainname_volume <- c(names(volume)[3:73])

volume_combat <- longCombat(idvar='src_subject_id',
                          timevar='eventname',  # could also be smri_visitid
                          batchvar='mri_info_deviceserialnumber', # should be a factor
                          features=brainname_volume,
                          formula='interview_age', # change formula
                          ranef=random,   # for random intercept or slope. Dani did a random intercept, (1|ID)
                          data=volume)


# get the harmonized data
volume_harmonized <-left_join(volume_combat$data_combat, demo, by = c("src_subject_id", "eventname"))
# save data
save(volume_harmonized, file="volume_harmonized.Rda")


# save combat feature names
#featurenames.combat <- names(volume_harmonized)[4:74] # change indices

# merge with original dataframe
#volume_final <- merge(volume, volume_harmonized[,c(1,2,4:41)], by=c('ID', 'TP'))  # change indices

#################################
# Model: subcortical volume
#################################
brainname_subvolume <- c(names(subvolume)[3:44])

subvolume_combat <- longCombat(idvar='src_subject_id',
                            timevar='eventname',  # could also be smri_visitid
                            batchvar='mri_info_deviceserialnumber', # should be a factor
                            features=brainname_subvolume,
                            formula='interview_age', # change formula
                            ranef=random,   # for random intercept or slope. Dani did a random intercept, (1|ID)
                            data=subvolume)


# get the harmonized data
subvolume_harmonized <-left_join(subvolume_combat$data_combat, demo, by = c("src_subject_id", "eventname"))
# save data
save(subvolume_harmonized, file="subvolume_harmonized.Rda")


#################################
# Model: thickness
#################################
brainname_thick <- c()

thick_combat <- longCombat(idvar='src_subject_id',
                          timevar='eventname',  # could also be smri_visitid
                          batchvar='batch', # should be a factor
                          features=brainname_thick,
                          formula='age + sex + TP + PDS', # change formula
                          ranef=random,   # for random intercept or slope. Dani did a random intercept, (1|ID)
                          data=thickness)


# get the harmonized data
thick_harmonized <- thick_combat$data_combat
# save combat feature names
featurenames.combat <- names(thick_harmonized)[5:45] # change indices

# merge with original dataframe
thick_final <- merge(thick, thick_harmonized[,c(1,2,4:41)], by=c('ID', 'TP'))  # change indices


########### Clean Up #########

# merge back into one dataframe that can be passed into model training
harmonized_mri <- left_join(area_harmonized, volume_harmonized, by = c('src_subject_id', 'eventname','mri_info_deviceserialnumber'))
harmonized_mri <- left_join(harmonized_mri, subvolume_harmonized, by = c('src_subject_id', 'eventname','mri_info_deviceserialnumber')) %>% 
  relocate(interview_age, .before=mri_info_deviceserialnumber)

# remove .combat from end of each feature
names(harmonized_mri) = gsub(pattern = ".combat", replacement = "", x = names(harmonized_mri))

save(harmonized_mri, file="harmonized_mri.Rda")


### DECIDE ABOUT SAVING BY WAVE ###

# split by wave, maybe remove scanner info and subject ids?
harmonized_mri_baseline <- harmonized_mri %>% 
  filter(eventname == "baseline_year_1_arm_1")

harmonized_mri_followup_2 <- harmonized_mri %>% 
  filter(eventname == "2_year_follow_up_y_arm_1")

harmonized_mri_followup_4 <- harmonized_mri %>% 
  filter(eventname == "4_year_follow_up_y_arm_1")

# save harmonized data
save(harmonized_mri_baseline, file="harmonized_mri_baseline.Rda")

save(harmonized_mri_followup_2, file="harmonized_mri_followup_2.Rda")

save(harmonized_mri_followup_4, file="harmonized_mri_followup_4.Rda")



# select just features for model
# split by wave
model1_features_baseline <- renamed_harmonized_mri %>% 
  filter(eventname == "baseline_year_1_arm_1")  %>% 
  select(-c(src_subject_id, eventname, mri_info_deviceserialnumber))

model1_features_followup <- renamed_harmonized_mri %>% 
  filter(eventname == "2_year_follow_up_y_arm_1")  %>% 
  select(-c(src_subject_id, eventname, mri_info_deviceserialnumber))


save(model1_features_baseline, file="model1_features_baseline.Rda")
save(model1_features_followup, file="model1_features_followup.Rda")



# at some point, split into training/testing (for model training)


# decide if we want names to be the same - probably
# load old training set - training_sample_baseline
# select subject IDs
# filter for only those subjects
# drop NA columns
load("/Volumes/devbrainlab/Lucy/BrainAGE/FYP/Sample1Baseline.Rda")

model2_features <- renamed_harmonized_mri %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% 
  filter(src_subject_id %in% sample_1_baseline$src_subject_id) %>% 
  select(-c(FS_BrainSeg_Vol ,FS_BrainSeg_Vol_No_Vent,FS_BrainSeg_Vol_No_Vent_Surf,
          FS_SupraTentorial_Vol_No_Vent,FS_SupraTentorial_No_Vent_Voxel_Count,FS_Mask_Vol,
          FS_BrainSegVol_eTIV_Ratio,FS_MaskVol_eTIV_Ratio,FS_L_Vessel_Vol, FS_L_ChoroidPlexus_Vol,
          FS_R_Vessel_Vol,FS_R_ChoroidPlexus_Vol,FS_OpticChiasm_Vol,FS_Total_GM_Vol, eventname, mri_info_deviceserialnumber))

save(model2_features, file="model2_features.Rda")
# resave


load("/Volumes/devbrainlab/Lucy/BrainAGE/FYP/Sample2Baseline.Rda")

model2_analysis <- renamed_harmonized_mri %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% 
  filter(src_subject_id %in% sample_2_baseline$src_subject_id) %>% 
  select(-c(FS_BrainSeg_Vol ,FS_BrainSeg_Vol_No_Vent,FS_BrainSeg_Vol_No_Vent_Surf,
            FS_SupraTentorial_Vol_No_Vent,FS_SupraTentorial_No_Vent_Voxel_Count,FS_Mask_Vol,
            FS_BrainSegVol_eTIV_Ratio,FS_MaskVol_eTIV_Ratio,FS_L_Vessel_Vol, FS_L_ChoroidPlexus_Vol,
            FS_R_Vessel_Vol,FS_R_ChoroidPlexus_Vol,FS_OpticChiasm_Vol,FS_Total_GM_Vol, eventname, mri_info_deviceserialnumber))

save(model2_analysis, file="model2_analysis.Rda")



# Model 2 - Follow-Up
model2_features_followup <- renamed_harmonized_mri %>% 
  filter(eventname == "2_year_follow_up_y_arm_1") %>% 
  filter(!(src_subject_id %in% sample_2_baseline$src_subject_id)) %>% 
  select(-c(FS_BrainSeg_Vol,FS_BrainSeg_Vol_No_Vent,FS_BrainSeg_Vol_No_Vent_Surf,
            FS_SupraTentorial_Vol_No_Vent,FS_SupraTentorial_No_Vent_Voxel_Count,FS_Mask_Vol,
            FS_BrainSegVol_eTIV_Ratio,FS_MaskVol_eTIV_Ratio,FS_L_Vessel_Vol, FS_L_ChoroidPlexus_Vol,
            FS_R_Vessel_Vol,FS_R_ChoroidPlexus_Vol,FS_OpticChiasm_Vol,FS_Total_GM_Vol, eventname, mri_info_deviceserialnumber))

save(model2_features_followup, file="model2_features_followup.Rda")


model2_analysis_followup <- renamed_harmonized_mri %>% 
  filter(eventname == "2_year_follow_up_y_arm_1") %>% 
  filter(src_subject_id %in% sample_2_baseline$src_subject_id) %>% 
  select(-c(FS_BrainSeg_Vol,FS_BrainSeg_Vol_No_Vent,FS_BrainSeg_Vol_No_Vent_Surf,
            FS_SupraTentorial_Vol_No_Vent,FS_SupraTentorial_No_Vent_Voxel_Count,FS_Mask_Vol,
            FS_BrainSegVol_eTIV_Ratio,FS_MaskVol_eTIV_Ratio,FS_L_Vessel_Vol, FS_L_ChoroidPlexus_Vol,
            FS_R_Vessel_Vol,FS_R_ChoroidPlexus_Vol,FS_OpticChiasm_Vol,FS_Total_GM_Vol, eventname, mri_info_deviceserialnumber))

save(model2_analysis_followup, file="model2_analysis_followup.Rda")
 