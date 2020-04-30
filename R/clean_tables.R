## This file provides a function to read/clean each section of the CRF tables.  If it has only one table, like demographics, then it reads just that one.  If it has multiple tables, like the treatment tables or endpoint tables then it works to merge them such that there is only a single entry per patient. This isn't always what you may want, and the logic is more likely to run into issues with changes on the Labkey side.

# Because there are new tables that split out columns that were in the old tables they need to be matched and merged.
match_cols<-function(mainTB, newTB) {
  mainTB = mainTB %>% ungroup %>% rename_all(list(~sub('^.+\\.', '',.))) 
  mainStudyCols = grep('Study', colnames(mainTB))

  newTB = newTB %>% ungroup %>% rename_all(list(~sub('^.+\\.', '',.))) 
  newStudyCols = grep('Study', colnames(newTB))
  
  x = colnames(mainTB %>% dplyr::select(-contains('Study')))
  y = colnames(newTB %>% dplyr::select(-contains('Study')))
  
  merged <- mainTB %>% dplyr::select(mainStudyCols,intersect(x,y)) %>% 
    dplyr::bind_rows( newTB %>% dplyr::select(newStudyCols, y) ) %>% 
    dplyr::mutate_at(vars(-matches('Study')), list(~str_to_upper(.)))
  
  dups = merged %>% dplyr::group_by(StudySubjectID) %>% filter(n()>1)
  
  uniqe_entry<-function(z) {
    if (length(na.omit(z)) <= 0) return(unique(z))
    z = na.omit(z)    
    while (length(z) > 1) {
      z = unique(z)
      
      if (class(z) %in% c('Date')) {
        z = min(z)
      } else if (is.numeric(z)) {
        z =  max(z)
      } else {
        z = paste(z,collapse=', ')
      }
    }
    return(z)
  }
  
  dups = dups %>% dplyr::group_by(StudySubjectID, StudySite) %>% dplyr::summarise_all(list(~uniqe_entry(.)))
  
  merged = merged %>% dplyr::group_by(StudySubjectID) %>% filter(n()==1) %>% bind_rows(dups) %>%
    dplyr::mutate_at(vars(-matches('Study|Protocol|Modality|Response')), list(~str_to_sentence(.))) %>%
    dplyr::mutate_at(vars(matches('Date')), list(~as.Date(format.Date(., '%Y-%m-%d'))))
  
  infoCols = colnames(merged %>% ungroup %>% dplyr::select(-matches('Study')))
  
  merged[which(apply(is.na(merged), 1, sum) < length(infoCols)),]
}

merge.patient.tables<-function(df1, df2, ...) {
  newDF <- merge.data.frame(df1, df2, by='row.names', ...)
  rownames(newDF) <- newDF[['Row.names']]
  newDF[['Row.names']] <- NULL
  return(newDF)
}

#' Read data from demographics CRF table.  
#' Strips out columns that include "Panther, CRF, age in months, and the "OpenClinica" informed consent columns.
#' @name read.demographics
#' @param ...
#'
#' @author
read.demographics<-function(ocs, table, occams_ids=NULL, rulesFile=NULL, ...) {
  table = get.tables.by.crf(ocs, 'demographics')
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  filter = make.id.filter(occams_ids,'DI_StudySubjectID')

  if (inherits(ocs, "OCCAMSLabkey")) {
    if (verbose) message(paste("Reading",table))
    di <- get.table.data(ocs=ocs, table=table, rulesFile=rulesFile, colFilter=filter)
  } else {
    di <- read.table.data(table,  rulesFile=rulesFile)
  }

  # Prefer to get DoD from FE tables, but Babs mentioned that people are supposed to go back and update the DI table with it as well so inconsistencies could be checked
  di <- di %>% dplyr::select(-DI.StudyNumber, -contains('Panther'), -contains('OpenClinica'),-DI.InformedConsentOCCAMSAndICGC, -DI.DateOfInformedConsentForCAMSAndICGC, -contains('CRF'), -DI.PatientEthnicityComment, -DI.ageAtDiagnosis.months)

  if (verbose) message(paste(nrow(di), "patients"))
  return(di)
}


#' Read data from expoosures CRF tables.  
#' This includes: family history for OGC, family history of other cancers, and the general exposures tables
#' Recalculate BMI, create EX.CurrentBMI.c and EX.5YearsAgoBMI.c columns
#' @name read.exposures
#' @param ...
#'
#' @author
read.exposures<-function(ocs, tables, occams_ids=NULL, rulesFiles=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  tables = get.tables.by.crf(ocs, 'exposures')

  if (length(tables) != 3)
    stop("Three tables expected for exposures")

  filter = make.id.filter(occams_ids,'EX_StudySubjectID')
  filter1 = make.id.filter(occams_ids,'EX1_StudySubjectID')
  filter2 = make.id.filter(occams_ids,'EX2_StudySubjectID')

  if (inherits(ocs, "OCCAMSLabkey")) {
    #stop("Create connection to OCCAMS Labkey first")
    if (verbose) message(paste("Reading",paste(tables, collapse=', ')))
    ex <- get.table.data(ocs=ocs, table=grep('^ex_',tables,value=T),  rulesFile=rulesFiles[1], colFilter=filter)
    ex_history_gastric <- get.table.data(ocs=ocs, table=grep('ex.*gastric',tables,value=T), rulesFile=rulesFiles[2], colFilter=filter1)
    ex_history_other <- get.table.data(ocs=ocs, table=grep('ex.*other_cancer',tables,value=T), rulesFile=rulesFiles[3], colFilter=filter2)

  } else {
    ex <- read.table.data(tables[1],  rulesFile=rulesFiles[1])
    ex_history_gastric <- read.table.data(tables[2], rulesFile=rulesFiles[2])
    ex_history_other <- read.table.data(tables[3], rulesFile=rulesFiles[3])
  }

  # -1 is often used for NA in the numeric columns
    ex <- ex %>% 
      dplyr::mutate_at(vars(grep('Height|Weight|BMI|Weeks|Months|Years|Days', colnames(.))), funs(as.numeric)) %>%
      dplyr::mutate_if( is.numeric, funs(ifelse(.<0,NA,.))) %>% 
      dplyr::group_by(EX.StudySubjectID) %>% 
      dplyr::mutate( #EX.CurrentWeightKG = ifelse(EX.CurrentWeightKG>140 | EX.CurrentWeightKG<22, NA, EX.CurrentWeightKG),
                     #EX.CurrentHeightCM = ifelse(EX.CurrentHeightCM>250 | EX.CurrentHeightCM<50, NA, EX.CurrentHeightCM),
                     EX.CurrentBMI.c = bmicalc(EX.CurrentWeightKG,EX.CurrentHeightCM),
                     EX.5YearsAgoBMI.c = bmicalc(EX.5YearsAgoWeightKG,EX.5YearsAgoHeightCM))

  # These two are always NA, seem to have moved to the other tables
  ex <- dplyr::select(ex, -EX.FamilyHistoryOfOtherCancerRelationship, -EX.FamilyHistoryOfOesophagoGastricCancerRelationship)
  
  if (verbose) message(paste(nrow(ex), "patients in",grep('^ex_',tables,value=T)))
  
  # Join family history tables by subject ID
  ex_history <- full_join(ex_history_gastric, ex_history_other, by=c("EX1.StudySubjectID"="EX2.StudySubjectID")) %>%
    dplyr::select(-matches('StudySite')) %>%
    dplyr::rename_all(list(~sub('^EX\\d','EX',.))) %>% ungroup %>%
    dplyr::mutate_if(is.character, list(~str_to_title(.))) %>%
    dplyr::mutate_at(vars(EX.FamilyHistoryOfOesophagoGastricCancerLocation, EX.FamilyHistoryOfOesophagoGastricCancerRelationship,EX.FamilyHistoryOfOtherCancerType,EX.FamilyHistoryOfOtherCancerRelationship), list(~factor(.)))

  if (verbose) message(paste(nrow(ex_history), "patients in", paste(grep('^ex.*[gastric|other]',tables,value=T), collapse=', ') ))

  return(list('ex'=ex, 'history'=ex_history))
}

#' Read data from endpoint CRF table.  
#' Strips out columns that include "OpenClinica(Start|Interview)Date
#' @name read.endpoint
#' @param ...
#'
#' @author
read.endpoint<-function(ocs, tables, occams_ids=NULL, rulesFile=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)
  
  table = get.tables.by.crf(ocs, c('endpoint'))
  filter = make.id.filter(occams_ids, 'EP_StudySubjectID')
  
  if (!is.null(rulesFile))
    rulesFile = grep('ep_', rulesFile, value=T)
  
  if (inherits(ocs, "OCCAMSLabkey")) {
    #stop("Create connection to OCCAMS Labkey first")
    if (verbose) message(paste("Reading",paste(table, collapse=",")))
    ep <- get.table.data(ocs,table, rulesFile=rulesFile, colFilter=filter)
  } else {
    if (verbose) message(paste("Reading",paste(basename(tables), collapse=", ")))
    ep <- read.table.data(table,  rulesFile=rulesFile)
  }
  
  ep <- dplyr::select(ep, -contains('OpenClinica'))
  
  if (nrow(ep) > 0) {
    ep <- ep %>% ungroup %>% dplyr::mutate_at(vars(EP.EndPoint,EP.ReasonForPatientDeath), list(~factor(.))) %>%
      dplyr::group_by(EP.StudySubjectID) 
  }

  if (verbose) message(paste0(nrow(ep), ' patients in ', table))

  return(ep)
}

#' Read data from followup, endpoint, and recurrence CRF tables.  
#' This requires a lot of cleaning around to get a single row per patient with the most recent follow up visit date, recurrence, and death. 
#' I have seen a few small inconsistencies in the dates, but they are days rather than years so for now this will work.  
#' TODO: add some date sanity checks
#' @name read.followup
#' @param ...
#'
#' @author
read.followup<-function(ocs, tables, occams_ids=NULL, rulesFiles=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  tables = get.tables.by.crf(ocs, c('endpoint','follow-up', 'recurrence'))

  if (length(tables) < 6)
    stop("Six tables expected for endpoint, follow-up, and recurrence.")

  if (inherits(ocs, "OCCAMSLabkey")) {
    #stop("Create connection to OCCAMS Labkey first")
    if (verbose) message(paste("Reading tables:\n\t",paste(tables, collapse="\n\t ")))

    # Not sure which direction this is going, but FE or EP must be superseded by the other at some point, except new values are continuing to appear in both so they'll need some merging
    ep <- read.endpoint(ocs,grep('ep',tables,value=T),occams_ids, rulesFile=rulesFiles) 
    if (!is.null(rulesFiles)) {
      fe <- get.table.data(ocs,grep('fe_',tables,value=T),  rulesFile=grep('fe_', rulesFiles,value=T), colFilter=make.id.filter(occams_ids, 'FE_StudySubjectID'))
      fe1 <- get.table.data(ocs,grep('fe1_',tables,value=T),  rulesFile=grep('fe1_', rulesFiles,value=T), colFilter=make.id.filter(occams_ids, 'FE1_StudySubjectID'))
    } else {
      fe <- get.table.data(ocs,grep('fe_',tables,value=T),  rulesFile=NULL, colFilter=make.id.filter(occams_ids, 'FE_StudySubjectID'))
      fe1 <- get.table.data(ocs,grep('fe1_',tables,value=T),  rulesFile=NULL, colFilter=make.id.filter(occams_ids, 'FE1_StudySubjectID'))
    }
    # followup and reoccurance - these are new tables
    #filterR = make.id.filter(occams_ids, 'FR_StudySubjectID')
    # This table currently only includes OpenClinica(Start|Interview)Date so there's no point in using it
    #fr <- get.table.data(ocs,grep('fr_',tables,value=T), rulesFile=NULL, colFilter=filterR)
    fr0 <- get.table.data(ocs,grep('fr0_',tables,value=T), rulesFile=NULL, colFilter=make.id.filter(occams_ids, 'FR0_StudySubjectID'))
    fr1 <- get.table.data(ocs,grep('fr1_',tables,value=T), rulesFile=NULL, colFilter=make.id.filter(occams_ids, 'FR1_StudySubjectID'))
  } else {
    if (verbose) message(paste("Reading",paste(basename(tables), collapse=", ")))
    ep <- read.endpoint(ocs, grep('ep',tables,value=T), occams_ids, rulesFile=grep('ep_', rulesFiles, value=T)) 
    fe <- read.table.data(ocs, grep('fe_',tables,value=T),  rulesFile=grep('fe_', rulesFiles,value=T))
    fr0 <- read.table.data(ocs,grep('fr0_',tables,value=T), rulesFile=NULL)
    fr1 <- read.table.data(ocs,grep('fr1_',tables,value=T), rulesFile=NULL)
  }
  
  if (nrow(ep) <= 0) {
    for(n in colnames(ep)) fe <- add_column(fe, !!n := NA_character_)
    feep <- fe %>% dplyr::select(-EP.StudySubjectID)
  } else if (nrow(fe) <= 0) {
    for(n in colnames(fe)) ep <- add_column(ep, !!n := NA_character_)
    feep <- ep %>% dplyr::rename(EP.StudySubjectID = FE.StudySubjectID) 
  } else {
    feep <- full_join(fe,ep,by=c('FE.StudySubjectID'='EP.StudySubjectID'))  
  }

  feep <- feep %>% dplyr::group_by(FE.StudySubjectID) %>%
    dplyr::mutate(
      FE.StudySite = na.omit(FE.StudySite, EP.StudySite)[1],
      FE.EndPoint = dplyr::case_when(
        (FE.EndPoint == 'Patient died' | EP.EndPoint == 'Patient died') ~ 'Patient died', 
        (FE.EndPoint == 'Follow-up ended' | EP.EndPoint == 'Follow-up ended') ~ 'Follow-up ended', 
        # this is because there was not a good option for "5-years follow up ended.  It did not mean the patient was withdrawing from the study
        (FE.EndPoint == 'Patient withdrawal' | EP.EndPoint == 'Patient withdrawal') ~ 'Follow-up ended', 
        TRUE ~ as.character(na.omit(FE.EndPoint, EP.EndPoint))[1] ),
      FE.DateOfPatientDeath = na.omit(FE.DateOfPatientDeath,EP.DateOfPatientDeath)[1],
      FE.ReasonForPatientDeath = na.omit(FE.ReasonForPatientDeath, EP.ReasonForPatientDeath)[1],
      FE.ReasonForPatientDeathOther = na.omit(FE.ReasonForPatientDeathOther, EP.ReasonForPatientDeathOther)[1] ) %>% 
    dplyr::select(starts_with('FE'), -contains('OpenClinica'), -matches('FE.ageAtDeath'))
  # get date of recurrence in a single row
  dis.r<-function(ID, hasR, dateR, dateU, i) {
    n = which(hasR == 'Yes')[1]
    if (length(n) <= 0 | is.na(n)) {
      return( c('No',NA)[i] )
    } else if (i == 1) {
      return(hasR[n])
    } else {
      date = dateR[n]
      if (is.na(date)) {
        date = dateU[n]
      }
      return(date)
    }
    return(NA_character_)            
  }
  
  fe1 <- fe1 %>% dplyr::group_by(FE1.StudySubjectID, FE1.StudySite) %>%
    dplyr::arrange(desc(FE1.DateOfUpdate)) %>% 
    dplyr::summarise(FE.LatestDateOfUpdate.c = FE1.DateOfUpdate[1],
                     FE.LatestReasonForFollowUp.c = FE1.ReasonForFollowUp[1],
                     FE.HasOriginalDiseaseReoccurred = dis.r(FE1.StudySubjectID,FE1.HasOriginalDiseaseReoccurred, FE1.DateOriginalDiseaseReoccurred, FE1.DateOfUpdate,1),
                     FE.DateOriginalDiseaseReoccurred = as.Date(dis.r(FE1.StudySubjectID,FE1.HasOriginalDiseaseReoccurred, FE1.DateOriginalDiseaseReoccurred, FE1.DateOfUpdate,2) )
                     ) 
  
  feep <- left_join(feep, fe1, by=c('FE.StudySubjectID'='FE1.StudySubjectID'))
  
  if (nrow(fr0) <= 0 & nrow(fr1) <= 0) {
    frr <- as_tibble(merge(fr0,fr1)) %>% dplyr::select(-FR1.StudySubjectID, -FR1.StudySite)
  } else if (nrow(fr0) <= 0) {
    for(n in colnames(fr0)) fr1 <- add_column(fr1, !!n := NA_character_)
    frr <- fr1 %>% dplyr::rename(FR0.StudySubjectID = FR1.StudySubjectID, FR0.StudySite = FR1.StudySite) 
  } else if (nrow(fr1) <= 0) {
    for(n in colnames(fr1)) fr0 <- add_column(fr0, !!n := NA_character_)
    frr <- fr0 %>% dplyr::select(-FR1.StudySubjectID, -FR1.StudySite) 
  } else {
    frr <- left_join(fr0,fr1, by=c('FR0.StudySubjectID' = 'FR1.StudySubjectID', 'FR0.StudySite' = 'FR1.StudySite'))
  }
  
  frr <- frr %>% dplyr::rename_at(vars(matches('^FE')), list(~sub('^FE', 'FR', .))) %>%
    dplyr::rename_at(vars(matches('Study(Subject|Site)')), list(~sub('FR[0-1]', 'FR',.))) %>%
    dplyr::group_by(FR.StudySubjectID, FR.StudySite) %>%
    dplyr::arrange(desc(FR.DateOfFollowUpEvent)) %>% 
    #filter(!grepl('death', FE1.ReasonForFollowUp)) %>%
    dplyr::summarise(FR.LatestDateOfUpdate.c = FR.DateOfFollowUpEvent[1],
                     FR.LatestReasonForFollowUp.c = FR.ReasonForFollowUp[1],
                     FR.HasOriginalDiseaseReoccurred = dis.r(FR.StudySubjectID,FR.HasOriginalDiseaseRecurred, FR.DateOriginalDiseaseRecurred, FR.DateOfFollowUpEvent,1),
                     FR.DateOriginalDiseaseReoccurred = as.Date(dis.r(FR.StudySubjectID,FR.HasOriginalDiseaseRecurred, FR.DateOriginalDiseaseRecurred, FR.DateOfFollowUpEvent,2)))
  
  # Merge and clean to get a single date of recurrence/endpoint/etc
  if (nrow(frr) <= 0) {
    for(n in colnames(frr)) feep <- add_column(feep, !!n := NA_character_)
    fer <- feep %>% dplyr::select(-FR.StudySubjectID, -FR.StudySite)
  } else {
    fer <- full_join(feep,frr,by=c('FE.StudySubjectID' = 'FR.StudySubjectID', 'FE.StudySite' = 'FR.StudySite')) 
  }
  
  fer <- fer %>% dplyr::mutate(FE.LatestReasonForFollowUp.c = sort(unique(na.omit(FE.LatestReasonForFollowUp.c), na.omit(FR.LatestReasonForFollowUp.c)))[1],
                  FE.DateOriginalDiseaseReoccurred = sort(unique(na.omit(c(FR.DateOriginalDiseaseReoccurred, FE.DateOriginalDiseaseReoccurred))))[1],
                  FE.LatestDateOfUpdate.c = sort(unique(na.omit(c(FR.LatestDateOfUpdate.c, FE.LatestDateOfUpdate.c))))[1],
                  FE.HasOriginalDiseaseReoccurred = unique(na.omit(FE.HasOriginalDiseaseReoccurred, FR.HasOriginalDiseaseReoccurred))[1]) %>%
    dplyr::select(-FR.LatestReasonForFollowUp.c, -FR.LatestDateOfUpdate.c, -FR.DateOriginalDiseaseReoccurred,-FR.HasOriginalDiseaseReoccurred)
  
  fer <- fer %>% group_by(FE.StudySubjectID) %>% 
    dplyr::mutate(Patient.Died.c = case_when(FE.EndPoint == 'Patient died' ~ 'yes', TRUE ~ 'no')) %>%
    ungroup %>% dplyr::mutate_at(vars(FE.EndPoint, matches('Died|Reason|HasOriginal')), list(~factor(.))) 
    
  return(list('fe'=fer))
}

#' Read data from prestaging CRF tables
#' There are 6 tables that appear to relate only to metastatic pre-staging. In these there are multiple entries per-patient with no
#' clear rationale for merging. Currently not needed for any analysis so leaving these for people to deal with on their own if they want them.
#' @name read.prestage
#' @param ...
#'
#' @author
read.prestage<-function(ocs, tables, occams_ids=NULL, rulesFiles=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  tables = sort(get.tables.by.crf(ocs, 'pretreatment staging'))

  if (length(tables) != 7)
    warning("Seven tables are expected for prestaging, currently only the final tables is read.")

  if (inherits(ocs, "OCCAMSLabkey")) {
    #stop("Create connection to OCCAMS Labkey first")
    #if (verbose) message(paste0("Reading pre-treatment tables:\n\t",paste(tables, collapse="\n\t")))
    
    ps <- get.table.data(ocs, grep('^ps_',tables, value=T), rulesFile=rulesFiles[1], colFilter=make.id.filter(occams_ids,'PS_StudySubjectID'))

    # These tables list different clinical modalities for staging (LAP, CT, PET) for metastasis only I think. For now not useful to us so leaving it off
    # eus_tbls <- grep('sampled|eus',tables, value=T)
    # ct_tbls <- grep('ct',tables, value=T)
    # lap_tbls <- grep('lap',tables, value=T)
    # pet_tbls <- grep('pet',tables, value=T)
    # fps_tbls <- grep('fps',tables, value=T)
    # 
    #     
    # ps$ID = ps$PS.StudySubjectID
    # staging = ps
    # for (i in 1:5) { #
    #   filter1 = make.id.filter(occams_ids, paste0('PS',i,'_StudySubjectID'))
    #   ps1 <- get.table.data(ocs,grep(paste0('^ps',i,'_'),tables, value=T), colFilter=filter1)
    #   ps1$ID = ps1[[group_vars(ps1)]]
    #   staging = full_join(staging, ps1, by='ID')
    # }
  } else {
    ps <- read.table.data(tables[grep('^ps_',basename(tables))], rulesFile=rulesFiles[1])
  }

  ps <- ps %>% dplyr::select(-matches('OpenClinica')) %>% ungroup %>%
    dplyr::mutate_at(vars(-matches('Study|PathologyReportNumber|Number')), list(~factor(.))) %>%
    dplyr::mutate_at(vars(matches('NumberOfNodes',ignore.case=T)), list(~as.numeric(.)))
  
  # There should not be duplicates in this table...
  dups <- which(with(ps, PS.StudySubjectID %in% names(which(table(PS.StudySubjectID) > 1))))
  duplicatePts <- ps[dups,]

  if (length(dups) > 0) ps <- ps[-dups, ]

  ps <- fix.tumor.factors(text.to.tstage(ps))

  ps <- ps %>% dplyr::group_by(PS.StudySubjectID) %>%
    # dplyr::mutate(PS.TNMStage.Final.c = tnmStage(PS.TStage.PrimaryTumour.FinalPretreatmentStaging, PS.NStage.PrimaryTumour.FinalPretreatmentStaging.TNM7, PS.MStage.PrimaryTumour.FinalPretreatmentStaging)) %>% 
    ungroup %>% mutate_at(vars(contains('NumberOf')), list(~as.integer(.)))
  
  if (verbose) message(paste0(nrow(ps), ' patients from ',grep('^ps_',tables, value=T)))
  
  return(list('ps' = ps, 'duplicates' = duplicatePts))
}

#' Read data from treatment plan CRF tables
#' @name read.treatment.plan
#' @param ...
#'
#' @author
read.treatment.plan<-function(ocs, tables, occams_ids=NULL, rulesFiles=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  tables = get.tables.by.crf(ocs, 'treatment plan')

  if (length(tables) < 2)
    stop("Two tables are expected for treatment plans, currently using only primary table.")

  if (inherits(ocs, "OCCAMSLabkey")) {
    if (verbose) message(paste("Reading",paste(tables, collapse=", ")))
    tp <- get.table.data(ocs=ocs, table=grep('tp_',tables,value=T), rulesFile=rulesFiles[1], colFilter=make.id.filter(occams_ids,'TP_StudySubjectID'))
    
    # This table just lists the hospital so leaving it out
    #tp1 <- get.table.data(ocs=ocs, table=grep('tp1_',tables,value=T), rulesFile=rulesFiles[2], colFilter=make.id.filter(occams_ids,'TP1_StudySubjectID'))
  } else {
    tp <- read.table.data( table=grep('tp_',tables,value=T), rulesFile=rulesFiles[1])
    #tp1 <- read.table.data(tables[2], rulesFile=rulesFiles[2])
  }
  
  # There are occasionally multiple entries typically when there's a treatment change. So not technically a duplicate but a separate entry that should be dealt with
  tp <- tp %>% dplyr::select(-matches('OpenClinica|SMDT|TP.DateFinalCarePlan')) %>% 
    group_by(TP.StudySubjectID, TP.StudySite) %>%
    dplyr::summarise_all( list(~paste(na.omit(.),collapse=' > ')) ) %>%
    dplyr::mutate_at(vars(-matches('Study')), list(~na_if(.,'')))
      
  if (verbose) message(paste(nrow(tp), "patients in ",grep('tp_',tables,value=T)))

  return(list('tp'=tp))
}

#' Read data from therapy plan CRF tables
#' Merged all of the various new therapy tables. Now if you want to look at whether patients recieved neoadj or endoscopic therapy you need to select the 
#' columns 'TR.NeoAdj' or 'TR.Endo'.  Bit confusing but currently easiest.
#' @name read.therapy
#' @param ...
#'
#' @author
read.therapy<-function(ocs, tables, occams_ids=NULL, rulesFiles=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  tables = get.tables.by.crf(ocs, 'therapy record')

  if (length(tables) < 7)
    warning("Seven tables are expected for therapy, currently only the primary table is read.")

  if (inherits(ocs, "OCCAMSLabkey")) {
    if (verbose) message(paste("Reading",paste(tables, collapse=",")))
    tr <- get.table.data(ocs,grep('^tr_', tables, value=T), colFilter=make.id.filter(occams_ids, 'TR_StudySubjectID'))

    read_and_exclude<-function(tb, exclude=NULL) {
      exclude = paste(c('OpenClinica',exclude),collapse = '|')
      tb %>% dplyr::select(-matches(exclude))
    }
    ## These were recently added so need to be merged with the older data in TR
    # neoadj chemo
    tr_neoadj <- get.table.data(ocs,grep('^trnc', tables, value=T), colFilter=make.id.filter(occams_ids, 'TRNC_StudySubjectID')) %>% read_and_exclude()
    # adj chemo
    tr_adj <- get.table.data(ocs,grep('^trac', tables, value=T), colFilter=make.id.filter(occams_ids, 'TRAC_StudySubjectID')) %>% read_and_exclude()
    # palliative chemo
    tr_pc <- get.table.data(ocs,grep('^trpc', tables, value=T), colFilter=make.id.filter(occams_ids, 'TRPC_StudySubjectID')) %>% read_and_exclude()
    # radiotherapy
    tr_rd <- get.table.data(ocs,grep('^trr', tables, value=T), colFilter=make.id.filter(occams_ids, 'TRR_StudySubjectID')) %>% read_and_exclude()
    # endoscopy
    tr_endo <- get.table.data(ocs,grep('^tre0_', tables, value=T), colFilter=make.id.filter(occams_ids, 'TRE0_StudySubjectID')) %>% read_and_exclude() %>%
      dplyr::rename_all(list(~sub('TRE0','TRE',.)))
    # other
    tr_other <- get.table.data(ocs,grep('^tro', tables, value=T), colFilter=make.id.filter(occams_ids, 'TRO_StudySubjectID')) %>% read_and_exclude()
  } else {
    tr <- read.table.data(grep('^tr_', tables, value=T), rulesFile=rulesFiles[1]) %>% read_and_exclude()
    tr_adj <- read.table.data(grep('^trac_', tables, value=T)) %>% read_and_exclude()
    tr_pc <- read.table.data(grep('^trpc_', tables, value=T)) %>% read_and_exclude()
    tr_rd <- read.table.data(grep('^trr_', tables, value=T)) %>% read_and_exclude()
    tr_endo <- read.table.data(grep('^tre0_', tables, value=T)) %>% read_and_exclude()
    tr_other <- read.table.data(grep('^tro_', tables, value=T)) %>% read_and_exclude()
  }

  # neoadj
  tr_neoadj = match_cols(tr,tr_neoadj) %>% dplyr::rename_at(vars(-matches('Study')), list(~paste0('TR.NeoAdj.',.))) 
  if (verbose) message(paste(nrow(tr_neoadj), 'neoadjuvant patients'))
  # endo
  tr_endo = match_cols(tr, tr_endo) %>% dplyr::rename_at(vars(-matches('Study')), list(~paste0('TR.Endo.',.)))
  if (verbose) message(paste(nrow(tr_endo), 'endoscopic patients'))
  # radio
  tr_rd = match_cols(tr, tr_rd) %>% dplyr::rename_at(vars(-matches('Study')), list(~paste0('TR.RT.',.)))
  if (verbose) message(paste(nrow(tr_rd), 'RFA patients'))
  # adjuvant -- the original table provided no way to enter both NeoAdj and Adj chemo so can't merge them
  tr_adj = tr_adj %>% dplyr::rename_at(vars(matches('Study')), list(~sub('^TR.+\\.','',.))) %>% 
    dplyr::rename_at(vars(-matches('Study')), list(~sub('^TR.+\\.', 'TR.Adj.',.)))
  if (verbose) message(paste(nrow(tr_adj), 'adjuvant patients'))
  # palliative (all therapies)
  tr_pc = match_cols(filter(tr, TR.TreatmentIntent == 'Palliative'), tr_pc) %>% dplyr::rename_at(vars(-matches('Study')), list(~paste0('TR.PalChemo.',.)))
  if (verbose) message(paste(nrow(tr_pc), 'palliative chemo patients'))
  # Other therapies
  tr_other = match_cols(tr, tr_other) %>% dplyr::rename_at(vars(-matches('Study')), list(~paste0('TR.Other.',.)))
  if (verbose) message(paste(nrow(tr_other), "'other' therapies patients"))

  tr <- tr %>% dplyr::rename_at(vars(matches('Study')), list(~sub('^TR\\.','',.))) %>% 
    dplyr::select(matches('Study|TreatmentIntent'))

  if (nrow(tr_neoadj) <= 0) {
    for (n in colnames(tr_neoadj)) tr <- add_column(tr, !!n := NA_character_)
    tr <- tr %>% dplyr::select(-TRNC.StudySubjectID, -TRNC.StudySite)
  } else {
    tr <- tr %>% dplyr::full_join(tr_neoadj, by=c('StudySubjectID','StudySite')) 
  }

  if (nrow(tr_endo) <= 0) {
    for (n in colnames(tr_endo)[-c(1:2)] ) tr <- add_column(tr, !!n := NA_character_)
    tr <- tr %>% dplyr::select(-TRNC.StudySubjectID, -TRNC.StudySite)
  } else {
    tr <- tr %>% dplyr::full_join(tr_endo, by=c('StudySubjectID','StudySite')) 
  }
  
  if (nrow(tr_rd) <= 0) {
    for (n in colnames(tr_rd)[-c(1:2)] ) tr <- add_column(tr, !!n := NA_character_)
  } else {
    tr <- tr %>% dplyr::full_join(tr_rd, by=c('StudySubjectID','StudySite')) 
  }  
    
  if (nrow(tr_adj) <= 0) {
    for (n in colnames(tr_adj)[-c(1:2)] ) tr <- add_column(tr, !!n := NA_character_)
  } else {
    tr <- tr %>% dplyr::full_join(tr_adj, by=c('StudySubjectID','StudySite'))
  }
  
  if (nrow(tr_pc) <= 0) {
    for (n in colnames(tr_pc)[-c(1:2)] ) tr <- add_column(tr, !!n := NA_character_)
  } else {
    tr <- tr %>% dplyr::full_join(tr_pc, by=c('StudySubjectID','StudySite')) 
  }

  if (nrow(tr_other) <= 0) {
    for (n in colnames(tr_other)[-c(1:2)] ) tr <- add_column(tr, !!n := NA_character_)
  } else {
    tr <- tr %>% dplyr::full_join(tr_other, by=c('StudySubjectID','StudySite')) 
  }

  if (verbose) message(paste(nrow(tr), "patients in treatment tables"))
  
  return(tr)
}

#' Read data from surgery CRF tables
#' @name read.surgery
#' @param ...
#'
#' @author
read.surgery<-function(ocs, table, occams_ids=NULL, rulesFile=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  table = get.tables.by.crf(ocs, 'surgical treatment')

  if (inherits(ocs, "OCCAMSLabkey")) {
    if (verbose) message(paste("Reading",table))
    st <- get.table.data(ocs,table, rulesFile=rulesFile, colFilter=make.id.filter(occams_ids,'ST_StudySubjectID'))
  } else {
    st <- read.table.data(table, rulesFile=rulesFile)
  }
  
  # There should not be duplicates in this table...
  dups <- which(duplicated(st[,'ST.StudySubjectID']))
  duplicatePts <- st[dups,]
  if (length(dups) > 0) st <- st[-dups, ]

  # A patient may be recorded as getting surgery (ST.MainSurgery) but ultimately not receive surgery. A patient may also undergo surgery but have no resection pathology as the surgeon will have not been able to complete it.
  st <- st %>% group_by(ST.StudySubjectID) %>% 
    dplyr::mutate(
    ST.SurgeryPerformed.c = factor(ifelse(ST.MainSurgery == 'yes' & (is.na(ST.ReasonIfNoSurgery) | is.na(ST.OtherReasonForNoSurgery)), 'yes', 'no'), levels=c('yes','no')),
    ST.PathologyReportGenerated.c = factor(ifelse(ST.SurgeryPerformed.c == 'yes' & !grepl("Open (and|&) Shut", ST.Procedure, ignore.case=T),'yes', 'no'), levels=c('yes','no'))) %>% ungroup %>% 
    dplyr::mutate_at(vars(ST.ASAGrade, matches('Percentage')), list(~as.numeric(.))) %>%
    dplyr::select(-matches('OpenClinica'))

  if(verbose) message(paste(nrow(st), "surgery patients"))

  return(st)
}

#' Read data from pathology CRF table
#' @name read.pathology
#' @param ...
#'
#' @author
read.pathology<-function(ocs, table, occams_ids=NULL, rulesFile=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  table = get.tables.by.crf(ocs, 'resection pathology')

  if (inherits(ocs, "OCCAMSLabkey")) {
    if (verbose) message(paste("Reading",table))
    rp <- get.table.data(ocs=ocs, table=table, uniqueID=NULL, rulesFile=rulesFile, colFilter=make.id.filter(occams_ids,'RP_StudySubjectID'))
  } else {
    rp <- read.table.data(table, uniqueID=NULL, rulesFile=rulesFile)
  }
  dups <- which(duplicated(rp$RP.StudySubjectID))
  duplicatePts <- rp[dups,]

  if (length(dups) > 0) rp <- rp[-dups, ]

  # The grades are too granular.
  diff.grade<-function(x) {
    if ( grepl('.*poor',x)  ) { # poor & moderate_to_poor
      return('poor')
    } else if ( grepl('.*well',x) ) { # well & moderate_to_well
      return('well')
    }
    return(x)
  }

  # Add a TNM column assessed based on the Tstage, Nstage, and Mstage. Add a differentiation column with fewer variables, recode the response column for easier reading. Finally, add a BE adjacent column that assess both Microscopic and Macroscopic columns for a binary response.
  if (nrow(rp) > 0) {
    rp <- fix.tumor.factors(text.to.tstage(rp))
  
    rp <- rp %>% ungroup %>% group_by(RP.StudySubjectID) %>% 
      dplyr::mutate(
        #RP.TNMStage.c = tnmStage(RP.TStage.PrimaryTumour, RP.Nstage.RP.TNM7, RP.MStage.DistantMetastasis),
        RP.TumourDifferentiation.c = diff.grade(RP.TumourGradingDifferentiationStatus),
        RP.TumourResponse = recode_factor(RP.TumourResponse,'0pc_remaining'='0%','less_than_20pc'='<20%', 'greater_than_or_equals_to_20pc'='≥20%', 'less_than_50pc'='<50%', 'greater_than_or_equals_to_50pc'='≥50%', .ordered=T),
        RP.Location = factor(RP.Location, levels=c('oesophageal','goj','gastric')),
        RP.MandardScoreForResponse = factor(RP.MandardScoreForResponse, levels=c('TRG1','TRG2','TRG3','TRG4','TRG5')),
        # If you see it macro, then it is there micro as well
        RP.BarrettsAdjacentToTumour.c = ifelse(RP.BarettsAdjacentToTumourMicroscopicIM == 'yes' | RP.BarettsAdjacentToTumourMacroscopic == 'yes','yes','no')
      ) %>% dplyr::select(-matches('OpenClinica'))

    # This was a check performed by skillcoyne and chughes in 2016 Sept to evaluate the entry of BE information into OpenClinica. Only AH patients were evaluated.  This file alters the BE adjacent information for only those patients.
    file = system.file("extdata", "be_updates_20160930.txt", package="openclinica.occams")
    if (!exists('be_updates') & (!is.null(file) & file.exists(file))) {
      if (verbose) message(paste("Reading", file))
      be_updates = readr::read_tsv(file, col_types='ccc')
  
      if (verbose) message("Updating RP.BarrettsAdjacent patients from Caitrona's worksheet.")
      be_updates = be_updates %>% dplyr::mutate( `Barret's Confirmed`= ifelse(`Barret's Confirmed` == '?', NA, `Barret's Confirmed`),
                                                 `Barret's Confirmed` = recode(`Barret's Confirmed`, 'N'='no','Y'='yes'))

      # Assume Caitrona's are the final say (they matched other than NA anyhow)
      rp = left_join(rp, be_updates, by=c('RP.StudySubjectID'='OCCAMS/ID')) %>% group_by(RP.StudySubjectID) %>% 
        dplyr::mutate(RP.BarrettsAdjacentToTumour.c = mr(`Barret's Confirmed`, RP.BarrettsAdjacentToTumour.c,1)) %>% 
        dplyr::select(-Source, -contains('Confirmed'))
    }

    rp = rp %>% ungroup %>% dplyr::mutate(RP.BarrettsAdjacentToTumour.c = factor(RP.BarrettsAdjacentToTumour.c)) %>% 
      dplyr::select(-contains('Nstage.RP.TNMSystem'),-contains('RP.TNM6'), -matches('PathologistName')) %>% group_by(RP.StudySubjectID)
  } else {
    rp <- rp %>% dplyr::mutate_at(vars(StudySite, StudySubjectID), list(~as.character(.))) -> rp
  }

  if(verbose) message(paste(nrow(rp), "pathology patients"))

  return(rp)
}

#' Read data from referral & diagnosis CRF tables
#' @name read.referral.diagnosis
#' @param ...
#'
#' @author
read.referral.diagnosis<-function(ocs, tables, occams_ids=NULL, rulesFile=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  tables = get.tables.by.crf(ocs, 'referral diagnosis')

  if (inherits(ocs, "OCCAMSLabkey")) {
    if (verbose) message(paste("Reading",tables))
    rd <- get.table.data(ocs, grep('^rd_', tables, value=T),  rulesFile=rulesFile, colFilter=make.id.filter(occams_ids,'RD_StudySubjectID'))
    rd2 <- get.table.data(ocs, grep('^rda_', tables, value=T),  rulesFile=NULL, colFilter=make.id.filter(occams_ids, 'RDA_StudySubjectID'))
  } else {
    rd <- read.table.data(tables[grep('^rd_', basename(tables))],  rulesFile=rulesFile)
  }

  rd <- rd %>% dplyr::select(-matches('OpenClinica'))
  rd2 <- rd2 %>% dplyr::select(-matches('OpenClinica'))
  
  merged <- match_cols(rd,rd2) #%>% dplyr::rename_all(list(~paste0('RD.',.)))

  rd <- rd %>% ungroup %>% 
    dplyr::select(RD.StudySite, RD.StudySubjectID, matches('Priority|Referral|Site|Status|Dysplasia'),RD.Histology,RD.OGD.Path.Number,matches('SignetRing'), matches('Date')) %>% dplyr::left_join(merged, by=c('RD.StudySite'='StudySite', 'RD.StudySubjectID'='StudySubjectID')) %>% dplyr::group_by(RD.StudySubjectID)
  
  if(verbose) message(paste(nrow(rd), "referral & diagnosis patients"))
  return(rd)
}

# Not a clinical table and I don't know how up to date this is but I pull it anyhow.  I should do the same for blood and ctdna but currently I don't.
read.tissue.collection<-function(ocs, tables, rulesFiles=NULL, occams_ids=NULL, verbose=F, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  tables = list.clinical.tables(ocs,'tc')

  filter = make.id.filter(occams_ids, 'TC_StudySubjectID')

  if (!inherits(ocs, "OCCAMSLabkey"))
    stop("Create connection to OCCAMS Labkey first")
  if (length(tables) < 2)
    stop("Two tables are expected for tissue collection")

  if (verbose) message(paste("Reading",paste(tables, collapse=",")))

  # all this table tells me is if there should be entries in the second table or not
  tc <- get.table.data(ocs=ocs, table=grep('^tc_',tables,value=T), colFilter=filter)

  # There should not be duplicates in this table...
  dups <- which(with(tc, TC.StudySubjectID %in% names(which(table(TC.StudySubjectID) > 1))))
  duplicatePts <- tc[dups,]

  warning(paste("Removing", nrow(duplicatePts), "duplicated patient rows"))
  if (length(dups) > 0) tc <- tc[-dups, ]

  samplesTaken <- tc %>% filter(TC.TissueSamplesTaken == "Yes")

  # useful table
  filter = make.id.filter(occams_ids, 'TC1_StudySubjectID')
  tc1 <- get.table.data(ocs=ocs, table=grep('^tc1_',tables,value=T), colFilter=filter)

  # But it seems the two tables may not actually match up! There are more patients in the second table then are listed in the first as having had tissue taken
  missingEntries = setdiff(samplesTaken$TC.StudySubjectID, tc1$TC1.StudySubjectID)

  # Recode for clarity
  tc1 = tc1 %>% dplyr::mutate(
    TC1.TissueType = recode(TC1.TissueType, 'E'="endoscopy", 'S'="surgical",'L'="laparoscopy",'R'="EMR"),
    TC1.TissueSource = recode(TC1.TissueSource, 'N'="normal oesophagus", 'B'="barretts", 'T'="tumour", 'G'="normal gastric", 'L'="lymph node", 'M'="metastasis")
  )

  tc1 = left_join(tc1, tc, by=c('TC1.StudySubjectID'='TC.StudySubjectID')) %>% dplyr::rename('StudySubjectID'='TC1.StudySubjectID') %>%
    dplyr::select(StudySubjectID, starts_with('TC.'), starts_with('TC1.')) %>% group_by(StudySubjectID)

  if (verbose) message(paste(nrow(tc1), "tissues"))

  return(tc1)
}

