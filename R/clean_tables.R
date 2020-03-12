## This file provides a function to read/clean each section of the CRF tables.  If it has only one table, like demographics, then it reads just that one.  If it has multiple tables, like the treatment tables or endpoint tables then it works to merge them such that there is only a single entry per patient. This isn't always what you may want, and the logic is more likely to run into issues with changes on the Labkey side.


merge.patient.tables<-function(df1, df2, ...) {
  newDF <- merge.data.frame(df1, df2, by='row.names', ...)
  rownames(newDF) <- newDF[['Row.names']]
  newDF[['Row.names']] <- NULL
  return(newDF)
}

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

  # Prefer to get this DoD from FE tables
  di <- di %>% select(-DI.StudyNumber, -contains('Panther'), -contains('OpenClinica'), -contains('CRF'),-contains('DateOfDeath'), -DI.PatientEthnicityComment, -DI.ageAtDiagnosis.months)

  if (verbose) message(paste(nrow(di), "patients"))
  return(di)
}

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
    ex = ex %>% 
      mutate_at(vars(grep('Height|Weight|BMI|Weeks|Months|Years|Days', colnames(.))), funs(as.numeric)) %>%
      mutate_if( is.numeric, funs(ifelse(.<0,NA,.))) %>% 
      group_by(EX.StudySubjectID) %>% 
      mutate(EX.CurrentWeightKG = ifelse(EX.CurrentWeightKG>140 | EX.CurrentWeightKG<22, NA, EX.CurrentWeightKG),
                     EX.CurrentHeightCM = ifelse(EX.CurrentHeightCM>215 | EX.CurrentHeightCM<60, NA, EX.CurrentHeightCM),
                     EX.CurrentBMI.c = bmicalc(EX.CurrentWeightKG,EX.CurrentHeightCM),
                     EX.5YearsAgoBMI.c = bmicalc(EX.5YearsAgoWeightKG,EX.5YearsAgoHeightCM))

  # These two are always NA, seem to have moved to the other tables
  ex = dplyr::select(ex, -EX.FamilyHistoryOfOtherCancerRelationship, -EX.FamilyHistoryOfOesophagoGastricCancerRelationship)
  
  if (verbose) message(paste(nrow(ex), "patients"))

  ex_history = full_join(ex_history_gastric, ex_history_other, by=c("EX1.StudySubjectID"="EX2.StudySubjectID")) %>%
    dplyr::select(-matches('StudySite')) %>%
    dplyr::rename_all(list(~sub('^EX\\d','EX',.))) %>% ungroup %>%
    dplyr::mutate_at(vars(EX.FamilyHistoryOfOesophagoGastricCancerLocation, EX.FamilyHistoryOfOesophagoGastricCancerRelationship,EX.FamilyHistoryOfOtherCancerType,EX.FamilyHistoryOfOtherCancerRelationship), list(~factor(.)))

  return(list('ex'=ex, 'history'=ex_history))
}


read.endpoint<-function(ocs, tables, occams_ids=NULL, rulesFile=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)
  
  table = get.tables.by.crf(ocs, c('endpoint'))
  filter = make.id.filter(occams_ids, 'EP_StudySubjectID')
  
  if (inherits(ocs, "OCCAMSLabkey")) {
    #stop("Create connection to OCCAMS Labkey first")
    if (verbose) message(paste("Reading",paste(tables, collapse=",")))
    ep <- get.table.data(ocs,grep('ep',tables,value=T), rulesFile=rulesFile, colFilter=filter)
  } else {
    if (verbose) message(paste("Reading",paste(basename(tables), collapse=", ")))
    ep <- read.table.data(table,  rulesFile=rulesFile)
  }
  ep %<>% ungroup %>% dplyr::mutate_at(vars(EP.EndPoint,EP.ReasonForPatientDeath), list(~factor(.)))
  
  ep %<>% group_by(EP.StudySubjectID) %>% 
    dplyr::mutate(Patient.Died = ifelse(EP.EndPoint == 'Patient died', 'yes','no') ) %>%
    ungroup %>% dplyr::mutate(Patient.Died = factor(Patient.Died))
  
  return(ep)
}


# TODO need to add the recurrence tables to this
read.followup<-function(ocs, tables, occams_ids=NULL, rulesFiles=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  tables = get.tables.by.crf(ocs, c('endpoint','follow-up', 'recurrence'))

  if (length(tables) < 6)
    stop("Six tables expected for endpoint, follow-up, and recurrence.")

  filterEP = make.id.filter(occams_ids, 'EP_StudySubjectID')
  filterFE = make.id.filter(occams_ids,'FE_StudySubjectID')
  filterFE1 = make.id.filter(occams_ids,'FE1_StudySubjectID')
  filterR = make.id.filter(occams_ids, 'FR_StudySubjectID')
  filterR0 = make.id.filter(occams_ids, 'FR0_StudySubjectID')
  filterR1 = make.id.filter(occams_ids, 'FR1_StudySubjectID')

  if (inherits(ocs, "OCCAMSLabkey")) {
    #stop("Create connection to OCCAMS Labkey first")
    if (verbose) message(paste("Reading tables:\n\t",paste(tables, collapse="\n\t ")))

    # EP appears to be set up to replace FE, except new values are continuing to appear in both so they'll need some merging
    ep <- get.table.data(ocs,grep('ep',tables,value=T), rulesFile=grep('ep_', rulesFiles, value=T), colFilter=filterEP)
    fe <- get.table.data(ocs,grep('fe_',tables,value=T),  rulesFile=grep('fe_', rulesFiles,value=T), colFilter=filterFE)

    feep = full_join(fe,ep,by=c('FE.StudySubjectID'='EP.StudySubjectID')) %>%
      dplyr::group_by(FE.StudySubjectID) %>%
      dplyr::mutate(
        FE.StudySite = mr(FE.StudySite, EP.StudySite),
        FE.EndPoint = case_when(
          (FE.EndPoint == 'Patient died' | EP.EndPoint == 'Patient died') ~ 'Patient died', 
          (FE.EndPoint == 'Follow-up ended' | EP.EndPoint == 'Follow-up ended') ~ 'Follow-up ended', 
          TRUE ~ as.character(mr(FE.EndPoint, EP.EndPoint)) ),
        FE.DateOfPatientDeath = mr(FE.DateOfPatientDeath,EP.DateOfPatientDeath),
        FE.ReasonForPatientDeath = mr(FE.ReasonForPatientDeath, EP.ReasonForPatientDeath),
        FE.ReasonForPatientDeathOther = mr(FE.ReasonForPatientDeathOther, EP.ReasonForPatientDeathOther) ) %>% 
      dplyr::select(starts_with('FE'), -contains('OpenClinica'), -matches('FE.ageAtDeath'))

    fe1 <- get.table.data(ocs,grep('fe1_',tables,value=T),  rulesFile=grep('fe1_', rulesFiles,value=T), colFilter=filterFE1)
    
    dis.r<-function(hasR, dateR, dateU, i) {
      n = which(hasR == 'Yes')[1]
      
      if (length(n) <= 0) {
        return( c('No',NA)[i] )
      } else if (i == 1) {
        return(hasR[n])
      } else {
        date = dateR[n]
        if (is.na(date)) date = dateU[n]
        return( as.Date(date) )
      }
    return(NA)            
    }
    
    fe1 %<>% dplyr::group_by(FE1.StudySubjectID) %>%
      arrange(desc(FE1.DateOfUpdate)) %>% 
      #filter(!grepl('death', FE1.ReasonForFollowUp)) %>%
      dplyr::summarise(FE.LatestDateOfUpdate.c=FE1.DateOfUpdate[1],
                       FE.LatestReasonForFollowUp.=FE1.ReasonForFollowUp[1],
                       #FE.HasOriginalDiseaseReoccurred.c=FE1.HasOriginalDiseaseReoccurred[1],
                       FE.HasOriginalDiseaseReoccurred = dis.r(FE1.HasOriginalDiseaseReoccurred,FE1.DateOriginalDiseaseReoccurred,FE1.DateOfUpdate,1),
                       FE.DateOriginalDiseaseReoccurred = dis.r(FE1.HasOriginalDiseaseReoccurred,FE1.DateOriginalDiseaseReoccurred,FE1.DateOfUpdate,2) ) 
                         
    feep = left_join(feep, fe1, by=c('FE.StudySubjectID'='FE1.StudySubjectID'))
    
    
    # Reoccurance
    
    filterR = make.id.filter(occams_ids, 'FR_StudySubjectID')
    filterR0 = make.id.filter(occams_ids, 'FR0_StudySubjectID')
    filterR1 = make.id.filter(occams_ids, 'FR1_StudySubjectID')
    
    
    fr <- get.table.data(ocs,grep('fr_',tables,value=T), rulesFile=NULL, colFilter=filterR)
    
    
    
  } else {
    #if (verbose) message(paste("Reading",paste(basename(tables), collapse=", ")))
    #fe <- read.table.data(tables[1],  rulesFile=rulesFiles[1])
    #fe2 <- read.table.data(tables[2], rulesFile=rulesFiles[2])
    #fe2 <- ddply( ddply(fe2, .(FE1.StudySubjectID), arrange, desc(FE1.DateOfUpdate)),
    #              .(FE1.StudySubjectID, FE1.StudySite), plyr::summarise, FE1.DateOfUpdate=FE1.DateOfUpdate[1], FE1.ReasonForFollowUp=FE1.ReasonForFollowUp[1], FE1.HasOriginalDiseaseReoccurred=FE1.HasOriginalDiseaseReoccurred[1])
  }
  
  feep %<>% ungroup %>% dplyr::mutate_at(vars(FE.EndPoint, FE.ReasonForPatientDeath), list(~factor(.)))

  
  
  

  return(list('fe'=feF, 'withdrawn'=withdrawn$FE.StudySubjectID))
}

read.prestage<-function(ocs, tables, occams_ids=NULL, rulesFiles=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  tables = get.tables.by.crf(ocs, 'pretreatment staging')

  if (length(tables) != 7)
    warning("Seven tables are expected for prestaging, currently only the final tables is read.")

  filter = make.id.filter(occams_ids,'PS_StudySubjectID')

  if (inherits(ocs, "OCCAMSLabkey")) {
    #stop("Create connection to OCCAMS Labkey first")
    if (verbose) message(paste("Reading",paste(tables, collapse=",")))
    ps <- get.table.data(ocs,grep('^ps_',tables, value=T), rulesFile=rulesFiles[1], colFilter=filter)
    # ps$ID = ps$PS.StudySubjectID
    #
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

  # There should not be duplicates in this table...
  dups <- which(with(ps, PS.StudySubjectID %in% names(which(table(PS.StudySubjectID) > 1))))
  duplicatePts <- ps[dups,]

  if (length(dups) > 0) ps <- ps[-dups, ]
  #ps <- rows.as.patient.id(ps, 2)

  ps <- fix.tumor.factors(text.to.tstage(ps))

  ps = ps %>% rowwise %>% dplyr::mutate(PS.TNMStage.Final.c = tnmStage(PS.TStage.PrimaryTumour.FinalPretreatmentStaging, PS.NStage.PrimaryTumour.FinalPretreatmentStaging.TNM7, PS.MStage.PrimaryTumour.FinalPretreatmentStaging)) %>% ungroup %>% mutate_at(vars(contains('NumberOf')), funs(as.integer))

  if (verbose) message(paste(nrow(ps), "patients"))

  return(list('ps' = ps, 'duplicates' = duplicatePts))
}


read.treatment.plan<-function(ocs, tables, occams_ids=NULL, rulesFiles=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  tables = get.tables.by.crf(ocs, 'treatment plan')

  if (length(tables) < 2)
    stop("Two tables are expected for treatment plans, currently using only primary table.")

  filter = make.id.filter(occams_ids,'TP_StudySubjectID')
  filter1 = make.id.filter(occams_ids,'TP1_StudySubjectID')

  if (inherits(ocs, "OCCAMSLabkey")) {
    if (verbose) message(paste("Reading",paste(tables, collapse=",")))
    tp <- get.table.data(ocs=ocs, table=grep('tp_',tables,value=T), rulesFile=rulesFiles[1], colFilter=filter)
    tp1 <- get.table.data(ocs=ocs, table=grep('tp1_',tables,value=T), rulesFile=rulesFiles[2], colFilter=filter1)
  } else {
    #tp <- read.table.data(tables[1], rulesFile=rulesFiles[1])
    #tp1 <- read.table.data(tables[2], rulesFile=rulesFiles[2])
  }
  # There should not be duplicates in this table...
  dups <- which(with(tp, TP.StudySubjectID %in% names(which(table(TP.StudySubjectID) > 1))))
  duplicatePts <- tp[dups,]

  if (length(dups) > 0) tp <- tp[-dups,]

  if (verbose) message(paste(nrow(tp), "patients"))

  return(list('tp'=tp,'duplicates'=duplicatePts))
}

# TODO, need to find a way to read/include the other tables
read.therapy<-function(ocs, tables, occams_ids=NULL, rulesFiles=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  tables = get.tables.by.crf(ocs, 'therapy record')

  if (length(tables) < 7)
    warning("Seven tables are expected for therapy, currently only the primary table is read.")

  filter = make.id.filter(occams_ids, 'TR_StudySubjectID')

  if (inherits(ocs, "OCCAMSLabkey")) {
    if (verbose) message(paste("Reading",paste(tables, collapse=",")))
    tr <- get.table.data(ocs,grep('^tr_', tables, value=T), rulesFile=rulesFiles[1], colFilter=filter)

    cl<-function(x) {
      if (is.na(unique(x))) return('')
      a = paste(na.omit(unique(x)),collapse=',')
      return(a)
    }

    # neoadj chemo
    trnc <- get.table.data(ocs,grep('^trnc', tables, value=T))
    trnc = full_join(tr, trnc, by=c('TR.StudySubjectID'='TRNC.StudySubjectID'))

    # Select one of the two records to report. This is because a new table was added that will record patients (starting Oct 2018 I think) but the old table still exists. No data was backfilled.
    trnc = trnc %>% dplyr::group_by(TR.StudySubjectID) %>% dplyr::summarise(
      TR.NeoAdj.Chemotherapy.Treatments = n(),
      TR.NeoAdj.ChemotherapyDateFirstCycleStarted = mr(TR.ChemotherapyDateFirstCycleStarted, TRNC.ChemotherapyDateFirstCycleStarted),
      TR.NeoAdj.ChemotherapyNumberOfCyclesPlanned = mr(TR.ChemotherapyNumberOfCyclesPlanned, TRNC.ChemotherapyNumberOfCyclesPlanned),
      TR.NeoAdj.ChemotherapyNumberOfCyclesGiven = mr(TR.ChemotherapyNumberOfCyclesGiven, TRNC.ChemotherapyNumberOfCyclesGiven),
      TR.NeoAdj.ChemotherapyTreatmentProtocol = paste(mr(TR.ChemotherapyTreatmentProtocol, TRNC.ChemotherapyTreatmentProtocol),collapse=','),
      TR.NeoAdj.ChemotherapyRegimenEpirubicin = mr(TR.ChemotherapyRegimenEpirubicin, TRNC.ChemotherapyRegimenEpirubicin,empty=''),
      TR.NeoAdj.ChemotherapyRegimenCisplatin = mr(TR.ChemotherapyRegimenCisplatin, TRNC.ChemotherapyRegimenCisplatin,empty=''),
      TR.NeoAdj.ChemotherapyRegimenOxaliplatin = mr(TR.ChemotherapyRegimenOxaliplatin, TRNC.ChemotherapyRegimenOxaliplatin,empty=''),
      TR.NeoAdj.ChemotherapyRegimen5FU = mr(TR.ChemotherapyRegimen5FU, TRNC.ChemotherapyRegimen5FU,empty=''),
      TR.NeoAdj.ChemotherapyRegimenCapecitabine = mr(TR.ChemotherapyRegimenCapecitabine, TRNC.ChemotherapyRegimenCapecitabine,empty=''),
      TR.NeoAdj.ChemotherapyRegimenLapatinib = mr(TR.ChemotherapyRegimenLapatinib, TRNC.ChemotherapyRegimenLapatinib,empty=''),
      TR.NeoAdj.ChemotherapyRegimenOther = mr(TR.ChemotherapyRegimenOther, TRNC.ChemotherapyRegimenOther,empty=''),
      TR.NeoAdj.ChemotherapyOutcomeOfTreatment = mr(TR.ChemotherapyOutcomeOfTreatment, TRNC.ChemotherapyOutcomeOfTreatment,empty=''),
      TR.NeoAdj.ChemotherapyOutcomeOfTreatmentNotes = mr(TR.ChemotherapyOutcomeOfTreatmentNotes, TRNC.ChemotherapyOutcomeOfTreatmentNotes,empty=''),
      TR.NeoAdj.Response = mr(TR.Response, TRNC.Response,empty=''),
      TR.NeoAdj.ResponseModality = TRNC.ResponseModality[1],
      TR.NeoAdj.ResponseDate = sort(na.omit(TR.ResponseDate, TRNC.ResponseDate))[1],
      TR.NeoAdj.TherapyType = unique(na.omit(TRNC.TherapyType))[1]
    )

    tr = left_join(trnc, tr %>% select(-contains('Chemotherapy')), by='TR.StudySubjectID') %>% group_by(TR.StudySubjectID) %>% select(-contains('NeoAdj'), everything(), -contains('OpenClinica')) %>% mutate_if(is.character,funs(ifelse(.=='', NA, .)))

    # endo therapy
    tre <- get.table.data(ocs,grep('^tre0_', tables, value=T))
    tr = full_join(tr, tre %>% select(-contains('StudySite')), by=c('TR.StudySubjectID'='TRE0.StudySubjectID')) %>% rename_all(funs(sub('^TRE0', 'TR',.))) %>% dplyr::group_by(TR.StudySubjectID) %>% select(-contains('StudySite')) %>% dplyr::summarise_all(funs(cl)) %>% ungroup %>% mutate_all(funs(ifelse(.=='', NA, .))) %>% group_by(TR.StudySubjectID)

    # TODO: These tabels were added only recently and they were not backfilled so there are not many patients in each, can deal with them later when/if needed
    # radiotherapy -- later
    trr <- get.table.data(ocs,grep('^trr_', tables, value=T))

    # adjuvant chemo -- later
    trac <- get.table.data(ocs,grep('^trac_', tables, value=T))

    # palliative chemo -- later
    trac <- get.table.data(ocs,grep('^trpc_', tables, value=T))

    # other therapy -- later
    tro <- get.table.data(ocs,grep('^tro_', tables, value=T))

  } else {
    #tr <- read.table.data(tables[grep('^tr_', basename(tables))],  rulesFile=rulesFiles[1])
  }
# TODO...this is missing patients or setting them as "no" when they are NA but have neoadj therapy info
  tr = tr %>% dplyr::mutate( TR.Chemotherapy = factor(ifelse((grepl('chemo', TR.CurativeTreatmentModality, ignore.case=T) |
                                                                grepl('chemo', TR.PalliativeTreatmentModality, ignore.case=T)), 'yes','no'),levels=c('yes','no') ),

                             TR.Radiotherapy = factor(ifelse((grepl('radio', TR.CurativeTreatmentModality, ignore.case=T) |
                                                                grepl('radio', TR.PalliativeTreatmentModality, ignore.case=T)), 'yes','no'),levels=c('yes','no') ),

                             TR.EndoscopicTherapy = factor(ifelse((grepl('endoscop[y|ic]', TR.CurativeTreatmentModality, ignore.case=T) |
                                                                     grepl('endoscop[y|ic]', TR.PalliativeTreatmentModality, ignore.case=T)), 'yes','no'),levels=c('yes','no') ) )


  neoadj.regimen<-function(fluro, cisp, epi, cape, oxal) {
    if (sum(fluro, cisp, epi, na.rm=T) == 3) return('ECF')
    if (sum(cape, cisp, epi, na.rm=T) == 3) return('ECX')
    if (sum(cape, oxal, epi, na.rm=T) == 3) return('EOX')

    if (sum(fluro, cisp,na.rm=T) == 2) return('CF')
    if (sum(cape, cisp,na.rm=T) == 2) return('CC')
    if (sum(cape, oxal,na.rm=T) == 2) return('CO')

    return(NA)
  }

  tr = tr %>% mutate_at(vars(contains('NeoAdj.ChemotherapyRegimen')), funs(recode), 'yes'=1, 'no'=0) %>% rowwise %>% dplyr::mutate(
    TR.NeoAdj.Regimen = (neoadj.regimen(TR.NeoAdj.ChemotherapyRegimen5FU, TR.NeoAdj.ChemotherapyRegimenCisplatin, TR.NeoAdj.ChemotherapyRegimenEpirubicin, TR.NeoAdj.ChemotherapyRegimenCapecitabine, TR.NeoAdj.ChemotherapyRegimenOxaliplatin)))

  if (verbose) message(paste(nrow(tr), "patients"))

  return(tr)
}


read.surgery<-function(ocs, table, occams_ids=NULL, rulesFile=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  table = get.tables.by.crf(ocs, 'surgical treatment')

  if (inherits(ocs, "OCCAMSLabkey")) {
    if (verbose) message(paste("Reading",table))

    filter = make.id.filter(occams_ids,'ST_StudySubjectID')
    st <- get.table.data(ocs,table, rulesFile=rulesFile,colFilter=filter)
  } else {
    st <- read.table.data(table, rulesFile=rulesFile)
  }
  # There should not be duplicates in this table...
  dups <- which(duplicated(st[,'ST.StudySubjectID']))
  duplicatePts <- st[dups,]
  if (length(dups) > 0) st <- st[-dups, ]

  # A patient may be recorded as getting surgery (ST.MainSurgery) but ultimately not receive surgery. A patient may also undergo surgery but have no resection pathology as the surgeon will have not been able to complete it.
  st = st %>% rowwise %>% dplyr::mutate(
    ST.SurgeryPerformed = factor(ifelse(ST.MainSurgery == 'yes' & (is.na(ST.ReasonIfNoSurgery) | is.na(ST.OtherReasonForNoSurgery)), 'yes', 'no'), levels=c('yes','no')),
    ST.PathologyReportGenerated = factor(ifelse(ST.SurgeryPerformed == 'yes' & !grepl("Open (and|&) Shut", ST.Procedure, ignore.case=T),'yes', 'no'), levels=c('yes','no')),
    ST.ASAGrade = as.numeric(ST.ASAGrade),
    ST.FEV1PercentagePredicted = as.numeric(ST.FEV1PercentagePredicted),
    ST.FVC1PercentagePredicted = as.numeric(ST.FVC1PercentagePredicted)
  )

  if(verbose) message(paste(nrow(st), "patients"))

  return(st)
}

read.pathology<-function(ocs, table, occams_ids=NULL, rulesFile=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  table = get.tables.by.crf(ocs, 'resection pathology')

  filter = make.id.filter(occams_ids,'RP_StudySubjectID')

  if (inherits(ocs, "OCCAMSLabkey")) {
    if (verbose) message(paste("Reading",table))
    rp <- get.table.data(ocs=ocs, table=table, uniqueID=NULL, rulesFile=rulesFile, colFilter=filter)
  } else {
    rp <- read.table.data(table, uniqueID=NULL, rulesFile=rulesFile)
  }

  dups <- which(duplicated(rp$RP.StudySubjectID))
  duplicatePts <- rp[dups,]

  if (length(dups) > 0) rp <- rp[-dups, ]

  rp <- fix.tumor.factors(text.to.tstage(rp))

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
  rp = rp %>% ungroup %>% rowwise %>% dplyr::mutate(
    RP.TNMStage.c = tnmStage(RP.TStage.PrimaryTumour, RP.Nstage.RP.TNM7, RP.MStage.DistantMetastasis),
    RP.TumourDifferentiation.c = diff.grade(RP.TumourGradingDifferentiationStatus)) %>%
    dplyr::mutate(
      RP.TumourResponse = recode_factor(RP.TumourResponse,'0pc_remaining'='0%','less_than_20pc'='<20%', 'greater_than_or_equals_to_20pc'='≥20%', 'less_than_50pc'='<50%', 'greater_than_or_equals_to_50pc'='≥50%', .ordered=T),
      RP.Location = factor(RP.Location, levels=c('oesophageal','goj','gastric')),
      RP.MandardScoreForResponse = factor(RP.MandardScoreForResponse, levels=c('TRG1','TRG2','TRG3','TRG4','TRG5')),
      # If you see it macro, then it is there micro as well
      RP.BarrettsAdjacentToTumour.c = ifelse(RP.BarettsAdjacentToTumourMicroscopicIM == 'yes' | RP.BarettsAdjacentToTumourMacroscopic == 'yes','yes','no')
    )

  # This was a check performed by skillcoyne and chughes in 2016 Sept to evaluate the entry of BE information into OpenClinica. Only AH patients were evaluated.  This file alters the BE adjacent information for only those patients.
  file = system.file("extdata", "be_updates_20160930.txt", package="openclinica.occams")
  if (!exists('be_updates') & (!is.null(file) & file.exists(file))) {
    if (verbose) message(paste("Reading", file))
    be_updates = readr::read_tsv(file, col_types='ccc')

    if (verbose) message("Updating RP.BarrettsAdjacent patients from Caitrona's worksheet.")
    be_updates = be_updates %>% dplyr::mutate( `Barret's Confirmed`= ifelse(`Barret's Confirmed` == '?', NA, `Barret's Confirmed`),
                                               `Barret's Confirmed` = recode(`Barret's Confirmed`, 'N'='no','Y'='yes'))

    # Assume Caitrona's are the final say (they matched other than NA anyhow)
    rp = left_join(rp, be_updates, by=c('RP.StudySubjectID'='OCCAMS/ID')) %>% rowwise %>% dplyr::mutate(
      RP.BarrettsAdjacentToTumour.c = mr(`Barret's Confirmed`, RP.BarrettsAdjacentToTumour.c,1) ) %>% select(-Source, -contains('Confirmed'))
  }

  rp = rp %>% ungroup %>% mutate(RP.BarrettsAdjacentToTumour.c = factor(RP.BarrettsAdjacentToTumour.c)) %>% select(-contains('Nstage.RP.TNMSystem'),-contains('RP.TNM6')) %>% group_by(RP.StudySubjectID)

  if(verbose) message(paste(nrow(rp), "patients"))

  return(rp)
}

read.referral.diagnosis<-function(ocs, tables, occams_ids=NULL, rulesFile=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  tables = get.tables.by.crf(ocs, 'referral diagnosis')

  filter = make.id.filter(occams_ids,'RD_StudySubjectID')
  filter2 = make.id.filter(occams_ids, 'RDA_StudySubjectID')


  if (inherits(ocs, "OCCAMSLabkey")) {
    if (verbose) message(paste("Reading",tables))
    rd <- get.table.data(ocs, grep('^rd_', tables, value=T),  rulesFile=rulesFile, colFilter=filter)
    rd2 <- get.table.data(ocs, grep('^rda_', tables, value=T),  rulesFile=NULL, colFilter=filter2)
    if (nrow(rd2) > 0)
      warning("RDA addendum table has data, update required to table read.")

  } else {
    rd <- read.table.data(tables[grep('^rd_', basename(tables))],  rulesFile=rulesFile)
  }

  if(verbose) message(paste(nrow(rd), "patients"))
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
    select(StudySubjectID, starts_with('TC.'), starts_with('TC1.')) %>% group_by(StudySubjectID)

  if (verbose) message(paste(nrow(tc1), "tissues"))

  return(tc1)
}

