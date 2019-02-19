
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

  samplesTaken <- rownames(subset(tc,TC.TissueSamplesTaken == "Yes"))

  # useful table
  filter = make.id.filter(occams_ids, 'TC1_StudySubjectID')
  tc1 <- get.table.data(ocs=ocs, table=grep('^tc1_',tables,value=T), colFilter=filter)

  # But it seems the two tables may not actually match up! There are more patients in the second table then are listed in the first as having had tissue taken
  missingEntries <- setdiff(unique(tc1$TC1.StudySubjectID), samplesTaken)

  # Recode for clarity
  types <- c("endoscopy", "surgical","laparoscopy","EMR")
  names(types) <- c('E','S','L','R')
  tc1$TC1.TissueType <- revalue(tc1$TC1.TissueType, types)

  sources <- c("normal oesophagus", "barretts", "tumour", "normal gastric", "lymph node", "metastasis")
  names(sources) <- c('N','B','T','G','L','M')
  tc1$TC1.TissueSource <- revalue(tc1$TC1.TissueSource, sources)

  if (verbose) message(paste(nrow(tc1), "samples"))

  return(tc1)
}

make.id.filter<-function(occams_ids, idCol) {
  if (!is.null(occams_ids)) {
    occams_ids = grep('^OCCAMS/[A-Z]{2}/[0-9]+', occams_ids, value=T)
    return( makeFilter(c(idCol, 'IN', paste(occams_ids, collapse=';'))) )
  }
  return(NULL)
}

read.demographics<-function(ocs, table, occams_ids=NULL, rulesFile=NULL, ...) {
  table = list.clinical.tables(ocs,'di')
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  filter = make.id.filter(occams_ids,'DI_StudySubjectID')

  if (inherits(ocs, "OCCAMSLabkey")) {
    #stop("Create connection to OCCAMS Labkey first")
    if (verbose) message(paste("Reading",table))
    di <- get.table.data(ocs=ocs, table=table,  rulesFile=rulesFile, colFilter=filter)
  } else {
    di <- read.table.data(table,  rulesFile=rulesFile)
  }

  #di$DI.ageAtDiagnosis <- as.numeric(di$DI.ageAtDiagnosis)
  # Prefer to get this from FE tables

  di <- di %>% select(-DI.StudyNumber, -contains('Panther'), -contains('OpenClinica'), -contains('CRF'),-contains('DateOfDeath'))

  if (verbose) message(paste(nrow(di), "patients"))
  return(di)
}

read.exposures<-function(ocs, tables, occams_ids=NULL, rulesFiles=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  tables = list.clinical.tables(ocs,'ex')

  if (length(tables) != 3)
    stop("Three tables expected for exposures")

  filter = make.id.filter(occams_ids,'EX_StudySubjectID')
  filter1 = make.id.filter(occams_ids,'EX1_StudySubjectID')
  filter2 = make.id.filter(occams_ids,'EX2_StudySubjectID')

  if (inherits(ocs, "OCCAMSLabkey")) {
    #stop("Create connection to OCCAMS Labkey first")
    if (verbose) message(paste("Reading",paste(tables, collapse=',')))
    ex <- get.table.data(ocs=ocs, table=grep('^ex_',tables,value=T),  rulesFile=rulesFiles[1], colFilter=filter)
    ex_history_gastric <- get.table.data(ocs=ocs, table=grep('ex.*gastric',tables,value=T), rulesFile=rulesFiles[2], colFilter=filter1)
    ex_history_other <- get.table.data(ocs=ocs, table=grep('ex.*other_cancer',tables,value=T), rulesFile=rulesFiles[3], colFilter=filter2)

  } else {
    ex <- read.table.data(tables[1],  rulesFile=rulesFiles[1])
    ex_history_gastric <- read.table.data(tables[2], rulesFile=rulesFiles[2])
    ex_history_other <- read.table.data(tables[3], rulesFile=rulesFiles[3])
  }

  ex = ex %>% mutate_at(vars(grep('Height|Weight|BMI|Weeks|Months|Years|Days', colnames(.))), funs(as.numeric))

  # -1 is often used for NA in the numeric columns
  ex = ex %>% mutate_if( is.numeric, funs(ifelse(.<0,NA,.)))

  ex = ex %>% mutate(EX.CurrentWeightKG = ifelse(EX.CurrentWeightKG>140 | EX.CurrentWeightKG<22, NA, EX.CurrentWeightKG),
                EX.CurrentHeightCM = ifelse(EX.CurrentHeightCM>215 | EX.CurrentHeightCM<60, NA, EX.CurrentHeightCM),
                EX.CurrentBMI = bmicalc(EX.CurrentWeightKG,EX.CurrentHeightCM))

  if (verbose) message(paste(nrow(ex), "patients"))

  ex_history = full_join(ex_history_gastric, ex_history_other, by=c("EX1.StudySubjectID"="EX2.StudySubjectID")) %>% ungroup

  return(list('ex'=ex, 'history'=ex_history))
}

read.endpoints<-function(ocs, tables, occams_ids=NULL, rulesFiles=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  tables = list.clinical.tables(ocs,c('fe','ep'))

  if (length(tables) != 3)
    stop("Three tables expected for endpoints")

  filter = make.id.filter(occams_ids,'FE_StudySubjectID')
  filter1 = make.id.filter(occams_ids,'FE1_StudySubjectID')
  filter2 = make.id.filter(occams_ids, 'EP_StudySubjectID')

  if (inherits(ocs, "OCCAMSLabkey")) {
    #stop("Create connection to OCCAMS Labkey first")
    if (verbose) message(paste("Reading",paste(tables, collapse=",")))

    fe <- get.table.data(ocs,grep('fe_',tables,value=T),  rulesFile=rulesFiles[1], colFilter=filter)
    # EP appears to be set up to replace FE
    ep <- get.table.data(ocs,grep('ep',tables,value=T), rulesFile=grep('ep_',rulesFiles,value=T), colFilter=filter2)

    feep = full_join(fe,ep,by=c('FE.StudySubjectID'='EP.StudySubjectID')) %>%
      dplyr::group_by(FE.StudySubjectID) %>%
      rowwise %>% dplyr::mutate(
      EP.StudySite = mr(FE.StudySite, EP.StudySite),
      EP.EndPoint = mr(FE.EndPoint, EP.EndPoint),
      EP.DateOfPatientDeath = mr(FE.DateOfPatientDeath,EP.DateOfPatientDeath),
      EP.ReasonForPatientDeath = mr(FE.ReasonForPatientDeath, EP.ReasonForPatientDeath),
      EP.ReasonForPatientDeathOther = mr(FE.ReasonForPatientDeathOther, EP.ReasonForPatientDeathOther)
    ) %>% dplyr::select(FE.StudySubjectID, starts_with('EP'), -contains('OpenClinica'))

    fe2 <- get.table.data(ocs,grep('fe1',tables,value=T), rulesFile=grep('fe1_',rulesFiles,value=T), colFilter=filter1)

    fe2 = fe2 %>% dplyr::group_by(FE1.StudySubjectID) %>% arrange(desc(FE1.DateOfUpdate)) %>% filter(!grepl('death', FE1.ReasonForFollowUp)) %>%
      dplyr::summarise(FE1.DateOfUpdate=FE1.DateOfUpdate[1],
                       FE1.ReasonForFollowUp=FE1.ReasonForFollowUp[1],
                       FE1.HasOriginalDiseaseReoccurred=FE1.HasOriginalDiseaseReoccurred[1])
  } else {
    if (verbose) message(paste("Reading",paste(basename(tables), collapse=", ")))

    #fe <- read.table.data(tables[1],  rulesFile=rulesFiles[1])
    #fe2 <- read.table.data(tables[2], rulesFile=rulesFiles[2])
    #fe2 <- ddply( ddply(fe2, .(FE1.StudySubjectID), arrange, desc(FE1.DateOfUpdate)),
    #              .(FE1.StudySubjectID, FE1.StudySite), plyr::summarise, FE1.DateOfUpdate=FE1.DateOfUpdate[1], FE1.ReasonForFollowUp=FE1.ReasonForFollowUp[1], FE1.HasOriginalDiseaseReoccurred=FE1.HasOriginalDiseaseReoccurred[1])
  }

  # Final endpoint
  feF = left_join(feep, fe2, by=c('FE.StudySubjectID'='FE1.StudySubjectID')) %>% group_by(FE.StudySubjectID) %>% dplyr::rename_all( funs(sub('^FE1','FE',.)) )

  # withdrawn from study -- turns out we can still use the data we have jsut cannot track them further
  withdrawn = feF %>% filter(grepl('withdraw', EP.EndPoint))

  #withdrawn <- feF[rows,]
  #if (length(rows) > 0) feF <- feF[-rows,]

  feF = feF %>% rowwise %>% dplyr::mutate(
    EP.EndPoint = ifelse(!is.na(EP.ReasonForPatientDeath),'Patient died', EP.EndPoint),
    FE.Patient.Died = factor(ifelse(!is.na(EP.EndPoint) & EP.EndPoint == 'Patient died', 'yes','no')))

    # this column disappeared when I checked 2019/02
    #feF$FE.CancerFree.Discharge <- as.factor(with(feF, grepl("discharge", FE.AdditionalDetailsForPatientEndPoint, ignore.case=T) & (FE.EndPoint != 'Patient died' | is.na(FE.EndPoint))))

  if (verbose) message(paste(nrow(feF), "patients"))

  return(list('fe'=feF, 'withdrawn'=withdrawn$FE.StudySubjectID))
}

read.prestage<-function(ocs, tables, occams_ids=NULL, rulesFiles=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  tables = list.clinical.tables(ocs,'ps')

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

  tables = list.clinical.tables(ocs,'tp')

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

read.therapy<-function(ocs, tables, occams_ids=NULL, rulesFiles=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  tables = list.clinical.tables(ocs,'tr')

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

      # These tabels were added only recently and they were not backfilled so there are not many patients in each, can deal with them later when/if needed
      # radiotherapy -- later
      trr <- get.table.data(ocs,grep('^trr_', tables, value=T))

      # adjuvant chemo -- later
      trac <- get.table.data(ocs,grep('^trac_', tables, value=T))

      # other therapy -- later
      tro <- get.table.data(ocs,grep('^tro_', tables, value=T))

    } else {
      #tr <- read.table.data(tables[grep('^tr_', basename(tables))],  rulesFile=rulesFiles[1])
    }



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

  table = list.clinical.tables(ocs,'st')

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
  require(plyr)

  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  table = list.clinical.tables(ocs,'rp')

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

  diff.grade<-function(x) {
    if ( grepl('.*poor',x)  ) { # poor & moderate_to_poor
      return('poor')
    } else if ( grepl('.*well',x) ) { # well & moderate_to_well
      return('well')
    }
    return(x)
  }


 rp = rp %>% dplyr::mutate(
   RP.TNMStage.c = tnmStage(RP.TStage.PrimaryTumour, RP.Nstage.RP.TNM7, RP.MStage.DistantMetastasis),
   RP.TumourDifferentiation.c = diff.grade(RP.TumourGradingDifferentiationStatus)) %>% ungroup %>% mutate(
      RP.TumourResponse = recode_factor(RP.TumourResponse,'0pc_remaining'='0%','less_than_20pc'='<20%', 'greater_than_or_equals_to_20pc'='≥20%', 'less_than_50pc'='<50%', 'greater_than_or_equals_to_50pc'='≥50%', .ordered=T),
      RP.Location = factor(RP.Location, levels=c('oesophageal','goj','gastric')),
      RP.MandardScoreForResponse = factor(RP.MandardScoreForResponse, levels=c('TRG1','TRG2','TRG3','TRG4','TRG5')),
      # If you see it macro, then it is there micro as well
      RP.BarrettsAdjacentToTumour = ifelse(RP.BarettsAdjacentToTumourMicroscopicIM == 'yes' | RP.BarettsAdjacentToTumourMacroscopic == 'yes','yes','no')
  )

  file = system.file("extdata", "be_updates_20160930.txt", package="openclinica.occams")
  if (!exists('be_updates') & (!is.null(file) & file.exists(file))) {
    if (verbose) message(paste("Reading", file))
    be_updates = readr::read_tsv(file, col_types='ccc')

    if (verbose) message("Updating RP.BarrettsAdjacent patients from Caitrona's worksheet.")
    be_updates = be_updates %>% dplyr::mutate( `Barret's Confirmed`=ifelse(`Barret's Confirmed` == '?', NA, `Barret's Confirmed`),
                                               `Barret's Confirmed` = recode(`Barret's Confirmed`, 'N'='no','Y'='yes'))

    # Assume Caitrona's are the final say (they matched other than NA anyhow)
    rp = left_join(rp, be_updates, by=c('RP.StudySubjectID'='OCCAMS/ID')) %>% rowwise %>% dplyr::mutate(
      RP.BarrettsAdjacentToTumour = mr(`Barret's Confirmed`, RP.BarrettsAdjacentToTumour,1) ) %>% select(-Source, -contains('Confirmed'))
  }

  rp = rp %>% ungroup %>% mutate(RP.BarrettsAdjacentToTumour = factor(RP.BarrettsAdjacentToTumour)) %>% select(-contains('Nstage.RP.TNMSystem'),-contains('RP.TNM6')) %>% group_by(RP.StudySubjectID)

  if(verbose) message(paste(nrow(rp), "patients"))

  return(rp)
}

read.referral.diagnosis<-function(ocs, tables, occams_ids=NULL, rulesFile=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  tables = list.clinical.tables(ocs,'rd')

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

merge.patient.tables<-function(df1, df2, ...) {
  newDF <- merge.data.frame(df1, df2, by='row.names', ...)
  rownames(newDF) <- newDF[['Row.names']]
  newDF[['Row.names']] <- NULL
  return(newDF)
}

dump.clinical.data<-function(ocs, occams_ids = NULL, prefixes=NULL, version=NULL, outdir="/tmp") {
  out <- paste(outdir, Sys.Date(), sep="/")

  if (is.null(prefixes))
    prefixes = get.prefixes(ocs, version)

  message(paste("Writing Labkey table output to", out))

  if (!dir.exists(out))
    dir.create(out)

  tables = list.clinical.tables(ocs, prefixes)

  for(table in tables) {
    message(paste("Reading table", table))
    data = get.table.data(ocs, table, uniqueID = 'StudySubjectID')

    if (!is.null(occams_ids))
      data = data %>% filter( !!rlang::sym(group_vars(data)) %in% occams_ids)

    message(paste("Writing", nrow(data), "rows"))
    write.table(data, quote=T, sep="\t",row.names=F, file=paste(out,"/",table, ".tsv", sep=""))
  }

  message(paste("Completed dump for", length(tables), "tables"))
}

#' Download in wide format the Labkey data for all or a selected group of patients.
#' @name download.all.tables
#' @param ocs
#' @param occams_ids Array of occams identifiers, this is optional. The default returns all patients, better to use get.patients(...) if you want this.
#' @param missing Replace missing values with this string, default is NA
#' @param verbose
#'
#' @author
#' @export
download.all.tables<-function(ocs, occams_ids=NULL, missing=NULL,versions='z1', verbose=T) {
  #prefixes=c('di','rd','ex','ps','tp','tr','st','rp','fe|ep','tc')
  prefixes = get.prefixes(ocs,versions)
  if (inherits(ocs, "OCCAMSLabkey")) {
    tables = list.clinical.tables(ocs, prefixes)
  } else {
    tables = list.files(ocs, full.names=T)
  }

  ordered_tables = sapply(prefixes, function(x) {
    tables[sort(grep(paste("^",x,sep=""), basename(tables)))]
  })

  di <- read.demographics(ocs, ordered_tables$di, rulesFile=paste(path.package('openclinica.occams')[1], "editrules/di_editrules.txt", sep='/'),occams_ids=occams_ids,verbose=verbose) %>% dplyr::mutate(StudySubjectID = DI.StudySubjectID)

  exes <- read.exposures(ocs, ordered_tables$ex, rulesFiles=NULL, occams_ids=occams_ids, verbose=verbose)
  ex = exes$ex %>% dplyr::mutate(StudySubjectID = EX.StudySubjectID)
  history = exes$history

  feep <- read.endpoints(ocs, c(ordered_tables$fe, ordered_tables$ep),
                       rulesFiles=paste(path.package('openclinica.occams'), c('editrules/fe_editrules.txt', 'editrules/fe1_editrules.txt', 'editrules/ep_editrules.txt'), sep='/'), occams_ids=occams_ids, verbose=verbose)
  fe = feep$fe %>% dplyr::mutate(StudySubjectID = FE.StudySubjectID)
  withdrawn <- feep$withdrawn

  rd <- read.referral.diagnosis(ocs, ordered_tables$rd, rulesFile=paste(path.package('openclinica.occams'),'editrules/rd_editrules.txt',sep='/'), occams_ids=occams_ids, verbose=verbose) %>% dplyr::mutate(StudySubjectID = RD.StudySubjectID)

  # these two have duplicates as well
  ps <-  read.prestage(ocs, ordered_tables$ps, rulesFiles=NULL, occams_ids=occams_ids, verbose=verbose)$ps %>% dplyr::mutate(StudySubjectID = PS.StudySubjectID)
  tp <-  read.treatment.plan(ocs, ordered_tables$tp, rulesFiles=NULL, occams_ids=occams_ids, verbose=verbose)$tp %>% dplyr::mutate(StudySubjectID = TP.StudySubjectID)

  tr <- read.therapy(ocs, ordered_tables$tr, rulesFiles=NULL, occams_ids=occams_ids, verbose=verbose) %>% dplyr::mutate(StudySubjectID = TR.StudySubjectID)
  st <- read.surgery(ocs, ordered_tables$st, rulesFile=NULL, occams_ids=occams_ids, verbose=verbose) %>% dplyr::mutate(StudySubjectID = ST.StudySubjectID)

  rp <- read.pathology(ocs, ordered_tables$rp, occams_ids=occams_ids, rulesFile=paste(path.package('openclinica.occams'),'editrules/rp_editrules.txt',sep='/'), verbose=verbose) %>% dplyr::mutate(StudySubjectID = RP.StudySubjectID)

  if (verbose) message("Creating final patient table")
  # Merge final table

  all = dplyr::full_join(di,ex,by='StudySubjectID') %>%
    dplyr::full_join(rd,by='StudySubjectID') %>%
    dplyr::full_join(tp,by='StudySubjectID') %>%
    dplyr::full_join(ps,by='StudySubjectID') %>%
    dplyr::full_join(tr,by='StudySubjectID') %>%
    dplyr::full_join(st,by='StudySubjectID') %>%
    dplyr::full_join(rp,by='StudySubjectID') %>%
    dplyr::full_join(fe,by='StudySubjectID') %>%
    dplyr::rename(ID=StudySubjectID, StudySite=DI.StudySite) %>% dplyr::group_by(ID) %>% dplyr::select(-contains('StudySubjectID'), -contains('CRF'), -contains('OpenClinica')) %>%
    ungroup %>% recode.siewert %>% recode.TNM %>% group_by(ID) %>% select(ID,StudySite, everything(),-matches('\\.StudySite'))

  # rm = which(rownames(all) %in% withdrawn)
  # if (length(rm) > 0) all = all[-rm,]

  if (verbose) message(paste("Final patient total:", nrow(all)))

  all <- prestage.to.path(all)

  # Calculated columns
  #all %>% select(matches('\\.c$'))

  # Dates after diagnosis that patients may have been seen
  clinDates <- grep('(EP|FE|RP|ST|TR).+Date', colnames(all), value=T)

  # In order
  clinDates = all %>% select(ID,EP.DateOfPatientDeath, FE.DateOfUpdate, matches('ST.+Date'), matches('TR.NeoAdj.+Date'), matches('TR.Radio.+Date'), matches('TR.Endo.+Date')) %>% mutate_all(as.Date) %>% group_by(ID)

  # can't work this one out in dplyr...
  clinDates$FE.LastSeenDate <- as.Date(apply(clinDates[-1], 1, function(x) {
    y <- sort(x, decreasing=T)
    lastDate = y[1]

    if (!is.na(y['FE.DateOfPatientDeath']) & y['FE.DateOfPatientDeath'] != y[1])
      lastDate = y['FE.DateOfPatientDeath']

    return(lastDate)
  }))

  all = dplyr::left_join(all, clinDates %>% dplyr::select(ID,FE.LastSeenDate), by='ID')

  # Dates after diagnosis that patients may have been seen
  diagDates <- grep('RD.DateOfOGCDiagnosis|DI.DateOfInformedConsentForCAMSAndICGC', colnames(all), value=T)
  all$RD.DiagnosisDate <- as.Date(apply(all[diagDates], 1, function(x) {
    y <- sort(x)
    firstDate <- y[1]

    if (!is.na(y['RD.DateOfOGCDiagnosis']))
      firstDate <- y['RD.DateOfOGCDiagnosis']

    return(firstDate)
  }))

  all = all %>% dplyr::mutate(Weeks.Survival = round(as.numeric(difftime(FE.LastSeenDate, RD.DiagnosisDate, units='weeks'))), 3) %>% select(ID, StudySite, Weeks.Survival, RD.DiagnosisDate, FE.LastSeenDate, everything())

  if ( length(which(is.na(all$EP.DateOfPatientDeath))) == nrow(all) ) {
    warning("Date of death has been removed, survival time not calculated. ")
    all$Weeks.Survival = NA
  }

  bad = with(all, which(Weeks.Survival < 1))
  message(paste(length(bad), "patients have a diagnosis date before their last seen (death/surgery/etc) date."))
  all[bad, c('RD.DiagnosisDate', 'FE.LastSeenDate', 'Weeks.Survival')] = NA

  tc <- read.tissue.collection(ocs, ordered_tables$tc, rulesFiles=NULL, occams_ids=occams_ids)

  return(list('patients'=all, 'tissues'=tc, 'withdrawn'=withdrawn))
}




# TODO rewrite transformation for NA instead of some random character string
transform.patient.table<-function(fd) {
  cols = grep('Survival',colnames(fd))
  fd = fd[-cols]

  newdf = as.data.frame(matrix(nrow=nrow(fd),ncol=0,dimnames=list(rownames(fd), c())))

  #num = grep('age(\\s|_)?At|Height|Width|Weight|Length|Weeks|Number|Percentage|BMI|UnitsPer|Years|Months|Weeks|Days',colnames(fd),ignore.case=T,value=T)
  #fd[num] = lapply(fd[num], as.numeric)
  #  fd[grep('SurgeryPerformed$|Chemotherapy$|Radiotherapy$|Endoscopy$', colnames(fd), value=T)] =
  #    lapply(fd[grep('SurgeryPerformed$|Chemotherapy$|Radiotherapy$|Endoscopy$', colnames(fd), value=T)], as.factor)

  for (col in names(which(unlist(lapply(fd[], is.factor))))) {
    if (  length(which(levels(fd[[col]]) == c("0","1"))) == 2  ) { # Binaries
      newdf[[ paste(col,'yes',sep=".") ]] = sapply(as.integer((fd[[col]] == 1 )), function(x) ifelse(is.na(x), 0, x))
      newdf[[ paste(col,'no',sep=".") ]] = sapply(as.integer((fd[[col]] == 0)), function(x) ifelse(is.na(x), 0, x))
      newdf[[ paste(col,'unknown',sep=".") ]] = as.integer((fd[[col]] == 'unknown' | is.na(fd[[col]])))
    } else {
      for (lev in levels(fd[[col]])) {
        if (grepl('unknown',lev,ignore.case=T)) {
          newdf[[ paste(col,lev,sep=".") ]] = as.integer((fd[[col]] == lev | is.na(fd[[col]])))
        } else {
          newdf[[ paste(col,lev,sep=".") ]] = as.integer((fd[[col]] == lev & !is.na(fd[[col]])))
        }
      }
    }
  }

  for (col in names(which(unlist(lapply(fd[], is.numeric))))) {
    #if (grepl('Log.Survival|Weeks.Survival|Chemotherapy|SurgeryPerformed|Endoscopy|Radiotherapy', col)) {
    #newdf[[col]] = as.numeric(fd[[col]])
    #next
    #}
    knowns = which(fd[[col]] > 0)
    if (length(knowns) <= length(fd[[col]])/2) { # Create a binary, unknown and > 1 maybe?
      newdf[[ paste(col, 'unknown',sep=".") ]] = as.numeric(fd[[col]] <= 0 | is.na(fd[[col]]))
      newdf[[ paste(col, '>=1',sep=".") ]] = sapply(as.numeric(fd[[col]] > 0), function(x) ifelse(is.na(x), 0,x))
    } else {
      p = 11
      quant = quantile(fd[knowns,col],probs=seq(0,1,length=p))
      while (length(which(quant == 0)) > 1) {
        p = p - 1
        quant = quantile(fd[knowns,col],probs=seq(0,1,length=p))
      }
      for (i in 2:length(quant))
        newdf[[ paste(col, names(quant[i]),sep=".") ]] = as.numeric(fd[[col]] <= quant[i] & fd[[col]] > quant[i-1] & !is.na(fd[[col]]) )
    }
  }

  if (length(which(apply(newdf,2,lenNA) > 0)) > 0)
    stop("NA's have been introduced in transformation")
  message(ncol(newdf))
  return(newdf)
}




