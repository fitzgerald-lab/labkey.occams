
read.tissue.collection<-function(ocs, tables, rulesFiles=NULL, ...) {
  if (!inherits(ocs, "OCCAMSLabkey"))
    stop("Create connection to OCCAMS Labkey first")
  if (length(tables) < 2)
    stop("Two tables are expected for tissue collection")
  message(paste("Reading",paste(tables, collapse=",")))

  # all this table tells me is if there should be entries in the second table or not
  tc <- get.table.data(ocs=ocs, table=grep('tc_',tables,value=T), ...)

  # There should not be duplicates in this table...
  dups <- which(with(tc, TC.StudySubjectID %in% names(which(table(TC.StudySubjectID) > 1))))
  duplicatePts <- tc[dups,]

  warning(paste("Removing", nrow(duplicatePts), "duplicated patient rows"))
  tc <- tc[-dups, ]
  tc <- rows.as.patient.id(tc, 2)

  samplesTaken <- rownames(subset(tc,TC.TissueSamplesTaken == "Yes"))

  # useful table
  tc1 <- get.table.data(ocs=ocs, table=grep('tc1_',tables,value=T))

  # But it seems the two tables may not actually match up! There are more patients in the second table then are listed in the first as having had tissue taken
  missingEntries <- setdiff(unique(tc1$TC1.StudySubjectID), samplesTaken)

  # Recode for clarity
  types <- c("endoscopy", "surgical","laparoscopy","EMR")
  names(types) <- c('E','S','L','R')
  tc1$TC1.TissueType <- revalue(tc1$TC1.TissueType, types)

  sources <- c("normal oesophagus", "barretts", "tumour", "normal gastric", "lymph node", "metastasis")
  names(sources) <- c('N','B','T','G','L','M')
  tc1$TC1.TissueSource <- revalue(tc1$TC1.TissueSource, sources)

  message(paste(nrow(tc1), "samples"))

  return(tc1)
}

read.demographics<-function(ocs, table, rulesFile=NULL, ...) {
  if (inherits(ocs, "OCCAMSLabkey")) {
    #stop("Create connection to OCCAMS Labkey first")
    message(paste("Reading",table))
    di <- get.table.data(ocs=ocs, table=table, uniqueID=2, rulesFile=rulesFile, ...)
  } else {
    di <- read.table.data(table, uniqueID=2, rulesFile=rulesFile)
  }

  #di$DI.ageAtDiagnosis <- as.numeric(di$DI.ageAtDiagnosis)
  # Prefer to get this from FE tables
  di$DI.PatientDateOfDeath = NULL

  message(paste(nrow(di), "patients"))
  return(di)
}

read.exposures<-function(ocs, tables, rulesFiles=NULL, ...) {
  if (length(tables) != 3)
    stop("Three tables expected for exposures")

  if (inherits(ocs, "OCCAMSLabkey")) {
    #stop("Create connection to OCCAMS Labkey first")
    message(paste("Reading",paste(tables, collapse=',')))
    ex <- get.table.data(ocs=ocs, table=grep('^ex_',tables,value=T), uniqueID=2, rulesFile=rulesFiles[1], ...)
    ex_history_gastric <- get.table.data(ocs=ocs, table=grep('ex.*gastric',tables,value=T), rulesFile=rulesFiles[2], ...)
    ex_history_other <- get.table.data(ocs=ocs, table=grep('ex.*other_cancer',tables,value=T), rulesFile=rulesFiles[3], ...)

  } else {
    ex <- read.table.data(tables[1], uniqueID=2, rulesFile=rulesFiles[1])
    ex_history_gastric <- read.table.data(tables[2], rulesFile=rulesFiles[2])
    ex_history_other <- read.table.data(tables[3], rulesFile=rulesFiles[3])
  }

  cols <- grep('Height|Weight|BMI|Weeks|Months|Years|Days', colnames(ex), value=T)
  ex[cols] <- lapply(ex[cols], as.numeric)

  # -1 is often used for NA in the numeric columns
  ex[cols] <- lapply(ex[cols], function(x) { ifelse(x < 0, NA, x) })

  ex$EX.CurrentBMI <- apply(ex[c('EX.CurrentBMI','EX.CurrentHeightCM', 'EX.CurrentWeightKG')], 1, function(x) {
    bmi <- x[1]
    if (is.na(x[1]) & !is.na(x[2]) & !is.na(x[3]))
      bmi <- bmicalc(x[2:3])
    return(bmi)
  })

  ex$EX.FamilyHistory.OGC.Relationship <- NA
  ex[unique(subset(ex_history_gastric, EX1.FamilyHistoryOfOesophagoGastricCancerRelationship %in% c('Brother', 'Sister'), select = EX1.StudySubjectID)[,1]), 'EX.FamilyHistory.OGC.Relationship'] <- "Sibling"

  ex[unique(subset(ex_history_gastric, EX1.FamilyHistoryOfOesophagoGastricCancerRelationship %in% c('Mother', 'Father'), select = EX1.StudySubjectID)[,1]), 'EX.FamilyHistory.OGC.Relationship'] <- "Parent"

  ex[unique( subset(ex_history_gastric, grepl('Uncle|Aunt|Grand', EX1.FamilyHistoryOfOesophagoGastricCancerRelationship), select = EX1.StudySubjectID)[,1]), 'EX.FamilyHistory.OGC.Relationship'] <- "Au.Unc.Grand"

  ex[unique(subset(ex_history_gastric, grepl('Cousin|Niece|Nephew', EX1.FamilyHistoryOfOesophagoGastricCancerRelationship), select = EX1.StudySubjectID)[,1]), 'EX.FamilyHistory.OGC.Relationship'] <- "Cousin"

  ex$EX.FamilyHistory.OGC.Relationship = ordered(ex$EX.FamilyHistory.OGC.Relationship, levels=c('Sibling','Parent','Au.Unc.Grand','Cousin'))

  message(paste(nrow(ex), "patients"))

  return(ex)
}

read.endpoints<-function(ocs, tables, rulesFiles=NULL, ...) {
  require(plyr)

  if (length(tables) != 2)
    stop("Two tables expected for endpoints")

  if (inherits(ocs, "OCCAMSLabkey")) {
    #stop("Create connection to OCCAMS Labkey first")
    message(paste("Reading",paste(tables, collapse=",")))

    fe <- get.table.data(ocs,grep('fe_',tables,value=T), uniqueID=2, rulesFile=rulesFiles[1], ...)
    fe2 <- get.table.data(ocs,grep('fe1',tables,value=T), rulesFile=rulesFiles[2], ...)
    fe2 <- ddply( ddply(fe2, .(FE1.StudySubjectID), arrange, desc(FE1.DateOfUpdate)),
                  .(FE1.StudySubjectID, FE1.StudySite), plyr::summarise, FE1.DateOfUpdate=FE1.DateOfUpdate[1], FE1.ReasonForFollowUp=FE1.ReasonForFollowUp[1], FE1.HasOriginalDiseaseReoccurred=FE1.HasOriginalDiseaseReoccurred[1])
  } else {
    message(paste("Reading",paste(basename(tables), collapse=", ")))

    #fe <- read.table.data(tables[1], uniqueID=2, rulesFile=rulesFiles[1])
    #fe2 <- read.table.data(tables[2], rulesFile=rulesFiles[2])
    #fe2 <- ddply( ddply(fe2, .(FE1.StudySubjectID), arrange, desc(FE1.DateOfUpdate)),
    #              .(FE1.StudySubjectID, FE1.StudySite), plyr::summarise, FE1.DateOfUpdate=FE1.DateOfUpdate[1], FE1.ReasonForFollowUp=FE1.ReasonForFollowUp[1], FE1.HasOriginalDiseaseReoccurred=FE1.HasOriginalDiseaseReoccurred[1])
  }

  # Final endpoint
  rownames(fe2) = fe2$FE1.StudySubjectID

  feF <- merge.patient.tables(fe,fe2, all=T)
  #feF <- merge(fe,fe2, by='row.names', all=T)

  # withdrawn from study
  rows <- which(with(feF, grepl("withdrawal", FE.EndPoint)))
  withdrawn <- feF[rows,]
  feF <- feF[-rows,]

  feF[which(with(feF, FE.EndPoint != 'Patient died' & !is.na(FE.ReasonForPatientDeath))), 'FE.EndPoint'] <- 'Patient died'

  feF$FE.CancerFree.Discharge <- as.factor(with(feF, grepl("discharge", FE.AdditionalDetailsForPatientEndPoint, ignore.case=T) & (FE.EndPoint != 'Patient died' | is.na(FE.EndPoint))))

  feF$FE.Patient.Died <- as.factor(ifelse(feF$FE.EndPoint == 'Patient died' & !is.na(feF$FE.EndPoint), "yes","no"))

  message(paste(nrow(feF), "patients"))

  return(list('fe'=feF, 'withdrawn'=rownames(withdrawn)))
}

read.prestage<-function(ocs, tables, rulesFiles=NULL, ...) {
  if (length(tables) != 7)
    warning("Seven tables are expected for prestaging, currently only the final tables is read.")

  if (inherits(ocs, "OCCAMSLabkey")) {
    #stop("Create connection to OCCAMS Labkey first")
    message(paste("Reading",paste(tables, collapse=",")))
    ps <- get.table.data(ocs,grep('^ps_',tables, value=T), rulesFile=rulesFiles[1])
  } else {
    ps <- read.table.data(tables[grep('^ps_',basename(tables))], rulesFile=rulesFiles[1])
  }

  # There should not be duplicates in this table...
  dups <- which(with(ps, PS.StudySubjectID %in% names(which(table(PS.StudySubjectID) > 1))))
  duplicatePts <- ps[dups,]

  ps <- ps[-dups, ]
  ps <- rows.as.patient.id(ps, 2)

  ps <- fix.tumor.factors(text.to.tstage(ps))

  message(paste(nrow(ps), "patients"))

  return(list('ps' = ps, 'duplicates' = duplicatePts))
}

read.treatment.plan<-function(ocs, tables, rulesFiles=NULL, ...) {
  if (length(tables) < 2)
    stop("Two tables are expected for treatment plans, currently using only primary table.")

  if (inherits(ocs, "OCCAMSLabkey")) {
    message(paste("Reading",paste(tables, collapse=",")))
    tp <- get.table.data(ocs=ocs, table=grep('tp_',tables,value=T), rulesFile=rulesFiles[1], ...)
    tp1 <- get.table.data(ocs=ocs, table=grep('tp1_',tables,value=T), rulesFile=rulesFiles[2], ...)
  } else {
    #tp <- read.table.data(tables[1], rulesFile=rulesFiles[1])
    #tp1 <- read.table.data(tables[2], rulesFile=rulesFiles[2])
  }
  # There should not be duplicates in this table...
  dups <- which(with(tp, TP.StudySubjectID %in% names(which(table(TP.StudySubjectID) > 1))))
  duplicatePts <- tp[dups,]

  tp <- tp[-dups,]
  tp <- rows.as.patient.id(tp,2)

  message(paste(nrow(tp), "patients"))

  return(list('tp'=tp,'duplicates'=duplicatePts))
}

read.therapy<-function(ocs, tables, rulesFiles=NULL, ...) {
  if (length(tables) < 7)
    warning("Seven tables are expected for therapy, currently only the primary table is read.")

    if (inherits(ocs, "OCCAMSLabkey")) {
      message(paste("Reading",paste(tables, collapse=",")))
      tr <- get.table.data(ocs,grep('^tr_', tables, value=T), uniqueID=2, rulesFile=rulesFiles[1])
    } else {
      tr <- read.table.data(tables[grep('^tr_', basename(tables))], uniqueID=2, rulesFile=rulesFiles[1])
    }

  dups <- which(with(tr, TR.StudySubjectID %in% names(which(table(TR.StudySubjectID) > 1))))

  # Currently there's nothing in the other tables
  #trc <- get.table.data(ocs,tables[3], rulesFile=rulesFiles[3])

  cols = c('TR.CurativeTreatmentModality','TR.PalliativeTreatmentModality')
  sapply(tr[cols], function(x) grep('chemo',x, value=T, ignore.case=T))

  tr$TR.Chemotherapy = 'no'
  tr[with(tr, grep('chemo', TR.CurativeTreatmentModality, ignore.case=T) ), 'TR.Chemotherapy'] =  'yes'
  tr[with(tr, grep('chemo', TR.PalliativeTreatmentModality, ignore.case=T) ), 'TR.Chemotherapy'] =  'yes'

  tr$TR.Radiotherapy = 'no'
  tr[with(tr, grep('radio', TR.CurativeTreatmentModality, ignore.case=T) ), 'TR.Radiotherapy'] =  'yes'
  tr[with(tr, grep('radio', TR.PalliativeTreatmentModality, ignore.case=T) ), 'TR.Radiotherapy'] =  'yes'

  tr$TR.EndoscopicTherapy = 'no'
  tr[with(tr, grep('endoscop[y|ic]', TR.CurativeTreatmentModality, ignore.case=T) ), 'TR.EndoscopicTherapy'] =  'yes'
  tr[with(tr, grep('endoscop[y|ic]', TR.PalliativeTreatmentModality, ignore.case=T) ), 'TR.EndoscopicTherapy'] =  'yes'


  tr[grep('Regimen',colnames(tr), value=T)] = lapply(tr[grep('Regimen',colnames(tr), value=T)], as.character)

  revalue.reg<-function(x) {
    reg=c(1,0)
    names(reg) = c('yes','no')
    as.integer(revalue(x, reg))
  }

  reg.char=c('yes','no')
  names(reg.char) = c(1,0)

  ecf = sapply(tr[grep('5FU|Cisplatin|Epirubicin', colnames(tr), value=T)], revalue.reg)
  tr$TR.Regimen.ECF = revalue(as.factor(as.integer(apply(ecf, 1, function(x) sum(x) == 3))), reg.char)

  cf = sapply(tr[grep('5FU','Cisplatin', colnames(tr), value=T)], revalue.reg)
  tr$TR.Regimen.CF = revalue(as.factor(as.integer(apply(ecf, 1, function(x) sum(x) == 2))), reg.char)

  ecx = sapply(tr[grep('Capecitabine|Cisplatin|Epirubicin', colnames(tr), value=T)], revalue.reg)
  tr$TR.Regimen.ECX = revalue(as.factor(as.integer(apply(ecf, 1, function(x) sum(x) == 3))), reg.char)

  eox = sapply(tr[grep('Capecitabine|Oxaliplatin|Epirubicin', colnames(tr), value=T)], revalue.reg)
  tr$TR.Regimen.EOX = revalue(as.factor(as.integer(apply(ecf, 1, function(x) sum(x) == 3))), reg.char)

  cc = sapply(tr[grep('Capecitabine|Cisplatin', colnames(tr), value=T)], revalue.reg)
  tr$TR.Regimen.CC = revalue(as.factor(as.integer(apply(ecf, 1, function(x) sum(x) == 2))), reg.char)

  co = sapply(tr[grep('Capecitabine|Oxaliplatin', colnames(tr), value=T)], revalue.reg)
  tr$TR.Regimen.CO = revalue(as.factor(as.integer(apply(ecf, 1, function(x) sum(x) == 2))), reg.char)

  message(paste(nrow(tr), "patients"))

  return(tr)
}

read.surgery<-function(ocs, table, rulesFile=NULL, ...) {
  if (inherits(ocs, "OCCAMSLabkey")) {
    message(paste("Reading",table))
    st <- get.table.data(ocs,table, rulesFile=rulesFile)
  } else {
    st <- read.table.data(table, rulesFile=rulesFile)
  }
  # There should not be duplicates in this table...
  dups <- which(duplicated(st[,2]))
  duplicatePts <- st[dups,]

  st <- st[-dups, ]
  st <- rows.as.patient.id(st, 2)

  st$ST.SurgeryPerformed <- with(st, ifelse (ST.MainSurgery == 'yes' & (is.na(ST.ReasonIfNoSurgery) | is.na(ST.OtherReasonForNoSurgery)), 'yes', 'no'))
  st$ST.PathologyReportGenerated <- ifelse(with(st,
                                               ST.SurgeryPerformed == 'yes' & !grepl("Open (and|&) Shut", ST.Procedure, ignore.case=T)),
                                          'yes', 'no')

  st <- st[grep("\\.[1-9]$",colnames(st), invert=T, value=T)]

  cols <- grep("ASAGrade|Percentage", colnames(st), value=T)
  st[cols] <- lapply(st[cols], as.numeric)

  message(paste(nrow(st), "patients"))

  return(st)
}

read.pathology<-function(ocs, table, rulesFile=NULL, ...) {
  require(plyr)

  if (inherits(ocs, "OCCAMSLabkey")) {
    message(paste("Reading",table))
    rp <- get.table.data(ocs=ocs, table=table, uniqueID=NULL, rulesFile=rulesFile, ...)
  } else {
    rp <- read.table.data(table, uniqueID=NULL, rulesFile=rulesFile)
  }

  dups <- which(duplicated(rp[,2]))
  duplicatePts <- rp[dups,]

  rp <- rp[-dups, ]

  cols <- grep("Siewert|NumberOf|width|length", colnames(rp), value=T)
  rp[cols] <- lapply(rp[cols], as.numeric)

  rp <- rows.as.patient.id(rp, 2)

  rp <- fix.tumor.factors(text.to.tstage(rp))

  # Order some of the factors
  rp$RP.TumourResponse <- ordered(as.factor(rp$RP.TumourResponse),
                                  levels=c('0pc_remaining','less_than_20pc','greater_than_or_equals_to_20pc', 'less_than_50pc','greater_than_or_equals_to_50pc'))

  rp$RP.SiewertClassification <- ordered(rp$RP.SiewertClassification, levels=c(1,2,3))
  rp$RP.Location <- ordered(rp$RP.Location, levels=c('oesophageal','goj','gastric'))
  rp$RP.MandardScoreForResponse <- ordered(rp$RP.MandardScoreForResponse, levels=names(table(rp$RP.MandardScoreForResponse)))

  rp$RP.TumourDifferentiation <- sapply(rp$RP.TumourGradingDifferentiationStatus, function(x) {
    if ( grepl('.*poor',x)  ) { # poor & moderate_to_poor
      return('poor')
    } else if ( grepl('.*well',x) ) { # well & moderate_to_well
      return('well')
    }
    return(x)
  })

  reg <- c('yes','no')
  names(reg) <- c(1,0)

  # If you see it macro, then it is there micro as well
  rp$RP.BarrettsAdjacent = 0
  rp[which(as.integer(with(rp, RP.BarettsAdjacentToTumourMicroscopicIM == 'yes' | RP.BarettsAdjacentToTumourMacroscopic == 'yes')) == 1), 'RP.BarrettsAdjacent'] = 1

  message("Updating Caitron's patients for RP.BarrettsAdjacent")
  be_updates = readr::read_tsv(system.file("extdata", "be_updates_20160930.txt", package="openclinica.occams"))
  be_updates = be_updates %>% dplyr::mutate( `Barret's Confirmed`=ifelse(`Barret's Confirmed` == '?', NA, `Barret's Confirmed`) )
  be_updates$`Barret's Confirmed` = as.integer(as.factor(be_updates$`Barret's Confirmed`))-1

  # Assume Caitron's are the final say (they matched other than NA anyhow)
  rp[be_updates$`OCCAMS/ID`,'RP.BarrettsAdjacent'] = be_updates$`Barret's Confirmed`

  rp$RP.BarrettsAdjacent <- revalue(as.factor(rp$RP.BarrettsAdjacent), reg)

  cols = c('RP.Nstage.RP.TNMSystem', 'RP.Nstage.RP.TNM6')
  rp[cols] = NULL

  message(paste(nrow(rp), "patients"))

  return(rp)
}

read.referral.diagnosis<-function(ocs, tables, rulesFile=NULL, ...) {
  if (inherits(ocs, "OCCAMSLabkey")) {
    message(paste("Reading",tables))
    rd <- get.table.data(ocs,grep('^rd_', tables, value=T), uniqueID=2, rulesFile=rulesFile, ...)
  } else {
    rd <- read.table.data(tables[grep('^rd_', basename(tables))], uniqueID=2, rulesFile=rulesFile)
  }
  rd$RD.SiewertClassification <- ordered(rd$RD.SiewertClassification, levels=c(1,2,3))

  message(paste(nrow(rd), "patients"))
  return(rd)
}

fix.tumor.factors<-function(fd) {
  # Need to fix the factors for tumor staging to ensure that they are essentially in size order 1 > 2 > 3 etc...
  orderedTstages = c('Tx','T0','Tis','T1','T1a','T1b','T2','T3','T4','T4a','T4b')
  orderedNstages = c('Nx','N0','N1','N2','N3')
  orderedMstages = c('Mx','M0','M1')

  tcols = grep('(Tstage|.+_T$)', colnames(fd), ignore.case=T,value=T)
  fd[tcols] = lapply(fd[tcols], ordered, levels=orderedTstages)

  ncols = grep('(Nstage.+TNM7|.+_N$)', colnames(fd), ignore.case=T,value=T)
  fd[ncols] = lapply(fd[ncols], ordered, levels=orderedNstages)

  mcols = grep('(Mstage|.+_M$)', colnames(fd), ignore.case=T,value=T)
  fd[mcols] = lapply(fd[mcols], ordered, levels=orderedMstages)

  return(fd)
}

text.to.tstage<-function(fd) {
  orderedTstages = c(NA,'Tx','T0','Tis','T1','T1a','T1b','T2','T3','T4','T4a','T4b')
  names(orderedTstages) = c('not_recorded',
                            'primary_tumour_cannot_be_assessed',
                            'no_evidence_of_primary_tumour',
                            'carcinoma_in_situ_high_grade_dysplasia',
                            'tumour_invades_lamina_propria_muscularis_mucosae_or_submucosa',
                            'tumour_invades_lamina_propria_or_muscularis_mucosae',
                            'tumour_invades_submucosa',
                            'tumour_invades_muscularis_propria',
                            'tumour_invades_adventitia',
                            'tumour_invades_adjacent_structures',
                            'tumour_invades_pleura_pericardium_or_diaphragm',
                            'tumour_invades_other_adjacent_structures_such_as_aorta_vertebral_body_or_trachea')

  tcols = grep('(Tstage|.+_T$)', colnames(fd), ignore.case=T,value=T)
  fd[tcols] = lapply(fd[tcols], ordered, levels=names(orderedTstages))

  fd[tcols] = lapply(fd[tcols], revalue, orderedTstages)

  return(fd)
}

merge.patient.tables<-function(df1, df2, ...) {
  newDF <- merge.data.frame(df1, df2, by='row.names', ...)
  rownames(newDF) <- newDF[['Row.names']]
  newDF[['Row.names']] <- NULL
  return(newDF)
}

bmicalc<-function(x) {
  wt <- x[1]
  ht <- x[2]

  if (is.na(wt) || is.na(ht) || wt <= 0 || ht <= 0)
    return(NA)
  return(wt/(ht/100))/(ht/100)
}

dump.clinical.data<-function(ocs, prefixes=NULL, version=NULL, outdir="/tmp") {
  out <- paste(outdir, Sys.Date(), sep="/")

  if (is.null(prefixes))
    prefixes = get.prefixes(ocs, version)

  message(paste("Writing Labkey table output to", out))

  if (!dir.exists(out))
    dir.create(out)

  tables = list.clinical.tables(ocs, prefixes)

  for(table in tables) {
    message(paste("Reading table", table))
    data = get.table.data(ocs, table)
    message(paste("Writing", nrow(data), "rows"))
    write.table(data, quote=T, sep="\t",row.names=F, file=paste(out,"/",table, ".tsv", sep=""))
  }

  message(paste("Completed dump for", length(tables), "tables"))
}

clean.missing<-function(df, missing) {
  if (!is.null(missing)) {
    cols = which(sapply(df[], is.character))
    #cols = which(sapply(df[], function(x) is.character(x) | is.factor(x) )) # maybe...
    dfT = df[cols]
    dfT[is.na(dfT)] = missing
    df[cols] = dfT
  }
  return(df)
}

download.all.tables<-function(ocs, prefixes=c('di','rd','ex','ps','tp','tr','st','rp','fe','tc'), missing=NULL) {
  require(plyr)

  if (inherits(ocs, "OCCAMSLabkey")) {
    tables = list.clinical.tables(ocs, prefixes)
  } else {
    tables = list.files(ocs, full.names=T)
  }

  ordered_tables = sapply(prefixes, function(x) {
    tables[sort(grep(paste("^",x,sep=""), basename(tables)))]
  })

  di <- clean.missing(
    read.demographics(ocs, ordered_tables$di, NULL),NULL )


  di <- clean.missing(
      read.demographics(ocs, ordered_tables$di, rulesFile=paste(path.package('openclinica.occams')[1], "editrules/di_editrules.txt", sep='/')),
    missing )

  ex <- clean.missing( read.exposures(ocs, ordered_tables$ex, rulesFiles=NULL), missing )
  fe <- read.endpoints(ocs, ordered_tables$fe, rulesFiles=paste(path.package('openclinica.occams')[1],c("editrules/fe_editrules.txt", "editrules/fe1_editrules.txt"), sep='/'))
  withdrawn <- fe$withdrawn
  fe <- clean.missing(fe$fe, missing)

  rd <- clean.missing( read.referral.diagnosis(ocs, ordered_tables$rd, rulesFile=NULL), missing )

  # these two have duplicates as well
  ps <- clean.missing( read.prestage(ocs, ordered_tables$ps, rulesFiles=NULL)$ps, missing )
  tp <- clean.missing( read.treatment.plan(ocs, ordered_tables$tp, rulesFiles=NULL)$tp, missing )

  tr <- clean.missing( read.therapy(ocs, ordered_tables$tr, rulesFiles=NULL), missing )
  st <- read.surgery(ocs, ordered_tables$st, rulesFile=NULL)

  rp <- clean.missing( read.pathology(ocs, ordered_tables$rp, rulesFile=NULL), missing)

  message("Creating final patient table")
  # Merge final table
  all <- merge.patient.tables(di, fe, all=T)
  all <- merge.patient.tables(all, tp, all=T)
  all <- merge.patient.tables(all, rd, all=T)
  all <- merge.patient.tables(all, ex, all=T)
  all <- merge.patient.tables(all, ps, all=T)
  all <- merge.patient.tables(all, tr, all=T)
  all <- merge.patient.tables(all, st, all=T)
  all <- merge.patient.tables(all, rp, all=T)

  all = all[-which(rownames(all) %in% withdrawn),]

  message(paste("Final patient total:", nrow(all)))

  all <- prestage.to.path(all)

  all[grep('CRF\\.(Name|Version)', colnames(all), value=T)] <- NULL
  all[grep('OpenClinicaStartDate|OpenClinicaInterviewDate', colnames(all), value=T)] <- NULL

  all$StudySubjectID = rownames(all)
  all$StudySite = all$DI.StudySite
  all[grep("[A-Z]{2}\\.StudySubjectID", colnames(all), value=T)] <- NULL
  all[grep("[A-Z]{2}\\.StudySite", colnames(all), value=T)] <- NULL

  # Dates after diagnosis that patients may have been seen
  clinDates <- grep('(FE|RP|ST|TR).+Date', colnames(all), value=T)
  all$FE.LastSeenDate <- as.Date(apply(all[clinDates], 1, function(x) {
    y <- sort(x, decreasing=T)
    lastDate = y[1]

    if (!is.na(y['FE.DateOfPatientDeath']) & y['FE.DateOfPatientDeath'] != y[1])
      lastDate = y['FE.DateOfPatientDeath']

    return(lastDate)
  }))

  # Dates after diagnosis that patients may have been seen
  diagDates <- grep('RD.DateOfOGCDiagnosis|DI.DateOfInformedConsentForCAMSAndICGC', colnames(all), value=T)
  all$RD.DiagnosisDate <- as.Date(apply(all[diagDates], 1, function(x) {
    y <- sort(x)
    firstDate <- y[1]

    if (!is.na(y['RD.DateOfOGCDiagnosis']))
      firstDate <- y['RD.DateOfOGCDiagnosis']

    return(firstDate)
  }))

  all$Weeks.Survival = round(as.numeric(with(all, difftime(FE.LastSeenDate, RD.DiagnosisDate, units='weeks'))), 3)

  bad = with(all, which(Weeks.Survival < 1))
  message(paste(length(bad), "patients have a diagnosis date before their last seen (death/surgery/etc) date."))
  all[bad, c('RD.DiagnosisDate', 'FE.LastSeenDate', 'Weeks.Survival')] = NA

  tc <- read.tissue.collection(ocs, ordered_tables$tc, rulesFiles=NULL)

  return(list('patients'=all, 'tissues'=tc))
}

prestage.to.path<-function(df) {
  # -1 to set 'unknown' to 0
  RP.TumourStage.Int <- as.integer(df$RP.TStage.PrimaryTumour)-1
  PS.TumourStage.Int <- as.integer(df$PS.TStage.PrimaryTumour.FinalPretreatmentStaging)-1

  # 0 == no change, negative #s indicate growth, positive indicate shrinkage
  df$RP.PS.Tumour.Diff <- PS.TumourStage.Int - RP.TumourStage.Int
  df$RP.Tumour.Growth <- sapply(df$RP.PS.Tumour.Diff, function(x) {
    if (is.na(x))
      return(NA)
    if (x > 0) {
      x <- 'shrink'
    } else if (x < 0) {
      x <- 'grow'
    } else if (x == 0) {
      x <- 'no change'
    }
    return(x)
  })
  return(df)
}

# Currently not in use due to possible introduced biases
remove.bad.dates<-function(df) {
  #Diagnosis
  fixRows = which(with(df, is.na(RD.DateOfOGCDiagnosis)))
  message(paste(length(fixRows), "patients missing diagnosis date, fixing to consent date"))
  df[fixRows, 'RD.DateOfOGCDiagnosis'] = df[fixRows, 'DI.DateOfInformedConsentForCAMSAndICGC']
  stillNa = which(with(df, is.na(RD.DateOfOGCDiagnosis)))

  #Chemo start
  noChemoStart = which(with(df, TR.Chemotherapy == 'yes' & is.na(TR.ChemotherapyDateFirstCycleStarted)))
  message(paste(length(noChemoStart), "patients given chemo have no treatment start date."))

  # Chemo date before diagnosis date, don't think I can fix these
  chemoBeforeDiag = which(with(df, TR.Chemotherapy == 'yes' & RD.DateOfOGCDiagnosis > TR.ChemotherapyDateFirstCycleStarted))
  message(paste(length(chemoBeforeDiag), "patients with diagnosis date after first chemo treatment"))

  # Operation
  opNoDate = which(with(df, ST.PathologyReportGenerated == 'yes' & is.na(ST.DateOfOperation)))
  message(paste(length(opNoDate), "patients had operation but no date entered"))

  opPreDiag = which(with(df, ST.PathologyReportGenerated == 'yes' & RD.DateOfOGCDiagnosis > ST.DateOfOperation))
  message(paste(length(opPreDiag), "patients with diagnosis date after operation date"))

  # Diagnosis after death
  diagPostDeath = which(with(df, FE.Patient.Died == 'yes' & RD.DateOfOGCDiagnosis > FE.DateOfPatientDeath))
  message(paste(length(diagPostDeath), "patients with diagnosis date after death date"))

  # Diagnosis after last seen date
  diagAfterLastSeen = which(with(df, RD.DateOfOGCDiagnosis > FE.LastSeenDate ))
  message(paste(length(diagAfterLastSeen), "patients with diagnosis date after last seen date"))

  #removeRows = unique(c(stillNa, noChemoStart, chemoBeforeDiag, opNoDate, opPreDiag, diagPostDeath, diagAfterLastSeen))
  removeRows = stillNa

  removed = df[removeRows, ]
  df = df[-removeRows,]

  return(list('data' = df, 'removed' = removed ))
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
