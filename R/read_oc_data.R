
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
  tc <- rows.as.patient.id(tc, 2)

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
    di <- get.table.data(ocs=ocs, table=table, uniqueID=2, rulesFile=rulesFile, colFilter=filter)
  } else {
    di <- read.table.data(table, uniqueID=2, rulesFile=rulesFile)
  }

  #di$DI.ageAtDiagnosis <- as.numeric(di$DI.ageAtDiagnosis)
  # Prefer to get this from FE tables
  di$DI.PatientDateOfDeath = NULL

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
    ex <- get.table.data(ocs=ocs, table=grep('^ex_',tables,value=T), uniqueID=2, rulesFile=rulesFiles[1], colFilter=filter)
    ex_history_gastric <- get.table.data(ocs=ocs, table=grep('ex.*gastric',tables,value=T), rulesFile=rulesFiles[2], colFilter=filter1)
    ex_history_other <- get.table.data(ocs=ocs, table=grep('ex.*other_cancer',tables,value=T), rulesFile=rulesFiles[3], colFilter=filter2)

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

  if (verbose) message(paste(nrow(ex), "patients"))

  return(ex)
}

read.endpoints<-function(ocs, tables, occams_ids=NULL, rulesFiles=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  require(plyr)

  tables = list.clinical.tables(ocs,'fe')

  if (length(tables) != 2)
    stop("Two tables expected for endpoints")

  filter = make.id.filter(occams_ids,'FE_StudySubjectID')
  filter1 = make.id.filter(occams_ids,'FE1_StudySubjectID')

  if (inherits(ocs, "OCCAMSLabkey")) {
    #stop("Create connection to OCCAMS Labkey first")
    if (verbose) message(paste("Reading",paste(tables, collapse=",")))

    fe <- get.table.data(ocs,grep('fe_',tables,value=T), uniqueID=2, rulesFile=rulesFiles[1], colFilter=filter)
    fe2 <- get.table.data(ocs,grep('fe1',tables,value=T), rulesFile=rulesFiles[2], colFilter=filter1)
    fe2 <- ddply( ddply(fe2, .(FE1.StudySubjectID), arrange, desc(FE1.DateOfUpdate)),
                  .(FE1.StudySubjectID, FE1.StudySite), plyr::summarise, FE1.DateOfUpdate=FE1.DateOfUpdate[1], FE1.ReasonForFollowUp=FE1.ReasonForFollowUp[1], FE1.HasOriginalDiseaseReoccurred=FE1.HasOriginalDiseaseReoccurred[1])
  } else {
    if (verbose) message(paste("Reading",paste(basename(tables), collapse=", ")))

    #fe <- read.table.data(tables[1], uniqueID=2, rulesFile=rulesFiles[1])
    #fe2 <- read.table.data(tables[2], rulesFile=rulesFiles[2])
    #fe2 <- ddply( ddply(fe2, .(FE1.StudySubjectID), arrange, desc(FE1.DateOfUpdate)),
    #              .(FE1.StudySubjectID, FE1.StudySite), plyr::summarise, FE1.DateOfUpdate=FE1.DateOfUpdate[1], FE1.ReasonForFollowUp=FE1.ReasonForFollowUp[1], FE1.HasOriginalDiseaseReoccurred=FE1.HasOriginalDiseaseReoccurred[1])
  }

  # Final endpoint
  rownames(fe2) = fe2$FE1.StudySubjectID

  feF <- merge.patient.tables(fe,fe2, all=T)

  # withdrawn from study
  rows <- which(with(feF, grepl("withdrawal", FE.EndPoint)))
  withdrawn <- feF[rows,]
  if (length(rows) > 0) feF <- feF[-rows,]

  feF[which(with(feF, FE.EndPoint != 'Patient died' & !is.na(FE.ReasonForPatientDeath))), 'FE.EndPoint'] <- 'Patient died'

  feF$FE.CancerFree.Discharge <- as.factor(with(feF, grepl("discharge", FE.AdditionalDetailsForPatientEndPoint, ignore.case=T) & (FE.EndPoint != 'Patient died' | is.na(FE.EndPoint))))

  feF$FE.Patient.Died <- as.factor(ifelse(feF$FE.EndPoint == 'Patient died' & !is.na(feF$FE.EndPoint), "yes","no"))

  if (verbose) message(paste(nrow(feF), "patients"))

  return(list('fe'=feF, 'withdrawn'=rownames(withdrawn)))
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
  } else {
    ps <- read.table.data(tables[grep('^ps_',basename(tables))], rulesFile=rulesFiles[1])
  }

  # There should not be duplicates in this table...
  dups <- which(with(ps, PS.StudySubjectID %in% names(which(table(PS.StudySubjectID) > 1))))
  duplicatePts <- ps[dups,]

  if (length(dups) > 0) ps <- ps[-dups, ]
  ps <- rows.as.patient.id(ps, 2)

  ps <- fix.tumor.factors(text.to.tstage(ps))

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
  tp <- rows.as.patient.id(tp,2)

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
      tr <- get.table.data(ocs,grep('^tr_', tables, value=T), uniqueID=2, rulesFile=rulesFiles[1], colFilter=filter)
    } else {
      tr <- read.table.data(tables[grep('^tr_', basename(tables))], uniqueID=2, rulesFile=rulesFiles[1])
    }

  #dups <- which(with(tr, TR.StudySubjectID %in% names(which(table(TR.StudySubjectID) > 1))))

  # Currently there's nothing in the other tables
  #trc <- get.table.data(ocs,tables[3], rulesFile=rulesFiles[3])

  #cols = c('TR.CurativeTreatmentModality','TR.PalliativeTreatmentModality')
  #sapply(tr[cols], function(x) grep('chemo',x, value=T, ignore.case=T))

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
  dups <- which(duplicated(st[,2]))
  duplicatePts <- st[dups,]

  if (length(dups) > 0) st <- st[-dups, ]
  st <- rows.as.patient.id(st, 2)

  st$ST.SurgeryPerformed <- with(st, ifelse (ST.MainSurgery == 'yes' & (is.na(ST.ReasonIfNoSurgery) | is.na(ST.OtherReasonForNoSurgery)), 'yes', 'no'))
  st$ST.PathologyReportGenerated <- ifelse(with(st,
                                               ST.SurgeryPerformed == 'yes' & !grepl("Open (and|&) Shut", ST.Procedure, ignore.case=T)),
                                          'yes', 'no')

  st <- st[grep("\\.[1-9]$",colnames(st), invert=T, value=T)]

  cols <- grep("ASAGrade|Percentage", colnames(st), value=T)
  st[cols] <- lapply(st[cols], as.numeric)

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

  dups <- which(duplicated(rp[,2]))
  duplicatePts <- rp[dups,]

  if (length(dups) > 0) rp <- rp[-dups, ]

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


  file = system.file("extdata", "be_updates_20160930.txt", package="openclinica.occams")
  if (!exists('be_updates') & (!is.null(file) & file.exists(file))) {
    if (verbose) message(paste("Reading", file))
    be_updates = readr::read_tsv(file)
  }

  if (verbose) message("Updating RP.BarrettsAdjacent patients from Caitrona's worksheet.")
  be_updates = be_updates %>% dplyr::mutate( `Barret's Confirmed`=ifelse(`Barret's Confirmed` == '?', NA, `Barret's Confirmed`) )
  be_updates$`Barret's Confirmed` = as.integer(as.factor(be_updates$`Barret's Confirmed`))-1
  # Assume Caitron's are the final say (they matched other than NA anyhow)
  be_updates = subset(be_updates, `OCCAMS/ID` %in% intersect(rownames(rp),be_updates$`OCCAMS/ID`))

  rp[ be_updates$`OCCAMS/ID`,'RP.BarrettsAdjacent'] = be_updates$`Barret's Confirmed`

  rp$RP.BarrettsAdjacent <- revalue(as.factor(rp$RP.BarrettsAdjacent), reg)

  cols = c('RP.Nstage.RP.TNMSystem', 'RP.Nstage.RP.TNM6')
  rp[cols] = NULL

  if(verbose) message(paste(nrow(rp), "patients"))

  return(rp)
}

read.referral.diagnosis<-function(ocs, tables, occams_ids=NULL, rulesFile=NULL, ...) {
  z = list(...)
  verbose = ifelse (!is.null(z$verbose), z$verbose, F)

  tables = list.clinical.tables(ocs,'rd')

  filter = make.id.filter(occams_ids,'RD_StudySubjectID')

  if (inherits(ocs, "OCCAMSLabkey")) {
    if (verbose) message(paste("Reading",tables))
    rd <- get.table.data(ocs,grep('^rd_', tables, value=T), uniqueID=2, rulesFile=rulesFile, colFilter=filter)
  } else {
    rd <- read.table.data(tables[grep('^rd_', basename(tables))], uniqueID=2, rulesFile=rulesFile)
  }
  rd$RD.SiewertClassification <- ordered(rd$RD.SiewertClassification, levels=c(1,2,3))

  if(verbose) message(paste(nrow(rd), "patients"))
  return(rd)
}

merge.patient.tables<-function(df1, df2, ...) {
  newDF <- merge.data.frame(df1, df2, by='row.names', ...)
  rownames(newDF) <- newDF[['Row.names']]
  newDF[['Row.names']] <- NULL
  return(newDF)
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

#' Download in wide format the Labkey data for all or a selected group of patients.
#' @name download.all.tables
#' @param ocs
#' @param occams_ids Array of occams identifiers, this is optional. The default returns all patients, better to use get.patients(...) if you want this.
#' @param missing Replace missing values with this string, default is NA
#' @param verbose
#'
#' @author
#' @export
download.all.tables<-function(ocs, occams_ids=NULL, missing=NULL, verbose=T) {
  require(plyr)
  prefixes=c('di','rd','ex','ps','tp','tr','st','rp','fe','tc')

  if (inherits(ocs, "OCCAMSLabkey")) {
    tables = list.clinical.tables(ocs, prefixes)
  } else {
    tables = list.files(ocs, full.names=T)
  }

  ordered_tables = sapply(prefixes, function(x) {
    tables[sort(grep(paste("^",x,sep=""), basename(tables)))]
  })

  di <- clean.missing(
      read.demographics(ocs, ordered_tables$di, rulesFile=paste(path.package('openclinica.occams')[1], "editrules/di_editrules.txt", sep='/'),occams_ids=occams_ids,verbose=verbose),
    missing )

  ex <- clean.missing( read.exposures(ocs, ordered_tables$ex, rulesFiles=NULL, occams_ids=occams_ids, verbose=verbose), missing )
  fe <- read.endpoints(ocs, ordered_tables$fe,
                       rulesFiles=paste(path.package('openclinica.occams'), c('editrules/fe_editrules.txt', 'editrules/fe1_editrules.txt'), sep='/'), occams_ids=occams_ids, verbose=verbose)
  withdrawn <- fe$withdrawn
  fe <- clean.missing(fe$fe, missing)

  rd <- clean.missing( read.referral.diagnosis(ocs, ordered_tables$rd, rulesFile=NULL, occams_ids=occams_ids, verbose=verbose), missing )

  # these two have duplicates as well
  ps <- clean.missing( read.prestage(ocs, ordered_tables$ps, rulesFiles=NULL, occams_ids=occams_ids, verbose=verbose)$ps, missing )
  tp <- clean.missing( read.treatment.plan(ocs, ordered_tables$tp, rulesFiles=NULL, occams_ids=occams_ids, verbose=verbose)$tp, missing )

  tr <- clean.missing( read.therapy(ocs, ordered_tables$tr, rulesFiles=NULL, occams_ids=occams_ids, verbose=verbose), missing )
  st <- read.surgery(ocs, ordered_tables$st, rulesFile=NULL, occams_ids=occams_ids, verbose=verbose)

  rp <- clean.missing( read.pathology(ocs, ordered_tables$rp, occams_ids=occams_ids, rulesFile=NULL, verbose=verbose), missing)

  if (verbose) message("Creating final patient table")
  # Merge final table
  all <- merge.patient.tables(di, fe, all=T)
  all <- merge.patient.tables(all, tp, all=T)
  all <- merge.patient.tables(all, rd, all=T)
  all <- merge.patient.tables(all, ex, all=T)
  all <- merge.patient.tables(all, ps, all=T)
  all <- merge.patient.tables(all, tr, all=T)
  all <- merge.patient.tables(all, st, all=T)
  all <- merge.patient.tables(all, rp, all=T)

  rm = which(rownames(all) %in% withdrawn)
  if (length(rm) > 0) all = all[-rm,]

  if (verbose) message(paste("Final patient total:", nrow(all)))

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

  if ( length(which(is.na(all$FE.DateOfPatientDeath))) == nrow(all) ) {
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




