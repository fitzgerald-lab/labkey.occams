#' Download all tables directly from Labkey for all or a subset of patients.
#' @name dump.clinical.data
#' @param ocs labkey connection object from connect.to.labkey()
#' @param occams_ids Array of occams identifiers, this is optional. The default returns all patients, better to use get.patients(...) if you want this.
#' @param prefixes Get only tables with the given prefixes. Call 'get.prefixes(ocs)' for a list of only CRF tables.
#' @param version Table version (default 'z1')
#' @outdir The path to output the files (default is <user home>/tmp)
#' @author skillcoyne
#' @export
dump.clinical.data<-function(ocs, occams_ids = NULL, prefixes=NULL, version='z1', outdir="~/tmp") {
  out <- path.expand(paste(outdir, Sys.Date(), sep="/"))

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
#' NOTE: Major changes from 2020 Apr 22 onwards. All of the therapy tables are merged.  There is no longer a TR.Chemotherapy column, to check if a patient has chemo
#' look at TR.NeoAdj or TR.Adj columns.  
#' All followup/endpoint/recurrence tables are now merged.  Only a single row per patient is provided with the most recent follow-up date, 
#' the date of recurrence (if there is one) and the date of death (if known).
#' @name download.wide.format
#' @param ocs Labkey connection object from connect.to.labkey()
#' @param occams_ids Array of occams identifiers, this is optional. The default returns all patients, better to use get.patients(...) if you want this.
#' @param missing Replace missing values with this string, default is NA
#' @param versions The table versions used from labkey, default is 'z1'
#' @param verbose
#' @return named list of tibbles: patients, family_history, tissues
#' @author skillcoyne
#' @export
download.wide.format<-function(ocs, occams_ids=NULL, missing=NULL, versions='z1', verbose=T) {
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

  di <- read.demographics(ocs, ordered_tables$di, rulesFile=paste(path.package('openclinica.occams')[1], "editrules/di_editrules.txt", sep='/'),occams_ids=occams_ids,verbose=verbose) %>% dplyr::rename_at(vars(matches('Study')), list(~sub('^DI\\.','',.)))

  exes <- read.exposures(ocs, ordered_tables$ex, rulesFiles=NULL, occams_ids=occams_ids, verbose=verbose)
  ex <- exes$ex %>%  dplyr::rename_at(vars(matches('Study')), list(~sub('^EX\\.','',.)))
  history <- exes$history

  feep <- read.followup(ocs, c(ordered_tables$fe,ordered_tables$ep,ordered_tables$fr), rulesFiles=paste(path.package('openclinica.occams'), c('editrules/fe_editrules.txt', 'editrules/fe1_editrules.txt', 'editrules/ep_editrules.txt'), sep='/'), occams_ids=occams_ids, verbose=verbose)
  fe <- feep$fe %>%  dplyr::rename_at(vars(matches('Study')), list(~sub('^FE\\.','',.)))
  
  rd <- read.referral.diagnosis(ocs, ordered_tables$rd, rulesFile=paste(path.package('openclinica.occams'),'editrules/rd_editrules.txt',sep='/'), occams_ids=occams_ids, verbose=verbose) %>%  dplyr::rename_at(vars(matches('Study')), list(~sub('^RD\\.','',.)))

  ps <- read.prestage(ocs, ordered_tables$ps, rulesFiles=NULL, occams_ids=occams_ids, verbose=verbose)$ps %>%  
    dplyr::rename_at(vars(matches('Study')), list(~sub('^PS\\.','',.)))
  
  tp <- read.treatment.plan(ocs, ordered_tables$tp, rulesFiles=NULL, occams_ids=occams_ids, verbose=verbose)$tp %>% 
    dplyr::rename_at(vars(matches('Study')), list(~sub('^TP\\.','',.)))
  
  tr <- read.therapy(ocs, ordered_tables$tr, rulesFiles=NULL, occams_ids=occams_ids, verbose=verbose) %>%  
    dplyr::rename_at(vars(matches('Study')), list(~sub('^TR\\.','',.)))

  st <- read.surgery(ocs, ordered_tables$st, rulesFile=NULL, occams_ids=occams_ids, verbose=verbose) %>% 
    dplyr::rename_at(vars(matches('Study')), list(~sub('^ST\\.','',.)))

  rp <- read.pathology(ocs, ordered_tables$rp, occams_ids=occams_ids, rulesFile=paste(path.package('openclinica.occams'),'editrules/rp_editrules.txt',sep='/'), verbose=verbose) %>% 
    dplyr::rename_at(vars(matches('Study')), list(~sub('^RP\\.','',.)))

  if (verbose) message("Creating final patient table")
  # Merge final table
  all <- dplyr::left_join(di,ex,by=c('StudySubjectID','StudySite')) %>%
    dplyr::left_join(rd,by=c('StudySubjectID','StudySite')) %>%
    dplyr::left_join(tp,by=c('StudySubjectID','StudySite')) %>%
    dplyr::left_join(ps,by=c('StudySubjectID','StudySite')) %>%
    dplyr::left_join(tr,by=c('StudySubjectID','StudySite')) %>%
    dplyr::left_join(st,by=c('StudySubjectID','StudySite')) %>%
    dplyr::left_join(rp,by=c('StudySubjectID','StudySite')) %>%
    dplyr::left_join(fe,by=c('StudySubjectID','StudySite')) %>%
    dplyr::rename(ID=StudySubjectID) %>% dplyr::group_by(ID) %>% 
    dplyr::select(-contains('CRF'), -contains('OpenClinica')) %>%
    dplyr::ungroup() %>% recode.siewert %>% 
    #recode.TNM %>%
    dplyr::group_by(ID) %>% dplyr::select(ID, StudySite, everything())

  if (verbose) message(paste("Final patient total:", nrow(all)))

  all <- prestage.to.path(all)

  # Dates after diagnosis that patients may have been seen
  clinDates <- grep('(FE|RP|ST|TR).+Date', colnames(all), value=T)

  # In order
  clinDates <- all %>% ungroup %>%
    dplyr::select(ID,FE.DateOfPatientDeath, FE.LatestDateOfUpdate.c, FE.DateOriginalDiseaseReoccurred, matches('ST.+Date'), matches('TR.NeoAdj.+Date'), matches('TR.Radio.+Date'), matches('TR.Endo.+Date')) %>%
    mutate_at(vars(-ID), list(~as.Date(format.Date(.,'%Y-%m-%d')))) 

  # can't work this one out in dplyr...
  clinDates$FE.LastSeenDate.c <- as.Date(apply(clinDates[-1], 1, function(x) {
    y <- sort(x, decreasing=T)
    lastDate = y[1]

    if (!is.na(y['FE.DateOfPatientDeath']) & y['FE.DateOfPatientDeath'] != y[1])
      lastDate = y['FE.DateOfPatientDeath']

    return(lastDate)
  })) 

  clinDates <- clinDates %>% dplyr::select(ID,FE.LastSeenDate.c) 
  
  all <- dplyr::inner_join(all, clinDates,  by='ID')

  # Dates after diagnosis that patients may have been seen
  diagDates <- grep('DI.DateOfDiagnosisOGC|RD.DateOfOGCDiagnosis|DI.DateOfInformedConsentForCAMSAndICGC', colnames(all), value=T)
  all$RD.DiagnosisDate.c <- as.Date(apply(all[diagDates], 1, function(x) {
    y <- sort(x, decreasing=F)
    firstDate <- y[1]

    if (!is.na(y['DI.DateOfDiagnosisOGC']))
      firstDate <- y['DI.DateOfDiagnosisOGC']

    return(firstDate)
  }))

  all <- all %>% dplyr::group_by(ID) %>%
    dplyr::mutate(Weeks.Survival.c = round(as.numeric(difftime(FE.LastSeenDate.c, RD.DiagnosisDate.c, units='weeks')))) %>% 
    dplyr::select(ID, StudySite, Weeks.Survival.c, RD.DiagnosisDate.c, FE.LastSeenDate.c, everything())

  bad = with(all, which(Weeks.Survival.c < 1))
  message(paste(length(bad), "patients have a diagnosis date before their last seen (death/surgery/etc) date."))
  all[bad, c('Weeks.Survival.c')] = NA

  tc <- read.tissue.collection(ocs, ordered_tables$tc, rulesFiles=NULL, occams_ids=occams_ids)
  return(list('patients'=all, 'family_history'=history, 'tissues'=tc))
}


#' Get selected patients in the wide table format
#' @name get.patients
#' @param ocs Connection from connect.to.labkey()
#' @param occams_ids Array of occams identifiers to select
#'
#' @author
#' @export
get.patients<-function(ocs, occams_ids, verbose=T) {
  require(Rlabkey)

  if (!inherits(ocs, "OCCAMSLabkey"))
    stop("Create connection to OCCAMS Labkey first")

  ids = grepl('^OCCAMS/[A-Z]{2}/[0-9]+', occams_ids)

  incorrect_ids = occams_ids[!ids]
  ids = occams_ids[ids]

  if (verbose) message(paste("Retrieving",length(ids),"patients."))

  if (length(incorrect_ids) > 0)
    warning(paste(length(incorrect_ids)), ' are not OCCAMS identifiers.')

  occams = download.wide.format(ocs=ocs, verbose=verbose, occams_ids=ids)

  return(list('patients'=occams$patients, 'tissues'=occams$tissues, 'incorrect_ids'=incorrect_ids, 'withdrawn'=occams$withdrawn))
}



make.id.filter<-function(occams_ids, idCol) {
  if (!is.null(occams_ids)) {
    occams_ids = grep('^OCCAMS/[A-Z]{2}/[0-9]+', occams_ids, value=T)
    return( makeFilter(c(idCol, 'IN', paste(occams_ids, collapse=';'))) )
  }
  return(NULL)
}





