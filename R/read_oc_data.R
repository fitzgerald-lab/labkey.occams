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


#' USE download.wide.format instead. This is the old method.
#' @name download.all.tables
#' @param ocs
#' @param occams_ids Array of occams identifiers, this is optional. The default returns all patients, better to use get.patients(...) if you want this.
#' @param missing Replace missing values with this string, default is NA
#' @param verbose
#'
#' @author
#' @export
download.all.tables<-function(ocs, occams_ids=NULL, missing=NULL,versions='z1', verbose=T) {
  warning('This method has been deprecated, please use download.wide.format(...) instead')
  download.wide.format(ocs,occams_ids, missing, versions, verbose)
}

#' Download in wide format the Labkey data for all or a selected group of patients.
#' @name download.wide.format
#' @param ocs Labkey connection object from connect.to.labkey()
#' @param occams_ids Array of occams identifiers, this is optional. The default returns all patients, better to use get.patients(...) if you want this.
#' @param missing Replace missing values with this string, default is NA
#' @param versions The table versions used from labkey, default is 'z1'
#' @param verbose
#'
#' @author skillcoyne
#' @export
download.wide.format<-function(ocs, occams_ids=NULL, missing=NULL,versions='z1', verbose=T) {
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
## TODO: CHECK THAT I"M LOOKINGF AT FR1 table
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





