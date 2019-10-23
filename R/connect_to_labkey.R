require(tidyverse)
require(Rlabkey)

#' Create connection to Labkey
#' @name connect.to.labkey
#' @param file labkey credential file path, default is ~/.labkey.cred
#'
#' @author
#' @export
connect.to.labkey<-function(file="~/.labkey.cred") {
  if (is.null(file))
    stop("File with url, path, user, pwd entries required for connection.")

    vars = readr::read_delim(file, delim="=", col_names = F, col_types = 'cc') %>% tidyr::spread(X1, X2)

  if (length(grep('path|pwd|url|user', colnames(vars) )) < 4) stop('Missing information in the labkey.cred file. path,pwd,url,user are required.')

  return(create.connection(url=vars$url, path=vars$path, user=vars$user, pwd=vars$pwd))
}

create.session<-function(url="https://occams.cs.ox.ac.uk/labkey", path="/ICGC/Cohorts/All Study Subjects", user, pwd) {
  require(Rlabkey)
  #  ssc <- Rlabkey::labkey.acceptSelfSignedCerts()
  #labkey.setCurlOptions(ssl_verifypeer=FALSE, ssl_verifyhost=FALSE)
  #labkey.setDefaults(apiKey="session|abcdef0123456789abcdef0123456789")
  session <- Rlabkey::getSession(baseUrl=url,
                                 folderPath=path,
                                 curlOptions= RCurl::curlOptions(userpwd=paste(user,pwd,sep=":")),
                                 lkOptions=NULL)
  return(session)
}

create.connection<-function(url="https://occams.cs.ox.ac.uk/labkey", path="/ICGC/Cohorts/All Study Subjects", user, pwd) {
  session <- create.session(url,path,user,pwd)

  if (is.null(session))
    stop("Failed to create session with provided credentials")

  oc <- tryCatch({
    Rlabkey::getSchema(session, "oc")
  }, error = function(e) {
    tryCatch({
      Rlabkey::getSchema(session, 'oc')
    }, error = function(e1){
      stop(paste("Failed to get schema from LabKey server:", e ))
    })
  })

  list(oc,session)->occams.session
  names(occams.session)<-c('schema','session')

  class(occams.session)<-c("OCCAMSLabkey","list")

  return(occams.session)
}


#' Get the table prefixes from the schema
#' @name get.prefixes
#' @param ocs Labkey connection
#' @param version version for tables, default 'z1'
#'
#' @author skillcoyne
#' @export
get.prefixes<-function(ocs, version='z1') {
  return(unique(sapply(grep(version, names(ocs$schema), value=T), substr, 1, 2)))
}


#' List the clinical tables available in Labkey
#' @name list.clinical.tables
#' @param ocs Labkey connection
#' @prefixes Only those with matches prefixes, default is the result of get.prefixes
#' @param version version for tables, default 'z1'
#'
#' @author skillcoyne
#' @export
list.clinical.tables<-function(ocs, prefixes=NULL, versions='z1' ) {
  if (is.null(prefixes))
    prefixes = get.prefixes(ocs,versions)

  if (!inherits(ocs, "OCCAMSLabkey"))
    stop("Create connection to OCCAMS Labkey first")

  tables <- grep( paste('^(', paste( paste(prefixes, '.+', versions,sep=''), collapse='|'), ')', sep=''), names(ocs$schema), value=T)
  if (length(tables) < length(prefixes))
    warning("Fewer tables found than prefixes provided.")

  return(tables)
}

# rows.as.patient.id<-function(rows, uniqueID=NULL) {
#   if (!is.null(uniqueID) & is.numeric(uniqueID))
#     rownames(rows) <- rows[,uniqueID]
#   return(rows)
# }

read.table.data<-function(f, columns=NULL, uniqueID=NULL, rulesFile=NULL, ...) {
  require(data.table)

  rows <- as.data.frame(data.table::fread(f, stringsAsFactors=F, sep="\t", header=T, quote='"'))
  rows <- clean.rows(rows, uniqueID, rulesFile)

  return(rows)
}

get.table.data<-function(ocs, table, columns=NULL, uniqueID='StudySubjectID', rulesFile=NULL, ...) {
  if (!inherits(ocs, "OCCAMSLabkey"))
    stop("Create connection to OCCAMS Labkey first")

  rows <- NULL
  if (is.null(columns)) {
    rows <- suppressWarnings(Rlabkey::getRows(ocs$session, ocs$schema[[table]], ...))
  } else {
    if (length(columns > 1)) columns <- paste(columns, collapse=",")
    rows <- suppressWarnings(Rlabkey::getRows(ocs$session, ocs$schema[[table]], colSelect=columns, ...))
  }

  rows <- suppressWarnings(clean.rows(rows, uniqueID, rulesFile))
  if (!is.null(uniqueID))
    rows = suppressWarnings(rows %>% dplyr::group_by_at(vars(contains(uniqueID))))

  rows <- rows %>% select(-contains('CRF'))

  return(rows)
}

clean.rows<-function(rows, uniqueID=NULL, rulesFile=NULL) {
  #message('clean rows')
  rows = rows %>% rename_all(funs(gsub('_','.',.)))

  # Trim whitespace & fix NAs
  #rows %>% rowwise %>% mutate_if(is.character, funs(clean.na))
  rows = rows %>% mutate_if(is.character, funs(trim))

  if (!is.null(rulesFile)) {
    failed_rules <- editrules(rulesFile, rows)
    rows <- failed_rules$df
  }

  rows = rows %>% rowwise %>% mutate_if(is.character, funs(clean.na))

  return(rows)
}

clean.na<-function(x, text=c('not known','not_recorded','not recorded', 'not assessed', 'unknown', '^unk$', '^n/a$','^n.a$','^na$')) {
    if(!is.na(x) & (grepl(paste(text, collapse="|"), x, ignore.case=T) | x == ''))
      return(NA)
    return(x)
}

editrules<-function(file, df, verbose=T) {
  if (!file.exists(file))
    stop(paste("File does not exist or is not readable: ", file))
  rules <- readr::read_delim(file, ':', col_types = 'ccc')

  if(verbose) message(paste("Applying rules found in", file))

  groups = group_vars(df)

  # Numeric
  cols = rules %>% filter(Def == 'numeric') %>% select(Column) %>% pull
  df = df %>% ungroup %>% mutate_at(cols, funs(as.numeric))

  # Date
  cols = rules %>% filter(Def == 'date') %>% select(Column) %>% pull
  df = df %>% ungroup %>% mutate_at(cols, funs(as.Date))

  # Character
  cols = rules %>% filter(Def == 'character') %>% select(Column) %>% pull
  df = df %>% ungroup %>% mutate_at(cols, funs(as.character))

  # Factor TODO -- not sure how to make specfic rules work with dplyr
  cols = rules %>% filter(Def == 'factor') %>% select(Column) %>% pull
  df = df %>% mutate_at(cols, funs(as.factor))

  failed <- apply(rules,1, function(x) {
    if (length(intersect(x[['Column']], colnames(df))) <= 0) {
      x[['Column']]
    } else {
      which(!with(df, eval(parse(text=paste(x[['Rule']], collapse="")))))
    }
  })
  names(failed) <- rules$Column

  for (n in names(failed)) {
    df[failed[[n]], n] <- NA
  }

  df = df %>% group_by_at(groups)

  return(list("f" = lapply(failed, length), "df"= df))
}


#' Get selected patients in the wide table format
#' @name get.patients
#' @param ocs Connection from connect.to.labkey()
#' @param occams_ids Array of occams identifiers to select
#'
#' @author
#' @export
get.patients<-function(ocs, occams_ids, verbose=F) {
  require(Rlabkey)

  if (!inherits(ocs, "OCCAMSLabkey"))
    stop("Create connection to OCCAMS Labkey first")

  ids = grepl('^OCCAMS/[A-Z]{2}/[0-9]+', occams_ids)

  incorrect_ids = occams_ids[!ids]
  ids = occams_ids[ids]

  if (verbose) message(paste("Retrieving",length(ids),"patients."))

  if (length(incorrect_ids) > 0)
    warning(paste(length(incorrect_ids)), ' are not OCCAMS identifiers.')

  occams = download.all.tables(ocs=ocs, verbose=verbose, occams_ids=ids)

  return(list('patients'=occams$patients, 'tissues'=occams$tissues, 'incorrect_ids'=incorrect_ids, 'withdrawn'=occams$withdrawn))
}
