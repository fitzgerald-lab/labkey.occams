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

  file = path.expand(file)

  vars = readr::read_delim(file, delim="=", col_names = F, col_types = 'cc') %>% tidyr::spread(X1, X2)

  if (length(grep('path|pwd|url|user', colnames(vars) )) < 4) stop('Missing information in the labkey.cred file. path,pwd,url,user are required.')

  return(create.connection(url=vars$url, path=vars$path, user=vars$user, pwd=vars$pwd, schema='oc'))
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

anonymized.ids<-function(file="~/.labkey.cred") {
  
  file = path.expand(file)
  vars = readr::read_delim(file, delim="=", col_names = F, col_types = 'cc') %>% tidyr::spread(X1, X2)
  
  ocs <- create.connection(vars$url, "/ICGC/Cohorts/OCCAMS SECURE CASE ID", vars$user, vars$pwd, schema='lists')

  table = names(ocs$schema)[1]
  
  if (!grepl('occams_case_ids_secure_case_ids', table))
    stop(paste0("Expected table name to include 'occams_case_ids_secure_case_ids', tables available:\n\t", paste(names(ocs$schema), collapse='\n\t')) )
  
  rows <- suppressWarnings(Rlabkey::getRows(ocs$session, ocs$schema[[table]]))
  
  as_tibble(rows) %>% set_names(c('OCCAMS.ID','SecureID'))
}

create.connection<-function(url="https://occams.cs.ox.ac.uk/labkey", path="/ICGC/Cohorts/All Study Subjects", user=NULL, pwd=NULL, schema='oc') {
  session <- create.session(url,path,user,pwd)

  if (is.null(session))
    stop("Failed to create session with provided credentials")

  oc <- tryCatch({
    Rlabkey::getSchema(session, schema)
  }, error = function(e) {
    tryCatch({
      Rlabkey::getSchema(session, schema)
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
#' @param type All or only CRF prefixes, default 'clinical'
#' @author skillcoyne
#' @export
get.prefixes<-function(ocs, version=NULL, type=c('clinical','all')) {
  rettype = match.arg(type)

  if (is.null(version))
    version = unique(crf.prefixes$version)

  if (rettype == 'clinical') return(crf.prefixes$prefix)

  return(unique(sapply(grep(version, names(ocs$schema), value=T), substr, 1, 2)))
}

#' Tibble for clinical table prefixes current 2019 Oct
#' @name crf.prefixes
#' @export
crf.prefixes<-tibble(
  table = c( 'demographics','exposures','referral diagnosis','pretreatment staging',
             'treatment plan','therapy record','surgical treatment','resection pathology',
             'recurrence','follow-up','endpoint'),
  prefix = c( 'di','ex','rd','ps','tp','tr','st','rp','fr','fe','ep'),
  version = 'z1'
)

#' List the clinical tables available in Labkey
#' @name list.clinical.tables
#' @param ocs Labkey connection
#' @prefixes Only those with matches prefixes, default is the result of get.prefixes
#' @param version version for tables, default 'z1'
#'
#' @author skillcoyne
#' @export
list.clinical.tables<-function(ocs, prefixes=NULL) {
  if (is.null(prefixes))
    prefixes = crf.prefixes$prefix

  if (!inherits(ocs, "OCCAMSLabkey"))
    stop("Create connection to OCCAMS Labkey first")

  tables <- grep( paste('^(', paste( paste(prefixes, '.+', unique(crf.prefixes$version),sep=''), collapse='|'), ')', sep=''), names(ocs$schema), value=T)
  if (length(tables) < length(prefixes))
    warning("Fewer tables found than prefixes provided.")

  return(tables)
}

#' Get tables for the given CRF term (see crf.prefixes) or prefix.
#' @name get.tables.by.crf
#' @param ocs Labkey connection
#' @term
#'
#' @author skillcoyne
#' @export
get.tables.by.crf<-function(ocs, crf) {
  crf = tolower(crf)
  tb = crf.prefixes %>% filter(table %in% crf | prefix %in% crf)
  list.clinical.tables(ocs, tb$prefix)
}



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
    if (nrow(rows) <= 0)
      stop("No rows found that matched the column filters.")
  }

    rows <- suppressWarnings(clean.rows(rows, uniqueID, rulesFile))

  if (!is.null(uniqueID))
    rows <- suppressWarnings(rows %>% dplyr::group_by_at(vars(contains(uniqueID))))

  rows <- rows %>% dplyr::select(-contains('CRF'))

  return(rows)
}

clean.rows<-function(rows, uniqueID=NULL, rulesFile=NULL) {
  #message('clean rows')
  rows = rows %>% rename_all(funs(gsub('_','.',.)))

  # Trim whitespace & fix NAs
  #rows %>% rowwise %>% mutate_if(is.character, funs(clean.na))
  rows = rows %>% dplyr::mutate_if(is.character, funs(trim))

  if (!is.null(rulesFile) & nrow(rows) > 0) {
    failed_rules <- editrules(rulesFile, rows)
    rows <- failed_rules$df
  }

  rows <- rows %>% rowwise %>% dplyr::mutate_if(is.character, list(clean.na))

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

  groups = dplyr::group_vars(df)

  # Numeric
  cols = rules %>% dplyr::filter(Def == 'numeric') %>% dplyr::select(Column) %>% pull
  if (length(cols) > 0) df = df %>% ungroup %>% dplyr::mutate_at(cols, funs(as.numeric))

  # Date
  cols = rules %>% dplyr::filter(Def == 'date') %>% dplyr::select(Column) %>% pull
  if (length(cols) > 0) df = df %>% ungroup %>% dplyr::mutate_at(cols, list(~as.Date(.,origin = "1960-10-01")))

  # Character
  cols = rules %>% dplyr::filter(Def == 'character') %>% dplyr::select(Column) %>% pull
  if (length(cols) > 0) df = df %>% ungroup %>% dplyr::mutate_at(cols, funs(as.character))

  # Factor TODO -- not sure how to make specfic rules work with dplyr
  cols = rules %>% dplyr::filter(Def == 'factor') %>% dplyr::select(Column) %>% pull
  if (length(cols) > 0) df = df %>% dplyr::mutate_at(cols, funs(as.factor))

  failed <- apply(rules,1, function(x) {
    if (length(intersect(x[['Column']], colnames(df))) <= 0) {
      x[['Column']]
    } else {
      which(!with(df, eval(parse(text=paste(x[['Rule']], collapse="")))))
    }
  })

  if (length(failed) > 0) {
    names(failed) <- rules$Column

    for (n in names(failed)) {
      df[failed[[n]], n] <- NA
    }
  }

  df = df %>% dplyr::group_by_at(groups)

  return(list("f" = lapply(failed, length), "df"= df))
}


