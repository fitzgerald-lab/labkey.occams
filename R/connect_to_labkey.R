
create.session<-function(url="https://occams.comlab.ox.ac.uk/labkey", path="/ICGC/Cohorts/All Study Subjects", user, pwd) {
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

connect.to.labkey<-function(file="~/.labkey.cred") {
  if (is.null(file))
    stop("File with url, path, user, pwd entries required for connection.")

  vars <- read.csv(file, header=F, sep="=", row.names=1, stringsAsFactors=F)

  return(create.connection(url=vars['url',], path=vars['path',], user=vars['user',], pwd=vars['pwd',]))
}

create.connection<-function(url="https://occams.comlab.ox.ac.uk/labkey", path="/ICGC/Cohorts/All Study Subjects", user, pwd) {
  session <- create.session(url,path,user,pwd)

  if (is.null(session))
    stop("Failed to create session with provided credentials")

  oc <- tryCatch({
    Rlabkey::getSchema(session, "oc")
  }, error = function(e) {
    tryCatch({
      Rlabkey::getSchema(session)
    }, error = function(e1){
      stop(paste("Failed to get schema from LabKey server:", e ))
    })
  })

  list(oc,session)->occams.session
  names(occams.session)<-c('schema','session')

  class(occams.session)<-c("OCCAMSLabkey","list")

  return(occams.session)
}

get.prefixes<-function(ocs, version='z1') {
  return(unique(sapply(grep(version, names(ocs$schema), value=T), substr, 1, 2)))
}

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

rows.as.patient.id<-function(rows, uniqueID=NULL) {
  if (!is.null(uniqueID) & is.numeric(uniqueID))
    rownames(rows) <- rows[,uniqueID]
  return(rows)
}

read.table.data<-function(f, columns=NULL, uniqueID=NULL, rulesFile=NULL, ...) {
  require(data.table)

  rows <- as.data.frame(data.table::fread(f, stringsAsFactors=F, sep="\t", header=T, quote='"'))
  rows <- clean.rows(rows, uniqueID, rulesFile)

  return(rows)
}

get.table.data<-function(ocs, table, columns=NULL, uniqueID=NULL, rulesFile=NULL, ...) {
  require(Rlabkey)

  if (!inherits(ocs, "OCCAMSLabkey"))
    stop("Create connection to OCCAMS Labkey first")

  rows <- NULL
  if (is.null(columns)) {
    rows <- suppressWarnings(Rlabkey::getRows(ocs$session, ocs$schema[[table]], ...))
  } else {
    if (length(columns > 1)) columns <- paste(columns, collapse=",")
    rows <- suppressWarnings(Rlabkey::getRows(ocs$session, ocs$schema[[table]], colSelect=columns, ...))
  }

  rows <- clean.rows(rows, uniqueID, rulesFile)

  return(rows)
}

clean.rows<-function(rows, uniqueID=NULL, rulesFile=NULL) {
  rows <- rows.as.patient.id(rows, uniqueID)

  colnames(rows) <- gsub("_", ".",colnames(rows))

  # Trim whitespace
  cols = which(sapply(rows[], is.character))
  rows[cols] = lapply(rows[cols],  trim)

  if (!is.null(rulesFile)) {
    failed_rules <- editrules(rulesFile, rows)
    rows <- failed_rules$df
  }

  # Fix the NAs
  cols <- which(sapply(rows, is.character))
  rows[cols] <- clean.na(rows[cols])

  return(rows)
}

clean.na<-function(df, text=c('not_recorded','not recorded', 'not assessed', 'unknown', '^unk$', 'n/a', '^na$')) {
  return( lapply(df[], function(x) { ifelse(grepl(paste(text, collapse="|"), x, ignore.case=T), NA, x) }) )
}


editrules<-function(file, df) {
  if (!file.exists(file))
    stop(paste("File does not exist or is not readable: ", file))
  rules <- read.table(file,header=T, stringsAsFactors=F, sep=":", quote="")
  message(paste("Applying rules found in", file))

  # Numeric
  cols = rules[which(rules$Def == 'numeric'), 'Column']
  df[cols] = lapply(df[cols], as.numeric)

  # Date
  cols = rules[which(rules$Def == 'date'), 'Column']
  df[cols] = lapply(df[cols], as.Date)

  # Character
  cols = rules[which(rules$Def == 'character'), 'Column']
  df[cols] = lapply(df[cols], as.character)

  # Factor
  cols = rules[which(rules$Def == 'factor'), 'Column']
  df[cols] = lapply(df[cols], as.factor)

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
  write.table(unlist(lapply(failed, length)), quote=F, col.names = F)

  return(list("f" = lapply(failed, length), "df"= df[rules$Column]))
}

