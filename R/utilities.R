# Select the filled in column, or the preferred one (i) if both are filled in
mr<-function(a,b,i=1,empty=NA) {
  c = unique(na.omit(c(a[i],b[i])))
  if (is_empty(c)) c = empty
  return(c[i])
}

# Order the tumor, nodal, and metastatic stages
fix.tumor.factors<-function(fd) {
  # Need to fix the factors for tumor staging to ensure that they are essentially in size order 1 > 2 > 3 etc...
  orderedTstages = c('Tx','T0','Tis','T1','T1a','T1b','T2','T3','T4','T4a','T4b')
  orderedNstages = c('Nx','N0','N1','N2','N3')
  orderedMstages = c('Mx','M0','M1')

  fd = fd %>% ungroup %>% dplyr::mutate_at(vars(matches('(Tstage|.+_T$)')), funs(factor), levels=orderedTstages, ordered=T) %>%
    dplyr::mutate_at(vars(matches('(Nstage.+TNM7|.+_N$)')), funs(factor), levels=orderedNstages, ordered=T) %>%
    dplyr::mutate_at(vars(matches('(Mstage|.+_M$)')), funs(factor), levels=orderedMstages, ordered=T)

  return(fd)
}

# translate the text entered in the database to the tumor stages
text.to.tstage<-function(fd) {
  orderedTstages = c(NA,'Tx','T0','Tis','T1','T1a','T1b','T2','T3','T4','T4a','T4b')

  fd = fd %>% ungroup %>% dplyr::mutate_at(vars(matches('(Tstage|.+_T$)')), funs(recode),
                          'not_recorded'=NA_character_,
                          'primary_tumour_cannot_be_assessed'='Tx',
                          'no_evidence_of_primary_tumour'='T0',
                          'carcinoma_in_situ_high_grade_dysplasia'='Tis',
                          'tumour_invades_lamina_propria_muscularis_mucosae_or_submucosa'='T1',
                          'tumour_invades_lamina_propria_or_muscularis_mucosae'='T1a',
                          'tumour_invades_submucosa'='T1b',
                          'tumour_invades_muscularis_propria'='T2',
                          'tumour_invades_adventitia'='T3',
                          'tumour_invades_adjacent_structures'='T4',
                          'tumour_invades_pleura_pericardium_or_diaphragm'='T4a',
                          'tumour_invades_other_adjacent_structures_such_as_aorta_vertebral_body_or_trachea'='T4b' ) %>%
    dplyr::mutate_at(vars(matches('(Tstage|.+_T$)')), funs(factor), levels=orderedTstages, ordered=T)

  return(fd)
}

# Calculate BMI as the one entered is often wrong
bmicalc<-function(wt,ht) {
  if (is.na(wt) || is.na(ht) || wt <= 0 || ht <= 0)
    return(NA)
  return (signif(wt/((ht/100)^2),3))
}

# Trim whitespace
trim<-function(x) {
  gsub("^\\s+|\\s+$", "", x)
}

# Recode and order siwert types
recode.siewert<-function(df) {
  grp = group_vars(df)
  df = df %>% ungroup %>% dplyr::mutate_at(vars(contains('SiewertClassification')), funs(recode_factor), '1'='Type I','2'='Type II','3'='Type III', .ordered=F, .missing=NA_character_ )
  if (length(grp) > 0) df = df %>% group_by(!!grp)
  return(df)
}

# order TNM stages
recode.TNM<-function(df) {
  grp = group_vars(df)
  df = df %>% ungroup %>% dplyr::mutate_at(vars(contains('TNMStage')), funs(factor), levels=c('I','II','III', 'IV'), ordered=T)
  if (length(grp) > 0) df = df %>% group_by(!!grp)
  return(df)
}


# Determine TNM stage (I,II,II,IV) from the individual tumor, node, and metastasis stages
tnmStage <- function(Ts,Ns,Ms) {
  if (is.na(Ts)) return(NA)
  if (!is.na(Ns) & Ns == 'unknown') Ns = NA
  if (!is.na(Ms) & Ms == 'unknown') Ms = NA

  stage = case_when(
    Ts %in% c('Tx','Tis') ~ '0',
    grepl('T1',Ts) & ( is.na(Ns) | Ns <= 'N0') & (is.na(Ms) | Ms <= 'M0') ~ 'I',
    Ts <= 'T3' & (is.na(Ns) | Ns <= 'N1') & (is.na(Ms) | Ms <= 'M0') ~ 'II',
    Ns > 'N0' & (is.na(Ms) | Ms <= 'M0') ~ 'III',
    Ms > 'M0' ~ 'IV'
  )
  return(stage)
}


# Calculate the difference between pathology tumor stage and pre-treatment tumor stage as a proxy for response to therapy
prestage.to.path<-function(df) {
  # -1 to set 'unknown' to 0
  RP.TumourStage.Int.c <- as.integer(df$RP.TStage.PrimaryTumour)-1
  PS.TumourStage.Int.c <- as.integer(df$PS.TStage.PrimaryTumour.FinalPretreatmentStaging)-1

  # 0 == no change, negative #s indicate growth, positive indicate shrinkage
  df$RP.PS.Tumour.Diff.c <- PS.TumourStage.Int.c - RP.TumourStage.Int.c
  df$RP.Tumour.Growth.c <- sapply(df$RP.PS.Tumour.Diff.c, function(x) {
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
