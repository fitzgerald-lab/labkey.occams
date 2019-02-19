mr<-function(a,b,i=1,empty=NA) {
  #print(paste('a',a))
  #print(paste('b',b))
  c = unique(na.omit(c(a[i],b[i])))
  #print(paste('c',is_empty(c)))
  if (is_empty(c)) c = empty
  return(c[i])
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

bmicalc<-function(wt,ht) {
  if (is.na(wt) || is.na(ht) || wt <= 0 || ht <= 0)
    return(NA)
  return (wt/(ht/100))/(ht/100)
}

trim<-function(x) {
  gsub("^\\s+|\\s+$", "", x)
}

recode.siewert<-function(df) {
  grp = group_vars(df)
  df = df %>% ungroup %>% mutate_at(vars(contains('SiewertClassification')), funs(recode_factor), '1'='Type I','2'='Type II','3'='Type III', .ordered=F, .missing=NA_character_ )
  if (length(grp) > 0) df = df %>% group_by(!!grp)
  return(df)
}

recode.TNM<-function(df) {
  grp = group_vars(df)
  df = df %>% ungroup %>% mutate_at(vars(contains('TNMStage')), funs(factor), levels=c('I','II','III', 'IV'), ordered=T)
  if (length(grp) > 0) df = df %>% group_by(!!grp)
  return(df)
}


tnmStage <- function(Ts,Ns,Ms) {
#  print(paste(Ts,Ns,Ms))

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


