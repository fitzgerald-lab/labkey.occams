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

bmicalc<-function(x) {
  wt <- x[1]
  ht <- x[2]

  if (is.na(wt) || is.na(ht) || wt <= 0 || ht <= 0)
    return(NA)
  return(wt/(ht/100))/(ht/100)
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

trim<-function(x) {
  gsub("^\\s+|\\s+$", "", x)
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

