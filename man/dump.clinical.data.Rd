% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_oc_data.R
\name{dump.clinical.data}
\alias{dump.clinical.data}
\title{Download all tables directly from Labkey for all or a subset of patients.}
\usage{
dump.clinical.data(
  ocs,
  occams_ids = NULL,
  prefixes = NULL,
  version = "z1",
  outdir = "~/tmp"
)
}
\arguments{
\item{ocs}{labkey connection object from connect.to.labkey()}

\item{occams_ids}{Array of occams identifiers, this is optional. The default returns all patients, better to use get.patients(...) if you want this.}

\item{prefixes}{Get only tables with the given prefixes. Call 'get.prefixes(ocs)' for a list of only CRF tables.}

\item{version}{Table version (default 'z1')}
}
\description{
Download all tables directly from Labkey for all or a subset of patients.
}
\author{
skillcoyne
}
