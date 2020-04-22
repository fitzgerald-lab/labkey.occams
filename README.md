# README #

This package helps to download LabKey information from the OCCAMS database into a readable format.

### How do I get set up? ###

Install R library from source, necessary packages are:

* plyr
* RCurl
* Rlabkey v. 2.1.136


_IMPORTANT NOTE:_ Until the OCCAMS LabKey server is updated to v.16 or above, the version of the Rlabkey package must NOT be higher than 2.2  .

Create a file with the following fields:

```R
 url=https://occams.cs.ox.ac.uk/labkey/
 path=/ICGC/Cohorts/All Study Subjects
 user=MyUser
 pwd=MyPassword
  ```


## Primary Use

These generate slightly different sets of data compared to what is found in the Labkey UI as the cleaning process merges and removes some fields, as well as adding others that are useful for survival analysis.

All of the labkey tables will be merged into a wide table with a single row per patient. The object returned contains a list with two names element.  'patient' is the semi-curated per-patient clinical information, 'tissues' is the tissue-related clinical information (uncurated currently).  


*Usage:*

Note that the default credential file path is ~/.labkey.cred so if you put it there, no need to provide the file path

```R
library(openclinica.occams)

ocs <- connect.to.labkey(file='labkey.creds') 
occams <- download.wide.format(ocs)

## OR if you have a list of OCCAMS identifiers
occams <- download.wide.format(ocs, occams_ids)

head(occams$patients)
```

# Version 0.3 Changes

- All of the followup/endpoint/recurrence tables are now being integrated.  A single row per-patient is created with the date of death (if know), date of recurrence, and last seen date is generated from the multiple tables.
- All of the new treatment tables are now being merged. Column names are changed, TR.NeoAdj indicates neoadjuvant therapies, TR.Adj indicates adjuvant etc. However, if a patient has 'NA' entries in the adjuvant columns, that does *not* mean they were not treated adjuvantly.  We were not previously recording that separately, so there is no way to look at that information prior to the 2015 CRFs.
- All columns that are created or inferred in the process of cleaning are indicated with a '.c' at the end of the column name: e.g. 'Weeks.Survival.c'
- Some previously created columns no longer exist, including TR.Chemotherapy. Users need to determine this information from the treatment (TR.*) columns as necessary for their analysis.


