# README #

This package helps to download LabKey information from the OCCAMS database into a readable format.

### How do I get set up? ###

Install R library from source, necessary packages are:

* plyr
* RCurl
* Rlabkey v. 2.1.136


_IMPORTANT NOTE:_ Until the OCCAMS LabKey server is updated to v.16 or above, the version of the Rlabkey package must NOT be higher than 2.2  .

Create a file with the following fields:

```
 url=https://occams.comlab.ox.ac.uk/labkey
 path=/ICGC/Cohorts/All Study Subjects
 user=MyUser
 pwd=MyPassword
 
 ```


## Primary Use

These generate slightly different sets of data compared to what is found in the Labkey UI as the cleaning process merges and removes some fields, as well as adding others that are useful for survival analysis.

All of the labkey tables will be merged into a wide table with a single row per patient. The object returned contains a list with two names element.  'patient' is the semi-curated per-patient clinical information, 'tissues' is the tissue-related clinical information (uncurated currently).  


*Usage:*

Note that the default credential file path is ~/.labkey.cred so if you put it there, no need to provide the file path

```

ocs <- connect.to.labkey(file='labkey.creds') 
occams <- download.all.tables(ocs)

head(occams$patients)
head(occams$tissues)

```



