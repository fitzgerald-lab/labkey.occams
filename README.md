# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

Downloading LabKey information from the OCCAMS database.

### How do I get set up? ###

IMPORTANT
Until the OCCAMS LabKey server is updated to >v.16 the version of Rlabkey must NOT be >=2.2  

Install R library from source
* plyr, RCurl, Rlabkey v. 2.1.136

Create a file with the following fields:

> url=https://occams.comlab.ox.ac.uk/labkey
> path=/ICGC/Cohorts/All Study Subjects
> user=< your user id >
> pwd=< your password >

### Two basic uses:

These generate slightly different sets of data as the cleaning process merges and removes some fields, as well as adding others necessary for survival analysis.

1. Create a dump of the data contained in all tables on LabKey
2. Generate a cleaned version of the data for analysis

## Usage:

default file is ~/.labkey.cred so if you put it there, no need to provide the file path

> ocs <- connect.to.labkey(file='labkey.creds') 

### For (1)

Each Labkey table will be output separately in the outdir, this is not recommended.

> dump.clinical.data(ocs, outdir="/tmp") 

### For (2)

All of the labkey tables will be merged into a wide table with a single row per patient for all tables. The object returned contains a list with two names element.  'patient' is the semi-curated per-patient clinical information, 'tissues' is the tissue-related clinical information (uncurated currently).  

> occams <- download.all.tables(ocs)
> occams$patients  # patient clinical data table
> occams$tissues  # tissue related clinical data table





