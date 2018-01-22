# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

Downloading LabKey information from the OCCAMS database.

### How do I get set up? ###

Install R library from source
* plyr, RCurl, Rlabkey

Create a file with the following fields:

url=https://occams.comlab.ox.ac.uk/labkey
path=/ICGC/Cohorts/All Study Subjects
user=< your user >
pwd=< your password >

### Two basic uses:

These generate slightly different sets of data as the cleaning process merges and removes some fields, as well as adding others necessary for survival analysis.

1. Create a dump of the data contained in all tables on LabKey
2. Generate a cleaned version of the data for analysis

Usage:

ocs <- connect.to.labkey(file='labkey.creds') # default is ~/.labkey.cred

# For (1)
#
dump.clinical.data(ocs, outdir="/tmp") # each Labkey table will be output separately in the outdir

# For (2)
occams <- download.all.tables(ocs)
occams$patients  # patient clinical data table
occams$tissues  # tissue related clinical data table





