mascotParser
============

A parser for Mascot search results ( both PMF and MSMS) in R package.


1. Example for parsing MSMS Mascot search result files:

library(MascotParser)
filename <- "MSMS_Mascot_search_result.dat";

#parse the dat file
myParse <- parse(new("MascotParser"),filename)

# Get MS1 peaks
myParse@qexp
myParse@qintensity

# how many query of MS2 in this data set?
myParse@MSMSnumQuery

# how many hits and hit scores from these query?
myParse@PIHnumHits
myParse@PIHhitScores

#Extract information of MS2 peaks from query 3

myms2=getMSMSqueryInfo(myParse,3)
myms2$mz
myms2$intensity

#plot mass spectrum of query 1
plotSpec(myParse,specTitle=NULL,queryId=1)

#draw a panel plot of MS2 peaks from query 1 and query 2

mz1=getMSMSqueryInfo(myParse,1)$mz
int1=getMSMSqueryInfo(myParse,1)$intensity
mz2=getMSMSqueryInfo(myParse,2)$mz
int2=getMSMSqueryInfo(myParse,2)$intensity

plotSpecPanel(mz1,int1,mz2,int2,xlimVal=500,plotTitle="test title")


2. Example for parsing PMF Mascot search result files:

library(MascotParser)
filename <- "PMF_Mascot_search_result.dat";

#parse the dat file
myParse <- parse(new("MascotParser"),filename)

#get unassigned peaks from the search
negPeaks <- myParse@qexp[getUnassignedPeaks(myParse, scoreThres=MascotThres)]

#plot PMF mass spetrum
plotSpec(myParse,specTitle=NULL,queryId=NULL)



