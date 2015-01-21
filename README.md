mascotParser
============

A parser for Mascot search results ( both PMF and MSMS) in R package.


##1. Example for parsing MSMS Mascot search result files:
```R
library(MascotParser)
filename <- "MSMS_Mascot_search_result.dat";
```
#####parse the dat file
```R
myParse <- parse(new("MascotParser"),filename)
```
##### Get MS1 peaks
```R
myParse@qexp
myParse@qintensity
```
##### how many query of MS2 in this data set?
```R
myParse@MSMSnumQuery
```
##### how many hits and hit scores from these query?
```R
myParse@PIHnumHits
myParse@PIHhitScores
```
##### Extract information of MS2 peaks from query 3
```R
myms2=getMSMSqueryInfo(myParse,3)
myms2$mz
myms2$intensity
```
###### plot MS1 mass spetrum
```R
plotSpec(myParse,specTitle=NULL,queryId=NULL)
```
##### plot mass spectrum of query 1
```R
plotSpec(myParse,specTitle=NULL,queryId=1)
```
##### draw a panel plot of MS2 peaks from query 1 and query 2
```R
mz1=getMSMSqueryInfo(myParse,1)$mz
int1=getMSMSqueryInfo(myParse,1)$intensity
mz2=getMSMSqueryInfo(myParse,2)$mz
int2=getMSMSqueryInfo(myParse,2)$intensity
plotSpecPanel(mz1,int1,mz2,int2,xlimVal=500,plotTitle="test title")
```
##2. Example for parsing PMF Mascot search result files:
```R
library(MascotParser)
filename <- "PMF_Mascot_search_result.dat";
```
###### parse the dat file
```
myParse <- parse(new("MascotParser"),filename)
```
###### get unassigned peaks from the search
```R
negPeaks <- myParse@qexp[getUnassignedPeaks(myParse, scoreThres=MascotThres)]
```
###### plot PMF mass spetrum
```R
plotSpec(myParse,specTitle=NULL,queryId=NULL)
```
