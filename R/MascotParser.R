
##############################################################################################################
# Package name: MascotParser packages in R language using S4 class
# Function: parse Mascot results file (*.dat) from both PMF and MSMS searches; support few functions to extract 
# useful information from the search results such as: positive hits, assigned (positive hit) peaks, etc..; 
# plot spectra in both MS1 and MS2 levels, etc..
# authors: Trung Nghia VU
# email: nghiavtr@gmail.com
# updated date: 14/08/2014
##############################################################################################################

##### Definition of MascotParser class
setClass("MascotParser",slots = c(filename = "character", remarks = "character", 
	boundaryCode="character",sourcefile="character",searchType="character",fastafile="character",database="character",
	enzyme="character",taxonomy="character",header="character",parameters="character",
	qmass="numeric", qexp="numeric",qcharge="character",qintensity="numeric",qmatch="numeric",qplughole="numeric",
	PIHaccessions="character",PIHhitScores = "numeric",PIHproteinMasses = "numeric",PIHnumHits = "numeric",PIHpeakInfo = "list",
	MSMSqueryHits="list",MSMSproteinInfo="character",MSMSqueryInfo="list",MSMSnumQuery="numeric")
	);

##############################################################################################################
##### Definition methods fo MascotParser

#parse dat file both PMF and MSMS

setGeneric("parse",function(object,filename){
	standardGeneric("parse")
	});

setMethod("parse","MascotParser", function(object,filename){

	fi  <- file(filename, open = "r");
	#general information
	boundaryCode=NULL;
	sourcefile=NULL;
	searchType=NULL;
	database=NULL;
	fastafile=NULL;
	enzyme=NULL;
	taxonomy=NULL;
	header=NULL;
	parameters=NULL;

	#general summary information
	qmass = NULL;
	qexp = NULL;
	qcharge = NULL;
	qintensity = NULL;
	qmatch = NULL;
	qplughole=NULL;
	
	#hits information (for PMF)
	PIHaccessions = NULL;
	PIHhitScores = NULL;
	PIHproteinMasses = NULL;
	PIHnumHits = NULL;
	PIHpeakInfo = NULL;

	# peptide hits and protein hits information (for MSMS)
	MSMSqueryHits = NULL;	
	MSMSproteinInfo = NULL;
	MSMSqueryInfo = NULL;
	MSMSnumQuery = NULL;


	#find boundary code of the file
	while (length(aLine <- readLines(fi, n = 1, warn = FALSE)) > 0) {
		if (length(grep("boundary=",aLine)) > 0){
			boundaryCode = strsplit(aLine,"=")[[1]][2];
			break;
		}
	} #end of while

	### start parsing the dat file

	#extract the information of parameters
	myflag = 0;
	while (length(aLine <- readLines(fi, n = 1, warn = FALSE)) > 0) {
		if (length(grep("Content-Type",aLine)) > 0 & length(grep("parameters",aLine)) > 0){
			while (length(aSubLine <- readLines(fi, n = 1, warn = FALSE)) > 0) {
				if (length(grep(boundaryCode,aSubLine)) > 0){ myflag = 1; break;}				
				if (nchar(aSubLine)) parameters<-paste(aSubLine,parameters,sep=";");
				if (length(grep("DB=",aSubLine)) > 0) database = strsplit(aSubLine,"=")[[1]][2];
				if (length(grep("FILE=",aSubLine)) > 0) sourcefile = strsplit(aSubLine,"=")[[1]][2];
				if (length(grep("SEARCH=",aSubLine)) > 0) searchType = strsplit(aSubLine,"=")[[1]][2];
			}
		}
	if (myflag == 1)break;
	}#end of while

	#extract the information of enzyme
	myflag = 0;
	while (length(aLine <- readLines(fi, n = 1, warn = FALSE)) > 0) {
		if (length(grep("Content-Type",aLine)) > 0 & length(grep("enzyme",aLine)) > 0){
			while (length(aSubLine <- readLines(fi, n = 1, warn = FALSE)) > 0) {
				if (length(grep(boundaryCode,aSubLine)) > 0){ myflag = 1; break;}
				if (nchar(aSubLine)) enzyme<-paste(aSubLine,enzyme,sep=";");
			}
		}
	if (myflag == 1)break;
	}#end of while


	#extract the information of taxonomy
	if (searchType != "PMF"){
		myflag = 0;
		while (length(aLine <- readLines(fi, n = 1, warn = FALSE)) > 0) {
			if (length(grep("Content-Type",aLine)) > 0 & length(grep("taxonomy",aLine)) > 0){
				while (length(aSubLine <- readLines(fi, n = 1, warn = FALSE)) > 0) {
					if (length(grep(boundaryCode,aSubLine)) > 0){ myflag = 1; break;}
					if (nchar(aSubLine)) taxonomy<-paste(aSubLine,taxonomy,sep=";");
				}
			}
			if (myflag == 1)break;
		}#end of while
	}else{
		taxonomy="";
	}
	#extract the information of header
	myflag = 0;
	while (length(aLine <- readLines(fi, n = 1, warn = FALSE)) > 0) {
		if (length(grep("Content-Type",aLine)) > 0 & length(grep("header",aLine)) > 0){
			while (length(aSubLine <- readLines(fi, n = 1, warn = FALSE)) > 0) {
				if (length(grep(boundaryCode,aSubLine)) > 0){ myflag = 1; break;}
				if (nchar(aSubLine)) header<-paste(aSubLine,header,sep=";");
				if (length(grep("fastafile=",aSubLine)) > 0) fastafile = strsplit(aSubLine,"=")[[1]][2];
			}
		}
		if (myflag == 1)break;
	}#end of while
	

	#extract masses
	myflag = 0;
	while (length(aLine <- readLines(fi, n = 1, warn = FALSE)) > 0) {
		if (length(grep("Content-Type",aLine)) > 0 & length(grep("summary",aLine)) > 0){
			while (length(aSubLine <- readLines(fi, n = 1, warn = FALSE)) > 0) {
				if (length(grep(boundaryCode,aSubLine)) > 0){ myflag = 1; break;}
				if (length(grep("num_hits",aSubLine)) > 0){ 
					PIHnumHits = as.integer(strsplit(aSubLine,"=")[[1]][2]);
					myflag = 1; 
					break;
				}
				if (length(grep("qmass",aSubLine)) > 0){
					qmass <- c(qmass, strsplit(aSubLine,"=")[[1]][2])
					aSubLine <- readLines(fi, n = 1, warn = FALSE);
					tmpStr <- strsplit(strsplit(aSubLine,"=")[[1]][2],",");
					qexp <- c(qexp, tmpStr[[1]][1]);
					qcharge <- c(qcharge, tmpStr[[1]][2]);
					aSubLine <- readLines(fi, n = 1, warn = FALSE);
					qintensity <- c(qintensity, strsplit(aSubLine,"=")[[1]][2])
					aSubLine <- readLines(fi, n = 1, warn = FALSE);
					qmatch <- c(qmatch, strsplit(aSubLine,"=")[[1]][2])				 
				}
			}#end of while
		if (myflag == 1) break;
		}
	}#end of while
	# finish reading masses


	# start reading protein identification hits (PIH) information
	if (PIHnumHits > 0){
	myflag = 0;
	numPeaks = length(qmass);
	hitCount = 1;
	while (length(aLine <- readLines(fi, n = 1, warn = FALSE)) > 0) {	
		keyStr = paste("h",hitCount,"=",sep = "");
		if (length(grep(keyStr,aLine)) > 0){
			infoArr = strsplit(aLine,",");
			PIHaccessions <- c(PIHaccessions,strsplit(infoArr[[1]][1],"=")[[1]][2]);
			PIHhitScores <- c(PIHhitScores,as.double(infoArr[[1]][2]));
			PIHproteinMasses <- c(PIHproteinMasses,as.double(infoArr[[1]][4]));
			
			#extract information of each peak
			PIHpeakInfoHit = NULL;
			peakCount = 1;
			while (length(aSubLine <- readLines(fi, n = 1, warn = FALSE)) > 0) {
				if (length(grep(boundaryCode,aSubLine)) > 0){ myflag = 1; break;}

				subkeyStr = paste("h",hitCount,"_q",peakCount,sep = "");
				if (length(grep(subkeyStr,aSubLine)) > 0){
					PIHpeakInfoHit <- c(PIHpeakInfoHit,strsplit(aSubLine,"=")[[1]][2]);
					peakCount = peakCount + 1;
				}
				if (peakCount > numPeaks) break;
			}			
			PIHpeakInfo[[hitCount]] = PIHpeakInfoHit;
			hitCount = hitCount + 1;
		}		
		
		if (hitCount > PIHnumHits) break;
		if (myflag == 1) break;
		}
	}		
	# finish reading protein identification hits (PIH) information
	
	# start to read peptide hits of each spectra (only for MSMS)
	if (searchType != "PMF"){
	myflag = 0;
#	MSMSqueryHits[[2]][[2]]<-NA;#initial 
	while (length(aLine <- readLines(fi, n = 1, warn = FALSE)) > 0) {	
		if (length(grep("Content-Type",aLine)) > 0 & length(grep("peptides",aLine)) > 0){
		aLine <- readLines(fi, n = 1, warn = FALSE)# remove a blank row			
			while (length(aSubLine <- readLines(fi, n = 1, warn = FALSE)) > 0) {
				if (length(grep(boundaryCode,aSubLine)) > 0){ myflag = 1; break;}

				tmpStr <- strsplit(aSubLine,"=");
				tmpStr2 <- strsplit(tmpStr[[1]][1],"_");
				qVal <- as.integer(substring(tmpStr2[[1]][1],2));
				pVal <- as.integer(substring(tmpStr2[[1]][2],2));
				
				if (length(MSMSqueryHits) < qVal){ # start a new query information
				
					MSMSqueryHits<- c(MSMSqueryHits,parseQueryhit(new("QueryHits"),tmpStr[[1]][2]));
				}
				else{
					if (pVal == length(MSMSqueryHits[[qVal]])){
						if (length(MSMSqueryHits[[qVal]])>1){				
							MSMSqueryHits[[qVal]][[pVal]]<-addQueryHitInfo(MSMSqueryHits[[qVal]][[pVal]],tmpStr[[1]][2])
						}
						else
							MSMSqueryHits[[qVal]]<-addQueryHitInfo(MSMSqueryHits[[qVal]],tmpStr[[1]][2])
					}
					else{ 
						MSMSqueryHits[[qVal]]<- c(MSMSqueryHits[[qVal]],parseQueryhit(new("QueryHits"),tmpStr[[1]][2]));
					}
				 } #end of else
				
			}	#end of while		
		}	#end of if
		if (myflag == 1) break;
	}#end of while
	}
	#finish reading hits of spectra
	
	######## ignore reading hits of decoy
	
	#start reading the list of proteins (only for MSMS)
	if (searchType != "PMF"){
	myflag = 0;
	while (length(aLine <- readLines(fi, n = 1, warn = FALSE)) > 0) {
		if (length(grep("Content-Type",aLine)) > 0 & length(grep("proteins",aLine)) > 0){
			aLine <- readLines(fi, n = 1, warn = FALSE)# remove a blank row
			#read protein hits
			while (length(aSubLine <- readLines(fi, n = 1, warn = FALSE)) > 0) {
				if (length(grep(boundaryCode,aSubLine)) > 0){ myflag = 1; break;}
				if (length(aLine) >0 ) MSMSproteinInfo <- c(MSMSproteinInfo,aSubLine);
			}#end of while
		}
		if (myflag == 1)break;
	}#end of while
	}
	#finish reading the list of proteins (only for MSMS)
	
	#start reading queries's information (only for MSMS)
	if (searchType != "PMF"){
	MSMSqueryInfo <- NULL;
	MSMSnumQuery = length(MSMSqueryHits);
	queryCount = 0;
	while (length(aLine <- readLines(fi, n = 1, warn = FALSE)) > 0) {
		if (length(grep("Content-Type",aLine)) > 0 & length(grep("query",aLine)) > 0){
			queryCount = queryCount+1;
#			cat(queryCount,"\n");
			aQuery <- NULL;
			while (length(aSubLine <- readLines(fi, n = 1, warn = FALSE)) > 0) {
				if (length(grep(boundaryCode,aSubLine)) > 0){ myflag = 1; break;}
				if (length(grep("title=",aSubLine)) > 0) aQuery$title = strsplit(aSubLine,"=")[[1]][2];
				if (length(grep("rtinseconds=",aSubLine)) > 0) aQuery$rtinseconds = as.double(strsplit(aSubLine,"=")[[1]][2]);
				if (length(grep("charge=",aSubLine)) > 0) aQuery$charge = strsplit(aSubLine,"=")[[1]][2];
				if (length(grep("mass_min=",aSubLine)) > 0) aQuery$mass_min = as.double(strsplit(aSubLine,"=")[[1]][2]);
				if (length(grep("mass_max=",aSubLine)) > 0) aQuery$mass_max = as.double(strsplit(aSubLine,"=")[[1]][2]);
				if (length(grep("int_min=",aSubLine)) > 0) aQuery$int_min = as.double(strsplit(aSubLine,"=")[[1]][2]);
				if (length(grep("int_max=",aSubLine)) > 0) aQuery$int_max = as.double(strsplit(aSubLine,"=")[[1]][2]);
				if (length(grep("num_vals=",aSubLine)) > 0) aQuery$num_vals = as.double(strsplit(aSubLine,"=")[[1]][2]);
				if (length(grep("num_used1=",aSubLine)) > 0) aQuery$num_used1 = as.double(strsplit(aSubLine,"=")[[1]][2]);
				if (length(grep("Ions1=",aSubLine)) > 0){
				 	aQuery$Ions1 = strsplit(aSubLine,"=")[[1]][2];
				 	ionsVec <- strsplit(aQuery$Ions1,",")[[1]];
				 	res <- sapply(ionsVec,function(x){
				 			as.double(strsplit(x,":")[[1]]);
				 		}
				 	)
				 	colnames(res) <- NULL;
				 	aQuery$mz <- res[1,];
				 	aQuery$intensity <- res[2,];
				 }				
			}#end of while			
			MSMSqueryInfo[[queryCount]] <- aQuery;
		}#end of if
		if (queryCount == MSMSnumQuery) break;
	}#end of while
	}
	#finish queries 's information (only for MSMS)
	
	close(fi)	

	#set values to the object
	object@filename=filename;
	object@boundaryCode=boundaryCode;
	object@sourcefile=sourcefile;
	object@searchType=searchType;
	object@fastafile=fastafile;
	object@database=database;
	object@enzyme=enzyme;
	object@taxonomy=taxonomy;
	object@header=header;
	object@parameters=parameters;


	object@qmass=as.numeric(qmass);
	object@qexp=as.numeric(qexp);
	object@qcharge=qcharge;
	object@qintensity=as.numeric(qintensity);
	object@qmatch=as.numeric(qmatch);
	object@qplughole=as.numeric(qplughole);

	if (!is.null(PIHaccessions)) object@PIHaccessions = PIHaccessions;
	if (!is.null(PIHhitScores)) object@PIHhitScores = PIHhitScores;
	if (!is.null(PIHproteinMasses)) object@PIHproteinMasses = PIHproteinMasses;
	if (!is.null(PIHnumHits)) object@PIHnumHits = PIHnumHits;
	if (!is.null(PIHpeakInfo)) object@PIHpeakInfo = PIHpeakInfo;

	if (!is.null(MSMSqueryHits)) object@MSMSqueryHits=MSMSqueryHits;
	if (!is.null(MSMSproteinInfo)) object@MSMSproteinInfo=MSMSproteinInfo;
	if (!is.null(MSMSqueryInfo)) object@MSMSqueryInfo=MSMSqueryInfo;
	if (!is.null(MSMSnumQuery)) object@MSMSnumQuery=MSMSnumQuery;
	
	
	return(object);
	
	});

#plot spectrum
setGeneric("plotSpec",function(object,specTitle=NULL,queryId=NULL){
	standardGeneric("plotSpec")
	});

setMethod("plotSpec","MascotParser", function(object,specTitle=NULL,queryId=NULL){
#	specTitle <- "";
	if (is.null(specTitle)) specTitle<-basename(object@filename);

	if (is.null(queryId)){
		#plot MS1 spectrum
		if (length(object@qmass) > 0){
			plot.default(object@qmass,object@qintensity,type='h',xlab="m/z",ylab="intensity",main=paste("MS1 spectrum of ",specTitle,sep=""))
		}else cat("Error, there is no data");
			
	}else{
		#plot MS2 spectrum of the queryId
		if (object@searchType!="PMF" && length(object@MSMSqueryInfo[[queryId]]$mz) > 0){
			plot.default(object@MSMSqueryInfo[[queryId]]$mz,object@MSMSqueryInfo[[queryId]]$intensity,type='h',xlab="m/z",ylab="intensity",main=paste("MS2 spectrum of ",specTitle, " with Id ",queryId,sep=""))
		}else cat("Error, there is no MS2 data of ",queryId);
	}
})

#get MS2 information
setGeneric("getMSMSqueryInfo",function(object,queryId=NULL){
	standardGeneric("getMSMSqueryInfo")
	});

setMethod("getMSMSqueryInfo","MascotParser", function(object,queryId=NULL){
#	specTitle <- "";
	return(object@MSMSqueryInfo[[queryId]]);
	
})


#get positive hits
setGeneric("getPosHits",function(object, scoreThres=41){
	standardGeneric("getPosHits")
	});

setMethod("getPosHits","MascotParser", function(object,scoreThres=41){
	if (object@searchType=="PMF"){
	return(which(object@PIHhitScores >=scoreThres));
	}else{
		cat("case of MSMS will be implemented soon \n");

	}
	
})


#get assigned peaks from possitive hits
setGeneric("getAssignedPeaks",function(object, scoreThres=41){
	standardGeneric("getAssignedPeaks")
	});

setMethod("getAssignedPeaks","MascotParser", function(object,scoreThres=41){
	if (object@searchType=="PMF"){
		posHits <- getPosHits(object,scoreThres);
		posPeakId <- NULL;
		if (length(posHits) > 0){
		for (i in 1:length(posHits)){
			posPeakId <- c(posPeakId,which(object@PIHpeakInfo[[posHits[i]]]!="-1"));
		}
		return(unique(posPeakId));
		}

	}else{
		cat("case of MSMS will be implemented soon \n");

	}
	
})


#get unassigned peaks from the spectrum
setGeneric("getUnassignedPeaks",function(object, scoreThres=41){
	standardGeneric("getUnassignedPeaks")
	});

setMethod("getUnassignedPeaks","MascotParser", function(object,scoreThres=41){
	if (object@searchType=="PMF"){
		posPeakId <- getAssignedPeaks(object,scoreThres);		
		return(c(1:length(myParse@qmass))[-posPeakId]);
	}else{
		cat("case of MSMS will be implemented soon \n");

	}
	
})

#get information by keywords
setGeneric("getKeywordsInfo",function(object,keywords=NULL){
	standardGeneric("getKeywordsInfo")
	});

setMethod("getKeywordsInfo","MascotParser", function(object,keywords=NULL){
	if (is.null(object@filename)){
		cat("There is no filename in the object");
		return;
	}

	if (length(keywords) < 1) {
		cat("The keyword list is empty");
		return;
	}
	# parse the file to get the information of the keywords
	keywordValues <- rep(NA,length(keywords));
	keywordIds <- c(1:length(keywords));
	fi  <- file(filename, open = "r");
	while (length(aLine <- readLines(fi, n = 1, warn = FALSE)) > 0) {
		if (length(keywords)<1)	break;
		for (i in 1:length(keywords)){
			if (length(grep(keywords[i],aLine)) > 0){
				keywordValues[keywordIds[i]] = strsplit(aLine,"=")[[1]][2];
				keywords <- keywords[-i];
				keywordIds <- keywordIds[-i];
				break;
			}
		}
	} #end of
	close(fi);
	return(keywordValues);
})


#get input filename
setGeneric("setInputfilename",function(object,filename=NULL){
	standardGeneric("setInputfilename")
	});

setMethod("setInputfilename","MascotParser", function(object,filename=NULL){
	object@filename=filename;
})


#plot two spectra in panel
setGeneric("plotSpecPanel",function(mz1,int1,mz2,int2,plotTitle=NULL,xlimVal=NULL){
	standardGeneric("plotSpecPanel")
	});

setMethod("plotSpecPanel","MascotParser", function(mz1,int1,mz2,int2,plotTitle=NULL,xlimVal=NULL){
	if (is.null(xlimVal)){
		xlimVal<-max(c(mz1,mz2));
		xlimVal<-xlimVal*1.05;
	}	

	def.par <- par(no.readonly = TRUE,mar = c(0,0,2,0), xaxs='i', yaxs='i')
 	layout(matrix(c(1,1,1,1,1,1,2,2,2, 2, 2, 2), 4, 3, byrow = TRUE))
  	matrix(c(1,1,1,1,1,1,2,2,2, 2, 2, 2), 4, 3, byrow = TRUE);
	# top-right 
	plot(mz1,int1, type="h", ylab = "", xlab = "", xlim = c(0,xlimVal), col="red",xaxt = "n",main=plotTitle) #
	# bottom-right
	par(mar = c(2,0,0,0))
	plot(mz2,int2, ylim = rev(range(int2)), type="h", xlim = c(0,xlimVal),col = "blue")
	par(def.par)
});
