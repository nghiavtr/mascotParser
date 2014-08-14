
##############################################################################################################
# Package name: QueryHits packages in R language using S4 class
# Function: parse Mascot results file (*.dat) from both PMF and MSMS searches; support few functions to extract 
# useful information from the search results such as: positive hits, assigned (positive hit) peaks, etc..; 
# plot spectra in both MS1 and MS2 levels, etc..
# authors: Trung Nghia VU
# email: nghiavtr@gmail.com
# updated date: 14/08/2014
##############################################################################################################

##############################################################################################################
##### Definition of QueryHits class
setClass("QueryHits",slots = c(queryHitText = "character", peptideSequence = "character", peptideMass="numeric",
	hitScore="numeric")
	);

##### Definition methods fo QueryHits

#parse a query hit string
setGeneric("parseQueryhit",function(object,queryHitText){
	standardGeneric("parseQueryhit")
	});

setMethod("parseQueryhit","QueryHits", function(object,queryHitText){

	object@queryHitText=queryHitText;
	
	tmpStr <- strsplit(queryHitText,",");
	object@peptideSequence=tmpStr[[1]][5];
	object@peptideMass=as.double(tmpStr[[1]][2]);
	object@hitScore=as.double(tmpStr[[1]][8]);
	
	return(object);	
	});

# add more information of the query hit 
setGeneric("addQueryHitInfo",function(object,queryHitText){
	standardGeneric("addQueryHitInfo")
	});

setMethod("addQueryHitInfo","QueryHits", function(object,queryHitText){
	object@queryHitText=paste(object@queryHitText,queryHitText, sep=";");	
	return(object);	
	});
