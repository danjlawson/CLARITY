## Comparing Linguistic phonetic data to linguistic lexical data
##
## Author: Igor Yanovich (igor.yanovich@gmail.com)
## Date: 2020
## Licence: GPLv3 


library("Clarity")


data_header <- read.table("header-NorthEuraLex.tsv", sep="\t", stringsAsFactors=FALSE)
lang_names <- read.table("nelex-0.9-langs-ie.txt", sep="\t", stringsAsFactors=FALSE)  


read_iwsa_data <- function(filename, lang_names, data_header, fix_diagonals=TRUE){
	iwsa_data <- read.table(filename, sep="\t", stringsAsFactors=FALSE)
	colnames(iwsa_data) <- data_header
	
	make_empty_diss_table <- function(entity_names){
		matrix(0, nrow=length(entity_names), ncol=length(entity_names), 
			dimnames=list(entity_names, entity_names))
	}

	phonetic = make_empty_diss_table(lang_names)
	lexical = make_empty_diss_table(lang_names)
	
	for (i in 1:dim(iwsa_data)[1]){
		phonetic[iwsa_data$L1[i], iwsa_data$L2[i]] <- phonetic[iwsa_data$L2[i], iwsa_data$L1[i]] <-  iwsa_data$average_normalized_IWDSC[i]
		lexical[iwsa_data$L1[i], iwsa_data$L2[i]] <- lexical[iwsa_data$L2[i], iwsa_data$L1[i]] <- 1 - iwsa_data$cognate_overlap[i]
	}
	
	if(fix_diagonals){
		phonetic <- c_fixdiagonals(phonetic)
		lexical <- c_fixdiagonals(lexical)
	}
		
	list(phonetic = phonetic, lexical = lexical)
}



actual_data <- read_iwsa_data("samples/northeuralex-0.9-cldf.tsv-mtx-originaldata.tsv", lang_names[,1], data_header)


samples.first.half <- lapply(1:100, function(i){
		filename_prefix = "samples/northeuralex-0.9-cldf.tsv-mtx-sample"
		filename = paste(filename_prefix, formatC(i, width=4, flag="0"), formatC(1, width=4, flag="0"), ".tsv", sep="")
		read_iwsa_data(filename, lang_names[,1], data_header)
	})
	
samples.second.half <- lapply(1:100, function(i){
		filename_prefix = "samples/northeuralex-0.9-cldf.tsv-mtx-sample"
		filename = paste(filename_prefix, formatC(i, width=4, flag="0"), formatC(2, width=4, flag="0"), ".tsv", sep="")
		read_iwsa_data(filename, lang_names[,1], data_header)
	})	
	
###
#PHONETIC FROM LEXICAL
###

lexical_scans <- Clarity_Scan(actual_data$lexical, kmax=dim(lang_names)[1]-1)
phonetic_predicts <- Clarity_Predict(actual_data$phonetic, lexical_scans)	

PhonFromLexList1 = lapply(1:100, function(i) 
	list(Lex1=samples.first.half[[i]]$lexical, 
	Lex2=samples.second.half[[i]]$lexical,
	Phon2=samples.second.half[[i]]$phonetic))
	
PhonFromLexList2 = lapply(1:100, function(i) 
	list(Lex1=samples.second.half[[i]]$lexical, 
	Lex2=samples.first.half[[i]]$lexical,
	Phon2=samples.first.half[[i]]$phonetic))	

PhonFromLexList = do.call(c, list(PhonFromLexList1, PhonFromLexList2))
	
PhonFromLexSignificance=Clarity_Compare(lexical_scans,Ynew=actual_data$phonetic,Ylist=PhonFromLexList,nbs=200,H0="scale") 
	# H0=“scale” applies Procrust transform, otherwise no transformation is done	

###
#LEXICAL FROM PHONETIC
###

phonetic_scans <- Clarity_Scan(actual_data$phonetic, kmax=dim(lang_names)[1]-1)
lexical_predicts <- Clarity_Predict(actual_data$lexical, phonetic_scans)	

LexFromPhonList1 = lapply(1:100, function(i) 
	list(Phon1=samples.first.half[[i]]$phonetic, 
	Phon2=samples.second.half[[i]]$phonetic,
	Lex2=samples.second.half[[i]]$lexical))
	
LexFromPhonList2 = lapply(1:100, function(i) 
	list(Phon1=samples.second.half[[i]]$phonetic, 
	Phon2=samples.first.half[[i]]$phonetic,
	Lex2=samples.first.half[[i]]$lexical))	

LexFromPhonList = do.call(c, list(LexFromPhonList1, LexFromPhonList2))
	
LexFromPhonSignificance=Clarity_Compare(phonetic_scans,Ynew=actual_data$lexical,Ylist=LexFromPhonList,nbs=200,H0="scale") 
	# H0=“scale” applies Procrust transform, otherwise no transformation is done	




###
#PLOTTING: ACTUAL DATA, THEN PERSISTENCE DIAGRAMS
###

pdf("Fig5-ling-example.pdf", height=12, width=12)

# layout(mat=matrix(1:4, nrow = 2, ncol = 2),
# 	heights = c(1, 2),
#     widths = c(2, 2)
#     )

par(mfrow=c(2,2))

Clarity_Chart(actual_data$phonetic,las=2,main="a) Phonetic dissimilarity", font.main = 1)
Clarity_Chart(actual_data$lexical,las=2,main="b) Lexical dissimilarity", font.main = 1)

plot(LexFromPhonSignificance, thresh=0.05, text=TRUE, cex.text=0.2, main="c) Persistences for Lexical predicted from Phonetic", xlab="Languages (ISO 639-3 codes)", ylab="K", las=2, font.main = 1)

plot(PhonFromLexSignificance, thresh=0.05, text=TRUE, cex.text=0.2, main="d) Persistences for Phonetic predicted from Lexical", xlab="Languages (ISO 639-3 codes)", ylab="K", las=2, font.main = 1)
	
dev.off()






###
#PLOTTING: PERSISTENCE DIAGRAMS WITH P=0.01 FOR SM
###

pdf("SFig1-ling-persistences-at-p-001.pdf", height=6, width=12)
par(mfrow=c(1,2))

plot(LexFromPhonSignificance, thresh=0.01, text=TRUE, cex.text=0.2, main="a) Persistences for Lexical predicted from Phonetic, p=0.01", xlab="Languages (ISO 639-3 codes)", ylab="K", las=2, font.main = 1)

plot(PhonFromLexSignificance, thresh=0.01, text=TRUE, cex.text=0.2, main="b) Persistences for Phonetic predicted from Lexical, p=0.01", xlab="Languages (ISO 639-3 codes)", ylab="K", las=2, font.main = 1)

dev.off()



###
#PLOTTING: RESAMPLED PHONETIC FROM LEXICAL PERSISTENCE BY LANGUAGE FOR SM 
###

PlotPersistenceCurvesByLanguage <- function(clarity_significance, label_orig, label_target){

	p0=Clarity_Persistence(clarity_significance$pred)
	p1=lapply(clarity_significance$splitlist,function(x)x$pdata)
	p2=lapply(clarity_significance$splitlist,function(x)x$ppred)

	for(language in rownames(p0)){
		tmpp1=t(sapply(p1,function(x)x[language,]))
		tmpp2=t(sapply(p2,function(x)x[language,]))
		tmppp=clarity_significance$pvals[language,]
		plot(p0[language,],log="y",type="n",
			main=paste(label_target, "from", label_orig, "persistence for",language),
			xlab="CLARITY complexity K", ylab="persistence P")
		for(i in 1:200) {
			lines(tmpp1[i,],col="burlywood")
			lines(tmpp2[i,],col="darkorchid1")
		}
		lines(p0[language,],lwd=2)
		points(p0[language,],cex=1.5,
			   pch=c(19,4,1)[1+(tmppp>.05)+(tmppp>.01)])
		legend("topright",legend=c(paste("P for full",  label_target, "(significant 0.01)"),
									paste("P for full",  label_target, "(significant 0.05)"),
								  	paste("P for full",  label_target, "(not significant)"),
								   	paste("P for an unseen resample from orig data (", label_orig,")", sep=""),
								   	paste("P for a resample from pred data (", label_target,")", sep="")),
			  						col=c(1,1,1,"burlywood","darkorchid1"),
			  						lty=c(1,1,1,1,1),pch=c(19,4,1,NA,NA))
	}

}



pdf("SFig2-persistence-curves-by-language.pdf", height=6, width=12)
PlotPersistenceCurvesByLanguage(LexFromPhonSignificance, label_orig="Phonetic", label_target="Lexical")
PlotPersistenceCurvesByLanguage(PhonFromLexSignificance, label_orig="Lexical", label_target="Phonetic")
dev.off()

