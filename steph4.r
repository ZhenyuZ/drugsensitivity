setwd("/Users/zhenyu/Downloads/working/Steph/ov")
library(httr)
library(data.table)
library(ggplot2)
library(GGally)
library(pRRophetic)
source("../TCGA_Module.r")

options(stringsAsFactors = FALSE)


###################################################
# Define OV metadata file in TCGA
###################################################

clin.drug.link <- "nationwidechildrens.org_clinical_drug_ov.txt"
clin.sample.link <- "nationwidechildrens.org_biospecimen_sample_ov.txt"
clin.radiation.link <- "nationwidechildrens.org_clinical_radiation_ov.txt"
clin.patient.link <- "nationwidechildrens.org_clinical_patient_ov.txt"
clin.omf.link <- "nationwidechildrens.org_clinical_omf_v4.0_ov.txt"
expr.4502A.sdrf.link <- "unc.edu_OV.AgilentG4502A_07_3.sdrf.txt"
expr.U133.sdrf.link <- "broad.mit.edu_OV.HT_HG-U133A.sdrf.txt"
expr.rnaseqV2.sdrf.link <- "unc.edu_OV.IlluminaHiSeq_RNASeqV2.1.1.0.sdrf.txt"


###################################################
# Read metadata
###################################################

drug.tbl <- GetTCGATable(clin.drug.link)
sample.tbl <- GetTCGATable(clin.sample.link)
patient.tbl <- GetTCGATable(clin.patient.link)
omf.tbl <- GetTCGATable(clin.omf.link)
expr.4502A.tbl <- GetTCGATable(expr.4502A.sdrf.link)
expr.U133.tbl <- GetTCGATable(expr.U133.sdrf.link)
expr.rnaseqV2.tbl <- GetTCGATable(expr.rnaseqV2.sdrf.link)


###################################################
# Process drug information and generate drug.use
###################################################

# In drug.use table, TRUE means the drug has been used; FALSE means the BCR record does not contain this drug; NA means either no report of this patient available, or it's not certain whether the drug has been used (drug name be [Not Available]  

# Define name alias/typo
taxol.alias = c("taxol", "paclitaxel", "abraxane", "xyotax", "paclitaxel; albumin-bount", "taxol/carboplatin", "paciltaxel", "paciltaxle", "paclitaxil", "pacliltaxel", "pacliatxel", "paciltaxal", "pacilatxel", "pacitaxel")
taxotere.alias = c("taxotere", "docetaxel", "doxetaxel", "taxoterecin")
cisplatin.alias = c("cisplatin", "cisplatinum", "cisplatin #2-7", "ccdp", "ciplastin")
carboplatin.alias = c("carboplatin", "paraplatin", "carboplatinum", "carboplatin 6th", "carboplatin #1", "taxol/carboplatin", "carbplatin", "carbopaltin", "carbobplatin", "carboplain", "carbo")
taxo.platin.alias = c(taxol.alias, taxotere.alias, cisplatin.alias, carboplatin.alias)
taxol.contained = taxol.alias
taxotere.contained = taxotere.alias 
cisplatin.contained = c(cisplatin.alias, "cisplatin/gemzar", "plfe") 
carboplatin.contained = c(carboplatin.alias, "topotecan/carboplatin")

# Initial drug.use
patients = unique(drug.tbl$bcr_patient_barcode)
drug.use = data.frame(matrix(F, length(patients), 5))
colnames(drug.use) = c("taxol", "taxotere", "carboplatin", "cisplatin", "other")
rownames(drug.use) = patients

# Generate drug.use
for (i in 1: nrow(drug.tbl)) {
	drug = tolower(drug.tbl[i,]$pharmaceutical_therapy_drug_name)
	patient = drug.tbl[i,]$bcr_patient_barcode
	if(drug %in% taxol.contained) drug.use[patient, "taxol"] = T
	if(drug %in% taxotere.contained) drug.use[patient, "taxotere"] = T
	if(drug %in% cisplatin.contained) drug.use[patient, "cisplatin"] = T
	if(drug %in% carboplatin.contained) drug.use[patient, "carboplatin"] = T
	if(!drug %in% taxo.platin.alias) {
		if(drug == tolower("[Not Available]")) {
			drug.use[patient, ] = drug.use[patient, ] | NA
		} else {
			drug.use[patient, "other"] = T
		}
	}
}

write.table(drug.use, "drug.use.txt", col.names=T, row.names=T, sep="\t", quote=F)
# drug.use = read.table("drug.use.txt", h=T)


###################################################
# Process patient information to generate surv
###################################################

	clin <- patient.tbl
	patients <- clin$bcr_patient_barcode
	death.days <- suppressWarnings(as.numeric(clin$death_days_to))
	contact.days <- suppressWarnings(as.numeric(clin$last_contact_days_to))
	surv <- data.frame(patients)
	surv$contact.days = contact.days
	surv$death.days = death.days
	# sanity check
	if(sum(contact.days > death.days, na.rm=T) != 0) 
		stop("Last contact after death?")
	death = rep(F, length(patients))
	death[which(!is.na(death.days))] = T
	days <- contact.days
	days[death] = death.days[death]
	surv$days = days
	surv$death = death
	surv$months = round(days/365.2425*12, 2)
	surv <- surv[!is.na(surv$days), ]
	write.table(surv, "surv.txt", col.names=T, row.names=F, sep="\t", quote=F)	
#	surv = read.table("surv.txt", h=T, sep="\t")


###################################################
# Combine drug use and survival data 
###################################################	
	m = match(surv$patients, rownames(drug.use))
	not.found = which(is.na(m))
	drug.use = drug.use[m, ]
	drug.use[not.found, ] = rep(F, 5)
	surv = cbind(surv, drug.use)
	rownames(surv) = surv$patients
	surv = surv[, -1]
	write.table(surv, "surv.4drugs.txt", col.names=T, row.names=T, sep="\t", quote=F)
#	surv = read.table("surv.4drugs.txt", h=T, sep="\t")


###################################################
# Load and process GDSC Taxol Drug Sensitivity Data
###################################################	
taxol.sensFile = "/Users/zhenyu/Downloads/working/Steph/ov/sensitivity_data_for_drug_11.csv"
taxotere.sensFile = "/Users/zhenyu/Downloads/working/Steph/ov/sensitivity_data_for_drug_1007.csv"
cisplatin.sensFile = "/Users/zhenyu/Downloads/working/Steph/ov/sensitivity_data_for_drug_1005.csv"

#### calculate taxol first
sensFile = taxol.sensFile


sens <- read.csv(sensFile, as.is=TRUE)
IC50s <- sens$IC.50
names(IC50s) <- sens$Cell.Line.Name
# tissues <- sens$Tissue
# names(tissues) <- sens$Cell.Line.Name


pData <- read.delim("/Users/zhenyu/Downloads/working/Steph/test/paper/Data/GdscPdata/E-MTAB-783.sdrf.txt", as.is=TRUE)
pDataUnique <- pData[pData$Source.Name %in% names(which(table(pData$Source.Name) == 1)), ]
rownames(pDataUnique) <- pDataUnique$Source.Name



commonCellLines <- rownames(pDataUnique)[rownames(pDataUnique) %in% names(IC50s)]
pDataUniqueOrd <- pDataUnique[commonCellLines, ]
IC50sOrd <- IC50s[commonCellLines]


###################################################
# Load and filter GDSC expression data
###################################################	

load(file="/Users/zhenyu/Downloads/working/Steph/test/paper/Data/GdscProcData/gdsc_brainarray_syms.RData")
trainDataOrd <- gdsc_brainarray_syms[, pDataUniqueOrd$Array.Data.File]


###################################################
# Load and process Bortezomib data
###################################################	
load("/Users/zhenyu/Downloads/working/Steph/test/paper/Data/bortezomibData/bortGeo.RData") # loads the geo data "bortezomib_mas5"
exprDataU133a <- cbind(exprs(bortezomib_mas5[[1]]), exprs(bortezomib_mas5[[2]]))
bortIndex <- c(which(pData(phenoData(bortezomib_mas5[[1]]))[,"characteristics_ch1.1"] == "treatment = PS341"), 255 + which(pData(phenoData(bortezomib_mas5[[2]]))[,"characteristics_ch1.1"] == "treatment = PS341"))
dexIndex <- c(which(pData(phenoData(bortezomib_mas5[[1]]))[,"characteristics_ch1.1"] == "treatment = Dex"), 255 + which(pData(phenoData(bortezomib_mas5[[2]]))[,"characteristics_ch1.1"] == "treatment = Dex"))
studyIndex <- c(as.character(pData(bortezomib_mas5[[1]])[, "characteristics_ch1"]), as.character(pData(bortezomib_mas5[[2]])[, "characteristics_ch1"]))


###################################################
# Change rownames of Bortezomib sample data exprDataU133a to gene names
###################################################	
library("hgu133a.db")
x <- hgu133aSYMBOL
mapped_probes <- mappedkeys(x)
names(mapped_probes) <- as.character(x[mapped_probes])
affy2sym <- as.character(x[mapped_probes])
sym2affy <- affy2sym
names(sym2affy) <- names(affy2sym)
dataU133aRowNames = sym2affy[rownames(exprDataU133a)]
save(dataU133aRowNames, file="/Users/zhenyu/Downloads/working/Steph/dataU133aRowNames.rda")
# load("/Users/zhenyu/Downloads/working/Steph/dataU133aRowNames.rda")
rownames(exprDataU133a) <- dataU133aRowNames


###################################################
# Filter aliquots by annotation database
###################################################	
file.tbl = data.frame(with(expr.U133.tbl, cbind(Comment..TCGA.Barcode., Comment..TCGA.Archive.Name., Array.Data.File, Comment..TCGA.Archive.Name..1, Derived.Array.Data.Matrix.File, Comment..TCGA.Archive.Name..2, Derived.Array.Data.Matrix.File.1)))
colnames(file.tbl) = c("aliquot", "lv1.url", "lv1.file", "lv2.url", "lv2.file", "lv3.url", "lv3.file")
tcga.url = "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/ov/cgcc/broad.mit.edu/ht_hg-u133a/transcriptome/"
file.tbl$lv1.url = with(file.tbl, paste0(tcga.url, lv1.url, "/", lv1.file))
file.tbl$lv2.url = with(file.tbl, paste0(tcga.url, lv2.url, "/", lv2.file))
file.tbl$lv3.url = with(file.tbl, paste0(tcga.url, lv3.url, "/", lv3.file))

# Query TCGA annotation database with disease OV
source("/Users/zhenyu/Downloads/working/Steph/ov/Module_Annotation.r")
disease = "OV"
annot.json <- GetAnnotationJSON(disease)
annot.tbl <- GetAnnotationTable(annot.json, outFile="/Users/zhenyu/Downloads/working/Steph/ov/OV.annotation.rda")

# Filter aliquots (and their corresponding analyte, protion, sample, patient) with any annotation items in the TCGA annotation database, except those with prior-malignancy and synchronized-malignancy
aliquots <- FilterByAnnotation(file.tbl$aliquot, annot.tbl, keep=c("6","204","6|204","204|6"))

# Only keep solid tumor aliquots
aliquots <- aliquots[substr(aliquots, 14, 15)=="01"]
file.tbl <- file.tbl[file.tbl$aliquot %in% aliquots, ]
aliquots <- file.tbl$aliquots
write.table(file.tbl, "file.tbl.txt", col.names=T, row.names=T, quote=F, sep="\t")



###################################################
# Process level 3 U133a expression data directly
###################################################

file = file.tbl[1, "lv3.url"]
test = read.table(textConnection(getURL(file)), sep="\t", skip=2, row.names=1)
expr = matrix(0, nrow(test), length(aliquots))
rownames(expr) = rownames(test)
colnames(expr) = aliquots

failedIndex = NULL
for(i in 1: length(aliquots)) {
	file = file.tbl$lv3.url[i]
	data = unlist(read.table(textConnection(getURL(file)), sep="\t", skip=2, row.names=1))
	if(length(data) == nrow(expr)) {
		expr[, i]  = data 
	} else {
		failedIndex = c(failedIndex, i)
	}
}
if(length(failedIndex) > 0) {
	cat("Expression data of some aliquots not retrieved.  Please try again! \n")
}
save(expr, file="raw.expr.rda")
# load("raw.expr.rda")

aliquots = colnames(expr)
patients = substr(aliquots, 1, 12)
# Combine expression value of the same patient, not quite tested 
if(sum(duplicated(patients)) >0) {
	expr = CollapseExpr(expr, patients)
} else {
	colnames(expr) = patients
}
save(expr, file="collapsed.expr.rda")
# load("collapsed.expr.rda")


###################################################
# Alternative to the step above, download the process OV expression data
###################################################	
library("affy")
celFilesLoc <- "/glusterfs/netapp/homes1/ZHANGZ18/steph/ov/raw/"
celFileName <- paste(celFilesLoc, dir(celFilesLoc, pattern=".CEL$"), sep="")
expr_raw <- ReadAffy(filenames=celFileName)
expr_brainarray_entrez <- rma(expr_raw)

library("hgu133plus2hsentrezg.db")
x <- hgu133plus2hsentrezgSYMBOL
mapped_probes_brain <- mappedkeys(x) 
names(mapped_probes_brain) <- as.character(x[mapped_probes_brain])
entrez2sym_brain <- as.character(x[mapped_probes_brain])
symBrain2entrez <- entrez2sym_brain
names(symBrain2entrez) <- names(entrez2sym_brain)
expr_brainarray_syms <- exprs(expr_brainarray_entrez)
sum(!is.na(match(rownames(exprs(expr_brainarray_entrez)), mapped_probes_brain)))
rownames(expr_brainarray_syms) <- symBrain2entrez[rownames(exprs(expr_brainarray_entrez))]
# file.tbl = read.table("/glusterfs/netapp/homes1/ZHANGZ18/steph/ov/file.tbl.txt", h=T, row.names=1, sep="\t")
m = match(colnames(expr_brainarray_syms), file.tbl$lv1.file)
colnames(expr_brainarray_syms) = file.tbl$aliquot[m]



###################################################
# Prediction
###################################################	

predictedSensitivity133a <- calcPhenotype(trainDataOrd, IC50sOrd, expr, selection=1)

write.table(predictedSensitivity133a, "predictedSensitivity133a.txt", col.names=F, row.names=T, quote=F, sep="\t")

common.patients = intersect(names(predictedSensitivity133a), rownames(surv))
# common.patients = intersect(rownames(predictedSensitivity133a), rownames(surv))

surv = surv[common.patients, ] 
surv$taxol.prediction = predictedSensitivity133a[common.patients]
# predictedSensitivity133a = predictedSensitivity133a[common.patients, ]
# surv = cbind(surv, predictedSensitivity133a)
save(surv, file="surv.lv3.pred.rda")


###################################################
# Compile final dataset for analysis
###################################################	

# 450 patients
data = surv
# 3 patients with ambiguous drug use information
patients.with.na.drug = unique(names(which(which(is.na(data), T)[, 2] >= 6)))
# 30 patients with no drug use information (no use or no info?)
patients.with.no.drug.info = rownames(data)[with(data, which(!(taxol | taxotere | cisplatin | carboplatin | other)))]
# 417 patients left
data = data[-which(rownames(data) %in% c(patients.with.na.drug, patients.with.no.drug.info)), ]

surv2 = surv[-which(rownames(surv) %in% patients.with.na.drug), ]
surv2 = rbind( 	cbind(surv[which(surv$taxol), ], cat = "Taxol"), 
				cbind(surv[which(surv$taxotere), ], cat = "Taxotere"),
				cbind(surv[which(surv$cisplatin), ], cat = "Cisplatin"),
				cbind(surv[which(surv$carboplatin), ], cat = "Carboplatin"),
				cbind(surv[which(surv$other), ], cat = "Other"),
				cbind(surv[patients.with.no.drug.info, ], cat = "None"))
ggsurv(survfit(with(surv2, Surv(months,death) ~ cat)))

surv3 = surv[-which(rownames(surv) %in% patients.with.na.drug), ]
surv3$cat = "Other"
surv3$cat[with(surv3, which(taxol|taxotere))] = "taxo_only"
surv3$cat[with(surv3, which(cisplatin|carboplatin))] = "platin_only"
surv3$cat[with(surv3, which((cisplatin|carboplatin) & (taxol|taxotere)))] = "taxo_and_platin"
surv3[patients.with.no.drug.info, ]$cat = "None"
ggsurv(survfit(with(surv3, Surv(months,death) ~ cat)))



# Plot Venn Diagram of drug distributions
library(VennDiagram)
x = with(data, list(which(taxol), which(taxotere), which(carboplatin), which(cisplatin)))
names(x) = c("Taxol", "Taxotere", "Carboplatin", "Cisplatin")
venn.diagram(x, filename = "Venn.tiff", col = "transparent", fill = c("cornflowerblue","green","yellow","darkorchid1"), alpha = 0.50, label.col = c("orange", "white", "darkorchid4", "white", "white",  "white", "white", "white", "darkblue", "white", "white", "white", "white",  "darkgreen", "white"), cex = 1.5, fontfamily = "serif", fontface = "bold", cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"), cat.cex = 1.5, cat.pos = 0, cat.dist = 0.07, cat.fontfamily = "serif", rotation.degree = 270, margin = 0.2)

####### taxol
drug = "taxol"
taxol.data = data[data$taxol, ]
taxol.data$pred = predictedSensitivity133a[rownames(taxol.data)]
pred = taxol.data$pred
quantile = rep("Q4", length(pred))
quantile[pred > quantile(pred, 0.25)] = "Q3"
quantile[pred > quantile(pred, 0.50)] = "Q2"
quantile[pred > quantile(pred, 0.75)] = "Q1"
taxol.data$quantile = quantile

qplot(death, pred, data=taxol.data, geom=c("boxplot", "jitter"), 
   fill=death, main="Predicted Taxol Sensitivity by Death",
   xlab="", ylab="Predicted Taxol Sensitivity")

ggsurv(survfit(with(taxol.data, Surv(months,death) ~ quantile)))

ggplot(taxol.data, aes(x=months, y=pred, color=death)) + geom_point(shape=1) + geom_smooth(method=lm)

data2 = taxol.data[taxol.data$death, ]
0.682 

ggsurv(survfit(with(surv, Surv(months,death) ~ taxol)))

# no NA in taxol
surv2 = surv1[which(!is.na(surv1$taxol)), ]
ggsurv(survfit(with(surv2, Surv(months,death) ~ taxol)))


# no other taxo.platin drugs
surv3 = surv2[with(surv2, which(!(taxotere|carboplatin|cisplatin))), ]
ggsurv(survfit(with(surv3, Surv(months,death) ~ taxol)))


# no no-info patient
surv4=surv2[-with(surv2, which(!(taxol|taxotere|carboplatin|cisplatin|other))),]
ggsurv(survfit(with(surv4, Surv(months,death) ~ taxol)))

surv.data = with(surv, Surv(months,death) ~ taxol)
ggsurv(survfit(surv.data))







surv.data = with(surv, Surv(months,death) ~ Taxotere)
surv.data = with(surv, Surv(months,death) ~ Carboplatin)
surv.data = with(surv, Surv(months,death) ~ Cisplatin)
surv.data = with(surv, Surv(months,death) ~ (carboplatin|cisplatin))
surv.data = with(surv, Surv(months,death) ~ (taxol|taxotere))
surv.data = with(surv, Surv(months,death) ~ ((carboplatin|cisplatin) & (taxol|taxotere)))

file = "bcgsc.ca_OV.IlluminaHiSeq_RNASeq.sdrf.txt"
tbl = GetTCGATable(file)
sample = unique(substr(tbl$Comment..TCGA.Barcode., 1, 15))
patient = unique(substr(sample[substr(sample, 14, 15)=="01"], 1, 12))
length(patient)





FilterDrugTable <- function(tbl){
  date.valid <- with(tbl, which(
                pharmaceutical_tx_started_days_to != "[Not Available]" 
                & pharmaceutical_tx_ended_days_to != "[Not Available]" ))  
  tbl <- tbl[date.valid,]
  date.valid <- with(tbl, which(as.numeric(pharmaceutical_tx_ended_days_to) 
                - as.numeric(pharmaceutical_tx_started_days_to) > 5))        
  tbl <- tbl[date.valid,]
  return(tbl)
}

  
  	
GetTCGATable <- function(url) {
  if(grepl("^http", url)) {
    s <- try(content(GET(url), as = "text"), silent = TRUE);
    if (class(s) == "try-error") 
      return(NULL);    
    tbl <- read.table(textConnection(s), h=T, sep="\t", quote="\"")	
  } else {
    tbl = data.frame(fread(url))
  }
  if(grepl("bcr_", tbl[1,1]) | grepl("CDE_ID", tbl[1,1]))
    tbl	<- tbl[-1,]
  if(grepl("bcr_", tbl[1,1]) | grepl("CDE_ID", tbl[1,1]))
    tbl	<- tbl[-1,]
  return(tbl);
}

# GetDrugAliasTable take 1 or 2 tables in CSV format, with the first column of TCGA drug names, and the third column of converted names, and return a conversion table with colnames "TCGA_drug_name", "real_drug_name" and "simplified_drug_name", with all drug names in lower cases
GetDrugAliasTable <- function(alias.csv, typo.csv = "") {
	options(stringsAsFactor = F)
	library(stringr)
	names = c("TCGA_drug_name", "real_drug_name", "simplified_drug_name")
	
	# Read first CSV
	alias = read.csv(alias.csv, h=T)
	names(alias) = names

	# Read second CSV if exists
	if (typo.csv != "") {
		typo = read.csv(typo.csv, h=T)
		names(typo) = names
		alias = rbind(alias, typo)
	}
	
	# Trim, tolower, and sanity check
	alias = apply(alias, 2, str_trim)
	alias = data.frame(apply(alias, 2, tolower))
	if(sum(alias$TCGA_drug_name %in% alias$simplified_drug_name) > 0 ) {
		cat("Warning! Some TCGA original names appear in the simplified names")
		print(alias$TCGA_drug_name[which(alias$TCGA_drug_name %in% alias$simplified_drug_name)])
	}
	return(alias)
}






