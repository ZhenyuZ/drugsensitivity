options(stringsAsFactors = F)
source("/Users/zhenyu/Downloads/working/Steph/TCGA_Module.r")

drug.alias.csv = "/Users/zhenyu/Downloads/working/Steph/drug_alias.csv"
drug.typo.csv = "/Users/zhenyu/Downloads/working/Steph/drug_typo.csv"

tbl = GetDrugAliasTable(drug.alias.csv, drug.typo.csv)
drugs = tolower(drugs)

output = character()
for (i in 1: length(drugs)) {
	drug = drugs[i]	
	w = which(tbl$TCGA_drug_name == drug)
	if(length(w) == 0) {
		output = c(output, drug)
	} else {
		output = c(output, tbl[w, "simplified_drug_name"])
	}
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


	lower.drugs = tolower(drugs)
	
	