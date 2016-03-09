#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. Marz 2016                        ###
#############################################################################
### SR top hits drugs names                                               ###
### From SR exection get top hits and annotate names using ontology info  ###
#############################################################################

##Getting HOME directory
home <- Sys.getenv("HOME")

# To use this script in ant first export this:
# export R_LIBS="/software/R/packages"

# library (ggplot2, lib.loc="/users/cn/jespinosa/R/library")
# library (plyr, lib.loc="/users/cn/jespinosa/R/library")
# library(reshape)
# library(xlsx)

args <- commandArgs ()
path_sr_tbl <- args[5]
path_ATC_names <- args[6]
name_out_tag <- args[7]

# sex <- args[7]
# age <- args[8]
# mode <- mode[9]

sr_df <- read.csv(file=path_sr_tbl, header=FALSE, sep="\t")
# sr_df <- read.csv(file="~/2015_viscMes/data/SR.csv", header=FALSE, sep="\t")

sr_df_ordered <- sr_df[with(sr_df, order(-V3)), ]

#######################################
# https://bioportal.bioontology.org/ontologies/ATC
atc_df <- read.csv(file=path_ATC_names, header=TRUE, sep=",")
# atc_df <- read.csv(file="/Users/jespinosa/2015_viscMes/data/SR_analysis/ATC.csv", header=TRUE, sep=",")
atc_vect <- gsub ("http://purl.bioontology.org/ontology/UATC/", "", atc_df$Class.ID)
atc_vector <- gsub ("http://purl.bioontology.org/ontology/STY/", "", atc_vect)
atc_df$code <- atc_vector
# head(atc_df)

top_hits <- head (sr_df_ordered, 30)

# head (sr_df)
# substring to 4 letter code 
sr_df$index_code <- substr(sr_df$V1, 1, 4)
sr_df$marker_code <- substr(sr_df$V2, 1, 4)
# head(sr_df)
# change colnames
colnames(sr_df) <- c("index_drug", "marker_drug", "SR", "index_code", "marker_code")

sr_df_index_lab <- merge (sr_df, atc_df[, c("code", "Preferred.Label")], by.x= "index_code", by.y = "code", sort=F)
sr_df_lab <- merge (sr_df_index_lab, atc_df[, c("code", "Preferred.Label")], by.x= "marker_code", by.y = "code", sort=F)
# head(sr_df_lab)

colnames(sr_df_lab) <- c(colnames(sr_df_lab)[1:5], "ATC_family_index", "ATC_family_marker")
# head(sr_df_lab)
sr_df_lab_index_name <- merge (sr_df_lab, atc_df[, c("code", "Preferred.Label")], by.x= "index_drug", by.y = "code", sort=F)
sr_df_lab_name <- merge (sr_df_lab_index_name, atc_df[, c("code", "Preferred.Label")], by.x= "marker_drug", by.y = "code", sort=F)

colnames(sr_df_lab_name) <- c(colnames(sr_df_lab_name)[1:7], "name_index", "name_marker")
head (sr_df_lab_name)

sr_df_lab_name_ordered <- sr_df_lab_name[with(sr_df_lab_name, order(-SR)), ]
# head(sr_df_lab_name_ordered,15)
# install.packages("xlsx")

tbl_file_ext <- ".csv"
name_tbl <-paste ("top_50_sr", "_" , name_out_tag, tbl_file_ext, sep="")

# write.xlsx(head(sr_df_lab_name_ordered, 50), name_tbl) 
write.table(head(sr_df_lab_name_ordered, 50), file = name_tbl, row.names=TRUE, sep="\t")


