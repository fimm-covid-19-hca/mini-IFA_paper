options(java.parameters = "-Xmx6g")
require(gridExtra)
library(reshape2)
library(xlsx)
library(plyr)
library(plotly)
library(shiny)
library(stringi)
library(htmlwidgets)
options(scipen = 999)


###
##Create directory structure
##Data\Images     For Storing the images
##Data\Normalisation_NegOnly           
##Results
##Code
##RDA
##Set working directories first
norm_name="Normalisation_NegOnly"

results_dir=paste0("..\\Results\\",norm_name)
#setwd(results_dir)
data_dir=paste0("..\\..\\Data\\",norm_name)
source("..\\..\\Code\\Misc_Functions.R")

mycol=c("N"="tomato","S"="dodgerblue","R"="seagreen1","M"="snow3",
        "CTRL1"="green1","CTRL2"="violet","EMPTY"="cyan","F5_100"="gray50","F5_25"="gray50","F6_100"="gray50","F6_25"="gray50","F7_100"="gray50","F7_25"="gray50","F8_100"="gray50","F8_25"="gray50",
        "NEG1"="firebrick","NEG2"="mediumorchid1","NEG3"="lightpink","POS1"="blue","POS2"="slateblue1")

myshape=c("CTRL1"=45,"CTRL2"=45,"EMPTY"=13,"F5_100"=1,"F5_25"=1,"F6_100"=1,"F6_25"=1,"F7_100"=1,"F7_25"=1,"F8_100"=1,"F8_25"=1,
          "NEG1"=25,"NEG2"=25,"NEG3"=25,"POS1"=43,"POS2"=25)
controls=c("CTRL1","CTRL2","EMPTY","F5_100","F5_25","F6_100","F6_25","F7_100","F7_25","F8_100","F8_25","NEG1","NEG2","NEG3","POS1","POS2")

controls_name=c("CTRL1","CTRL2","EMPTY","F5_100","F5_25","F6_100","F6_25","F7_100","F7_25","F8_100","F8_25","NEG1","NEG2","NEG3","POS1","POS2")
proteins=c("S","R","N","M")

file_names <- list.files(data_dir,pattern  = "ratios_per_well_plate*", all.files = TRUE,full.names = TRUE, recursive = TRUE,ignore.case = TRUE, include.dirs = TRUE)
screen_table=do.call("rbind.fill",lapply(as.character(file_names),function(file_name)
{
  tbl_dataframe= read.csv(file_name, header=TRUE, blank.lines.skip=TRUE,skip=0,sep=",",colClasses=NA,stringsAsFactors =FALSE)
  tbl_dataframe$Sample_Type=ifelse(tbl_dataframe$Sample.Name %in% c(controls,controls_name),tbl_dataframe$Sample.Name,"Sample")
  protein=split_give_last_n_vect(tbl_dataframe$Destination.Plate.Barcode[1],"-",1)
  protein=unname(substr(protein,1,1))
  tbl_dataframe$Protein=protein
  tbl_dataframe$Plate_set=unname(split_remove_last_n_vect(tbl_dataframe$Destination.Plate.Barcode[1],"-",1))
  names(tbl_dataframe)[names(tbl_dataframe) == "well"] <- "Well"
  tbl_dataframe
}))

screen_table$DRow=substr(screen_table$Well,1,1)
screen_table$DCol=as.numeric(substr(screen_table$Well,2,3))

screen_table=screen_table[,c("Plate_set","Protein","Destination.Plate.Barcode","Well","DRow","DCol","Transferred","Sample_Type","Content","Sample","Sample.Name","Dilution","Cell.Count","Raw.Mean.NUCLEUS.INTENSITY.MEAN.DAPI",
                             "Ratio.positive.IgM","Ratio.negative.IgM","Ratio.atypical.IgM","Ratio.small.bright.IgM","Ratio.trash.IgM", "Raw.Mean.CELL.INTENSITY.MEAN.Alexa.IgM","Raw.Mean.DONUT.INTENSITY.MEAN.Alexa.IgM","Raw.Mean.Positives.CELL.INTENSITY.MEAN.Alexa.IgM","Raw.Mean.Positives.DONUT.INTENSITY.MEAN.Alexa.IgM",
                             "Ratio.positive.IgA","Ratio.negative.IgA","Ratio.atypical.IgA","Ratio.small.bright.IgA","Ratio.trash.IgA", "Raw.Mean.CELL.INTENSITY.MEAN.Alexa.IgA","Raw.Mean.DONUT.INTENSITY.MEAN.Alexa.IgA","Raw.Mean.Positives.CELL.INTENSITY.MEAN.Alexa.IgA","Raw.Mean.Positives.DONUT.INTENSITY.MEAN.Alexa.IgA",
                             "Ratio.positive.IgG","Ratio.negative.IgG","Ratio.atypical.IgG","Ratio.small.bright.IgG","Ratio.trash.IgG", "Raw.Mean.CELL.INTENSITY.MEAN.Alexa.IgG","Raw.Mean.DONUT.INTENSITY.MEAN.Alexa.IgG","Raw.Mean.Positives.CELL.INTENSITY.MEAN.Alexa.IgG","Raw.Mean.Positives.DONUT.INTENSITY.MEAN.Alexa.IgG")]                          


screen_table=screen_table[with(screen_table, order(Plate_set,Protein,Destination.Plate.Barcode,DRow,DCol)),]
save(screen_table,file="..\\..\\RDA\\screen_table.rda")
write.xlsx2(screen_table,paste0("Combined_Data_",norm_name,".xlsx"),sheetName = "Data",append=FALSE)