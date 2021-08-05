source("..\\..\\Code\\Misc_Functions.R")
options(java.parameters = "-Xmx15g")
require(gridExtra)
library(reshape2)
library(xlsx)
library(plyr)
library(plotly)
library(shiny)
library(stringi)
library(htmlwidgets)
library(Hmisc)
library(plotrix)
options(scipen = 999)

#norm_name="Normalisation_NegOnly"
#results_dir=paste0(".\\Results\\",norm_name,"\\Plots")
#setwd(results_dir)

data_dir=paste0("..\\..\\Data\\",norm_name)
images_dir=paste0("..\\..\\Data\\Images\\resample")
load(file="..\\..\\RDA\\screen_table.rda")

############################################################################################
############################################################################################
stderr <- compiler::cmpfun(function(x){sqrt(var(x,na.rm=T)/length(na.omit(x)))})
lowsd <- compiler::cmpfun(function(x){return(mean(x)-stderr(x))})
highsd <- compiler::cmpfun(function(x){return(mean(x)+stderr(x))})
pop.sd <- compiler::cmpfun(function(x)(sqrt(var(x)*(length(x)-1)/length(x)))); 
ssmd <- compiler::cmpfun(function(x,y)round((mean(x)-mean(y))/sqrt(var(x)+var(y))));
zfactor <- compiler::cmpfun(function(x,y)round((1-(3 * (pop.sd(x)+pop.sd(y)))/(abs(mean(x)-mean(y)))),2));
robustzfactor <- compiler::cmpfun(function(x,y)round((1-(3 * (mad(x)+mad(y)))/(abs(median(x)-median(y)))),2));


controls=c(controls,"CTRL1","CTRL2","EMPTY","F5_100","F5_25","F6_100","F6_25","F7_100","F7_25","F8_100","F8_25","NEG1","NEG2","NEG3","POS1","POS2")
controls_name=c("CTRL1","CTRL2","EMPTY","F5_100","F5_25","F6_100","F6_25","F7_100","F7_25","F8_100","F8_25","NEG1","NEG2","NEG3","POS1","POS2")
proteins=c("S","R","N","M")
mycol=c("Ratio.positive"="skyblue1","Ratio.negative"="tomato","Ratio.atypical"="slategray","Ratio.small.bright"="yellow","Ratio.trash"="springgreen",
        "Negative"="salmon","COVID-19"="deepskyblue","NA"="salmon",">2Weeks"="navy","<2Weeks"="deepskyblue","Transfection_%"="#999999",
        "CTRL1"="#CC6666","CTRL2"="#66CC99","EMPTY"="#9999CC",
        "NEG1"="firebrick","NEG2"="mediumorchid1","NEG3"="mediumvioletred","POS1"="blue","POS2"="slateblue1",
        "N"="#E69F00","M"="#999999","R"="#009E73","S"= "#56B4E9","IgA"="#0072B2","IgG"="#D55E00","IgM"="#CC79A7",
        "Sample"="gray30")

myshape=c("CTRL1"=8,"CTRL2"=13,"EMPTY"=10,"Sample"=21,
          "NEG1"=95,"NEG2"=95,"NEG3"=18,"POS1"=43,"POS2"=43,
          "M"=21,"N"=22,"R"=23,"S"=24)
          
############################################################################################
#####Convert images to base 64

file_names <- list.files(images_dir,pattern  = "^[^~]*.jpg", all.files = TRUE,full.names = TRUE, recursive = TRUE,ignore.case = TRUE, include.dirs = TRUE)
img_table=do.call("rbind.fill",lapply(as.character(file_names),function(file_name)
{
  print(file_name)
  fn=remove_directory_path(file_name)
  file_info=unlist(strsplit(file_name,"\\/"))[2]
  Channel=unname(split_give_last_n_vect(file_info,"_",1))
  Plate=unlist(strsplit(file_info,"_"))[1]
  Protein=unname(split_give_last_n_vect(Plate,"-",1))
  Protein=unname(substr(Protein,1,1))
  Sample_Name_Full=remove_file_extension(fn)
  
  Image_base64=base64enc::dataURI(file = file_name)
  
  img_64 <- data.frame(Destination.Plate.Barcode=Plate,Channel=Channel,Protein=Protein,Sample_Name_Full=Sample_Name_Full,Image_base64=Image_base64,stringsAsFactors=FALSE)  
  img_64
}))

save(img_table,file="..\\..\\RDA\\img_table.rda")

############################################################################################

f_samples=c("F5_100","F5_25","F6_100","F6_25","F7_100","F7_25","F8_100","F8_25")

screen_table=screen_table[!screen_table$Sample.Name %in% f_samples,]

dt_long=gather(data=screen_table, key=Feature, value=Score,c(15:41), na.rm = FALSE, convert = FALSE)
dt_long=dt_long[!is.na(dt_long$Score),]
dt_long$Channel=stri_sub(dt_long$Feature,-3,-1)
#############################################
clinical_list= read.xlsx2("..\\..\\Annotations\\20210420_mini_IFA_CoV_Samples_Clinical_data_MP_tarkistettu.xlsx",1,colIndex=1:36, as.data.frame=TRUE, header=TRUE, colClasses="character",stringsAsFactors =FALSE)
clinical_list=clinical_list[,c("Sample.code","FIG_2A","FIG_2C","sample","date","type","PCR.Sample.group","Timw.bw.sampling.and.PCR..days.","Severity","Sample.type","serial")]
clinical_list$Sample.Date=format(as.POSIXct(as.numeric(clinical_list$date)*86400, origin="1899-12-30",tz="GMT"),'%d.%m.%Y')
names(clinical_list)[names(clinical_list) == "sample"] <- "Sample.Name"
names(clinical_list)[names(clinical_list) == "Sample.type"] <- "Sample.Origin"
names(clinical_list)[names(clinical_list) == "Timw.bw.sampling.and.PCR..days."] <- "Time_Bw_Sampling_PCR"
clinical_list$SN=clinical_list$Sample.Name
clinical_list$Sample.Name=gsub("CoV-","",clinical_list$Sample.Name)
clinical_list$Sample.Name=str_pad(clinical_list$Sample.Name, width=3, pad="0")
clinical_list$Sample.Name=paste0("S2020_",clinical_list$Sample.Name)

dt_long_Clinical=merge(dt_long,clinical_list,by.x=c("Sample.Name","Sample"),by.y=c("Sample.Name","serial"),all.x=TRUE)

################################3
dt_fig2a=dt_long_Clinical[dt_long_Clinical$FIG_2A %in% c("x") | grepl("^S2017",dt_long_Clinical$Sample.Name),]
dt_fig2a$Sample_Group=ifelse(grepl("^S2017",dt_fig2a$Sample.Name),"Negative","COVID-19")
dt_fig2a=dt_fig2a[grepl("Ratio.positive",dt_fig2a$Feature),]
dt_fig2a$PCR_Sampling_Time=ifelse(is.na(dt_fig2a$Time_Bw_Sampling_PCR),"NA",ifelse(dt_fig2a$Time_Bw_Sampling_PCR >= 14,">2Weeks","<2Weeks"))


########Just generate unique set of values and store in rda file. Same values will be used for all the plotting
sample_negs=unique(dt_fig2a$Sample.Name[dt_fig2a$Sample_Group %in% c("Negative")])
sample_positives=unique(dt_fig2a$Sample.Name[dt_fig2a$Sample_Group %in% c("COVID-19")])
sample_positives_unique=unique(dt_fig2a$Sample.code[dt_fig2a$Sample_Group %in% c("COVID-19")])

 sample_negatives=unique(sample(sample_negs,size=length(sample_positives_unique),replace=FALSE))
 save(sample_negatives,file="..\\..\\RDA\\sample_negatives.rda")
 sample_negatives_fig2c=unique(sample(sample_negatives,size=5,replace=FALSE))
 save(sample_negatives_fig2c,file="..\\..\\RDA\\sample_negatives_fig2c.rda")

load("..\\..\\RDA\\sample_negatives.rda")
load("..\\..\\RDA\\sample_negatives_fig2c.rda")

dt_fig2a=dt_fig2a[dt_fig2a$Sample.Name %in% c(sample_negatives,sample_positives),]
###########
###Table of ratios, positive, negatives,

corr_table_dotplots=data.frame(
                               Antibody=character(),Protein=character(),
                               P.Value=numeric(),
                               Mean_Negative=numeric(),SD_Negative=numeric(),StdErr_Negatives=numeric(),
                               Mean_COVID_19=numeric(),SD_COVID_19=numeric(),StdErr_COVID_19=numeric(),
                               stringsAsFactors=FALSE)


for(ab in unique(sort(dt_fig2a$Channel)))
{
  ab_table=dt_fig2a[dt_fig2a$Channel==ab,]
  
  for(prot in unique(sort(ab_table$Protein)))
  {
    prot_ab_table=ab_table[ab_table$Protein==prot,]
    
    sample_positives_scores=prot_ab_table$Score[prot_ab_table$Sample.Name %in% c(sample_positives)]
    sample_negs_scores=prot_ab_table$Score[prot_ab_table$Sample.Name %in% c(sample_negs)]
    P.Value=wilcox.test(sample_negs_scores,sample_positives_scores,paired=FALSE)
    
    corr_record=data.frame(Antibody=ab,Protein=prot,
                           P.Value=P.Value$p.value,
                           Mean_Negative=round(mean(sample_negs_scores),5),
                           Mean_COVID_19=round(mean(sample_positives_scores),5),
                           SD_Negative=round(sd(sample_negs_scores),5),
                           SD_COVID_19=round(sd(sample_positives_scores),5),
                           StdErr_Negatives=round(std.error(sample_negs_scores,na.rm = TRUE),5),
                           StdErr_COVID_19=round(std.error(sample_positives_scores,na.rm = TRUE),5),
                           stringsAsFactors=FALSE)
                           
    corr_table_dotplots=rbind(corr_table_dotplots,corr_record)
  }
}  

write.xlsx2(corr_table_dotplots,"P_Value_Stats_etc.xlsx",sheetName = "Pvalues",append=FALSE)

dt_fig2a_IgA_IgG_IgM=dt_fig2a

dt_fig2a=dt_fig2a[dt_fig2a$Channel %in% c("IgA","IgG"),]

p1 <-ggplot(data=dt_fig2a, aes(x=Sample_Group,y=Score)) +
  geom_boxplot(lwd=0.2,fatten=0.5,width = 0.5, position = position_dodge(0),outlier.colour = NA) +
  #geom_violin(trim = FALSE) +
  geom_dotplot(aes(color = PCR_Sampling_Time,fill = PCR_Sampling_Time),
               #binwidth = 0.1,
               #count=5,
               width=1,
               binpositions="all",
               method = "histodot",#dotdensity,stackratio=0.5,
               binaxis='y', stackdir='center',
               dotsize = 0.3,position = position_dodge(0.2))+
  scale_color_manual(values = mycol)+
  scale_fill_manual(values = mycol,labels=c("<2Weeks",">2Weeks","Negative Control"))+
  scale_y_continuous(limits = c(-0.01, 0.35), breaks = seq(0, 0.35, by = 0.05))+
  facet_grid(Channel ~fct_relevel(Protein,'S','R','N','M'), scales = 'free')+
  ggtitle(paste0("Predicted Positive Ratios")) +
  labs(y = "Positive Ratio",x="Samples",fill="Time Between Sampling and PCR")+
  theme_classic()+
  guides(color = FALSE)+
  theme(
    aspect.ratio=1,
    legend.position = "bottom",
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white", colour = "grey50"),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.001, "lines"),
    #legend.title = element_blank(),
    strip.text.x = element_text(color="black", size=7, face="bold"),
    strip.text.y = element_text(color="black", size=7, face="bold"),
    plot.title = element_text(color="grey20", size=10, face="bold"),
    axis.title.y = element_text(color="grey50", size=9, face="bold"),
    axis.title.x = element_text(color="grey50", size=9, face="bold"),
    axis.text.x = element_text(size=5,angle = 0, vjust = 1, hjust=0.5),
    axis.text.y = element_text(size=5,angle = 0, vjust = 1, hjust=0.5),
    legend.title=element_text(size=5), 
    legend.text=element_text(size=4))
print(p1)
ggsave(paste0("Fig2A_Dotplots","_",norm_name,".tiff"),plot=p1, width = 4, height = 4, dpi=600, compression = "lzw")
####################################3
#####
###Fig2A With IgM
#dt_fig2a_IgA_IgG_IgM

p2 <-ggplot(data=dt_fig2a_IgA_IgG_IgM, aes(x=Sample_Group,y=Score)) +
  geom_boxplot(lwd=0.2,fatten=0.5,width = 0.5, position = position_dodge(0),outlier.colour = NA) +
  #geom_violin(trim = FALSE) +
  geom_dotplot(aes(color = PCR_Sampling_Time,fill = PCR_Sampling_Time),
               #binwidth = 0.1,
               #count=5,
               width=1,
               binpositions="all",
               method = "histodot",#dotdensity,stackratio=0.5,
               binaxis='y', stackdir='center',
               dotsize = 0.3,position = position_dodge(0.2))+
  #stat_summary(fun.data=mean_sdl, fun.args = list(mult=3), geom="pointrange", size=0.1,color="black")+
  scale_color_manual(values = mycol)+
  scale_fill_manual(values = mycol,labels=c("<2Weeks",">2Weeks","Negative Control"))+
  scale_y_continuous(limits = c(-0.01, 0.8), breaks = seq(0, 0.8, by = 0.1))+
  facet_grid(Channel ~fct_relevel(Protein,'S','R','N','M'), scales = 'free')+
  ggtitle(paste0("Predicted Positive Ratios")) +
  labs(y = "Positive Ratio",x="Samples",fill="Time Between Sampling and PCR")+
  theme_classic()+
  guides(color = FALSE)+
  theme(
    aspect.ratio=1,
    legend.position = "bottom",
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white", colour = "grey50"),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.001, "lines"),
    #legend.title = element_blank(),
    strip.text.x = element_text(color="black", size=7, face="bold"),
    strip.text.y = element_text(color="black", size=7, face="bold"),
    plot.title = element_text(color="grey20", size=10, face="bold"),
    axis.title.y = element_text(color="grey50", size=9, face="bold"),
    axis.title.x = element_text(color="grey50", size=9, face="bold"),
    axis.text.x = element_text(size=5,angle = 0, vjust = 1, hjust=0.5),
    axis.text.y = element_text(size=5,angle = 0, vjust = 1, hjust=0.5),
    legend.title=element_text(size=5), 
    legend.text=element_text(size=4))
print(p2)
ggsave(paste0("Fig2A_Dotplots_With_IgM","_",norm_name,".tiff"),plot=p2, width = 4, height = 4, dpi=600, compression = "lzw")



##################################################3
######Figure 2C dt_fig2C
#######################################3

#dt_fig2C=dt_long_Clinical[dt_long_Clinical$FIG_2C %in% c("x (hospital)","x (home)","x (ICU)") | dt_long_Clinical$Sample.Name %in% sample_negatives_fig2c,]
dt_fig2C=dt_long_Clinical[dt_long_Clinical$FIG_2C %in% c("x (hospital)","x (home)","x (ICU)"),]

dt_long_median_neg=dt_long_Clinical[dt_long_Clinical$Sample.Name %in% sample_negatives,]
dt_long_median_neg=dt_long_median_neg[,c("Channel","Protein","Cell.Count","Raw.Mean.NUCLEUS.INTENSITY.MEAN.DAPI","Feature","Score")]

dt_long_median_neg=dt_long_median_neg[grepl("Ratio.positive",dt_long_median_neg$Feature),]
#write.xlsx(dt_long_median_neg,"dt_long_median_neg.xlsx",sheetName = "Fig2C",append = FALSE)
dt_long_median_neg <- dt_long_median_neg %>% group_by(Feature,Channel,Protein) %>% summarise_all(median) %>% as.data.frame()
dt_long_median_neg$Sample.Name="Median_Neg"

dt_fig2C=rbind.fill(dt_fig2C,dt_long_median_neg)

dt_fig2C=dt_fig2C[grepl("Ratio.positive",dt_fig2C$Feature),]
dt_fig2C$Severity=gsub("Hosp., non-ICU","Hospital(Non-ICU)",dt_fig2C$Severity)
dt_fig2C$Severity=gsub("home","Home",dt_fig2C$Severity)
dt_fig2C$Severity=ifelse(is.na(dt_fig2C$Severity),"Negative Control",dt_fig2C$Severity)


dt_fig2C_IgA_IgG=dt_fig2C[dt_fig2C$Channel %in% c("IgG","IgA"),]

p2 <-ggplot(data=dt_fig2C_IgA_IgG, aes(x=Sample.Name,y=Score,fill=Channel,shape=Protein)) +
  geom_point(position = position_jitter(w=0.1),size=1,alpha=0.8)+
  scale_fill_manual(values = mycol)+
  scale_shape_manual(values=myshape)+
  facet_grid(Channel ~ Severity, scales = 'free_x')+
  ggtitle(paste0("Predicted Positive Ratios")) +
  labs(x = "COVID-19 Patients",y="Positivity Ratio")+
  scale_y_continuous(limits = c(-0.01, 0.35), breaks = seq(0, 0.35, by = 0.05))+
  guides(fill=guide_legend(override.aes=list(colour=c("IgA"="#0072B2","IgG"="#D55E00"))))+
  theme_classic()+
  theme(
    aspect.ratio=1,
    legend.position = "bottom",
    panel.border = element_rect(fill=NA, size=0,colour="black"),
    panel.background = element_rect(fill = "white", colour = "grey50"),
    panel.grid.minor.x = element_blank(),
    panel.spacing.x = unit(0,"line"),
    #strip.placement = 'outside',
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(color="black", size=7, face="bold"),
    strip.text.y = element_text(color="black", size=7, face="bold"),
    plot.title = element_text(color="grey20", size=10, face="bold"),
    axis.title.y = element_text(color="grey50", size=9, face="bold"),
    axis.title.x = element_text(color="grey50", size=9, face="bold"),
    axis.text.x = element_text(size=5,angle = 60, vjust = 1, hjust=1),
    axis.text.y = element_text(size=5,angle = 0, vjust = 1, hjust=0.5),
    legend.title=element_text(size=5), 
    legend.text=element_text(size=4))
print(p2)
ggsave(paste0("Fig2D_Scatterplots","_",norm_name,".tiff"),plot=p2,width = 4, height = 4, dpi=600, compression = "lzw")

dt_fig2C_IgM=dt_fig2C[dt_fig2C$Channel %in% c("IgG","IgA","IgM"),]

#write.xlsx2(dt_fig2C,"Fig2C_list.xlsx",sheetName = "Fig2c",append=FALSE)

p2 <-ggplot(data=dt_fig2C_IgM, aes(x=Sample.Name,y=Score,fill=Channel,shape=Protein)) +
  geom_point(position = position_jitter(w=0.1),size=1,alpha=0.8)+
  scale_fill_manual(values = mycol)+
  scale_shape_manual(values=myshape)+
  facet_grid(Channel ~ Severity, scales = 'free_x')+
  ggtitle(paste0("Predicted Positive Ratios")) +
  labs(x = "COVID-19 Patients",y="Positivity Ratio")+
  #scale_x_continuous(expand = c(0, 0),breaks=seq(0.1:08,by=0.1),position = "top")+
  scale_y_continuous(limits = c(-0.01, 0.5), breaks = seq(0, 0.5, by = 0.05))+
  guides(fill=guide_legend(override.aes=list(colour=c("IgA"="#0072B2","IgG"="#D55E00","IgM"="#CC79A7"))))+
  #guides(fill=guide_legend(override.aes=list(colour=c("IgA"="#0072B2","IgG"="#D55E00","IgM"="#CC79A7"))))+
  theme_classic()+
  theme(
    aspect.ratio=1,
    legend.position = "bottom",
    panel.border = element_rect(fill=NA, size=0,colour="black"),
    panel.background = element_rect(fill = "white", colour = "grey50"),
    panel.grid.minor.x = element_blank(),
    panel.spacing.x = unit(0,"line"),
    #strip.placement = 'outside',
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(color="black", size=7, face="bold"),
    strip.text.y = element_text(color="black", size=7, face="bold"),
    plot.title = element_text(color="grey20", size=10, face="bold"),
    axis.title.y = element_text(color="grey50", size=9, face="bold"),
    axis.title.x = element_text(color="grey50", size=9, face="bold"),
    axis.text.x = element_text(size=5,angle = 60, vjust = 1, hjust=1),
    axis.text.y = element_text(size=5,angle = 0, vjust = 1, hjust=0.5),
    legend.title=element_text(size=5), 
    legend.text=element_text(size=4))
print(p2)
ggsave(paste0("Fig2D_Scatterplots_With_IgM","_",norm_name,".tiff"),plot=p2,width = 4, height = 4, dpi=600, compression = "lzw")

######################################################################################
######QC Plots
#######################################################################################
st=screen_table
st$Sample.Category=st$Sample_Type
st$Sample.Category=ifelse(grepl("^S2017",st$Sample.Name),"Negative",st$Sample.Category)
st$Sample.Category=ifelse(grepl("^S2020",st$Sample.Name),"COVID-19",st$Sample.Category)


pfname <- paste("All_QC",norm_name,".pdf",sep="")
pdf(pfname,width=15, height=8,family="Helvetica-Narrow")
old.par <- par(no.readonly = TRUE)

qcstat <- data.frame(Plate=character(),MEAN_NUCLEI_COUNT=numeric(),SD_NUCLEI_COUNT=numeric(),MEAN_DAPI_INTENSITY=numeric(),SD_DAPI_INTENSITY=numeric(),
                     Z_Prime_IgA=numeric(),Z_Prime_IgG=numeric(),Z_Prime_IgM=numeric(),
                     Mean_NEG1_IgA=numeric(),SD_NEG1_IgA=numeric(),Mean_NEG1_IgG=numeric(),SD_NEG1_IgG=numeric(),Mean_NEG1_IgM=numeric(),SD_NEG1_IgM=numeric(),
                     Mean_POS1_IgA=numeric(),SD_POS1_IgA=numeric(),Mean_POS1_IgG=numeric(),SD_POS1_IgG=numeric(),Mean_POS1_IgM=numeric(),SD_POS1_IgM=numeric(),
                     stringsAsFactors=FALSE)
barcodes=c("HX133-09-S01","HX133-09-S02","HX133-09-S03","HX133-09-S04","HX133-09-N01","HX133-09-N02","HX133-09-N03","HX133-09-N04",
           "HX133-09-R01","HX133-09-R02","HX133-09-R03","HX133-09-R04","HX133-09-M01","HX133-09-M02","HX133-09-M03","HX133-09-M04")
#unique(screen_table$Destination.Plate.Barcode)
for(bc in barcodes)
{
  plate_table=dt_long[dt_long$Destination.Plate.Barcode==bc,]
  plate_table=plate_table[grepl("Ratio.positive",plate_table$Feature),]
  plate_table=plate_table[!plate_table$Sample.Name %in% c("EMPTY","BLANK"),]
    plate_name=bc
    MEAN_NUCLEI_COUNT=mean(plate_table$Cell.Count)
    SD_NUCLEI_COUNT=sd(plate_table$Cell.Count)
    
    MEAN_DAPI_INTENSITY=mean(plate_table$Raw.Mean.NUCLEUS.INTENSITY.MEAN.DAPI)
    SD_DAPI_INTENSITY=sd(plate_table$Raw.Mean.NUCLEUS.INTENSITY.MEAN.DAPI)
    
    channel_table_IgA=plate_table[plate_table$Channel=="IgA",]
    channel_table_IgG=plate_table[plate_table$Channel=="IgG",]
    channel_table_IgM=plate_table[plate_table$Channel=="IgM",]

    pos_IgA= channel_table_IgA$Score[channel_table_IgA$Sample.Name %in% c("pos","POS","POS1")]
    neg_IgA= channel_table_IgA$Score[channel_table_IgA$Sample.Name %in% c("neg","NEG1")]

    pos_IgG= channel_table_IgG$Score[channel_table_IgG$Sample.Name %in% c("pos","POS","POS1")]
    neg_IgG= channel_table_IgG$Score[channel_table_IgG$Sample.Name %in% c("neg","NEG1")]
    
    pos_IgM= channel_table_IgM$Score[channel_table_IgM$Sample.Name %in% c("pos","POS","POS1")]
    neg_IgM= channel_table_IgM$Score[channel_table_IgM$Sample.Name %in% c("neg","NEG1")]
    
    plate_z_IgA	<-	zfactor(neg_IgA,pos_IgA)
    plate_z_IgG	<-	zfactor(neg_IgG,pos_IgG)
    plate_z_IgM	<-	zfactor(neg_IgM,pos_IgM)
    Mean_NEG1_IgA=round(mean(neg_IgA),4)
    Mean_NEG1_IgG=round(mean(neg_IgG),4)
    Mean_NEG1_IgM=round(mean(neg_IgM),4)
    SD_NEG1_IgA=round(pop.sd(neg_IgA),4)
    SD_NEG1_IgG=round(pop.sd(neg_IgG),4)
    SD_NEG1_IgM=round(pop.sd(neg_IgM),4)
    mean_POS1_IgA=round(mean(pos_IgA),4)
    mean_POS1_IgG=round(mean(pos_IgG),4)
    mean_POS1_IgM=round(mean(pos_IgM),4)
    sd_POS1_IgA=round(pop.sd(pos_IgA),4)
    sd_POS1_IgG=round(pop.sd(pos_IgG),4)
    sd_POS1_IgM=round(pop.sd(pos_IgM),4)
    qcplate <- data.frame(Plate=plate_name,MEAN_NUCLEI_COUNT=MEAN_NUCLEI_COUNT,SD_NUCLEI_COUNT=SD_NUCLEI_COUNT,
                          MEAN_DAPI_INTENSITY=MEAN_DAPI_INTENSITY,SD_DAPI_INTENSITY=SD_DAPI_INTENSITY,
                          Z_Prime_IgA=plate_z_IgA,Z_Prime_IgG=plate_z_IgG,Z_Prime_IgM=plate_z_IgM,
                          Mean_NEG1_IgA=Mean_NEG1_IgA,SD_NEG1_IgA=SD_NEG1_IgA,Mean_NEG1_IgG=Mean_NEG1_IgG,SD_NEG1_IgG=SD_NEG1_IgG,Mean_NEG1_IgM=Mean_NEG1_IgM,SD_NEG1_IgM=SD_NEG1_IgM,
                          Mean_POS1_IgA=mean_POS1_IgA,SD_POS1_IgA=sd_POS1_IgA,Mean_POS1_IgG=mean_POS1_IgG,SD_POS1_IgG=sd_POS1_IgG,Mean_POS1_IgM=mean_POS1_IgM,SD_POS1_IgM=sd_POS1_IgM,
                          stringsAsFactors=FALSE)
    qcstat=rbind(qcstat,qcplate)
  }

save(qcstat,file ="..\\..\\RDA\\QC.rda")
fname <- "QC_All_Plates.xlsx"
write.xlsx2(qcstat,fname,sheetName="Data",row.names=FALSE, append=FALSE)
par(mfrow=c(1,1),mar=c(0,0,2,0))

maxrow = 28;
npages = ceiling(nrow(qcstat)/maxrow);
for (i in 1:npages)
{
  idx = seq(1+((i-1)*maxrow), i*maxrow)
  qcs=qcstat[idx,]
  qcs=qcs[complete.cases(qcs),]
  colcs=qcs
  colcs$Plate<-"black"
  colcs$Z_Prime_IgA <- ifelse(colcs$Z_Prime_IgA < 0.5,"red","blue")
  colcs$Z_Prime_IgG <- ifelse(colcs$Z_Prime_IgG < 0.5,"red","blue")
  colcs$Z_Prime_IgM <- ifelse(colcs$Z_Prime_IgM < 0.5,"red","blue")
  colcs$MEAN_NUCLEI_COUNT <-"blue"
  colcs$SD_NUCLEI_COUNT <-"blue"
  colcs$MEAN_DAPI_INTENSITY<-"blue"
  colcs$SD_DAPI_INTENSITY <-"blue"
  colcs$Mean_NEG1_IgA <-"blue"
  colcs$SD_NEG1_IgA <-"blue"
  colcs$Mean_POS1_IgA <-"blue"
  colcs$SD_POS1_IgA <-"blue"
  colcs$Mean_NEG1_IgG <-"blue"
  colcs$SD_NEG1_IgG <-"blue"
  colcs$Mean_POS1_IgG <-"blue"
  colcs$SD_POS1_IgG <-"blue"
  colcs$Mean_NEG1_IgM <-"blue"
  colcs$SD_NEG1_IgM <-"blue"
  colcs$Mean_POS1_IgM <-"blue"
  colcs$SD_POS1_IgM <-"blue"
  textplot(qcs,col.data=as.matrix(colcs),show.rownames = FALSE,show.colnames = TRUE, col.colnames="red", valign="top",halign="left",hadj=0,cex=0.4,rmar=1.1)
  title("QC Summary:",cex.main=2,col.main="red",adj=0)
}
dev.off()

#############################
##########Plate Heatmaps############

plotlist = list()
plot_list <- lapply(unique(st$Destination.Plate.Barcode), function(bc)
{
  nanmat <- matrix(NaN,nrow=16,ncol=24,byrow=FALSE)
  rownames(nanmat) <- LETTERS[1:16]
  colnames(nanmat) <- 1:24
  nan_tbl <- melt(nanmat)
  colnames(nan_tbl) <- c("DRow","DCol")
  nan_tbl=nan_tbl[,c(1:2)]
  nan_tbl$DRow=as.character(nan_tbl$DRow)
  
  plate_table=st[st$Destination.Plate.Barcode==bc,]
  plate_table=merge(plate_table,nan_tbl,by=c("DRow","DCol"),all=TRUE)
  plate_table=plate_table[with(plate_table, order(DCol,DRow)), ]
  plate_table=plate_table[,c("DRow","DCol","Well","Cell.Count","Sample.Name","Sample_Type","Sample.Category")]
  plate_table=unique(plate_table)
  
  neg=plate_table$Cell.Count[plate_table$Sample.Category %in% c("NEG","Neg","NEG2","NEG1","NEG3")]
  cutoff=mean(neg)/10
  plate_table$Cell.Count_Samples=ifelse(plate_table$Sample_Type %in% c("EMPTY","BLANK","empty"),NA,plate_table$Cell.Count)
  plate_table$Cell.Count_Z=scale(plate_table$Cell.Count_Samples)
  plate=paste0("Cell Count : ",bc)
   
  collevels =LETTERS[1:16]
  p <- ggplot(plate_table,aes(x=DCol,y=ordered(DRow, levels=rev(collevels))))+
    geom_tile(aes(fill = Cell.Count_Z),colour="grey50")+
    geom_text(aes(label = Cell.Count_Samples))+
    scale_x_continuous(expand = c(0, 0),breaks=seq(1:24),position = "top")+
    scale_y_discrete(expand = c(0, 0),breaks=LETTERS[1:16])+
    scale_fill_gradient2(low="blue",high = "red", na.value="grey50" )+
    labs( title=plate,x="Columns",y="Rows")+
    theme_bw()+
    theme(legend.position = "none"  
          #panel.ontop = TRUE,
          #panel.background = element_rect(fill = "transparent")
    )
  print(p)
})

ml=marrangeGrob(grobs = plot_list, ncol=2,nrow=2,as.table = TRUE)
ggsave(paste0("Cell_Count_Plate_Heatmap","_",norm_name,".pdf"), ml, width = 24, height = 12,dpi=300, units = "in")

#######################################################################################################
##############################################################
###########Mean.NUCLEUS.INTENSITY.MEAN.DAPI_Z
plotlist = list()
unique(dt_long$Plate_set)
plot_list <- lapply(unique(screen_table$Destination.Plate.Barcode), function(bc)
{
  nanmat <- matrix(NaN,nrow=16,ncol=24,byrow=FALSE)
  rownames(nanmat) <- LETTERS[1:16]
  colnames(nanmat) <- 1:24
  nan_tbl <- melt(nanmat)
  colnames(nan_tbl) <- c("DRow","DCol")
  nan_tbl=nan_tbl[,c(1:2)]
  nan_tbl$DRow=as.character(nan_tbl$DRow)
  
  plate_table=st[st$Destination.Plate.Barcode==bc,]
  plate_table=merge(plate_table,nan_tbl,by=c("DRow","DCol"),all=TRUE)
  plate_table=plate_table[with(plate_table, order(DCol,DRow)),]
  plate_table=plate_table[,c("DRow","DCol","Well","Raw.Mean.NUCLEUS.INTENSITY.MEAN.DAPI","Cell.Count","Sample_Type")]
  plate_table=unique(plate_table)
  
  plate_table$Raw.Mean.NUCLEUS.INTENSITY.MEAN.DAPI_samples=ifelse(plate_table$Sample_Type %in% c("EMPTY","BLANK","empty"),NA,plate_table$Raw.Mean.NUCLEUS.INTENSITY.MEAN.DAPI)
  plate_table$Raw.Mean.NUCLEUS.INTENSITY.MEAN.DAPI_Z=scale(plate_table$Raw.Mean.NUCLEUS.INTENSITY.MEAN.DAPI_samples)
  
  plate=paste0("Raw.Mean.NUCLEUS.INTENSITY.MEAN.DAPI : ",bc)
  collevels =LETTERS[1:16]
  p <- ggplot(plate_table,aes(x=DCol,y=ordered(DRow, levels=rev(collevels))))+
    geom_tile(aes(fill = Raw.Mean.NUCLEUS.INTENSITY.MEAN.DAPI_Z),colour="grey50")+
    #geom_text(aes(label = round(Raw.Mean.NUCLEUS.INTENSITY.MEAN.DAPI_samples,1)))+
    geom_text(aes(label = round(Raw.Mean.NUCLEUS.INTENSITY.MEAN.DAPI_Z,2)))+
    scale_x_continuous(expand = c(0, 0),breaks=seq(1:24),position = "top")+
    scale_y_discrete(expand = c(0, 0),breaks=LETTERS[1:16])+
    scale_fill_gradient2(low="blue", high = "red", na.value="gray50" )+
    labs( title=plate,x="Columns",y="Rows")+
    theme_bw()+
    theme(legend.position = "none"  
          #panel.ontop = TRUE,
          #panel.background = element_rect(fill = "transparent")
    )
  print(p)
})

ml=marrangeGrob(grobs = plot_list, ncol=2,nrow=2,as.table = TRUE)
ggsave(paste0("Raw.MEAN_DAPI_Intensity_Plate_Heatmap","_",norm_name,".pdf"), ml, width = 24, height = 12,dpi=300, units = "in")

###########################################################################################################
######################################################
#############Control Barplots########################
plot_list=list()
for(plate_set in unique(dt_long$Plate_set))
{ 
  pst=dt_long[dt_long$Plate_set %in% c(plate_set),]
  for(prot in sort(unique(pst$Protein)))
  {  
    prot_table=pst[pst$Protein %in% c(prot),]
    prot_table=prot_table[!prot_table$Sample_Type %in% c("Sample","sample","CTRL2","CTRL1","empty","EMPTY","BLANK","anti-HIS/HA antibody","rabbit Anti-RBD and Anti-NP serum"),]
    selected_features=unique(grep(pattern ="Ratio", x = prot_table$Feature,value=TRUE))
    ###Only Positive
    selected_features=unique(grep(pattern ="Ratio.positive", x = prot_table$Feature,value=TRUE))
    prot_table=prot_table[prot_table$Feature %in% c(selected_features),]
    prot_table$Feature=gsub(".Ig[A|G|M]","",prot_table$Feature,perl=TRUE)
    p <- NULL
    p <-ggplot(prot_table,aes(Sample_Type,Score,fill=Feature))+
      geom_bar(position="Dodge",stat = "summary", fun= "mean")+
      stat_summary(fun.y=mean, fun.min=lowsd, fun.max=highsd, geom="errorbar", position=position_dodge(.9),color = 'black', size=.5, width=0.2) +
      facet_grid(Channel~Destination.Plate.Barcode,scales='free')+
      scale_fill_manual(values = mycol)+
      labs( title="Positive Ratios",x="Sample Type",y="Positive Ratio")+
      theme_bw() +
      theme(
        legend.position = "none",
        panel.border = element_rect(fill=NA, size=0,colour="black"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0,"line"),
        #strip.placement = 'outside',
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        strip.text.x = element_text(color="black", size=7, face="bold"),
        strip.text.y = element_text(color="black", size=7, face="bold"),
        plot.title = element_text(color="grey20", size=10, face="bold"),
        axis.title.y = element_text(color="grey50", size=9, face="bold"),
        axis.title.x = element_text(color="grey50", size=9, face="bold"),
        axis.text.x = element_text(size=5,angle = 0, vjust = 1, hjust=0.5),
        axis.text.y = element_text(size=5,angle = 0, vjust = 1, hjust=0.5))
    print(p)
    plot_list[[prot]]=p
  }
}

graphics.off()

layout_matrix=matrix(1:3, 1, 1, TRUE)

ml=marrangeGrob(grobs = plot_list, ncol=5,nrow=4,top=NULL,layout_matrix=layout_matrix)
ggsave(paste0("Control_Barplots","_",norm_name,".pdf"), ml, dpi=600, width = 9, height = 8, units = "in")
################################################
###Positive Ratios Heatmaps#################
################################################

dtl=dt_long
dtl$AgAb = paste0(dtl$Channel,"_",dtl$Protein)
dtl=dtl[!dtl$Sample.Name %in% c("CTRL1","CTRL2","EMPTY","BLANK","empty","anti-HIS/HA antibody","rabbit Anti-RBD and Anti-NP serum"),]
selected_features=unique(grep(pattern ="Ratio.positive", x = dtl$Feature,value=TRUE))
dtl=dtl[dtl$Feature %in% c(selected_features),]
dtl=dtl[,c("Sample.Name","Feature","Protein","Channel","AgAb","Score")]
dtl = dtl %>% 
  group_by(Sample.Name,AgAb) %>% 
  mutate(
    Score =mean(Score)
  ) %>% as.data.frame(dtl)

dtl = dtl %>% 
  group_by(Sample.Name) %>% 
  mutate(
    Avg_Score =mean(Score)
  ) %>% as.data.frame(dtl)

dtl=dtl[dtl$Channel %in% c("IgG","IgA"),]
names(dtl)[names(dtl) == "Score"] <- "Ratio.positive"
dtl=dtl[with(dtl,order(Avg_Score)),]
dtl_wider=unique(dtl)

p=ggplot(dtl_wider, aes(AgAb,Sample.Name)) +
  geom_tile(aes(fill = Ratio.positive)) +
  geom_text(aes(label = format(round(Ratio.positive,4),nsmall=2)))+
  scale_fill_gradient2( low = "blue", high = "red", na.value="gray",name = "")+
  ggtitle("Positive Ratios Heatmap") +
  labs(x = "Features", y = "Positive.Ratio") +scale_x_discrete(position = "top") +
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(hjust = 0.5))
pdf(paste0("Positive_Ratios_Heatmap_All_",norm_name,".pdf"), height=100, width=18)
print(p)
dev.off()
pl=ggplotly(p,width = 1400, height = 5000) 
cat(repr::repr_html(pl), file = paste0("Positive_Ratios_Heatmap_All","_",norm_name,".html"))

####################
#####Positive_Ratios New
require(gridExtra)
require(dplyr)

gg=dtl_wider %>% group_by(Channel) %>% 
  do(gg = {ggplot(., aes(AgAb,Sample.Name, fill = Ratio.positive)) + 
      geom_tile() + facet_grid(~Channel) + 
      geom_text(aes(label = format(round(Ratio.positive,4))))+
      scale_fill_gradient2( low = "blue", high = "red", na.value="gray",name = "")+
      #labs(title = "Positive Ratio Heatmap",x="Samples")+
      scale_x_discrete(position = "top") +
      guides(fill = guide_colourbar(title.position = "top")) +
      theme(
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            plot.title = element_text(color="grey90", face="bold",size=25,angle = 0, vjust = 1, hjust=0.5),
            axis.text.x = element_text(color="grey50", face="bold",size=15,angle = 0, vjust = 1, hjust=0.5),
            axis.title.x = element_blank(),
            legend.position = "top")}) %>% 
  .$gg %>% arrangeGrob(grobs = ., nrow = 1) %>% grid.arrange()

ggsave(paste0("Positive_Ratios_Heatmap","_",norm_name,".pdf"), gg, width = 24, height = 100, units = "in",limitsize = FALSE)
