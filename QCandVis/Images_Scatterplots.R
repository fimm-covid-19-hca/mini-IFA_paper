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

# norm_name="Normalisation_NegOnly"
# results_dir=paste0(".\\Results\\",norm_name,"\\Plots")
#setwd(results_dir)
load(file="..\\..\\RDA\\screen_table.rda")

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

load(file="..\\..\\RDA\\img_table.rda")
st=screen_table
st$Sample.Category=st$Sample_Type
st$Sample.Category=ifelse(grepl("^S2017",st$Sample.Name),"Negative",st$Sample.Category)
st$Sample.Category=ifelse(grepl("^S2020",st$Sample.Name),"COVID-19",st$Sample.Category)
f_samples=c("F5_100","F5_25","F6_100","F6_25","F7_100","F7_25","F8_100","F8_25")
st=st[!st$Sample.Name %in% f_samples,]
##################################################################################
######
###Plate wise interactive Plots with images 
######
st$Sample_Name_Full=paste0(st$Sample.Name,"_",st$Sample)
st=st[!st$Sample_Type %in% c("EMPTY","BLANK","empty"),]
ig_plate_table=do.call("rbind.fill",lapply(as.character(unique(st$Destination.Plate.Barcode)),function(bc){
  print(bc)
  plate_table=st[st$Destination.Plate.Barcode %in% c(bc),]
  plate_table=plate_table[,c("Destination.Plate.Barcode","Sample_Name_Full","Sample.Name","Sample","Sample_Type","Protein","Ratio.positive.IgM","Ratio.positive.IgA","Ratio.positive.IgG")]
  
  st_long=gather(data=plate_table, key=Feature, value=Score,c(7:9), na.rm = FALSE, convert = FALSE)
  st_long=st_long[!is.na(st_long$Score),]
  
  st_long=ddply(st_long,.(Destination.Plate.Barcode,Sample_Name_Full,Sample.Name,Sample,Sample_Type,Protein,Feature),numcolwise(mean))
  st_wide=spread(st_long, key = Feature, value = Score)
  
  
  img_table_plate=img_table[img_table$Destination.Plate.Barcode==bc,]  
  img_table_plate=img_table_plate[,c("Sample_Name_Full","Image_base64" )]
  st_wide=merge(st_wide,img_table_plate,by=c("Sample_Name_Full"),all.x=TRUE)
  st_wide
}))

ig_table=ig_plate_table
ig_table=ig_table %>% mutate_if(is.numeric, ~round(., 5))%>%as.data.frame()

colnames(ig_table)=gsub("Ratio.positive.","",colnames(ig_table))

for(bc in sort(unique(ig_table$Destination.Plate.Barcode)))
{  
  plate_table=ig_table[ig_table$Destination.Plate.Barcode==bc,]
p1 <-ggplot(data=plate_table, aes(x=reorder(Sample_Name_Full,IgG), y=IgG,customdata = Image_base64,color=Sample_Type,shape=Sample_Type)) +
  geom_point(position = position_jitter(h=0.0001),size=5,alpha=0.6,aes(text = paste(
    "\n\t   Predicted Positive Ratios","",
    "\nSample Name : ", Sample_Name_Full,
    "\n\n  \t  \t   IgA   \t      IgG \t       IgM","\t ",
    "\n \t",Protein[1],"  \t  ",format(round(IgA,4),nsmall=2),"\t",format(round(IgG,4),nsmall=2),"\t",format(round(IgM,4),nsmall=2)
    )),show.legend = TRUE)+
  scale_shape_manual(values=myshape)+
  scale_colour_manual(values=mycol) +
  scale_x_discrete(expand=c(0,10))+
  ggtitle(paste0("Predicted Positive Ratios ",bc)) +
  labs(x = "Sample Names")+
  guides(fill=guide_legend(override.aes=list(size=5)))+
  guides(colour=guide_legend(override.aes=list(size=5)))+
  #guides(fill=guide_legend(override.aes=list(size=10)))+
  #guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
  theme_classic()+
  theme(
    panel.border = element_rect(fill=NA, size=2,colour="black"),
    panel.background = element_rect(fill = "white", colour = "grey50"),
    plot.title = element_text(color="grey20", size=25, face="bold"),
    axis.title.y = element_text(color="grey50", size=20, face="bold"),
    axis.title.x = element_text(color="grey50", size=20, face="bold"),
    #legend.title = element_blank(),
    #legend.position = c(0.7, 0.2),
    axis.text.x = element_text(size=4,angle = 45, vjust = 1, hjust=1))
print(p1)


p=ggplotly(p1,tooltip = c("text"), dynamicTicks = FALSE,
           layerData = 1, width = 1200, height = 900, originalData = FALSE)%>%
style(hoverlabel = list(align = "left")) %>% layout(legend = list(orientation = "v", x = 0.6, y = 0.95,font = list(size = 20))) %>%
htmlwidgets::onRender(readLines("..\\..\\code/tooltip-image.js"))
cat(repr::repr_html(p), file = paste0("IgG_",bc,"_",norm_name,"","_Images.html"))
}
################################################################################

