options(java.home="C:\\Program Files\\Java\\jdk-11.0.2\\")
options(java.parameters = "-Xmx8g")
library(gdata)
library(MASS)
library(raster)
library(stats)
library(ggplot2)
library(gplots)
library(graphics)
library(xlsx)
library(rJava)
require(graphics)
require(biwt)
require(gridExtra)
library(reshape2)
library(plyr)
library(dplyr)
library(scales)
library(stringr)
library(extrafont)
library(tidyverse)

split_remove_last_n <- function(str,split_char,n){
  ###Usage: Splits the string with specified character and removes last n elements
  ###plate_name=split_remove_last_n(plate_id,"_",1)
  items <- unlist(strsplit(str,split_char))
  nlen=length(items)-n
  new_str <- paste(items[c(1:nlen)],collapse=split_char)
  return(new_str)
}
split_remove_last_n_vect <- Vectorize(split_remove_last_n, SIMPLIFY = TRUE)

split_remove_first_n <- function(str,split_char,n){
  ###Usage: Splits the string with specified character and removes first n elements
  ###plate_name=split_remove_last_n(plate_id,"_",1)
  items <- unlist(strsplit(str,split_char))
  nlen=n+1
  new_str <- paste(items[c(nlen:length(items))],collapse=split_char)
  return(new_str)
}
split_remove_first_n_vect <- Vectorize(split_remove_first_n, SIMPLIFY = TRUE)

split_give_last_n <- function(str,split_char,n){
  ###Usage: Splits the string with specified character and returns last n elements
  ###plate_name=split_remove_last_n(plate_id,"_",1)
  items <- unlist(strsplit(str,split_char))
  nlen=(length(items)-n) + 1
  new_str <- paste(items[c(nlen:length(items))],collapse=split_char)
  return(new_str)
}
split_give_last_n_vect <- Vectorize(split_give_last_n, SIMPLIFY = TRUE)

split_give_first_n <- function(str,split_char,n){
  ###Usage: Splits the string with specified character and returns last n elements
  ###plate_name=split_remove_last_n(plate_id,"_",1)
  items <- unlist(strsplit(str,split_char))
  new_str <- paste(items[c(1:n)],collapse=split_char)
  return(new_str)
}
split_give_first_n_vect <- Vectorize(split_give_first_n, SIMPLIFY = TRUE)


split_give_nth_starting_back <- function(str,split_char,n){
  ###Usage: Splits the string with specified character and returns last n elements
  ###plate_name=split_remove_last_n(plate_id,"_",1)
  items <- unlist(strsplit(str,split_char))
  if(length(items)>1)
  {
    nlen=(length(items)-n) + 1
    new_str <- items[nlen]
  }else{new_str=""}
  return(new_str)
}
split_give_nth_starting_back_vect <- Vectorize(split_give_nth_starting_back, SIMPLIFY = TRUE)

pi_return <- function(mean_pos,mean_neg,well_signal){
  result <- round (100* ((mean_neg- well_signal) / (mean_neg - mean_pos)),1)
  #PI	Input
  #pi_return(852,14137,10491)
  #pos1	852		10491
  #neg1	14137
  #result 27.44448626
  return(result)
}

split_give_nth_starting_front <- function(str,split_char,n){
  ###Usage: Splits the string with specified character and returns last n elements
  ###plate_name=split_give_nth_starting_front(p,"_",1)
  items <- unlist(strsplit(str,split_char))
  new_str <- items[n]
  return(new_str)
}
sn <- Vectorize(split_give_nth_starting_front, SIMPLIFY = TRUE)
split_give_nth_starting_front_vect <- Vectorize(split_give_nth_starting_front, SIMPLIFY = TRUE)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

ncharRight <- function(x,n){
  substr(x, nchar(x)-n+1,nchar(x)-n+1)
}


file_extension <- function(str){
  ###Usage: Splits the string with . and returns the last element which means file extension
  ###file_extension(a)
  items <- unlist(strsplit(str,"\\."))
  new_str=items[length(items)]
  return(new_str)
}

remove_file_extension <- function(str){
  ###Usage: Splits the string with . and returns the last element which means file extension
  ###file_extension(a)
  new_str=paste0("\\.",file_extension(str))
  result_str=gsub(new_str,"",str)
  return(result_str)
}

remove_directory_path <- function(file_name)
{
  fhead=unlist(strsplit(file_name,"/",fixed = T))
  fname=strsplit(fhead[length(fhead)],"/",fixed=T)[[1]]
  return(fname)
  
  #file_name="..\\Data\\910913-OS.xlsx"
  #file_name<-split_give_last_n(file_name,"\\\\",1)
  
}