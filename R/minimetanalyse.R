### Started 19 August 2019 ###
### By Lizzie ###

## Reviewing mini meta-analysis of tracking x traits papers ##

######################
### To do items!!! ###
######################
# 
######################

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## wd
setwd("~/Documents/git/projects/temporalvar/R")

## get part of the data
d <- read.csv("minimeta/accept_or_reject.csv", header=TRUE)

table(d$Accept.Reject) # this has 68 a, but below shows 69 -- WHY?
table(d$why)
# From the above I combined:
# model + theoretical model + theory/model - no data
# no phenology change measured + no trait or phenological change measured + no trat or phenological change measured
# meta-analysis + review - no data + no data
table(d$accept.reject_2) # this round is single species studies
table(d$why_3) # these are added to above....
# Traits not associated with phenological shifts - focus on herbivory ... need more info from Kelley
46+45+32+10+8+6+2 # studies removed

