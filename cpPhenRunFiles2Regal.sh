#!/bin/bash

#Script to copy runtime files from STORAGE (git repo) to regal run NODE
#  STORAGE location, which is git repo:  /n/wolkovich_lab/temporalvar/
#  run NODE where output is written:  /n/regal/wolkovich_lab/temporalvar/
#  Note that this overwrites any changes in the NODE versions

STORAGE="/n/wolkovich_lab/temporalvar"
NODE="/n/regal/wolkovich_lab/temporalvar"

cp ${STORAGE}/R/sourcefiles/*.R ${NODE}/R/sourcefiles
cp ${STORAGE}/R/getInputParms.txt ${NODE}/R
cp ${STORAGE}/R/PhenologyModel.R ${NODE}/R
