#!/bin/bash

SAMPLEFILE=$1

OUTDIR="/eos/user/c/clange/HGCal/ScaleResolution/skip_24_11_11/"
PREFIX="root://eoscms.cern.ch//eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/"
GUNTYPE="pt"
POSTFIX="NTUP"
REFNAME="genpart"
# OBJNAMES="pfcluster megacluster"
OBJNAMES="megacluster"
QUEUE="8nh"

for SAMPLE in `cat $SAMPLEFILE`; do
  SAMPLEDIR="${PREFIX}${SAMPLE}/${POSTFIX}/"
  PARTICLE=`echo ${SAMPLE} | gawk 'match($0, /.*Single(.*)Pt.*/, arr) { print arr[1] }'`
  PTVAL=`echo ${SAMPLE} | gawk -v part=".*Single${PARTICLE}Pt(.*)Eta.*" 'match($0, part, arr) { print arr[1] }'`
  TAG=`echo ${SAMPLE} | gawk 'match($0, /.*Fall17DR-(.*)FEVT.*/, arr) { print arr[1] }'`
  PID=0
  if [ "$PARTICLE" == "Pi" ]; then
    PID=211
  elif [ "$PARTICLE" == "Gamma" ]; then
    PID=22
  fi;
  for OBJNAME in $OBJNAMES; do
    echo python resolutionScaleFiller_batchWrapper.py --inputdir ${SAMPLEDIR} --outdir ${OUTDIR} --gunType $GUNTYPE --pid $PID --genValue $PTVAL --tag $TAG --ref $REFNAME --obj $OBJNAME --queue $QUEUE --skipLayers 3,10,17,24,34,46
    # --skipLayers 2,4,7,10,13,16,19,22,25,27,30,34,38,42,46,50
    python resolutionScaleFiller_batchWrapper.py --inputdir ${SAMPLEDIR} --outdir ${OUTDIR} --gunType $GUNTYPE --pid $PID --genValue $PTVAL --tag $TAG --ref $REFNAME --obj $OBJNAME --queue $QUEUE --skipLayers 3,10,17,24,34,46
    # --skipLayers 2,4,7,10,13,16,19,22,25,27,30,34,38,42,46,50
  done
done
