#!/bin/bash

SAMPLEFILE=$1
# new filenames: FlatRandomPtGunProducer_PDGid211_nPart1_Pt*_Eta1p6To2p8_noPU_predragm_cmssw1040pre1_20181119
# HGC TDR names: _SinglePiPt*Eta1p6_2p8_PhaseIITDRFall17DR-noPUFEVT_93X_upgrade2023_realistic_v2-v1_GEN-SIM-RECO
# Ilya's names:  FlatRandomPtGunProducer_gorbunov_PiPt*_PUP200_eta16to28_20171101

OUTDIR="/eos/user/p/predragm/Physics/HGCAL/ScaleResolution/resScaleInputs_20181120"
PREFIX="root://eoscms.cern.ch//eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/"
GUNTYPE="pt"
POSTFIX="NTUP"
REFNAME="genpart"
OBJNAMES="pfcluster"
#OBJNAMES="pfcluster megacluster"
QUEUE="8nh"

for SAMPLE in `cat $SAMPLEFILE`; do
  SAMPLEDIR="${PREFIX}${SAMPLE}/${POSTFIX}/"
  PARTICLE=`echo ${SAMPLE} | gawk 'match($0, /.*PDGid(.*)_nPart.*/, arr) { print arr[1] }'`
  PTVAL=`echo ${SAMPLE} | gawk -v part=".*PDGid${PARTICLE}_nPart1_Pt(.*)_Eta.*" 'match($0, part, arr) { print arr[1] }'`
  TAG=`echo ${SAMPLE} | gawk 'match($0, /.*Eta1p6To2p8_(.*)_predragm.*/, arr) { print arr[1] }'`
  PID=-1
  if [ "$PARTICLE" == "211" ]; then
    PID=211
  elif [ "$PARTICLE" == "22" ]; then
    PID=22
  fi;
  for OBJNAME in $OBJNAMES; do
    echo python resolutionScaleFiller_batchWrapper.py --inputdir ${SAMPLEDIR} --outdir ${OUTDIR} --gunType $GUNTYPE --pid $PID --genValue $PTVAL --tag $TAG --ref $REFNAME --obj $OBJNAME --queue $QUEUE
    python resolutionScaleFiller_batchWrapper.py --inputdir ${SAMPLEDIR} --outdir ${OUTDIR} --gunType $GUNTYPE --pid $PID --genValue $PTVAL --tag $TAG --ref $REFNAME --obj $OBJNAME --queue $QUEUE
  done
done
