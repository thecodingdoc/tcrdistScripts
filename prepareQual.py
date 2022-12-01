#!/usr/bin/python

import sys

######################################################################
# prepareQual.py                                                     #
# Author:  Dario Ghersi                                              #
# Version: 20171127                                                  #
# Goal:    process fasta and quality files and prepares an input     #
#          file compatible with TCR-DIST                             #
# Usage:   prepareQual.py FASTA QUAL OUTPUT                          #
######################################################################

######################################################################
# CONSTANTS                                                          #
######################################################################

EPITOPES = ("NAI", "GLC", "YVL")

######################################################################
# FUNCTIONS                                                          #
######################################################################

def storeFasta(fastaFileName):
  """
  """

  fastaSeqs = {}
  fastaFile = open(fastaFileName, "r")
  seq = ""
  for line in fastaFile:
    if line[0] == ">":
      if seq != "":
        fastaSeqs[header] = seq
      header = line[1:-1].split("_")[0]
      seq = ""
    else:
      seq += line[:-1]

  fastaSeqs[header] = seq
  fastaFile.close()

  return fastaSeqs

######################################################################

def storeQual(qualFileName):
  """
  """

  qualScores = {}
  qualFile = open(qualFileName, "r")
  scores = []
  for line in qualFile:
    if line[0] == ">":
      if len(scores) > 0:
        qualScores[header] = ".".join(scores)
      header = line[1:-1].split()[0].split("_")[0]
      scores = []
    else:
      temp = map(lambda x: int(x), line[:-1].split())
      temp = map(lambda x: str(x), temp)
      scores.extend(temp)

  qualFile.close()
  qualScores[header] = ".".join(scores)

  return qualScores

######################################################################

def printResults(fastaSeqs, qualScores, outFileName, subject):
  """
  print the results by matching alpha and beta
  """

  outFile = open(outFileName, "w")
  outFile.write("\t".join(("id", "epitope", "subject", "a_nucseq",\
                           "b_nucseq", "a_quals", "b_quals")) + "\n")
  
  ## match up the receptors
  for rec in fastaSeqs:
    if rec[:3] == "tra":
        
      # check whether a beta chain exists
      beta = rec.replace("tra", "trb")
      if fastaSeqs.has_key(beta):
          
        # identify the epitope
        currEp = ""
        for epitope in EPITOPES:
          if rec.upper().find(epitope) != -1:
            currEp = epitope
            break

        # write the pair
        toWrite = "\t".join((rec[4:], currEp, subject, fastaSeqs[rec],\
                            fastaSeqs[beta], qualScores[rec],\
                             qualScores[beta]))
        outFile.write(toWrite + "\n")
      else:
        print "Missing beta:" + rec      
  
  outFile.close()

######################################################################
# MAIN PROGRAM                                                       #
######################################################################

## parse the parameters
if len(sys.argv) != 5:
  print "Usage: prepareQual.py FASTA QUAL OUTPUT SUBJECT"
  sys.exit(1)
fastaFileName, qualFileName, outFileName, subject = sys.argv[1:]

## store the fasta sequences
fastaSeqs = storeFasta(fastaFileName)

## store the quality scores
qualScores = storeQual(qualFileName)

## print the results
printResults(fastaSeqs, qualScores, outFileName, subject)
