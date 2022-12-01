#!/usr/bin/python

######################################################################
# filterThomasParsed.py                                              #
# Author:  Dario Ghersi                                              #
# Version: 20180930                                                  #
# Goal:    process the .tsv parsed output from Paul Thomas' program  #
#          and retain the top clonotypes                             #
# Usage:   filterThomasParsed.py TSV NUM_TIMES CHAIN                 #
######################################################################

import sys

######################################################################
# FUNCTIONS                                                          #
######################################################################

def findPosChain(header, chain):
  """
  find the position in the header of the chain of interest
  """

  pos = -1
  fields = header[:-1].split("\t")
  for i in range(len(fields)):
    if (chain == "alpha" and fields[i] == "cdr3a") or\
       (chain == "beta" and fields[i] == "cdr3b"):
        pos = i

  return pos

######################################################################

def storeData(tsvFile, posChain):
  """
  count the frequency of each clonotype
  """

  freqID = {}
  for line in tsvFile:
    fields = line[:-1].split("\t")
    cdr = fields[posChain]

    # consider only productive clonotypes
    if cdr.find("#") == -1 and cdr.find("*") == -1:
      if freqID.has_key(cdr):
        freqID[cdr] += 1
      else:
        freqID[cdr] = 1
        
  return freqID
  
######################################################################
# MAIN PROGRAM                                                       #
######################################################################

## parse the parameters
if len(sys.argv) != 4:
  print "Usage: filterThomasParsed.py TSV NUM_TIMES CHAIN"
  sys.exit(1)
tsvFileName = sys.argv[1]
numTimesThresh = int(sys.argv[2])
chain = sys.argv[3]

## extract the position of the chain of interest
tsvFile = open(tsvFileName, "r")
header = tsvFile.readline()
posChain = findPosChain(header, chain)

## store the data
freqID = storeData(tsvFile, posChain)

## rewind the file
tsvFile.seek(0)

## print the filtered output
header = tsvFile.readline()
print header,
for line in tsvFile:
  fields = line[:-1].split("\t")
  cdr = fields[posChain]
  if len(cdr) > 1 and cdr.find("#") == -1 and cdr.find("*") == -1:
    if freqID[cdr] >= numTimesThresh:
      print line,

tsvFile.close()
