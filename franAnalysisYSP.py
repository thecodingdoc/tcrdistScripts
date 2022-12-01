#!/usr/bin/python

######################################################################
# franAnalysisYSP.py                                                 #
# Author:  Dario Ghersi                                              #
# Goal:    Analyze TSV files and return the following:               #
#                                                                    #
#          1. # of nucleotide seqs per aa sequence                   #
#          2. # of nucleotide insertions                             #
#          3. # of glycines in the aa sequence                       #
#                                                                    #
#          The analysis is done at the epitope level, classifying    #
#          sequences as private or public if they are unique or      #
#          shared by more than one donor, respectively               #
# Version: 20190821                                                  #
# Usage:   franAnalysisYSP.py FOLDER                                 #
######################################################################

import glob
import sys

######################################################################
# CONSTANTS                                                          #
######################################################################

DNA_POS = 0
AA_POS = 1
FRAME_POS = 38
N1_INS_POS = 27
N2_INS_POS = 30

######################################################################
# FUNCTIONS                                                          #
######################################################################

def countOccurrenceClonotypes(folder):
  """
  store the occurrence of the amino acid sequence encoded by a
  particular nucleotide sequence. Note that a DNA sequence may
  only occur once, but receive a value > 1 if its aa sequence occurs
  multiple times
  """

  ## process the sequence files
  files = glob.glob(folder + "/*.tsv")
  aaSeqs = {}
  for item in files:
    infile = open(item, "r")
    header = infile.readline()
    for line in infile:
      fields = line[:-1].split("\t")
      dna = fields[DNA_POS]
      aa = fields[AA_POS]
      if fields[FRAME_POS] == "In":
        if aaSeqs.has_key(aa):
          if not item in aaSeqs[aa]:
            aaSeqs[aa].append(item)
        else:
          aaSeqs[aa] = [item]
    infile.close()

  ## count clonotype occurrence
  for item in aaSeqs:
    aaSeqs[item] = len(aaSeqs[item])

  return aaSeqs

######################################################################
# MAIN PROGRAM                                                       #
######################################################################

## parse the parameters
if len(sys.argv) != 2:
  print "Usage: franAnalysisYSP.py FOLDER"
  sys.exit(1)
folder = sys.argv[1]

## count the occurrence of each clonotype
clonotypeOcc = countOccurrenceClonotypes(folder)

## print the statistics for each clonotype
outputPrivateNT = open(folder + "_private_nt.txt", "w")
outputPublicNT = open(folder + "_public_nt.txt", "w")
outputPrivateNT.write("AA_SEQ\tNUM_NT\n")
outputPublicNT.write("AA_SEQ\tNUM_NT\n")
outputPublicIns = open(folder + "_public_ins.txt", "w")
outputPrivateIns = open(folder + "_private_ins.txt", "w")
outputPrivateIns.write("DNA_SEQ\tAA_SEQ\tNUM_INS\n")
outputPublicIns.write("DNA_SEQ\tAA_SEQ\tNUM_INS\n")
outputPrivateGly = open(folder + "_private_gly.txt", "w")
outputPublicGly = open(folder + "_public_gly.txt", "w")
outputPublicGly.write("AA_SEQ\tNUM_GLY\n")
outputPrivateGly.write("AA_SEQ\tNUM_GLY\n")

## process all files
files = glob.glob(folder + "/*.tsv")
nt = {}
insertions = {}
for item in files:
  infile = open(item, "r")
  header = infile.readline()
  for line in infile:
    fields = line[:-1].split("\t")
    dna = fields[DNA_POS]
    aa = fields[AA_POS]
    if fields[FRAME_POS] == "In":
      if not insertions.has_key(dna):
        insertions[dna] = [aa, str(int(fields[N1_INS_POS]) + int(fields[N1_INS_POS]))]
        
      if nt.has_key(aa):
        nt[aa] += 1
      else:
        nt[aa] = 1
  infile.close()

## print the gly and nt results
for aa in nt:
  if clonotypeOcc[aa] > 1:
    outputPublicNT.write(aa + "\t" + str(nt[aa]) + "\n")
    outputPublicGly.write(aa + "\t" + str(aa.count("G")) + "\n")
  else:
    outputPrivateNT.write(aa + "\t" + str(nt[aa]) + "\n")
    outputPrivateGly.write(aa + "\t" + str(aa.count("G")) + "\n")

## print the ins results
for dna in insertions:
  if clonotypeOcc[insertions[dna][0]] > 1:
    outputPublicIns.write(dna + "\t" + insertions[dna][0] + "\t" + insertions[dna][1] + "\n")
  else:
    outputPrivateIns.write(dna + "\t" + insertions[dna][0] + "\t" + insertions[dna][1] + "\n")

## close all files
outputPrivateNT.close()
outputPublicNT.close()
outputPrivateGly.close()
outputPublicGly.close()
outputPublicIns.close()
outputPrivateIns.close()
