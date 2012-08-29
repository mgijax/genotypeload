#!/usr/local/bin/python

#
# Program: genotypeload.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To load new Genotypes into MGI
#
# Requirements Satisfied by This Program:
#
# Usage:
#	genotypeload.py
#
# Envvars:
#
# Inputs:
#
#	A tab-delimited file in the format:
#
#	field 1:  Strain ID
#	field 2:  Strain Name
#	field 3:  MGI Marker ID
#	field 4:  Allele 1 ID
#	field 5:  Allele 2 ID
#	field 6:  Conditional (yes/no)
#	field 7:  Exist As Term (ex. 'Mouse Line', 'Cell Line', 'Chimeric')
#	field 8:  Exist As Term (ex. 'Mouse Line', 'Cell Line', 'Chimeric')
#	field 9:  General Notes (ex. 1027)
#	field 10: Private Notes (ex. 1028)
#	field 11: Pair State (ex. 'Homozygous', 'Hemizygous X-linked', etc.)
#	field 12: Compound (ex. 'Top', 'Bottom', 'Not Applicable'
#	field 13: Creation Date
#
# This assumes only one allele pair
# Once more than one allele pair is needed, revisions will need to be made
# GXD_Genotype.note field exist but has never been used and is always null.
# _MGIType_key = 12
# To add at a later date: Images (use imageload)
#
# Postprocessing:
#	exec GXD_orderAllelePairs 
#	exec GXD_orderGenotypesAll (GXD_AlleleGenotype cache)
#	allele combination ('allelecombination.csh')
#	OMIM cache
#
# Outputs:
#
#       BCP files:
#
#       GXD_Genotype.bcp                master Genotype records
#       GXD_AllelePair.bcp              master Genotype records
#       GXD_AlleleGenotype.bcp          master Genotype records
#
#       MGI_Note/MGI_NoteChunk          notes
#
#       ACC_Accession.bcp               Accession records
#
#       Diagnostics file of all input parameters and SQL commands
#       Error file
#
# Exit Codes:
#
# Assumes:
#
#	That no one else is adding such records to the database.
#
# Bugs:
#
# Implementation:
#
# History
#
# 08/28/2012	lec
#	- TR10273/sanger MP annotations
#

import sys
import os
import accessionlib
import db
import mgi_utils
import loadlib

#globals

#
# from configuration file
#
user = os.environ['MGD_DBUSER']
passwordFileName = os.environ['MGD_DBPASSWORDFILE']
mode = os.environ['GENOTYPEMODE']
inputFileName = os.environ['GENOTYPEDATAFILE']
outputDir = os.environ['GENOTYPEDATADIR']

DEBUG = 0		# if 0, not in debug mode
TAB = '\t'		# tab
CRT = '\n'		# carriage return/newline

bcpon = 1		# can the bcp files be bcp-ed into the database?  default is yes.

diagFile = ''		# diagnostic file descriptor
errorFile = ''		# error file descriptor
inputFile = ''		# file descriptor
genotypeFile = ''       # file descriptor
allelepairFile = ''	# file descriptor
accFile = ''            # file descriptor
noteFile = ''		# file descriptor
noteChunkFile = ''	# file descriptor

genotypeTable = 'GXD_Genotype'
allelepairTable = 'ALL_Marker_Assoc'
accTable = 'ACC_Accession'
noteTable = 'MGI_Note'
noteChunkTable = 'MGI_NoteChunk'
newGenotypeFile = 'newAllele.txt'

genotypeFileName = outputDir + '/' + genotypeTable + '.bcp'
allelepairFileName = outputDir + '/' + allelepairTable + '.bcp'
accFileName = outputDir + '/' + accTable + '.bcp'
noteFileName = outputDir + '/' + noteTable + '.bcp'
noteChunkFileName = outputDir + '/' + noteChunkTable + '.bcp'
newGenotypeFileName = outputDir + '/' + newGenotypeFile

diagFileName = ''	# diagnostic file name
errorFileName = ''	# error file name

genotypeKey = 0         # GXD_Genotype._Genotype_key
allelepairKey =  0	# GXD_AllelePair._AllelePair_key
accKey = 0              # ACC_Accession._Accession_key
noteKey = 0		# MGI_Note._Note_key
mgiKey = 0              # ACC_AccessionMax.maxNumericPart
mgiNoteObjectKey = 12   # MGI_Note._MGIType_key
mgiNoteSeqNum = 1       # MGI_NoteChunk.sequenceNum
mgiGeneralNoteTypeKey = 1027   # MGI_Note._NoteType_key
mgiPrivateNoteTypeKey = 1028   # MGI_Note._NoteType_key

mgiTypeKey = 12		# ACC_MGIType._MGIType_key for Genotype
strainTypeKey = 10      # ACC_MGIType._MGIType_key for Strain
alleleTypeKey = 11      # ACC_MGIType._MGIType_key for Allele
mgiPrefix = "MGI:"

loaddate = loadlib.loaddate

# Purpose: prints error message and exits
# Returns: nothing
# Assumes: nothing
# Effects: exits with exit status
# Throws: nothing

def exit(
    status,          # numeric exit status (integer)
    message = None   # exit message (string)
    ):

    if message is not None:
        sys.stderr.write('\n' + str(message) + '\n')
 
    try:
        diagFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
        errorFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
        diagFile.close()
        errorFile.close()
	inputFile.close()
    except:
        pass

    db.useOneConnection(0)
    sys.exit(status)
 
# Purpose: process command line options
# Returns: nothing
# Assumes: nothing
# Effects: initializes global variables
#          exits if files cannot be opened
# Throws: nothing

def init():
    global diagFile, errorFile, inputFile, errorFileName, diagFileName
    global genotypeFile, allelepairFile
    global accFile, noteFile, noteChunkFile
    global newGenotypeFile
 
    db.useOneConnection(1)
    db.set_sqlUser(user)
    db.set_sqlPasswordFromFile(passwordFileName)
 
    head, tail = os.path.split(inputFileName) 

    diagFileName = outputDir + '/' + tail + '.diagnostics'
    errorFileName = outputDir + '/' + tail + '.error'

    try:
        diagFile = open(diagFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % diagFileName)
		
    try:
        errorFile = open(errorFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % errorFileName)
		
    try:
        inputFile = open(inputFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inputFileName)

    try:
        genotypeFile = open(genotypeFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % genotypeFileName)

    try:
        allelepairFile = open(allelepairFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % allelepairFileName)

    try:
        accFile = open(accFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % accFileName)

    try:
        noteFile = open(noteFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % noteFileName)

    try:
        noteChunkFile = open(noteChunkFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % noteChunkFileName)

    try:
        newGenotypeFile = open(newGenotypeFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % newGenotypeFileName)

    # Log all SQL
    db.set_sqlLogFunction(db.sqlLogAll)

    # Set Log File Descriptor
    db.set_sqlLogFD(diagFile)

    diagFile.write('Start Date/Time: %s\n' % (mgi_utils.date()))
    diagFile.write('Server: %s\n' % (db.get_sqlServer()))
    diagFile.write('Database: %s\n' % (db.get_sqlDatabase()))

    errorFile.write('Start Date/Time: %s\n\n' % (mgi_utils.date()))

    return

# Purpose: verify processing mode
# Returns: nothing
# Assumes: nothing
# Effects: if the processing mode is not valid, exits.
#	   else, sets global variables
# Throws:  nothing

def verifyMode():

    global DEBUG

    if mode == 'preview':
        DEBUG = 1
        bcpon = 0
    elif mode != 'load':
        exit(1, 'Invalid Processing Mode:  %s\n' % (mode))

# Purpose:  sets global primary key variables
# Returns:  nothing
# Assumes:  nothing
# Effects:  sets global primary key variables
# Throws:   nothing

def setPrimaryKeys():

    global genotypeKey, allelepairKey, accKey, noteKey, mgiKey

    results = db.sql('select maxKey = max(_Genotype_key) + 1 from GXD_Genotype', 'auto')
    genotypeKey = results[0]['maxKey']

    results = db.sql('select maxKey = max(_AllelePair_key) + 1 from GXD_AllelePair', 'auto')
    allelepairKey = results[0]['maxKey']

    results = db.sql('select maxKey = max(_Accession_key) + 1 from ACC_Accession', 'auto')
    accKey = results[0]['maxKey']

    results = db.sql('select maxKey = max(_Note_key) + 1 from MGI_Note', 'auto')
    noteKey = results[0]['maxKey']

    results = db.sql('select maxKey = maxNumericPart + 1 from ACC_AccessionMax ' + \
        'where prefixPart = "%s"' % (mgiPrefix), 'auto')
    mgiKey = results[0]['maxKey']

# Purpose:  BCPs the data into the database
# Returns:  nothing
# Assumes:  nothing
# Effects:  BCPs the data into the database
# Throws:   nothing

def bcpFiles():

    bcpdelim = "|"

    if DEBUG or not bcpon:
        return

    genotypeFile.close()
    allelepairFile.close()
    accFile.close()
    noteFile.close()
    noteChunkFile.close()

    bcpI = 'cat %s | bcp %s..' % (passwordFileName, db.get_sqlDatabase())
    bcpII = '-c -t\"|" -S%s -U%s' % (db.get_sqlServer(), db.get_sqlUser())

    bcp1 = '%s%s in %s %s' % (bcpI, genotypeTable, genotypeFileName, bcpII)
    bcp2 = '%s%s in %s %s' % (bcpI, allelepairTable, allelepairFileName, bcpII)
    bcp3 = '%s%s in %s %s' % (bcpI, accTable, accFileName, bcpII)
    bcp4 = '%s%s in %s %s' % (bcpI, noteTable, noteFileName, bcpII)
    bcp5 = '%s%s in %s %s' % (bcpI, noteChunkTable, noteChunkFileName, bcpII)

    for bcpCmd in [bcp1, bcp2, bcp3, bcp4, bcp5]:
	diagFile.write('%s\n' % bcpCmd)
	os.system(bcpCmd)

    return

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  verifies and processes each line in the input file
# Throws:   nothing

def processFile():

    global genotypeKey, allelepairKey, accKey, noteKey, mgiKey

    lineNum = 0
    # For each line in the input file

    for line in inputFile.readlines():

        error = 0
        lineNum = lineNum + 1

        # Split the line into tokens
        tokens = line[:-1].split('\t')

#	field 1:  Strain ID
#	field 2:  Strain Name
#	field 3:  MGI Marker ID
#	field 4:  Allele 1 ID
#	field 5:  Allele 2 ID
#	field 6:  Conditional (yes/no)
#	field 7:  Exist As Term (ex. 'Mouse Line', 'Cell Line', 'Chimeric')
#	field 8:  General Notes (ex. 1027)
#	field 9:  Private Notes (ex. 1028)
#	field 10: Pair State (ex. 'Homozygous', 'Hemizygous X-linked', etc.)
#	field 11: Compound (ex. 'Top', 'Bottom', 'Not Applicable'
#	field 12: Creation Date

        try:
	    strainID = tokens[0]
	    strainName = tokens[1]
	    markerID = tokens[2]
	    allele1ID = tokens[3]
	    allele2ID = tokens[4]
	    conditional = tokens[5]
	    existsAs = tokens[6]
	    generalNote = tokens[7]
	    privateNote = tokens[8]
	    pairState = tokens[9]
	    pairCompound = tokens[10]
	    createdBy = tokens[11]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

	# strain key
	strainKey = loadlib.verifyObject(strainID, strainTypeKey, None, lineNum, errorFile)

	# marker key
	markerKey = loadlib.verifyMarker(markerID, lineNum, errorFile)

	# allele1 key
	allele1Key = loadlib.verifyObject(allele1ID, alleleTypeKey, None, lineNum, errorFile)

	# allele2 key
	allele2Key = loadlib.verifyObject(allele2ID, alleleTypeKey, None, lineNum, errorFile)

	if conditional == 'yes':
	    conditionalKey = 1
	else:
	    conditionalKey = 0

	# _vocab_key = 60 (Genotype Exists As)
	existsAsKey = loadlib.verifyTerm('', 60, existsAs, lineNum, errorFile)

	# _vocab_key = 39 (Allele Pair State)
	pairStateKey = loadlib.verifyTerm('', 39, pairState, lineNum, errorFile)

	# _vocab_key = 42 (Allele Compound)
	pairCompundKey = loadlib.verifyTerm('', 42, pairCompound, lineNum, errorFile)

	createdByKey = loadlib.verifyUser(createdBy, lineNum, errorFile)

	if strainKey == 0 \
		or markerKey == 0 \
		or allele1Key == 0 \
		or allele2Key == 0 \
		or existsAsKey == 0 \
		or pairStateKey == 0 \
		or pairCompundKey == 0 \
		or createdByKey == 0:
	    error = 1

	# this tag only supports one allele pair per genotype
	sequenceNum = 1

        # if errors, continue to next record
        if error:
            continue

        # if no errors, process the allele

        genotypeFile.write('%s|%s|%s||%s|%s|%s|%s|%s\n' \
            % (genotypeKey, strainKey, conditionalKey, existsAsKey, \
	    createdByKey, createdByKey, loaddate, loaddate))

        allelepairFile.write('%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
	    % (allelepairKey, genotypeKey, allele1Key, allele2Key, markerKey, \
	       pairStateKey, pairCompundKey, sequenceNum, \
	       createdByKey, createdByKey, loaddate, loaddate))

        # MGI Accession ID for the new genotype

        accFile.write('%s|%s%d|%s|%s|1|%d|%d|0|1|%s|%s|%s|%s\n' \
            % (accKey, mgiPrefix, mgiKey, mgiPrefix, mgiKey, genotypeKey, mgiTypeKey, \
	       createdByKey, createdByKey, loaddate, loaddate))

	# Print out a new text file and attach the new MGI Allele IDs as the last field

        newGenotypeFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s\n' \
	    % (markerID, \
	    symbol, \
	    name, \
	    alleleStatus, \
	    alleleType, \
	    mgi_utils.prvalue(germLine), \
	    references, \
	    mgi_utils.prvalue(strainOfOrigin), \
	    mgi_utils.prvalue(molecularNotes), \
	    mgi_utils.prvalue(driverNotes), \
	    createdBy, mgiPrefix, mgiKey))

        accKey = accKey + 1
        mgiKey = mgiKey + 1
	allelepairKey = allelepairKey + 1
        genotypeKey = genotypeKey + 1

    #	end of "for line in inputFile.readlines():"

    #
    # Update the AccessionMax value
    #

    if not DEBUG:
        db.sql('exec ACC_setMax %d' % (lineNum), None)

#
# Main
#

init()
verifyMode()
setPrimaryKeys()
processFile()
bcpFiles()
exit(0)

