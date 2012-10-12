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
#	MGD_DBUSER
#	MGD_DBPASSWORDFILE
#	GENOTYPE_INPUT_FILE
#	GENOTYPELOAD_OUTPUT
#	GENOTYPELOAD_MODE
#	OUTPUTDIR
#
# Inputs:
#
#	A tab-delimited file in the format (GENOTYPE_INPUT_FILE):
#
#	field 1:  Genotype Order #
#	field 2:  Genotype ID
#	field 3:  Strain ID
#	field 4:  Strain Name
#	field 5:  MGI Marker ID
#	field 6:  Allele 1 ID
#	field 7:  Mutant Cell Line 1
#	field 8:  Allele 2 ID
#	field 9:  Mutant Cell Line 2
#	field 10: Conditional (yes/no)
#	field 11: Exist As Term (ex. 'Mouse Line', 'Cell Line', 'Chimeric')
#	field 12: General Notes (ex. 1027)
#	field 13: Private Notes (ex. 1028)
#	field 14: Pair State (ex. 'Homozygous', 'Hemizygous X-linked', etc.)
#	field 15: Compound (ex. 'Top', 'Bottom', 'Not Applicable'
#	field 16: Created By User
#
# This assumes only one allele pair
# Once more than one allele pair is needed, revisions will need to be made
# GXD_Genotype.note field exist but has never been used and is always null.
# _MGIType_key = 12
# To add at a later date: Images (use imageload)
#
# Postprocessing:
#	exec GXD_orderAllelePairs : no necessary since we will order as part of this load
#	exec GXD_orderGenotypesAll (GXD_AlleleGenotype cache): done
#	allele combination ('allelecombination.csh'): done
#	OMIM cache: not done
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
import db
import mgi_utils
import loadlib
import alleleloadlib

#globals

#
# from configuration file
#
user = os.environ['MGD_DBUSER']
passwordFileName = os.environ['MGD_DBPASSWORDFILE']
mode = os.environ['GENOTYPELOAD_MODE']
inputFileName = os.environ['GENOTYPE_INPUT_FILE']
outputDir = os.environ['OUTPUTDIR']
genotypeOutput = os.environ['GENOTYPELOAD_OUTPUT']

accSetMax = 'exec ACC_setMax %d'

# this stored-procedure sets the GXD_AlleleGenotype information of a given genotype (delete/reload)
orderGenotypes = 'exec GXD_orderGenotypesAll %d'

# this product sets the MGI_Note/MGI_NoteChunk information
alleleCombination = os.environ['ALLCACHELOAD'] + '/allelecombinationByGenotype.py' + \
		' -S' + os.environ['MGD_DBSERVER'] + \
		' -D' + os.environ['MGD_DBNAME'] + \
		' -U' + os.environ['MGD_DBUSER'] + \
		' -P' + os.environ['MGD_DBPASSWORDFILE'] + \
		' -K%d\n'

DEBUG = 0		# if 0, not in debug mode
TAB = '\t'		# tab
CRT = '\n'		# carriage return/newline

bcpon = 1		# can the bcp files be bcp-ed into the database?  default is yes.

diagFile = ''		# diagnostic file descriptor
errorFile = ''		# error file descriptor
inputFile = ''		# file descriptor
genotypeFile = ''       # file descriptor
genotypeCacheFile = ''	# file descriptor
allelepairFile = ''	# file descriptor
accFile = ''            # file descriptor
noteFile = ''		# file descriptor
noteChunkFile = ''	# file descriptor

genotypeTable = ''	# table name
genotypeCacheTable = ''	# table name
allelepairTable = ''	# table name
accTable = ''		# table name
noteTable = ''		# table name
noteChunkTable = ''	# table name

genotypeFileName = ''	# bcp file
genotypeCacheFileName = ''	# bcp file
allelepairFileName = ''	# bcp file
accFileName = ''	# bcp file
noteFileName = ''	# bcp file
noteChunkFileName = ''	# bcp file
genotypeOutputName = ''	# output file name

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

runAlleleCombination = []

loaddate = loadlib.loaddate

skipBCP = 1		# does the BCP file load need to be skipped?

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

def initialize():
    global diagFile, errorFile, inputFile, errorFileName, diagFileName
    global genotypeFile, genotypeCacheFile, allelepairFile
    global accFile, noteFile, noteChunkFile
    global genotypeOutput

    global genotypeTable, genotypeCacheTable, allelepairTable
    global accTable, noteTable, noteChunkTable

    global genotypeFileName, genotypeCacheFileName, allelepairFileName
    global accFileName, noteFileName, noteChunkFileName, genotypeOutputName
 
    db.useOneConnection(1)
    db.set_sqlUser(user)
    db.set_sqlPasswordFromFile(passwordFileName)
 
    head, tail = os.path.split(inputFileName) 
    tail = tail + '.'

    genotypeTable = 'GXD_Genotype'
    genotypeCacheTable = 'GXD_AlleleGenotype'
    allelepairTable = 'GXD_AllelePair'
    accTable = 'ACC_Accession'
    noteTable = 'MGI_Note'
    noteChunkTable = 'MGI_NoteChunk'

    diagFileName = outputDir + '/' + tail + 'diagnostics'
    errorFileName = outputDir + '/' + tail + 'error'
    genotypeFileName = outputDir + '/' + tail + genotypeTable + '.bcp'
    genotypeCacheFileName = outputDir + '/' + tail + genotypeCacheTable + '.bcp'
    allelepairFileName = outputDir + '/' + tail + allelepairTable + '.bcp'
    accFileName = outputDir + '/' + tail + accTable + '.bcp'
    noteFileName = outputDir + '/' + tail + noteTable + '.bcp'
    noteChunkFileName = outputDir + '/' + tail + noteChunkTable + '.bcp'
    genotypeOutputName = genotypeOutput

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
        genotypeCacheFile = open(genotypeCacheFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % genotypeCacheFileName)

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
        genotypeOutput = open(genotypeOutputName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % genotypeOutputName)

    # Log all SQL
    db.set_sqlLogFunction(db.sqlLogAll)

    # Set Log File Descriptor
    db.set_sqlLogFD(diagFile)

    diagFile.write('Start Date/Time: %s\n' % (mgi_utils.date()))
    diagFile.write('Server: %s\n' % (db.get_sqlServer()))
    diagFile.write('Database: %s\n' % (db.get_sqlDatabase()))

    errorFile.write('Start Date/Time: %s\n\n' % (mgi_utils.date()))

    return 0

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

    return 0

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

    return 0

# Purpose:  BCPs the data into the database
# Returns:  nothing
# Assumes:  nothing
# Effects:  BCPs the data into the database
# Throws:   nothing

def bcpFiles():

    bcpdelim = "|"

    if DEBUG or not bcpon:
        return 0

    if skipBCP:
	return 0

    genotypeFile.close()
    genotypeCacheFile.close()
    allelepairFile.close()
    accFile.close()
    noteFile.close()
    noteChunkFile.close()

    bcpI = 'cat %s | bcp %s..' % (passwordFileName, db.get_sqlDatabase())
    bcpII = '-c -t\"|" -S%s -U%s' % (db.get_sqlServer(), db.get_sqlUser())

    bcp1 = '%s%s in %s %s' % (bcpI, genotypeTable, genotypeFileName, bcpII)
    bcp2 = '%s%s in %s %s' % (bcpI, genotypeCacheTable, genotypeCacheFileName, bcpII)
    bcp3 = '%s%s in %s %s' % (bcpI, allelepairTable, allelepairFileName, bcpII)
    bcp4 = '%s%s in %s %s' % (bcpI, accTable, accFileName, bcpII)
    bcp5 = '%s%s in %s %s' % (bcpI, noteTable, noteFileName, bcpII)
    bcp6 = '%s%s in %s %s' % (bcpI, noteChunkTable, noteChunkFileName, bcpII)

    for bcpCmd in [bcp1, bcp2, bcp3, bcp4, bcp5, bcp6]:
	diagFile.write('%s\n' % bcpCmd)
	os.system(bcpCmd)

    # run alleleCombination for each genotype added
    diagFile.write('%s\n' % runAlleleCombination)
    os.system(''.join(runAlleleCombination))

    return 0

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  verifies and processes each line in the input file
# Throws:   nothing

def processFile():

    global genotypeKey, allelepairKey, accKey, noteKey, mgiKey
    global runAlleleCombination
    global skipBCP

    lineNum = 0
    # For each line in the input file

    for line in inputFile.readlines():

        error = 0
        lineNum = lineNum + 1

        # Split the line into tokens
        tokens = line[:-1].split('\t')

        try:
	    genotypeOrder = tokens[0]
	    genotypeID = tokens[1]
	    strainID = tokens[2]
	    strainName = tokens[3]
	    markerID = tokens[4]
	    allele1ID = tokens[5]
	    mutant1ID = tokens[6]
	    allele2ID = tokens[7]
	    mutant2ID = tokens[8]
	    conditional = tokens[9]
	    existsAs = tokens[10]
	    generalNote = tokens[11]
	    privateNote = tokens[12]
	    pairState = tokens[13]
	    pairCompound = tokens[14]
	    createdBy = tokens[15]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

	# if the genotypeID is specified, then skip this row
	# this signifies that this genotype already exists in the system

	if len(genotypeID) > 0:
            genotypeOutput.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (\
	        genotypeOrder, genotypeID, strainID, strainName, markerID, allele1ID, allele2ID, conditional, \
	        existsAs, generalNote, privateNote, pairState, pairCompound, createdBy))
	    continue

	# strain key
	strainKey = loadlib.verifyObject(strainID, strainTypeKey, None, lineNum, errorFile)

	# marker key
	markerKey = loadlib.verifyMarker(markerID, lineNum, errorFile)

	# allele1 key
	allele1Key = loadlib.verifyObject(allele1ID, alleleTypeKey, None, lineNum, errorFile)

	# mutant1 key
	if len(mutant1ID) > 0:
	    mutant1Key = alleleloadlib.verifyMutnatCellLine(mutant1ID, lineNum, errorFile)
	else:
	    mutant1Key = ''

	# allele2 key
	if len(allele2ID) > 0:
	    allele2Key = loadlib.verifyObject(allele2ID, alleleTypeKey, None, lineNum, errorFile)
        else:
	    allele2Key = ''

	# mutant2 key
	if len(mutant2ID) > 0:
	    mutant2Key = alleleloadlib.verifyMutnatCellLine(mutant2ID, lineNum, errorFile)
	else:
	    mutant2Key = ''

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

        genotypeCacheFile.write('%s|%s|%s|1|%s|%s|%s|%s\n' \
            % (genotypeKey, markerKey, allele1Key, \
	    createdByKey, createdByKey, loaddate, loaddate))

        allelepairFile.write('%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
	    % (allelepairKey, genotypeKey, allele1Key, allele2Key, markerKey, \
	       mutant1Key, mutant2Key, \
	       pairStateKey, pairCompundKey, sequenceNum, \
	       createdByKey, createdByKey, loaddate, loaddate))

        # MGI Accession ID for the new genotype

        accFile.write('%s|%s%d|%s|%s|1|%d|%d|0|1|%s|%s|%s|%s\n' \
            % (accKey, mgiPrefix, mgiKey, mgiPrefix, mgiKey, genotypeKey, mgiTypeKey, \
	       createdByKey, createdByKey, loaddate, loaddate))

	# Print out a new text file and attach the new MGI Allele IDs as the last field

	genotypeID = mgiPrefix + str(mgiKey)

        genotypeOutput.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (\
	    genotypeOrder, genotypeID, strainID, strainName, markerID, allele1ID, allele2ID, conditional, \
	    existsAs, generalNote, privateNote, pairState, pairCompound, createdBy))

	# call allele-combinatin re-fresh for this genotype
	runAlleleCombination.append(alleleCombination % (genotypeKey))

        accKey = accKey + 1
        mgiKey = mgiKey + 1
	allelepairKey = allelepairKey + 1
        genotypeKey = genotypeKey + 1

	# don't skip the bcp file loading...data exists that needs to be loaded
	skipBCP = 0	

    #	end of "for line in inputFile.readlines():"

    #
    # Update the AccessionMax value
    #

    if not DEBUG and not skipBCP:
        db.sql(accSetMax % (lineNum), None)

    return 0

#
# Main
#

if initialize() != 0:
    exit(1)

if verifyMode() != 0:
    exit(1)

if setPrimaryKeys() != 0:
    exit(1)

if processFile() != 0:
    exit(1)

if bcpFiles() != 0:
    exit(1)

exit(0)

