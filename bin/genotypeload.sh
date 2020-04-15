#!/bin/sh
#
#  genotypeload.sh
###########################################################################
#
#  Purpose:
#
#      This script is a wrapper around the entire Allele load process.
#
#  Usage:
#
#      genotypeload.sh
#
#  Env Vars:
#
#      See the configuration file (genotypeload.config)
#
#  Inputs:  None
#
#  Outputs:
#
#      - Log file (${LOG_DIAG})
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  Fatal error occurred
#
#  Assumes:  Nothing
#
#  Implementation:
#
#      This script will perform following steps:
#
#      1) Source the configuration file to establish the environment.
#      2) Establish the log file.
#      3) Call Genotype load
#
#  Notes:  None
#
###########################################################################

cd `dirname $0`

CONFIG=$1
echo $CONFIG

#
# Make sure the configuration file exists and source it.
#
if [ -f ${CONFIG} ]
then
    . ${CONFIG}
else
    echo "Missing configuration file: ${CONFIG}"
    exit 1
fi

#
#  Source the DLA library functions.
#

if [ "${DLAJOBSTREAMFUNC}" != "" ]
then
    if [ -r ${DLAJOBSTREAMFUNC} ]
    then
        . ${DLAJOBSTREAMFUNC}
    else
        echo "Cannot source DLA functions script: ${DLAJOBSTREAMFUNC}" | tee -a ${LOG}
        exit 1
    fi
else
    echo "Environment variable DLAJOBSTREAMFUNC has not been defined." | tee -a ${LOG}
    exit 1
fi

#
# createArchive
#
if [ ${GENOTYPELOAD_STANDALONE} = 1 ]
then
    preload ${OUTPUTDIR}
fi

#
# Establish the log file.
#
LOG=${LOG_DIAG}

if [ ${GENOTYPELOAD_STANDALONE} = 1 ]
then
    rm -rf ${LOG}
    touch ${LOG}
fi

#
# Load Genotype file
#
echo "" >> ${LOG}
date >> ${LOG}
echo "Call genotypeload.py" | tee -a ${LOG}
${PYTHON} ./genotypeload.py 2>&1 >> ${LOG}
STAT=$?
checkStatus ${STAT} "Call genotypeload.py"

#
# these cache loads have been intentionally commented out - since this cache is run later on as part of weekly tasks
# while testing, turn these cahe loads back on
#
# Run Allele Combination load
#
#echo "" >> ${LOG}
#date >> ${LOG}
#echo "Call allelecombination.csh" | tee -a ${LOG}
#${ALLCACHELOAD}/allelecombination.csh
#STAT=$?
#checkStatus ${STAT} "Call genotypeload.py/allelecombination.csh"
#
# Run OMIM load
#
#echo "" >> ${LOG}
#date >> ${LOG}
#echo "Call mrkomim.csh" | tee -a ${LOG}
#${MRKCACHELOAD}/mrkomim.csh 2>&1 >> ${LOG}
#STAT=$?
#checkStatus ${STAT} "Call genotypeload.py/mrkomim.csh"

#
# run postload cleanup and email logs
#
if [ ${GENOTYPELOAD_STANDALONE} = 1 ]
then
    shutDown
fi

exit 0
