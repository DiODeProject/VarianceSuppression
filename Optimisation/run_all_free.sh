#!/bin/bash

###################################
# Synopsis of the script
###################################
EXPECTED_ARGS=0
if [ $# -lt ${EXPECTED_ARGS} ]; then
echo "This script generates N experiment files (.argos files)."
echo "Usage $0 <param> "
echo $'\t'"[MANDATORY] <param> "
exit


else

    #Experiment name           
    RHO_VALUES=('free' 100 1)
    #TIME_VALUES=(10 100)
    TIME_VALUES=( 1000 )
        
    #Looping variables
    echo RHO_VALUES=${RHO_VALUES[*]}
    echo TIME_VALUES=${TIME_VALUES[*]}

#	if [[ ${TIME_VALUES} -eq 50 ]]; then
#		HRS="00"
#		MIN="10"
#	else
#		HRS="00"
#		MIN="30"
#	fi
	HRS="168"
	#HRS="50"
	MIN="00"
	RUNJOB=runjob_l.sh
	
	EXP_NAME=free-slope
	#Set job name to distanguish between jobs on the queu
	sed -e "s|jobname|${EXP_NAME}|"   \
	    -e "s|min|${MIN}|"   \
	    -e "s|hrs|${HRS}|"   \
	    runjob_template.sh > ${RUNJOB}
	
	#Path variables   
	RESULT_DIR=${HOME}/IFD/${EXP_NAME} 
	mkdir -p ${RESULT_DIR}

		
	for TIME_VALUE in ${TIME_VALUES[*]}
	do
		for RHO in ${RHO_VALUES[*]}
		do

			COMMAND="qsub ${RUNJOB} ${TIME_VALUE} 0 ${RHO} ${RESULT_DIR}"
			${COMMAND}
			COMMAND="qsub ${RUNJOB} ${TIME_VALUE} free ${RHO} ${RESULT_DIR}"
			${COMMAND}
		done
	done
fi
