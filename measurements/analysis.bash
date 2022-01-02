#!/bin/bash 

#
# Course: High Performance Computing 2021/2022
#
# Lecturer: Francesco Moscato	fmoscato@unisa.it
#
# Group:
# Marseglia	Mattia		0622701697	    m.marseglia1@studenti.unisa.it
# Spingola     Camilla		0622701698  	c.spingola@studenti.unisa.it
# Turi		    Vito		0622701795  	v.turi3@studenti.unisa.it
# Sica 		Ferdinando	0622701794	    f.sica24@studenti.unisa.it
#
# Copyright (C) 2021 - All Rights Reserved
#
# This file is part of Contest-OMP: RadixSort.
#
# Contest-OMP: RadixSort is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Contest-OMP: RadixSort is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Contest-OMP: RadixSort.  If not, see <http://www.gnu.org/licenses/>.
#

#specific the formats of time
TIME_STAMP=$(date +%s.%N)
TIMEFORMAT='%3U;%3E;%3S;%P'

#definitions of some variables used in this script:

#number of measurements to be made for each combination 
NUM_MEASURES=100

#dimension of item in program vector
VECT_DIMENSIONS=(5000000 20000000)

#number of threads used in our analysis to evaluate the performance variations
#with the different types of parallelized and non-parallelized algorithms.
#N.B. 0 is used for considerate serial execution 
NUM_PROCESS=(0 1 2 4 8 16)

#different options for compiler optimizations in back-end
COMP_OPT=(1 2 3)

#reference to programs 0 for radix sort based on counting sort, 1 for radix based on brutal algorithms
ALGORITHMS=(0 1)

# QUA DEVESCRIVERE CAMILLAAAAAreference to programs 0 for radix sort based on counting sort, 1 for radix based on brutal algorithms
INIT_MODE=(0 1)

#MAX_DIGIT saved all length of max digit that we want to try in measurements in loops operations
MAX_DIGIT=(9999 99999999)

#the path in which this script is placed
START_PATH=$(  cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P)

#bash directive to detecting the and of a script and reacting to it in a pre-calculated manner
trap "exit" INT

#function used to execute the program with different input value passed
execute(){


          # $1     2     3          4       5      6            7       8    
#           $dim $c_opt $init_mode $num_p $dest  $path prog $name_prog $algo  $max digit
    for ((i=0; i<$NUM_MEASURES; i++)); do
        if [[ $4 -eq 0 ]]; then
            program=$7_seq_O$2
            (time $6/$program $1 $8 $9) 2>&1 | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/;/g' -e 's/,/./g' -e 's/;/,/g' >> $5
        else
            program=$7_O$2
            (time mpirun -np $4 $6/$program $1 $3 $8 $9) 2>&1 | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/;/g' -e 's/,/./g' -e 's/;/,/g' >> $5
        fi
        

        printf "\r> %d/%d %3.1d%% " $(expr $i + 1) $NUM_MEASURES $(expr \( \( $i + 1 \) \* 100 \) / $NUM_MEASURES)
        printf "#%.0s" $(seq -s " " 1 $(expr \( $i \* 40 \) / $NUM_MEASURES))
		
   
    done
}

#function that execute code so many times which are the different number of max digit, 
#it establisced the path for every file changing throw parameters of size, optimization, 
#schedule and chunck-size passed
gen_execute(){
    for max_digit in ${MAX_DIGIT[@]}; do
           # $1     2   3           4       5               6           7           8
#           $dim $c_opt $init_mode $num_p $NAME_FILE_DEST $path_prog $name_prog $algo 
        DEST=$START_PATH/measures/ALGORITHM-$8/INIT_MODE-$3/SIZE-$1/MAX_DIGIT-$max_digit/OPTIMIZATION-O$2/$5
        
        #creation of necessary folders
        mkdir -p $(dirname $DEST) 2> /dev/null

        #for each configuration of number of threads, dimensions of rows and columns and for each compiler's optimization
        #the program will execute num_measures measurements and the results will be write on the correct destination file.
        
        #print the name of the file now being processed 
        echo -e "\n$DEST"
        #print on the file the parameters obtained from measurement, to be analyzed 
        echo "algo,init_mode,size,processes,init,funct,user,elapsed,sys,pCPU" >$DEST

        execute $1 $2 $3 $4 $DEST $6 $7 $8 $max_digit
    done

}



#function to create the file in format csv where the different values obtained 
#during the analysis have to be saved.
#The purpose is to obtain values of different times for each dimesion of vector, 
#for each number of threads and for each option of compiler's optimization.
#used to set how many execution we have to do changing parameters descibed up
generate(){

    path_prog=$1
    name_prog=$2
    
    for max_digit in ${MAX_DIGIT[@]}; do
        for dim in ${VECT_DIMENSIONS[@]}; do
            #create vector 
            $1/WriteVectOnFile $dim $max_digit

            for algo in ${ALGORITHMS[@]}; do
                for c_opt in ${COMP_OPT[@]}; do
                    for init_mode in ${INIT_MODE[@]}; do
                        for num_p in ${NUM_PROCESS[@]}; do
                            #definition of destination file
                            if [[ $num_p -eq 0 ]]; then
                                NAME_FILE_DEST=N_PROCESS-$num_p-SERIAL.csv
                            else
                                NAME_FILE_DEST=N_PROCESS-$num_p.csv
                            fi
                            
                            DEST=$START_PATH/measures/ALGORITHM-$algo/INIT_MODE-$init_mode/SIZE-$dim/MAX_DIGIT-$max_digit/OPTIMIZATION-O$c_opt/$NAME_FILE_DEST
                            
                            #creation of necessary folders
                            mkdir -p $(dirname $DEST) 2> /dev/null

                            #for each configuration of number of threads, dimensions of rows and columns and for each compiler's optimization
                            #the program will execute num_measures measurements and the results will be write on the correct destination file.
                            
                            #print the name of the file now being processed 
                            echo -e "\n$DEST"
                            #print on the file the parameters obtained from measurement, to be analyzed 
                            echo "algo,init_mode,size,processes,init,funct,user,elapsed,sys,pCPU" >$DEST

                            execute $dim $c_opt $init_mode $num_p $DEST $path_prog $name_prog $algo $max_digit
                        done
                    done
                done
            done
        done
    done
}





#if the script has been called with first_command command, is called
#the function generate

#$1 current binary dir, $2 project name
generate $1 $2
