/*
 * Course: High Performance Computing 2021/2022
 *
 * Lecturer: Francesco Moscato	fmoscato@unisa.it
 *
 * Group:
 * Marseglia	Mattia		0622701697	    m.marseglia1@studenti.unisa.it
 * Spingola     Camilla		0622701698  	c.spingola@studenti.unisa.it
 * Turi		    Vito		0622701795  	v.turi3@studenti.unisa.it
 * Sica 		Ferdinando	0622701794	    f.sica24@studenti.unisa.it
 *
 * Copyright (C) 2021 - All Rights Reserved
 *
 * This file is part of Contest-MPI: RadixSort.
 *
 * Contest-MPI: RadixSort is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Contest-MPI: RadixSort is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Contest-MPI: RadixSort.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
        @file main.c

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "RadixSortSeq.h"

int main(int argc, char **argv) {
    int *array;
    double algo_end_time = 0.0;
    double read_end_time = 0.0;
    int length, algorithm, max_digit;
    if (argc < 4) {
        printf("ERROR! Usage: ./main length algorithm max_digit");
        exit(1);
    }

    length = atoi(argv[1]);
    algorithm = atoi(argv[2]);
    max_digit = atoi(argv[3]);

    STARTTIME(1);
    init_structures(&array, length, "VectGruppo12");  // due tipi di input da file
    ENDTIME(1, read_end_time);

    if (algorithm == 0) {
        STARTTIME(2);
        myRadixsort(array, length);
        ENDTIME(2, algo_end_time);
    } else {
        STARTTIME(2);
        radix_sort(array, length);
        ENDTIME(2, algo_end_time);
    }
    printf("%d;0;%d;0;%.5f;%.5f\n", algorithm, length, read_end_time, algo_end_time);

    return 0;
}