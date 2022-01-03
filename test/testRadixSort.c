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
        @file testRadixSort.c
*/

#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "RadixSort.h"

typedef void (*decorableMM)(int *, int, int, int);

/**
 * @brief This function tests the dot product between the matrix 'a' and the array 'b'.
 * @param expec      expected dot product result.
 * @param size       size of the result.
 * @param a          pointer to the matrix used in the dot product.
 * @param b          pointer to the array used in the dot product.
 * @param result     pointer to array used to store the result of dot product.
 * @param rows       number of rows.
 * @param columns    number of columns.
 * @param threads    number of threads.
 * @param dot        decorated dot product function.
 */
void test_radixsort(int *vect, int length, int num_process, int rank, decorableMM mm) {
    mm(vect, length, num_process, rank);
    printf("Starting raxisorttest\n");
    for (int i = 1; i < length; i++) {
        assert(vect[i - 1] <= vect[i]);
    }
}
void test_correct_elements(int *vect1, int *vect2, int length) {
    int i, j;
    int outOfRange = -4;
    for (i = 0; i < length; i++) {
        if (vect1[i] > outOfRange)
            outOfRange = vect1[i] + 1;
    }
    printf("Starting test_correct_elemetns\n");
    for (i = 0; i < length; i++) {
        for (j = 0; j < length; j++) {
            if (vect1[i] == vect2[j]) {
                vect2[j] = outOfRange;
                break;
            }
        }
        if (j == length)
            assert(0);
    }
}

int main(int argc, char *argv[]) {
    printf("Starting\n");
    int rank, num_process;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);
    int length = 27;
    printf("Initialized test variables\n");
    int a1[] = {15, 40, 65, 90, 0, 115, 30, 80, 130, 180, 230, 45, 120, 195, 270, 345, 0, 60, 160, 260, 360, 460, 75, 200, 325, 450, 575};
    int a2[length];
    int a3[length];
    for (int i = 0; i < length; i++) {
        if (i % 2)
            a3[i] = -a1[1];
        else
            a3[i] = a1[1];
        a2[i] = -a1[i];
    }

    int result[length];

    printf("Initialized test countings B\n");

    memcpy(result, a1, length * sizeof(int));
    test_radixsort(result, length, num_process, rank, myRadixsort);
    test_correct_elements(a1, result, length);

    memcpy(result, a2, length * sizeof(int));
    test_radixsort(result, length, num_process, rank, myRadixsort);
    test_correct_elements(a2, result, length);

    memcpy(result, a3, length * sizeof(int));
    test_radixsort(result, length, num_process, rank, myRadixsort);
    test_correct_elements(a3, result, length);

    memcpy(result, a1, length * sizeof(int));
    test_radixsort(result, length, num_process, rank, radix_sort);

    test_correct_elements(a1, result, length);

    memcpy(result, a2, length * sizeof(int));
    test_radixsort(result, length, num_process, rank, radix_sort);

    test_correct_elements(a2, result, length);

    memcpy(result, a3, length * sizeof(int));
    test_radixsort(result, length, num_process, rank, radix_sort);

    test_correct_elements(a3, result, length);

    printf("Tested mm ... Done");
    exit(EXIT_SUCCESS);
}
