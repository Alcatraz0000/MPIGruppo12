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
        @file testReadWriteFile.c
*/

#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "RadixSort.h"

#define FILE_A "TestVectGruppo12"

/**
 * @brief This function tests the lectures from file of the first algorithm.
 * @param vect1      is the vector were isert element from the file.
 * @param length       size of the array "vect1".
 * @param mode          mode; 0 for sequential lecture, 1 for parallel.
 * @param rank          rank of the process.
 * @param num_process          total number of processes.
 */
void test_init_structures_algorithm0(int *vect1, int length, int mode, int rank, int num_process) {
    FILE *file = fopen(FILE_A, "w");
    fwrite(vect1, sizeof(int), length, file);
    fclose(file);

    int *result;
    init_structures(&result, length, mode, rank, num_process, FILE_A);
    if (rank == 0)
        for (int i = 0; i < length; i++)
            if (vect1[i] != result[i])
                assert(0);
}

/**
 * @brief This function tests the lectures from file of the first algorithm.
 * @param vect1      is the vector were isert element from the file.
 * @param length       size of the array "vect1".
 * @param mode          mode; 0 for sequential lecture, 1 for parallel.
 * @param rank          rank of the process.
 * @param num_process          total number of processes.
 */
void test_init_structures_algorithm1(int *vect1, int length, int mode, int rank, int num_process) {
    FILE *file = fopen(FILE_A, "w");
    fwrite(vect1, sizeof(int), length, file);
    fclose(file);

    int *result;
    init_structuresAlgo1(&result, length, rank, num_process, FILE_A);
    if (rank == 0)
        for (int i = 0; i < length; i++)
            if (vect1[i] != result[i])
                assert(0);
}

/**
 * @brief This main function call function for test the correctness of lectures from file functions.
 * @param argc      classic value not used.
 * @param argv       classic value not used.
 */
int main(int argc, char *argv[]) {
    
    printf("\nPREPARING VARIABLES FOR LECTURES FROM FILE'S TESTS\n");
    int rank, num_process;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);
    printf("I'm using %d process\n", num_process);
    int length = 27;
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

    int *result;
    printf("\nSTARTING TESTS FOR LECTURE FROM FILE FUNCTIONS:\n");

    printf("Starting tests for sequential lecture used in first algorithm\n");
    test_init_structures_algorithm0(a1, length, 0, rank, num_process);
    printf("...PASSED 1/3\n");
    test_init_structures_algorithm0(a2, length, 0, rank, num_process);
    printf("...PASSED 2/3\n");
    test_init_structures_algorithm0(a3, length, 0, rank, num_process);
    printf("...PASSED 3/3\n");
    

    printf("Starting tests for parallel lecture used in first algorithm\n");
    test_init_structures_algorithm0(a1, length, 1, rank, num_process);
    printf("...PASSED 1/3\n");
    test_init_structures_algorithm0(a2, length, 1, rank, num_process);
    printf("...PASSED 2/3\n");
    test_init_structures_algorithm0(a3, length, 1, rank, num_process);
    printf("...PASSED 3/3\n");
    

    printf("Starting tests for parallel lecture used in second algorithm\n");
    test_init_structures_algorithm1(a1, length, 1, rank, num_process);
    printf("...PASSED 1/3\n");
    test_init_structures_algorithm1(a2, length, 1, rank, num_process);
    printf("...PASSED 2/3\n");
    test_init_structures_algorithm1(a3, length, 1, rank, num_process);
    printf("...PASSED 3/3\n");
    

    printf("Tested ... Done");
    exit(EXIT_SUCCESS);
}
