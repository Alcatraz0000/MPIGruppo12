/*
 * Course: High Performance Computing 2020/2021
 *
 * Lecturer: Francesco Moscato	fmoscato@unisa.it
 *
 * Group:
 *
 * Copyright (C) 2021 - All Rights Reserved
 *
 * This file is part of CommonAssignmentMPI01.
 *
 * CommonAssignmentMPI01 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CommonAssignmentMPI01 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CommonAssignmentMPI01.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
        @file test.c
*/

#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "RadixSort.h"

#define FILE_A "TestVectGruppo12"

void test_init_structures(int *vect1, int length, int mode, int rank, int num_process) {
    FILE *file = fopen(FILE_A, "w");
    fwrite(vect1, sizeof(int), length, file);
    fclose(file);

    int *result;
    init_structures(&result, length, mode, rank, num_process, FILE_A);
    if (rank == 0)
        for (int i = 0; i < length; i++) {
            if (vect1[i] != result[i])
                assert(0);
        }
}

int main(int argc, char *argv[]) {
    printf("Starting\n");
    int rank, num_process;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);
    printf("I'm using %d process", num_process);
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

    printf("Initialized test countings B\n");

    test_init_structures(a1, length, 0, rank, num_process);
    test_init_structures(a2, length, 0, rank, num_process);
    test_init_structures(a3, length, 0, rank, num_process);

    test_init_structures(a1, length, 1, rank, num_process);
    test_init_structures(a2, length, 1, rank, num_process);
    test_init_structures(a3, length, 1, rank, num_process);

    printf("Tested mm ... Done");
    exit(EXIT_SUCCESS);
}
