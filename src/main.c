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
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "RadixSort.h"

int main(int argc, char **argv) {
    int rank, num_process;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);
    int *array;
    int *tmp;
    double algo_start_time, algo_end_time, read_start_time, read_end_time;
    int length, algorithm, max_digit, init_mode;

    if (argc < 5) {
        printf("ERROR! Usage: ./main length init_mode algorithm max_digit");
        exit(1);
    }

    length = atoi(argv[1]);
    init_mode = atoi(argv[2]);
    algorithm = atoi(argv[3]);
    max_digit = atoi(argv[4]);
    read_start_time = MPI_Wtime();
    if (algorithm == 0)
        init_structures(&array, length, init_mode, rank, num_process, "VectGruppo12");
    else {
        init_structuresAlgo1(&tmp, length, rank, num_process, "VectGruppo12");
    }
    read_end_time = MPI_Wtime();

    if (algorithm == 0) {
        algo_start_time = MPI_Wtime();
        myRadixsort(array, length, num_process, rank);
        algo_end_time = MPI_Wtime();
    } else {
        algo_start_time = MPI_Wtime();
        radix_sort(&array, tmp, length, num_process, rank);
        algo_end_time = MPI_Wtime();
    }
    double read_time = read_end_time - read_start_time;
    double algo_time = algo_end_time - algo_start_time;

    double global_read_time, global_algo_time;
    MPI_Reduce(&read_time, &global_read_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&algo_time, &global_algo_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) printf("%d;%d;%d;%d;%.5f;%.5f\n", algorithm, init_mode, length, num_process, global_read_time, global_algo_time);

    MPI_Finalize();
    return 0;
}