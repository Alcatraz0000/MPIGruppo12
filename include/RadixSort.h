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

#ifndef RADIXSORT_H_ 
#define RADIXSORT_H_


#include <mpi.h>
#include <time.h>

void init_structures(int **array, int length, int mode, int rank, int num_process, char *FILE_A);
void init_structuresAlgo1(int **array, int length, int rank, int num_process, char *FILE_A);
void write_on_File(int size, int max_digit, char *FILE_A);
void getMaxandMin(int *arr, int n, int *min, int *max);
void countingSortAlgo1(int *rec_buf, int digit, int rank, int dim, int min, int *count);
void radix_sort(int **glob_array, int *tmp_array, int n, int num_process, int rank);
void countingSortAlgo0(int array[], int base, int size, int raw_index, int *matrix);
void myRadixsort(int *array, int length, int num_process, int rank);
void serviceRadixsort(int array[], int size, int max, int rank, int num_process, MPI_Comm comm);
void getMaxDigitSeq(int array[], int length, int *array_pos, int *array_neg, int *max_pos, int *max_neg, int *pos, int *neg);
#endif
