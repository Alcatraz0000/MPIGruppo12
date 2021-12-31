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
 * This file is part of Contest-OMP: RadixSort.
 *
 * Contest-OMP: RadixSort is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Contest-OMP: RadixSort is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Contest-OMP: RadixSort.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
        @file RadixSort.c
*/

#include "RadixSort.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int init_structures(int **array, int length, int mode, int rank, int num_process, char *FILE_A) {  // implementazione I/O doppio da file

    int *tmp_array;

    if (mode == 0) {
        tmp_array = (int *)malloc(length * sizeof(int));
        if (tmp_array == NULL)
            perror("Memory Allocation - tmp_array");
        if (rank == 0) {
            FILE *file = fopen(FILE_A, "r");
            if (fread(tmp_array, sizeof(int), length, file) != length)
                perror("error during lecture from file");
            fclose(file);
            *array = tmp_array;
        }
    } else {
        int dim = length / num_process + 1;
        if (rank == 0)
            *array = (int *)malloc(length * sizeof(int));
        int *tmp_array = (int *)calloc(dim, sizeof(int));
        if (tmp_array == NULL)
            perror("Memory Allocation - tmp_array");
        MPI_File fh_a;
        MPI_Datatype dt_row_a;
        MPI_Type_contiguous(dim, MPI_INT, &dt_row_a);
        MPI_Type_commit(&dt_row_a);
        MPI_File_open(MPI_COMM_WORLD, FILE_A, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh_a);
        int displacement = rank * (dim) * sizeof(int);

        MPI_File_set_view(fh_a, displacement, MPI_INT, dt_row_a, "native", MPI_INFO_NULL);
        if (MPI_File_read(fh_a, tmp_array, 1, dt_row_a, MPI_STATUS_IGNORE) != MPI_SUCCESS)
            perror("error during lecture from file with MPI");
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(tmp_array, dim, MPI_INT, *array, dim, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
    }

    return 1;
}

void getMaxDigitSeq(int array[], int size, int *array_pos, int *array_neg, int *max_pos, int *max_neg, int *pos, int *neg) {
    *max_pos = 0;
    *max_neg = 0;
    *neg = 0;
    *pos = 0;
    int tmp_neg = 0, tmp_pos = 0, i = 0;
    for (i = 0; i < size; i++) {
        if (array[i] < 0) {
            array_neg[*neg] = array[i];
            (*neg)++;
        }
        if (array[i] >= 0) {
            array_pos[*pos] = array[i];
            (*pos)++;
        }
    }

    for (i = 0; i < size; i++) {
        if (array[i] < tmp_neg)
            tmp_neg = array[i];
        if (array[i] > tmp_pos)
            tmp_pos = array[i];
    }

    while (tmp_neg < 0) {
        tmp_neg /= 10;
        (*max_neg)++;
    }
    while (tmp_pos > 0) {
        tmp_pos /= 10;
        (*max_pos)++;
    }
}

void serviceRadixsort(int array[], int size, int max, int rank, int num_process, MPI_Comm comm) {
    // Get maximum element
    int base = 10;
    int *vect = (int *)calloc(((base + 1) * max), sizeof(int));
    int *vect_local = (int *)calloc((base + 1), sizeof(int));
    int i = 0, j = 0;
    MPI_Request request;
    int dest = rank % 2;
    // Apply counting sort to sort elements based on place value.
    for (int i = rank; i < max; i += num_process) {
        for (int k = 0; k < base + 1; k++)
            vect_local[k] = 0;
        countingSortAlgo0(array, base, size, i, vect_local);

        if (rank == 0) {
            for (int k = 0; k < base + 1; k++) {
                vect[i * (base + 1) + k] = vect_local[k];
                vect_local[k] = 0;
            }
            for (j = i + 1; j < max && j < i + num_process; j++) {
                MPI_Recv(vect_local, base + 1, MPI_INT, j % num_process, 0, comm, MPI_STATUS_IGNORE);
                for (int k = 0; k < base + 1; k++) {
                    vect[j * (base + 1) + k] = vect_local[k];
                }
            }
        } else {
            MPI_Ssend(vect_local, base + 1, MPI_INT, 0, 0, comm);
        }
    }

    MPI_Barrier(comm);

    int place = 1;
    if (rank == 0) {
        int *output = (int *)calloc(size, sizeof(int));

        for (i = 0; i < max; i++) {
            for (j = size - 1; j >= 0; j--) {
                output[vect[i * (base + 1) + (((array[j] / place) % base) - vect[i * (base + 1) + base])] - 1] = array[j];
                vect[i * (base + 1) + (((array[j] / place) % base) - vect[i * (base + 1) + base])]--;
            }

            for (int j = 0; j < size; j++) {
                array[j] = output[j];
            }
            place = place * 10;
        }
    }
    MPI_Barrier(comm);
    MPI_Bcast(array, size, MPI_INT, 0, comm);

    MPI_Barrier(comm);
}

/**

 * @brief This function allows to find the maximum in an array.

 * @param arr      array.

 * @param n        array size.

 */

int getMaxandMin(int *arr, int n, int *min, int *max) {
    *min = arr[0];

    *max = arr[0];

    for (int i = 1; i < n; i++) {
        if (arr[i] > *max)

            *max = arr[i];

        if (arr[i] < *min)

            *min = arr[i];
    }
}

void countingSortAlgo0(int array[], int base, int size, int raw_index, int *matrix) {
    int place = 1;

    for (int i = 0; i < raw_index; i++)
        place = place * 10;

    int max = (array[0] / place) % base;
    int min = (array[0] / place) % base;
    int i;
    for (int i = 1; i < size; i++) {
        if (((array[i] / place) % base) > max)
            max = ((array[i] / place) % base);

        if (((array[i] / place) % base) < min)
            min = ((array[i] / place) % base);
    }

    matrix[base] = min;

    int length = max - min + 1;
    for (int j = 0; j < size; j++)
        matrix[(((array[j] / place) % base) - min)]++;

    // Calculate cumulative count
    for (int j = 1; j < length; j++)
        matrix[j] += matrix[(j - 1)];
}

/**
 * @brief The main function that sorts the array of size n.
 * @param array          array.
 * @param rec_buf        sub-array of each process.
 * @param n              array size.
 * @param digit          number which represent ciphers.
 * @param num_process    number of processes.
 * @param rank           rank of the current process.
 * @param dim            dimension of the rec_buf.
 */
void countingSortAlgo1(int *array, int *rec_buf, int n, int digit, int num_process, int rank, int dim, int min, int *count) {
    // Compute local count for each processes
    int i, position, local_count[10] = {0};
    for (i = 0; i < dim; i++) {
        local_count[((rec_buf[i] - min) / digit) % 10]++;
    }
    // Reduce all the sub counts to root process
    if (rank == 0) {
        MPI_Reduce(local_count, count, 10, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        for (i = 1; i < 10; i++) {
            count[i] += count[i - 1];
        }

    } else {
        MPI_Reduce(local_count, 0, 10, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
}
/**
 * @brief The function that separates the array in subarray and starts the sorting process.
 * @param array          array.
 * @param n              array size.
 * @param num_process    number of processes.
 * @param rank           rank of the current process.
 */
void radix_sort(int *array, int n, int num_process, int rank) {
    int rem = n % num_process;  // elements remaining after division among processes
    int dim, displacement;

    if (rank < rem) {
        dim = n / num_process + 1;  // ne metto uno in piu a tutti
        displacement = rank * dim;  // a che displacement scrivere
    } else {
        dim = n / num_process;
        displacement = rank * dim + rem;
    }

    int *rec_buf = (int *)malloc(sizeof(int) * dim);
    int *sendcounts = NULL;
    int *displs = NULL;
    if (rank == 0) {
        sendcounts = malloc(sizeof(int) * num_process);
        displs = malloc(sizeof(int) * num_process);
    }
    MPI_Gather(&dim, 1, MPI_INT, sendcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);       // tutti mandano allo 0 il numero degli elementi da loro ordinati
    MPI_Gather(&displacement, 1, MPI_INT, displs, 1, MPI_INT, 0, MPI_COMM_WORLD);  // tutti mandano allo 0 il proprio displacment

    MPI_Scatterv(array, sendcounts, displs, MPI_INT, rec_buf, dim, MPI_INT, 0, MPI_COMM_WORLD);  // scatters a buffer in parts to all processes in a comunicator
    // ora ogni processo avrà il suo insieme di elementi da ordinare
    if (rank == 0) {
        free(sendcounts);
        free(displs);
    }
    // ogi processo calcolerà il massimo tra i suoi elementi
    int local_max;
    int local_min;

    getMaxandMin(array, n, &local_min, &local_max);

    int global_max, global_min;
    // ora bisogna calcolare un massimo globale tra tutti i processi
    MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_min, &global_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    // ciascun processo chiama la counting sort per quante sono le cifre decimali del max elemento globale
    int max_pos = 0;
    int tmp_pos = global_max - global_min;
    while (tmp_pos > 0) {
        tmp_pos /= 10;
        max_pos++;
    }

    int *frequencies[max_pos];
    if (rank == 0) {
        for (int i = 0; i < max_pos; i++) {
            frequencies[i] = (int *)calloc(10, sizeof(int));
        }
    }

    int decimal_digit = 0;
    for (int digit = 1; (global_max - global_min) / digit > 0; digit *= 10) {
        countingSortAlgo1(array, rec_buf, n, digit, num_process, rank, dim, global_min, frequencies[decimal_digit]);
        decimal_digit++;
    }

    if (rank == 0) {
        int *temp_array = (int *)malloc(sizeof(int) * n);
        int val = 1;
        for (int j = 0; j < max_pos; j++) {
            for (int i = n - 1; i >= 0; i--) {
                temp_array[frequencies[j][((array[i] - global_min) / val) % 10] - 1] = array[i];
                frequencies[j][((array[i] - global_min) / val) % 10]--;
            }
            val *= 10;
            memcpy(array, temp_array, sizeof(int) * n);
        }
        free(temp_array);
    }
    free(rec_buf);
}
void myRadixsort(int *array, int length, int num_process, int rank) {
    int *array_pos = (int *)malloc(length * sizeof(int));
    if (array_pos == NULL)
        perror("Memory Allocation - a");
    int *array_neg = (int *)malloc(length * sizeof(int));
    if (array_neg == NULL)
        perror("Memory Allocation - a");
    int max_pos, max_neg;
    int size_pos, size_neg;
    int old_rank = rank;
    int old_num_process = num_process;
    MPI_Comm pari;
    MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank, &pari);
    MPI_Comm_rank(pari, &rank);
    MPI_Comm_size(pari, &num_process);
    if (old_rank == 0)
        getMaxDigitSeq(array, length, array_pos, array_neg, &max_pos, &max_neg, &size_pos, &size_neg);
    MPI_Barrier(MPI_COMM_WORLD);

    if (old_num_process > 1) {
        MPI_Bcast(&max_neg, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&max_pos, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&size_neg, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&size_pos, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(array_neg, size_neg, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(array_pos, size_pos, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (old_num_process <= 1) {
        serviceRadixsort(array_neg, size_neg, max_neg, rank, old_num_process, MPI_COMM_WORLD);

        memcpy(array, array_neg, size_neg * sizeof(int));

        serviceRadixsort(array_pos, size_pos, max_pos, rank, old_num_process, MPI_COMM_WORLD);
        memcpy(array + size_neg, array_pos, size_pos * sizeof(int));
    } else {
        // ho sicuramente >2 processi --> need to bcast array neg e array pos

        if ((old_rank % 2) == 0) {
            // rank 0 master negativi sicur
            serviceRadixsort(array_neg, size_neg, max_neg, rank, num_process, pari);
            MPI_Barrier(pari);
            if (old_rank == 0) {
                memcpy(array, array_neg, size_neg * sizeof(int));
                MPI_Recv(array_pos, size_pos, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                memcpy(array + size_neg, array_pos, size_pos * sizeof(int));
            }
        } else {
            // rank 1 master poisitivi sicuro
            serviceRadixsort(array_pos, size_pos, max_pos, rank, num_process, pari);
            if (old_rank == 1) {
                MPI_Send(array_pos, size_pos, MPI_INT, 0, 0, MPI_COMM_WORLD);
            }
        }
    }
}
