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

void init_structures(int **array, int length, int mode, int rank, int num_process, char *FILE_A) {
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
        MPI_Gather(tmp_array, dim, MPI_INT, *array, dim, MPI_INT, 0, MPI_COMM_WORLD);
    }
}

/**
 * @brief This function initializes the array structure needed in the program. It uses two modalities to take inputs depending on the value of mode parameter.
 * With mode 0 the inputs are read from file sequentially otherwise they are read in a parallel mode with MPI file view.
 * @param array         pointer to the array to sort in the sorting algorithm.
 * @param length        dimension of array, so the number of elements contained in the array.
 * @param mode          reading from file modality.
 * @param rank          rank of the process that is executing the function.
 * @param num_process   the number of processes that are executing in parallel this function.
 * @param FILE_A        the name of the file to read on.
 */
void init_structuresAlgo0(int **array, int length, int rank, int num_process, char *FILE_A) {  // implementazione I/O doppio da file

    int *tmp_array;

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
}

void init_structuresAlgo1(int **array, int length, int rank, int num_process, char *FILE_A) {  // implementazione I/O doppio da file
    int dim;
    if (rank == 0) {
        dim = length / num_process + length % num_process;
    } else {
        dim = length / num_process;
    }
    *array = (int *)calloc(dim, sizeof(int));
    if (*array == NULL)
        perror("Memory Allocation - tmp_array");
    MPI_File fh_a;
    MPI_Datatype dt_row_a;

    MPI_Type_contiguous(dim, MPI_INT, &dt_row_a);
    MPI_Type_commit(&dt_row_a);
    MPI_File_open(MPI_COMM_WORLD, FILE_A, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh_a);
    int displacement;
    if (rank != 0) displacement = (rank * dim + length % num_process) * sizeof(int);
    if (rank == 0) displacement = (rank * dim) * sizeof(int);

    MPI_File_set_view(fh_a, displacement, MPI_INT, dt_row_a, "native", MPI_INFO_NULL);
    if (MPI_File_read(fh_a, *array, 1, dt_row_a, MPI_STATUS_IGNORE) != MPI_SUCCESS)
        perror("error during lecture from file with MPI");

    // MPI_Gather(tmp_array, dim, MPI_INT, *array, dim, MPI_INT, 0, MPI_COMM_WORLD);
}

/**
 * @brief This function calculates the number of positive and negative elements in the array and puts it in the corresponding arrays.
 * It also calculates the maximum positive element and the minimum negative element.
 * @param array         pointer to the array with total elements to sort.
 * @param size          dimension of the array containing all elements.
 * @param array_pos     pointer to array used to store the positive elements.
 * @param array_neg     pointer to array used to store the negative elements.
 * @param max_pos       pointer to the number of digits of the maximum positive element.
 * @param max_neg       pointer to the number of digits of the minimum negative element.
 * @param pos           pointer to the number of positive elements.
 * @param neg           pointer to the number of negative elements.
 */
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

/**
 * @brief This function implements the radix sort algorithm for algorithm 0 for a only positive or only negative array, based on counting sort algorithm,
 * storing during the computation the frequencies of each element in a vector. Every rank executes a part of the for loop with the corresponding parameters.
 * The rank 0 process, after each call of the counting sort, saves the results in a vector, restores 0 values in its local vector and waits for other processes sends of the results.
 * @param array       pointer to the array with the elements to sort, that are only positive or negative.
 * @param size        dimension of the subarray containing all positive or negative elements.
 * @param max         maximum number of digits contained in the longer element of this subarray
 * @param rank        rank of the process that is executing the function.
 * @param num_process the number of processes that are executing in parallel this function.
 * @param comm        the Communicator of the MPI processes that are executing the function.
 */
void serviceRadixsort(int array[], int size, int max, int rank, int num_process, MPI_Comm comm) {
    // Get maximum element
    int base = 10;
    int *vect;
    if (rank == 0) vect = (int *)calloc(((base + 1) * max), sizeof(int));
    int *vect_local = (int *)calloc((base + 1), sizeof(int));
    int i = 0, j = 0;

    // Apply counting sort to sort elements based on place value.
    for (int i = rank; i < max; i += num_process) {
        memset(vect_local, 0, base * sizeof(int));
        countingSortAlgo0(array, base, size, i, vect_local);
        if (rank != 0) {
            MPI_Send(vect_local, base + 1, MPI_INT, 0, 0, comm);
        } else {
            for (int k = 0; k < base + 1; k++) {
                vect[i * (base + 1) + k] = vect_local[k];
                vect_local[k] = 0;
            }
            for (j = i + 1; j < max && j < i + num_process; j++) {
                MPI_Recv(vect + (j) * (base + 1), base + 1, MPI_INT, j % num_process, 0, comm, MPI_STATUS_IGNORE);
            }
        }
    }

    if (rank == 0) {
        int *output = (int *)calloc(size, sizeof(int));
        int place = 1;
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
}

/**
 * @brief This function allows to find the maximum and the minimum in an array.
 * @param arr      pointer to the array that has to be sorted.
 * @param n        array size.
 * @param min      pointer to the minimum of the array.
 * @param max      pointer to the maximum of the array.
 */
void getMaxandMin(int *arr, int n, int *min, int *max) {
    *min = arr[0];
    *max = arr[0];
    for (int i = 1; i < n; i++) {
        if (arr[i] > *max)
            *max = arr[i];
        if (arr[i] < *min)
            *min = arr[i];
    }
}

/**
 * @brief This function implements the counting sort algorithm based on significant places, for the algorithm 0.
 * @param array         pointer to the array with the elements to sort, that are only positive or negative.
 * @param base          the base of the elements analized.
 * @param size          dimension of the array to order.
 * @param raw_index     the digit on which the counting has to be executed.
 * @param vect          the pointer to the vector in which are stored the frequencies of each positional element.
 */
void countingSortAlgo0(int array[], int base, int size, int raw_index, int *vect) {
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

    vect[base] = min;

    int length = max - min + 1;
    for (int j = 0; j < size; j++)
        vect[(((array[j] / place) % base) - min)]++;

    // Calculate cumulative count
    for (int j = 1; j < length; j++)
        vect[j] += vect[(j - 1)];
}

/**
 * @brief This function implements the counting sort for algorithm 1. Every process with rank different from 0 sends its frequencies vector to rank 0 with a MPI reduce,
 * then the rank 0 calculates the total vector of frequencies with another MPI reduce that sums the local frequencies.
 * @param rec_buf        pointer to sub-array of each process.
 * @param digit          number which represents cipher on which the counting sort is executed.
 * @param rank           rank of the current process.
 * @param dim            dimension of the sub-array rec_buf.
 * @param count          pointer to the output array with the total frequencies of the ciphers.
 */
void countingSortAlgo1(int *rec_buf, int digit, int rank, int dim, int min, int *count) {
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

void radix_sort(int **glob_array, int *tmp_array, int n, int num_process, int rank) {
    int dim;
    int other;
    int *array;
    if (rank == 0) {
        array = (int *)calloc(n, sizeof(int));
        if (array == NULL)
            perror("Memory Allocation - tmp_array");
        dim = n / num_process + n % num_process;
        other = n / num_process;
    } else {
        dim = n / num_process;
    }

    // ogi processo calcolerà il massimo tra i suoi elementi
    int local_max, local_min;
    getMaxandMin(tmp_array, dim, &local_min, &local_max);

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
        countingSortAlgo1(tmp_array, digit, rank, dim, global_min, frequencies[decimal_digit]);
        decimal_digit++;
    }
    /*
    if (rank != 0)
        MPI_Send(tmp_array, dim, MPI_INT, 0, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        memcpy(array, tmp_array, sizeof(int) * dim);
        for (int i = 1; i < num_process; i++)
            MPI_Recv(array + dim + (i - 1) * other, other, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }*/
    int *dime;
    dime = (int *)malloc(sizeof(int) * num_process);
    dime[0] = n / num_process + n % num_process;
    for (int i = 1; i < num_process; i++) {
        dime[i] = n / num_process;
    }
    int *displa;
    displa = (int *)malloc(sizeof(int) * num_process);
    displa[0] = 0;
    for (int i = 1; i < num_process; i++) {
        displa[i] = i * n / num_process + n % num_process;
    }
    MPI_Gatherv(tmp_array, dim, MPI_INT, array, dime, displa, MPI_INT, 0, MPI_COMM_WORLD);

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
        *glob_array = array;
    }
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

    if (old_num_process > 1) {
        MPI_Bcast(&max_neg, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&max_pos, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&size_neg, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&size_pos, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (old_rank == 1) {
            MPI_Recv(array_pos, size_pos, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (old_rank == 0) {
            MPI_Send(array_pos, size_pos, MPI_INT, 1, 0, MPI_COMM_WORLD);
        }
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
            MPI_Bcast(array_neg, size_neg, MPI_INT, 0, pari);
            serviceRadixsort(array_neg, size_neg, max_neg, rank, num_process, pari);
            if (old_rank == 0) {
                memcpy(array, array_neg, size_neg * sizeof(int));
                MPI_Recv(array + size_neg, size_pos, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        } else {
            // rank 1 master poisitivi sicuro

            MPI_Bcast(array_pos, size_pos, MPI_INT, 0, pari);
            serviceRadixsort(array_pos, size_pos, max_pos, rank, num_process, pari);
            if (old_rank == 1) {
                MPI_Send(array_pos, size_pos, MPI_INT, 0, 0, MPI_COMM_WORLD);
            }
        }
    }
}
