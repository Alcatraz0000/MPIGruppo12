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

#include "RadixSortSeq.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void init_structures(int **array, int length, char *FILE_A) {  // implementazione I/O doppio da file
    int *tmp_array = (int *)malloc(length * sizeof(int));
    if (tmp_array == NULL)
        perror("Memory Allocation - tmp_array");
    FILE *file = fopen(FILE_A, "r");
    if (fread(tmp_array, sizeof(int), length, file) != length)
        perror("error during lecture from file");
    fclose(file);
    *array = tmp_array;
}

void write_on_File(int length, int max_digit, char *FILE_A) {
    int *tmp_array = (int *)malloc(length * sizeof(int));
    if (tmp_array == NULL)
        perror("Memory Allocation - tmp_array");
    srand(time(NULL));
    for (int i = 0; i < length; i++) {
        if (i % 2)
            tmp_array[i] = -(rand() % max_digit);
        else
            tmp_array[i] = (rand() % max_digit);
    }
    FILE *file = fopen(FILE_A, "w");
    fwrite(tmp_array, sizeof(int), length, file);
    fclose(file);
}

/**
 * @brief The function that separates the array in subarray and starts the sorting process.
 * @param array          array.
 * @param n              array size.
 * @param num_process    number of processes.
 * @param rank           rank of the current process.
 */
void radix_sort(int *array, int n) {
    int max;
    int min;
    getMaxandMin(array, n, &min, &max);

    for (int digit = 1; (max - min) / digit > 0; digit *= 10) {
        countingSortAlgo1(array, max, min, n, digit);
    }
}
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

void countingSortAlgo1(int *vet, int max, int min, int n, int dig) {
    int Count[10] = {0};

    for (int i = 0; i < n; i++) {
        Count[((vet[i] - min) / dig) % 10]++;
    }

    for (int i = 1; i < 10; i++) {
        Count[i] += Count[i - 1];
    }

    int *Output = (int *)malloc(sizeof(int) * n);
    for (int i = n - 1; i >= 0; i--) {
        Output[Count[((vet[i] - min) / dig) % 10] - 1] = vet[i];
        Count[((vet[i] - min) / dig) % 10]--;
    }

    memcpy(vet, Output, sizeof(int) * n);
    free(Output);
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

void serviceRadixsort(int array[], int size, int max) {
    // Get maximum element
    int base = 10;
    int *vect[max];
    for (int i = 0; i < max; i++) {
        vect[i] = (int *)calloc((base + 1), sizeof(int));
    }

    int *output = (int *)calloc(size, sizeof(int));
    int i = 0, j = 0;

    for (i = 0; i < max; i++) {
        countingSortAlgo0(array, base, size, i, vect[i]);
    }

    int place = 1;

    for (i = 0; i < max; i++) {
        for (j = size - 1; j >= 0; j--) {
            output[vect[i][(((array[j] / place) % base) - vect[i][base])] - 1] = array[j];
            vect[i][(((array[j] / place) % base) - vect[i][base])]--;
        }

        for (int j = 0; j < size; j += 2) {
            array[j] = output[j];
            array[j + 1] = output[j + 1];
        }
        place = place * 10;
    }
}

void myRadixsort(int *array, int length) {
    int *array_pos = (int *)malloc(length * sizeof(int));
    if (array_pos == NULL)
        perror("Memory Allocation - a");
    int *array_neg = (int *)malloc(length * sizeof(int));
    if (array_neg == NULL)
        perror("Memory Allocation - a");
    int max_pos, max_neg;
    int size_pos, size_neg;
    getMaxDigitSeq(array, length, array_pos, array_neg, &max_pos, &max_neg, &size_pos, &size_neg);

    serviceRadixsort(array_neg, size_neg, max_neg);

    memcpy(array, array_neg, size_neg * sizeof(int));

    serviceRadixsort(array_pos, size_pos, max_pos);
    memcpy(array + size_neg, array_pos, size_pos * sizeof(int));
}
