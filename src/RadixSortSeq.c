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
        @file RadixSort.c
*/

#include "RadixSortSeq.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief This function initializes the array structure needed in the program.
 * The inputs are read from file sequentially.
 * @param array         pointer to the portion of array to sort in the sorting algorithm.
 * @param length        dimension of array, so the number of elements contained in the array.
 * @param FILE_A        the name of the file to read on.
 */
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
 * @brief The function that starts the sorting process.
 * First we search max and min in the array. We calculate the max number of digit (max_pos) and call countingSort on all of this.
 * @param array          array.
 * @param n              array size.
 */
void radix_sort(int *array, int n) {
    int max;
    int min;
    getMaxandMin(array, n, &min, &max);
    int max_pos = 0;
    int tmp_pos = max - min;
    while (tmp_pos > 0) {
        tmp_pos /= 10;
        max_pos++;
    }
    int *frequencies[max_pos];

    for (int i = 0; i < max_pos; i++) {
        frequencies[i] = (int *)calloc(10, sizeof(int));
    }
    int decimal_digit = 0;
    for (int digit = 1; (max - min) / digit > 0; digit *= 10) {
        countingSortAlgo1(array, min, n, frequencies[decimal_digit], digit);
        decimal_digit++;
    }
    int *temp_array = (int *)malloc(sizeof(int) * n);
    int val = 1;
    for (int j = 0; j < max_pos; j++) {
        for (int i = n - 1; i >= 0; i--) {
            temp_array[frequencies[j][((array[i] - min) / val) % 10] - 1] = array[i];
            frequencies[j][((array[i] - min) / val) % 10]--;
        }
        val *= 10;
        memcpy(array, temp_array, sizeof(int) * n);
    }
    free(temp_array);
}

/**
 * @brief This function implements the counting sort for algorithm 1, calculates the total vector of frequencies.
 * @param vet        pointer to array.
 * @param min           the minimum value find in the array.
 * @param n            dimension of the array.
 * @param Count          pointer to the output array with the total frequencies of the ciphers.
 * @param dig          the digit we are computing.
 */

void countingSortAlgo1(int *vet, int min, int n, int *Count, int dig) {
    for (int i = 0; i < n; i++) {
        Count[((vet[i] - min) / dig) % 10]++;
    }

    for (int i = 1; i < 10; i++) {
        Count[i] += Count[i - 1];
    }
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
 * @brief This function implements the counting sort algorithm based on significant places, for the algorithm 0.
 * @param array         pointer to the array with the elements to sort, that are only positive or negative.
 * @param base          the base of the elements analized.
 * @param size          dimension of the array to order.
 * @param raw_index     the digit on which the counting has to be executed.
 * @param vect          the pointer to the vector in which are stored the frequencies of each positional element.
 */

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
 * @brief This function implements the radix sort algorithm for algorithm 0 for a only positive or only negative array, based on counting sort algorithm,
 * storing during the computation the frequencies of each element in a vector.
 * @param array       pointer to the array with the elements to sort, that are only positive or negative.
 * @param size        dimension of the subarray containing all positive or negative elements.
 * @param max         maximum number of digits contained in the longer element of this subarray
 */

void serviceRadixsort(int array[], int size, int max) {
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

/**
 * @brief This function instantiates two different arrays for store separately positive and negative values, in order to separate the operation of sorting.
 * @param array      pointer to the array with total elements to sort.
 * @param length       dimension of the array containing all elements.
 */
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
