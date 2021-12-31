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
        @file main.c

*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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

int main(int argc, char **argv) {
    int *array;
    int length, max_digit;

    if (argc < 3) {
        printf("ERROR! Usage: ./main rows columns threads");
        exit(1);
    }

    length = atoi(argv[1]);
    max_digit = atoi(argv[2]);

    write_on_File(length, max_digit, "VectGruppo12");
}