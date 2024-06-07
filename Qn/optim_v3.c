#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define N 26
#define MASK ((((int64_t)1) << (N-1)) - 1)
#define p 16

#define EARLY_STOPPING 100

// Custom data-types
typedef struct {
    int64_t* arr;
    int size;
    int capacity;
    size_t element_size;
} Array64Bit;

typedef struct {
    Array64Bit** arr;
    int size;
} PrefixArr;

void generate_all_rows(int** result_arr, int* temp_arr, int at, int temp_size, int* size);
int* generate_complement(int* arr);
void solve(size_t idx, unsigned int *cols, int64_t *diags, int64_t *anti_diags, int *rows, Array64Bit *result_arr, bool left);
double central_binomial_coefficient(int n);
int compare_ints(const void *a, const void *b);
int64_t prefixer_combine(int64_t* A, int64_t* L, int64_t* prefix2countA, int64_t* prefix2countL, int sizeA, int sizeL);
int64_t naive_combine(int64_t* A, int64_t* L, int leftA, int rightA, int leftL, int rightL);
void appendToArray64Bit(Array64Bit* array, int64_t to_append);
void print_array64(int64_t* arr, int size);
void shuffleArray(int** array, int n);

// Main program
int main() {
    char filename[50]; 
    sprintf(filename, "./outputs/output_v2_%d_%d_%d.txt", N, p, EARLY_STOPPING);
    FILE *outputFile = fopen(filename, "w");

    printf("N = %d\n", N);
    fprintf(outputFile, "N = %d\n", N);
    double _number_of_rows = central_binomial_coefficient(N/2);
    int number_of_rows = round(_number_of_rows);
    printf("Number of rows: %d\n", number_of_rows);
    fprintf(outputFile, "Number of rows: %d\n", number_of_rows);
    int** all_rows = malloc(number_of_rows * sizeof(int*));
    int size = 0;
    int temp_arr[N/2];
    generate_all_rows(all_rows, temp_arr, 0, 0, &size);
    printf("All %d rows generated\n", number_of_rows);
    fprintf(outputFile, "All %d rows generated\n", number_of_rows);

    shuffleArray(all_rows, number_of_rows); // Needed for testing
    printf("Shuffled\n");

    float total_combine_time_interval = 0;
    int total_combines_interval = 0;
    int64_t total_count = 0;
    // int64_t naive_count = 0;

    for (size_t row_idx = 0; row_idx < number_of_rows; row_idx++) {
        if (total_combines_interval >= EARLY_STOPPING) break;

        int* rowsA = all_rows[row_idx];
        int* rowsL = generate_complement(rowsA);

        Array64Bit A = {NULL, 0, 0, sizeof(int64_t)};
        Array64Bit L = {NULL, 0, 0, sizeof(int64_t)};

        unsigned int A_cols = 0, L_cols = 0;
        int64_t A_diags = 0, A_anti_diags = 0, L_diags = 0, L_anti_diags = 0;

        // Solve half-boards
        solve(0, &A_cols, &A_diags, &A_anti_diags, rowsA, &A, false);
        solve(0, &L_cols, &L_diags, &L_anti_diags, rowsL, &L, true);

        // Sort
        qsort(A.arr, A.size, sizeof(int64_t), compare_ints);
        qsort(L.arr, L.size, sizeof(int64_t), compare_ints);

        // prefix2countX[prefix] is when prefix ends (exclusive)
        int64_t* prefix2countA = calloc(1 << p,  sizeof(int64_t));
        int64_t* prefix2countL = calloc(1 << p,  sizeof(int64_t));

        // O(2**p)
        for (size_t A_idx = 0; A_idx < A.size; A_idx++) prefix2countA[A.arr[A_idx] >> (2*N - 2 - p)]++;
        for (size_t L_idx = 0; L_idx < L.size; L_idx++) prefix2countL[L.arr[L_idx] >> (2*N - 2 - p)]++;
        
        clock_t start_time = clock();
        total_count += prefixer_combine(A.arr, L.arr, prefix2countA, prefix2countL, A.size, L.size);
        clock_t end_time = clock();

        float duration = (float)(end_time - start_time)/CLOCKS_PER_SEC;
        printf("A.size: %d, L.size: %d, Time: %f\n", A.size, L.size, duration);
        fprintf(outputFile, "A.size: %d, L.size: %d, Time: %f\n", A.size, L.size, duration);

        total_combine_time_interval += duration;
        total_combines_interval++;
    
        printf("Total count: %llu\t---\tAvg time: %f\n", total_count, total_combine_time_interval/total_combines_interval);
        fprintf(outputFile, "Total count: %llu\t---\tAvg time: %f\n", total_count, total_combine_time_interval/total_combines_interval);

        free(A.arr); free(L.arr);
    }

    printf("\n\n\n---------------\n\n\n");
    fprintf(outputFile, "\n\n\n---------------\n\n\n");
    printf("Total count: %llu\t---\tAvg time: %f\n", total_count, total_combine_time_interval/total_combines_interval);
    fprintf(outputFile, "Total count: %llu\t---\tAvg time: %f\n", total_count, total_combine_time_interval/total_combines_interval);

    for (int i = 0; i < number_of_rows; i++) {
        free(all_rows[i]); // Since all_rows.arr is an array of int pointers
    }
    free(all_rows);
    
    fclose(outputFile);

    return 0;
}

// Utility functions
// For qsort
int compare_ints(const void *a, const void *b)  {
    // Casting the void pointers to int pointers and dereferencing them
    int64_t intA = *(const int64_t*)a;
    int64_t intB = *(const int64_t*)b;

    // Comparison logic
    if (intA < intB) {
        return -1;
    } else if (intA > intB) {
        return 1;
    } else {
        return 0;
    }
}

// For counting # of rows to generate
double central_binomial_coefficient(int n) {
    if (n == 0) return 1;

    double result = 1;
    for (int i = 1; i <= n; ++i) {
        result *= (double)(n + i) / i;
    }

    return result;
}

void generate_all_rows(int** result_arr, int* temp_arr, int at, int temp_size, int* size) {
    if (temp_size == N/2) {
        int* combination = malloc(sizeof(int) * (N/2));
        memcpy(combination, temp_arr, sizeof(int) * (N/2));
        result_arr[*size] = combination;
        (*size)++;
    } else {
        for (int i = at; i < N; ++i) {
            temp_arr[temp_size] = i;
            generate_all_rows(result_arr, temp_arr, i + 1, temp_size + 1, size);
        }
    }
}

int* generate_complement(int* arr) {
    int* complement = malloc(sizeof(int) * (N/2));

    int pointer = 0;
    for (size_t i = 0; i < N; i++)
    {
        bool in_arr = false;
        for (size_t j = 0; j < N/2; j++)
        {
            if (arr[j] == i) {
                in_arr = true;
                break;
            }
        }
        if (!in_arr) {
            complement[pointer++] = i;
        }
    }

    return complement;
}


void solve(size_t idx, unsigned int *cols, int64_t *diags, int64_t *anti_diags, int *rows, Array64Bit *result_arr, bool left) {
    if (idx >= (N/2)) {
        if (left)
        {
            int64_t result = ((*diags >> ((N/2))) & MASK) | (((*anti_diags >> ((N/2))) & MASK) << (N - 1));
            appendToArray64Bit(result_arr, result);
        } else {
            int64_t result = ((*diags) & MASK) | (((*anti_diags) & MASK) << (N - 1));
            appendToArray64Bit(result_arr, result);
        }
        return;
    }

    unsigned int row = rows[idx];
    for (size_t col = 0; col < (N/2); col++) {
        if (((((int64_t)1) << (N - 1 + col - row)) & *diags) | ((((int64_t)1) << (col + row)) & *anti_diags) | ((((int64_t)1) << col) & *cols)) continue;

        *diags |= (((int64_t)1) << (N - 1 + col - row));
        *anti_diags |= (((int64_t)1) << (col + row));
        *cols |= (1 << col);

        solve(idx + 1, cols, diags, anti_diags, rows, result_arr, left);

        *diags ^= (((int64_t)1) << (N - 1 + col - row));
        *anti_diags ^= (((int64_t)1) << (col + row));
        *cols ^= (1 << col);
    }
}

int64_t prefixer_combine(int64_t* A, int64_t* L, int64_t* prefix2countA, int64_t* prefix2countL, int sizeA, int sizeL) {
    int64_t combine_count = 0;
    int64_t pA_count = 0, pL_count = 0;

    for (size_t prefixA = 0; prefixA < (1 << p); prefixA++)
    {
        for (size_t prefixL = 0; prefixL < (1 << p); prefixL++)
        {
            if ((prefixA & prefixL) == 0 && prefix2countA[prefixA] && prefix2countL[prefixL]) {
                combine_count += naive_combine(A, L, pA_count, pA_count + prefix2countA[prefixA], pL_count, pL_count + prefix2countL[prefixL]);
            }

            pL_count += prefix2countL[prefixL];
        }
        pL_count = 0;
        pA_count += prefix2countA[prefixA];        
    }
    
    return combine_count;
}

int64_t naive_combine(int64_t* A, int64_t* L, int leftA, int rightA, int leftL, int rightL) {
    int64_t count = 0;
    // printf("leftA: %d, rightA: %d, leftL: %d, rightL: %d\n", leftA, rightA, leftL, rightL);

    for (size_t i = leftA; i < rightA; i++)
    {
        for (size_t j = leftL; j < rightL; j++)
        {
            if ((A[i] & L[j]) == 0)
            {
                count++;
            }   
        }
    }
    
    return count;
}

void appendToArray64Bit(Array64Bit* array, int64_t to_append) {
    if (array->size == array->capacity) {
        // Resize array if necessary
        int new_capacity = array->capacity > 0 ? array->capacity * 2 : 1;
        int64_t* new_arr = realloc(array->arr, new_capacity * sizeof(int64_t));
        if (!new_arr) {
            fprintf(stderr, "Failed to allocate memory\n");
            exit(EXIT_FAILURE);
        }
        array->arr = new_arr;
        array->capacity = new_capacity;
    }

    // Append the value and increase the size
    array->arr[array->size] = to_append;
    array->size++;
}

void print_array64(int64_t* arr, int size) {
    printf("[");
    for (int i = 0; i < size; i++) {
        printf("%llu", arr[i]);
        if (i < size - 1) {
            printf(", ");
        }
    }
    printf("]\n");
}

// Purely for testing purposes
void shuffleArray(int** array, int n) {
    // Seed the random number generator
    srand(time(NULL));

    for (int i = n - 1; i > 0; i--) {
        // Generate a random index between 0 and i (note, the squared is to due to RAND_MAX being to small to shuffle properly for N large)
        int j = ((int64_t)rand() * (int64_t)rand()) % (i + 1);

        // Swap array[i] with array[j]
        int* temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
}