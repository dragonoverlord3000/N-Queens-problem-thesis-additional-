#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

// Macros
#define N 26
#define MASK ((1 << (N-1)) - 1)

// Custom data-types
typedef struct {
    int64_t* arr;
    int size;
    int capacity;
    size_t element_size;
} Array64Bit;

// Count is how many has that prefix, the rest is standard binary tree
typedef struct BTNode {
    int count;
    struct BTNode* left;
    struct BTNode* right;
} BTNode;

// Function prototypes
void solve(size_t idx, unsigned int *cols, int64_t *diags, int64_t *anti_diags, int *rows, Array64Bit *result_arr, bool left);
void appendToArray64Bit(Array64Bit* array, int64_t to_append);
void generate_all_rows(int** result_arr, int* temp_arr, int at, int temp_size, int* size);
int* generate_complement(int* arr);
void print_array64(int64_t* arr, int size);
void shuffleArray(int** array, int n);
int binary_search(int64_t* arr, int lower, int upper, int digit);
int compare_ints(const void *a, const void *b);
double central_binomial_coefficient(int n);
int naive_combine(int64_t* A, int64_t* L, int sizeA, int sizeL);

void construct_BT(BTNode* node, Array64Bit* array, int digit, int lower, int upper);
void freeBT(BTNode* node);
int64_t BT_combiner(BTNode* ABT, BTNode* LBT, int digit);

void printBT(const char* prefix, const BTNode* node, int isLeft);
void printBTNode(const BTNode* node);
BTNode* newNode(int count);


int main() {
    // Define the outputfile
    char filename[50]; 
    sprintf(filename, "./outputs/output_v1_%d.txt", N);
    FILE *outputFile = fopen(filename, "w");

    // Give basic info about the run
    printf("N = %d\n", N);
    fprintf(outputFile, "N = %d\n", N);

    double _number_of_rows = central_binomial_coefficient(N/2);
    printf("Number of rows: %f\n", _number_of_rows);
    fprintf(outputFile, "Number of rows: %f\n", _number_of_rows);
    int number_of_rows = round(_number_of_rows);
    printf("Number of rows: %d\n", number_of_rows);
    fprintf(outputFile, "Number of rows: %d\n", number_of_rows);
    int** all_rows = malloc(number_of_rows * sizeof(int*));
    int size = 0;
    int temp_arr[N/2];
    generate_all_rows(all_rows, temp_arr, 0, 0, &size);
    printf("All %d rows generated\n", number_of_rows);
    fprintf(outputFile, "All %d rows generated\n", number_of_rows);

    // TESTING
    shuffleArray(all_rows, number_of_rows); // Needed for testing
    printf("Shuffled\n");
    
    // Variable for keeping the result
    int total_combines = 0;
    float total_combine_time = 0;
    int64_t total_count = 0;

    for (int row_idx = 0; row_idx < number_of_rows; row_idx++) {
        int* rowsA = all_rows[row_idx];
        int* rowsL = generate_complement(rowsA);

        // Solve right board and left board
        Array64Bit A = {NULL, 0, 0, sizeof(int64_t)};
        Array64Bit L = {NULL, 0, 0, sizeof(int64_t)};

        unsigned int A_cols = 0, L_cols = 0;
        int64_t A_diags = 0, A_anti_diags = 0, L_diags = 0, L_anti_diags = 0;

        solve(0, &A_cols, &A_diags, &A_anti_diags, rowsA, &A, false);
        solve(0, &L_cols, &L_diags, &L_anti_diags, rowsL, &L, true);
        // printf("Solved\n");

        qsort(A.arr, A.size, sizeof(int64_t), compare_ints);
        qsort(L.arr, L.size, sizeof(int64_t), compare_ints);
        // printf("Sorted\n");


        printf("counting\n");
        clock_t start_time = clock();
        //
        clock_t end_time = clock();
        // total_count_naive += naive_combine(A.arr, L.arr, A.size, L.size);

        total_combine_time += (float)(end_time - start_time)/CLOCKS_PER_SEC;
        total_combines++;
        printf("Total count %lld\t---\tAverage combine time: %f\n", total_count, total_combine_time/(float)total_combines);
        fprintf(outputFile, "Total count %lld\t---\tAverage combine time: %f\n", total_count, total_combine_time/(float)total_combines);
        
        free(A.arr); free(L.arr);
        free(rowsA); free(rowsL);
    }
    printf("Total count: %lld\t---\tAverage combine time: %f\t---\tTotal combine time: %f\n", total_count, total_combine_time/(float)total_combines, total_combine_time);
    fprintf(outputFile, "Total count: %lld\t---\tAverage combine time: %f\t---\tTotal combine time: %f\n", total_count, total_combine_time/(float)total_combines, total_combine_time);

    return 0;
}

double central_binomial_coefficient(int n) {
    if (n == 0) return 1;

    double result = 1;
    for (int i = 1; i <= n; ++i) {
        result *= (double)(n + i) / i;
    }

    return result;
}

// The half-board solver
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

// Generate all the N choose N/2 rows - simple recursion
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

// If N = 10 and arr = [0,1,4,6,9], then it's complement would be [2,3,5,7,8]
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

// lower is inclusive, upper is exclusive - searches on most important idx
int binary_search(int64_t* arr, int lower, int upper, int digit) {
    // Special case - all of `digit` is 1
    if ((((int64_t)1) << digit) & arr[lower]) return lower;
    
    while (upper - lower > 1) {
        int new = (lower + upper)/2;
        if (((((int64_t)1) << digit) & arr[new]) == 0) lower = new;
        else upper = new;
    }
    return upper;
}

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

