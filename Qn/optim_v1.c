#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <limits.h>

#define N 28
#define MASK ((((int64_t)1) << (N-1)) - 1)

#define TESTING false
#define EARLY_STOPPING 11 // For testing purposes, INT_MAX is equivalent to no early stopping

const int TREE_DEPTH = 18; //(N < 18) ? N : 18; // The ternary is just for making testing easier

// Custom data-types
typedef struct {
    int64_t* arr;
    int size;
    int capacity;
    size_t element_size;
} Array64Bit;

typedef struct {
    int A_left;
    int A_right;
    int L_left;
    int L_right;
} Interval;

typedef struct {
    int count;
    int index;
} CountIndexPair;

/*Prototypes*/
void appendToArray64Bit(Array64Bit* array, int64_t to_append);
void generate_all_rows(int** result_arr, int* temp_arr, int at, int temp_size, int* size);
void solve(size_t idx, unsigned int *cols, int64_t *diags, int64_t *anti_diags, int *rows, Array64Bit *result_arr, bool left);
int naive_combine(int64_t* A, int64_t* L, int sizeA, int sizeL);
int* generate_complement(int* arr);
void print_array(int* arr, int size);
void print_array64(int64_t* arr, int size);
double central_binomial_coefficient(int n);
int binary_search(int64_t* arr, int lower, int upper, int digit);
void construct_intervals(Interval** al_intervals, int* interval_capacity, int* interval_size, int64_t* A, int64_t* L, int A_left, int A_right, int L_left, int L_right, int lvl);
int compare_ints(const void *a, const void *b);
int64_t interval_combine(Interval* intervals, int num_intervals, int64_t* A, int64_t* L, int sizeA, int sizeL);
void shuffleArray(int** array, int n);
int64_t reorder(int64_t num, int* permutation, int len);
void digit_counter(int64_t* arr, int arr_len, int num_bits, int* digits);
int compareCountIndexPairs(const void* a, const void* b);
void createPermutation(int* digit_counts, int num_bits, int* permutation);
int countUnique(int64_t* arr, int size);
void write_array_to_file(const char* filename, int64_t* array, size_t size);


int main() {
    char filename[50]; 
    sprintf(filename, "./outputs/output_v1_%d_%d_%d.txt", N, TREE_DEPTH, EARLY_STOPPING);
    FILE *outputFile = fopen(filename, "w");

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

    shuffleArray(all_rows, number_of_rows); // Needed for testing
    printf("Shuffled\n");

    // interval combine
    float total_combine_time_interval = 0;
    int total_combines_interval = 0;
    int64_t total_count = 0;

    for (size_t row_idx = 0; row_idx < number_of_rows; row_idx++)
    {
        if (total_combines_interval >= EARLY_STOPPING && TESTING) break;

        int* rowsA = all_rows[row_idx];
        int* rowsL = generate_complement(rowsA);

        print_array(rowsA, N/2);
        print_array(rowsL, N/2);

        Array64Bit A = {NULL, 0, 0, sizeof(int64_t)};
        Array64Bit L = {NULL, 0, 0, sizeof(int64_t)};

        unsigned int A_cols = 0, L_cols = 0;
        int64_t A_diags = 0, A_anti_diags = 0, L_diags = 0, L_anti_diags = 0;

        solve(0, &A_cols, &A_diags, &A_anti_diags, rowsA, &A, false);
        solve(0, &L_cols, &L_diags, &L_anti_diags, rowsL, &L, true);
        
        printf("Solved\n");
        printf("A.size: %d, L.size: %d\n", A.size, L.size);

        // Reorder optimization (before sorting)
        int* digit_counts = calloc(2*N-2, sizeof(int));
        int* digit_counts_A = calloc(2*N-2, sizeof(int));
        int* digit_counts_L = calloc(2*N-2, sizeof(int));
        
        digit_counter(A.arr, A.size, 2*N-2, digit_counts_A);
        digit_counter(L.arr, L.size, 2*N-2, digit_counts_L);
        for (size_t i = 0; i < 2*N-2; i++) digit_counts[i] = digit_counts_A[i] + digit_counts_L[i];
        int permutation[2*N-2];
        createPermutation(digit_counts, 2*N-2, permutation);
        for (size_t i = 0; i < A.size; i++) A.arr[i] = reorder(A.arr[i], permutation, 2*N-2);
        for (size_t i = 0; i < L.size; i++) L.arr[i] = reorder(L.arr[i], permutation, 2*N-2);

        printf("Reordered\n");


        write_array_to_file("./data/A5.bin", A.arr, A.size);
        write_array_to_file("./data/L5.bin", L.arr, L.size);

        // Sort / construct implicit trees - note that in principle we only need to perform radix sort on the first TREE_DEPTH ints
        qsort(A.arr, A.size, sizeof(int64_t), compare_ints);
        qsort(L.arr, L.size, sizeof(int64_t), compare_ints);

        

        if (TESTING) {
            printf("sizeA: %d, uniqueA: %d, sizeL: %d, uniqueL: %d\n", A.size, countUnique(A.arr, A.size), L.size, countUnique(L.arr, L.size));
            print_array64(A.arr + A.size/4, 100);
            break;
        }

        Interval* al_intervals = NULL;
        int interval_capacity = 0;
        int interval_size = 0;
        construct_intervals(&al_intervals, &interval_capacity, &interval_size, A.arr, L.arr, 0, A.size, 0, L.size, 0);

        printf("Intervals constructed\n");

        clock_t start_time = clock();
        total_count += interval_combine(al_intervals, interval_size, A.arr, L.arr, A.size, L.size);
        clock_t end_time = clock();

        total_combine_time_interval += (float)(end_time - start_time)/CLOCKS_PER_SEC;
        total_combines_interval++;
        printf("Total count %llu\t---\tAverage combine time interval: %f\n", total_count, total_combine_time_interval/total_combines_interval);
        fprintf(outputFile, "Total count %llu\t---\tAverage combine time interval: %f\n", total_count, total_combine_time_interval/total_combines_interval);
        free(A.arr); free(L.arr); free(al_intervals); free(digit_counts); free(digit_counts_A); free(digit_counts_L);
    }
    
    printf("Total count %llu\t---\tAverage combine time interval: %f\n", total_count, total_combine_time_interval/total_combines_interval);
    fprintf(outputFile, "Total count %llu\t---\tAverage combine time interval: %f\n", total_count, total_combine_time_interval/total_combines_interval);


    for (int i = 0; i < number_of_rows; i++) {
        free(all_rows[i]); // Since all_rows.arr is an array of int pointers
    }
    free(all_rows);
    
    fclose(outputFile);
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

// Just a utility function
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

// This method is just horrible in terms of speed, but a good reference for other methods
int naive_combine(int64_t* A, int64_t* L, int sizeA, int sizeL) {
    int count = 0;
    for (size_t i = 0; i < sizeA; i++)
    {
        for (size_t j = 0; j < sizeL; j++)
        {
            if ((A[i] & L[j]) == 0) count++;
        }
    }
    
    return count;
}

// This method acheives many x the speedup (above >100x for N big enough and TREE_DEPTH big enough)
// But does take up a lot of memory for TREE_DEPTH big enough - TREE_DEPTH = 18 is about 5 GB, every additional lvl multiplies by a factor of ~3x
int64_t interval_combine(Interval* intervals, int num_intervals, int64_t* A, int64_t* L, int sizeA, int sizeL) {
    int64_t count = 0; // int64_t is needed when N is > 24
    
    for (size_t i = 0; i < num_intervals; i++)
    {
        int A_left = intervals[i].A_left, A_right = intervals[i].A_right, L_left = intervals[i].L_left, L_right = intervals[i].L_right;
        for (size_t j = A_left; j < A_right; j++)
        {
            for (size_t k = L_left; k < L_right; k++)
            {
                if ((A[j] & L[k]) == 0) count++;
            }   
        }
    }

    return count;
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

// Simple utility function for adding to the interval array
void addInterval(Interval** al_intervals, int* interval_capacity, int* interval_size, int A_left, int A_right, int L_left, int L_right) {
    if (*interval_size == *interval_capacity) {
        // Need to resize
        *interval_capacity = (*interval_capacity > 0) ? (*interval_capacity * 2) : 1;
        Interval* new_intervals = realloc(*al_intervals, *interval_capacity * sizeof(Interval));
        if (!new_intervals) {
            fprintf(stderr, "Failed to allocate memory for intervals\n");
            exit(EXIT_FAILURE);
        }
        *al_intervals = new_intervals;
    }

    // Add the new interval
    Interval new_interval = {A_left, A_right, L_left, L_right};
    (*al_intervals)[*interval_size] = new_interval;
    (*interval_size)++;
}


// Right's are exclusive, Left's are inclusive
void construct_intervals(Interval** al_intervals, int* interval_capacity, int* interval_size, int64_t* A, int64_t* L, int A_left, int A_right, int L_left, int L_right, int lvl) {
    if (A_right == A_left || L_right == L_left) return;
    if (lvl == TREE_DEPTH) {
        addInterval(al_intervals, interval_capacity, interval_size, A_left, A_right, L_left, L_right);
        return;
    }
    // Digit switch idxs - result of sorting
    int switch_idx_A = binary_search(A, A_left, A_right, 2*N-3-lvl);
    int switch_idx_L = binary_search(L, L_left, L_right, 2*N-3-lvl);

    // A0L0
    construct_intervals(al_intervals, interval_capacity, interval_size, A, L, A_left, switch_idx_A, L_left, switch_idx_L, lvl+1);
    // A1L0
    construct_intervals(al_intervals, interval_capacity, interval_size, A, L, switch_idx_A, A_right, L_left, switch_idx_L, lvl+1);
    // A0L1
    construct_intervals(al_intervals, interval_capacity, interval_size, A, L, A_left, switch_idx_A, switch_idx_L, L_right, lvl+1);
    return;
}

void print_array(int* arr, int size) {
    printf("[");
    for (int i = 0; i < size; i++) {
        printf("%d", arr[i]);
        if (i < size - 1) {
            printf(", ");
        }
    }
    printf("]\n");
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

// Reorders the bits of a 64-bit integer based on a given permutation
int64_t reorder(int64_t num, int* permutation, int len) {
    int64_t new_num = 0;
    for (int i = 0; i < len; i++) {
        int64_t digit = (1 & (num >> permutation[i]));
        new_num |= digit << i;
    }
    return new_num;
}

// Counts the number of 1s in each bit position of numbers in a 64-bit integer array
void digit_counter(int64_t* arr, int arr_len, int num_bits, int* digits) {
    for (int i = 0; i < num_bits; i++) {
        digits[i] = 0;
    }

    for (int i = 0; i < arr_len; i++) {
        for (int j = 0; j < num_bits; j++) {
            digits[j] += (arr[i] >> j) & 1;
        }
    }
}

// Comparison function for qsort (for ascending order)
int compareCountIndexPairs(const void* a, const void* b) {
    CountIndexPair* pair1 = (CountIndexPair*) a;
    CountIndexPair* pair2 = (CountIndexPair*) b;
    return pair1->count - pair2->count; // Ascending order
}

// Function to create a permutation array based on the counts
void createPermutation(int* digit_counts, int num_bits, int* permutation) {
    // Create an array of CountIndexPair
    CountIndexPair pairs[num_bits];
    for (int i = 0; i < num_bits; i++) {
        pairs[i].count = digit_counts[i];
        pairs[i].index = i;
    }

    // Sort the pairs in descending order of counts
    qsort(pairs, num_bits, sizeof(CountIndexPair), compareCountIndexPairs);

    // Extract the sorted indices to form the permutation array
    for (int i = 0; i < num_bits; i++) {
        permutation[i] = pairs[i].index;
    }
}

// This function is just for calculating stats - expects sorted input
int countUnique(int64_t* arr, int size) {
    if (size == 0) return 0;
    int count = 1;
    for (size_t i = 1; i < size; i++) {
        if (arr[i] != arr[i-1]) {
            count += 1;
        }
    }
    return count;
}

void write_array_to_file(const char* filename, int64_t* array, size_t size) {
    FILE* file = fopen(filename, "wb");
    if (file == NULL) {
        perror("Error opening file for writing");
        return;
    }

    // Write the size of the array first
    fwrite(&size, sizeof(size), 1, file);

    // Write the array data
    fwrite(array, sizeof(int64_t), size, file);

    fclose(file);
}