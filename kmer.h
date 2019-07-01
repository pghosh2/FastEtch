#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "countmin.h"
#include "prng.h"

#define DEFAULT_SIGMA_SIZE              26
#define ADJUST_SEQ_SIZE                 110
#define BKT_CAPACITY_THRESH             20
#define DEF_NUM_CONTIGS                 100
#define DEF_CONTIG_SIZE                 100

void omp_set_lock(omp_lock_t *lock);
void omp_unset_lock(omp_lock_t *lock);

//omp_lock_t hash_lock;
//FILE * fp_outfile;

//data structure for storing each Read
typedef struct Rd_Sequence
{
	// current sequence size
	int rlen;
	// data
	char *read_data;
} Sequence_r;

/* Each key-value pair entry */
/*
typedef struct htable_each_entry
{
        char kmer_plus_one[WIN_SIZE_PLUS];
        int flag_acgt[4];
        int kmer_count;
	struct htable_each_entry *hh_next;
        struct htable_each_entry *succ_ptr;
        int linklist_id;
        int my_rank_in_ll;
        bool pred_flag;
        long int my_kmer_id;

} htable_entry;
*/

/* structure for each bucket entry */
/*
typedef struct htable_each_bucket
{
	htable_entry *series;

	// Lock for each bucket entry
	omp_lock_t writelock;

	int count;

} htable_bucket;
*/

/* structure for the Hash table */  
/*
typedef struct hash_table_t {
	size_t size;
	htable_bucket *buckets;

	// Lock for the entire table
	//omp_lock_t table_lock;

} hash_table_t; 
*/

typedef struct each_contig_entry
{
       long int my_contig_id;
       char *contig_data;

} contig_entry;

typedef struct contig_series
{
       int contig_count;
       contig_entry *c_series;

} contig_thrd_list;

//hash_table_t *create_hash_table (int size);
int hash_fcn (const char* word, unsigned M);
void parse_alphabet_file (FILE * fp);
void Initialize_input_read_seq (Sequence_r **sequence, int num_reads);
void free_reads (int num_reads);
void free_contigs (contig_thrd_list *global_clist, int num_threads, int *malloc_contig_count);
bool is_sigma (char c);
void parse_input_reads_file (FILE *fp, Sequence_r *rd_sequence, int num_reads, size_t size);
//unsigned generate_hash (hash_table_t *hashtab, char *key);
unsigned hash_to_bucket (char *key, unsigned num_buckets);
//htable_entry *Find_key_exists (hash_table_t *hashtab, char *key, unsigned bucket_id, htable_entry ** last_entry);
//htable_entry *Find_succ_exists (hash_table_t *hashtab, char *key, unsigned bucket_id);
//void Expand_hash_table (hash_table_t *hashtab);
//void Insert_key (hash_table_t *hashtab, htable_entry *tag, unsigned bucket_id, htable_entry *last_entry);
int count_hashtab_entries (CM_type_Heap *cm);
//void print_bucket_entries (hash_table_t *hashtab);
void print_hashtab_entries (CM_type *cm, CM_type_Heap *cm_h, double collison_rate);
//void total_all_kmer_count(hash_table_t *hashtab);
//void free_table (hash_table_t *hashtable);
void change_to_num (char *pred_kmer, int *num);
void change_to_char (char *pred_kmer, int num);
int find_max (int array[]);
int find_min (int array[], int *min);
int compare (const void * a, const void * b);
int FindIndex( const int a[], int size, int value);
//void calculate_heuristic_cutoff (hash_table_t *hashtable, int num_entries, int num_kmers);
void Sliding_window (int idx, CM_type *cm, CM_type_Heap *cm_h);
//void enumerate_contigs_to_file (hash_table_t *my_hash_table, int bucket_id, FILE* f);
void initialize_contig_array_per_thread (contig_entry **c_seq, size_t def_num_contigs, size_t def_contig_size);
void enumerate_contigs_to_array (CM_type *cm_c, CM_type_Heap *cm, contig_entry **c_seq, int *tot_contigs, int *alloc_contig_num, 
                                 size_t max_contig_length, double collison_rate);
//CM_type * CM_Init(int width, int depth, int seed);
CM_type * CM_Init(int width, int depth, int seed, bool is_hash_present, char *hashes_filename);
//CM_type_Heap * CM_Heap_Init(int width, int seed);
CM_type_Heap * CM_Heap_Init(int width, int depth, int seed, bool is_hash_present, char *hashes_filename);
int convert_hash_fcn (const char* word);
unsigned int hash_str(const char *str);
void CM_Update(CM_type * cm, unsigned int item, int diff);
void CM_Update_str(CM_type *cm, CM_type_Heap *cm_h, char *key, int val, int Heap_Threshold);
int CM_PointEst(CM_type * cm, unsigned int query);
int CM_Estimate_str(CM_type *cm, char *key);
void CM_Destroy(CM_type * cm);
void CM_Heap_Destroy(CM_type_Heap * cm);

long int CMS_Size (CM_type *cm);
long int CM_Size_Heap (CM_type_Heap *cm);
long int Reads_size (int num_reads, Sequence_r *rd_seq);
long int Contigs_size (int num_threads, contig_thrd_list *global_clist, int *malloc_contig_count);

unsigned createMask(unsigned a, unsigned b);
void convert_char_to_binary (char *kmer_name, unsigned char *result);
void convert_binary_to_char (unsigned char *result, char output_str[WIN_SIZE_PLUS]);

