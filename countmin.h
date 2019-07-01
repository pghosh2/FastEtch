// Two different structures: 
//   1 -- The basic CM Sketch
//   2 -- The hierarchical CM Sketch: with log n levels, for range sums etc. 

#define min(x,y)	((x) < (y) ? (x) : (y))
#define max(x,y)	((x) > (y) ? (x) : (y))
#include <iostream> 
#include <algorithm>
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define K_SIZE          (WIN_SIZE+1)/4

extern unsigned int *hashes;

typedef struct CM_type{
  long long count;
  int depth;
  int width;
  int **counts; // total aggregate count of all kmer's hashed to this cell
  unsigned int *hasha, *hashb;

} CM_type;

typedef struct cm_cell_llist{
  unsigned char kmer_plus_one[K_SIZE];
  bool kmer_pred_flag;
  int kmer_partial_count;
  //struct cm_cell_llist *kl_next;

} cm_cell_list;

typedef struct cm_cell_entry{
  //cm_cell_list *kmer_entries; // List of kmer's that hash to this cell
  std::vector <cm_cell_list> kmer_entries;
  int num_coll_kmers; // length of cm_cell_list
  omp_lock_t writelock; // Lock for each CM cell entry 

} CM_cell_entry;

typedef struct CM_type_Heap{
  int width;
  CM_cell_entry *tab_entries;
  unsigned int hasha, hashb;

} CM_type_Heap;

typedef struct CMF_type{ // shadow of above stucture with floats
  double count;
  int depth;
  int width;
  double ** counts;
  unsigned int *hasha, *hashb;

} CMF_type;

extern CM_type * CM_Init(int, int, int, bool, char*);
extern CM_type_Heap * CM_Heap_Init(int, int, int, bool, char*);
extern CM_type * CM_Copy(CM_type *);
extern void CM_Destroy(CM_type *);
extern void CM_Heap_Destroy(CM_type_Heap *);
extern int CM_Size(CM_type *);
cm_cell_list *find_kmer_exists(CM_type_Heap *cm, unsigned char *key, int loc);
void CM_Update_insert (CM_type_Heap *cm, char *key, int loc, int val);

extern void CM_Update(CM_type *, unsigned int, int); 
extern int CM_PointEst(CM_type *, unsigned int);
extern int CM_PointMed(CM_type *, unsigned int);
extern int CM_InnerProd(CM_type *, CM_type *);
extern int CM_Residue(CM_type *, unsigned int *);

extern CMF_type * CMF_Init(int, int, int);
extern CMF_type * CMF_Copy(CMF_type *);
extern void CMF_Destroy(CMF_type *);
extern int CMF_Size(CMF_type *);
extern void CMF_Update(CMF_type *, unsigned int, double); 
extern double CMF_InnerProd(CMF_type *, CMF_type *);
extern double CMF_PointProd(CMF_type *, CMF_type *, unsigned int);

typedef struct CMH_type{
  long long count;
  int U; // size of the universe in bits
  int gran; // granularity: eg 1, 4 or 8 bits
  int levels; // function of U and gran
  int freelim; // up to which level to keep exact counts
  int depth;
  int width;
  int ** counts;
  unsigned int **hasha, * *hashb;
} CMH_type;

extern CMH_type * CMH_Init(int, int, int, int);
extern CMH_type * CMH_Copy(CMH_type *);
extern void CMH_Destroy(CMH_type *);
extern int CMH_Size(CMH_type *);

extern void CMH_Update(CMH_type *, unsigned int, int);
extern int * CMH_FindHH(CMH_type *, int);
extern int CMH_Rangesum(CMH_type *, int, int);

extern int CMH_FindRange(CMH_type * cmh, int);
extern int CMH_Quantile(CMH_type *cmh,float);
extern long long CMH_F2Est(CMH_type *);
