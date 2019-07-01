#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "prng.h"
#include "massdal.h"
//#include "countmin.h"
#include "kmer.h"
#include <omp.h>

long int CMS_Size (CM_type *cm)
{ // return the size of the sketch in bytes
  long int counts, hashes, admin;

  if (!cm) return 0;
  admin=sizeof(CM_type);
  counts=cm->width*cm->depth*sizeof(int);
  hashes=cm->depth*2*sizeof(unsigned int);
  return(admin + hashes + counts);
}

long int CM_Size_Heap (CM_type_Heap *cm)
{
    long int admin=0, tab_entries=0, kmer_cells=0, kmer_entries=0;
    int i=0;

    if (!cm) return 0;
    admin=sizeof(CM_type_Heap);
    tab_entries=sizeof(CM_cell_entry)*cm->width;

    for(i=0; i<cm->width; i++) {
        kmer_entries += cm->tab_entries[i].num_coll_kmers; 
        kmer_cells += ((cm->tab_entries[i].num_coll_kmers) * sizeof(cm_cell_list));
    }
 
    printf("Total number of kmer entries in Heap: %ld, size of each kmer_entry: %ld \n", kmer_entries, sizeof(cm_cell_list));
    printf("memory used by admin: %f, tab_entries: %f, kmer_cells: %f \n",
           ((double)admin/(double)1000000), ((double)tab_entries/(double)1000000), ((double)kmer_cells/(double)1000000));

    return (admin + tab_entries + kmer_cells);
}

long int Reads_size (int num_reads, Sequence_r *rd_seq)
{
    int i=0;
    long int all_reads_size=0, admin=0;

    admin=sizeof(Sequence_r) * num_reads;
    for (i=0; i<num_reads; i++)
	 all_reads_size += strlen(rd_seq[i].read_data);

    return (admin + all_reads_size);
}

long int Contigs_size (int num_threads, contig_thrd_list *global_clist, int *malloc_contig_count)
{
    int i=0,j=0;
    long int all_contigs_size=0, tot_contig_length=0, admin=0;

    admin=sizeof(contig_thrd_list)*num_threads;
    for(i=0; i<num_threads; i++) {
        all_contigs_size += (sizeof(contig_entry) * malloc_contig_count[i]);
        for(j=0; j<malloc_contig_count[i]; j++)
           tot_contig_length += strlen(global_clist[i].c_series[j].contig_data);
    }
   
    return (admin + all_contigs_size + tot_contig_length);
}

