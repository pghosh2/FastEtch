#include "kmer.h"
#include <assert.h>
#include <omp.h>

int num_expand_ops = 0;
long int total_delta = 0;
long int this_kmer_id = 0;
long int this_contig_id = 0;
long long int add_all_kmer_ids = 0;
int sigma_size = 0;
int COVERAGE = 0;
int delta = 0;
int Heap_Threshold = 0;

Sequence_r *rd_seq = NULL;
char alphabet_array[27];
unsigned int *hashes;

int hash_fcn (const char* word, unsigned M)
{
    unsigned int hash = 0;
    for (int i = 0 ; word[i] != '\0' ; i++)
    {
        hash = 31*hash + word[i];
    }
    return hash % M;
}

void parse_alphabet_file (FILE * fp)
{
	char ch;
	int num_sigma = 0;

	alphabet_array[num_sigma] = '$';
	num_sigma++;

	while ((ch = fgetc (fp)) != EOF)
	{
		if (ch != ' ')
		{
			alphabet_array[num_sigma] = ch;
			num_sigma++;
		}
	}
	alphabet_array[num_sigma] = '\0';
}

void Initialize_input_read_seq (Sequence_r **sequence, int num_reads)
{
	int j=0;

	// allocate sequence
	*sequence = (Sequence_r *) malloc (sizeof (Sequence_r) * num_reads);

	for(j=0; j<num_reads; j++)
	{
		(*sequence)[j].read_data = (char *) malloc (sizeof (char) * ADJUST_SEQ_SIZE);

	}
}

void free_reads (int num_reads)
{
	int i=0;

        /* print the amount of memory consumed by the reads */
        printf("\nAmount of memory consumed by all the Reads: %f MB\n", ((double)Reads_size(num_reads, rd_seq)/(double)1000000));
          
	for (i=0; i<num_reads; i++)
		free(rd_seq[i].read_data);

	free(rd_seq);

}

void free_contigs (contig_thrd_list *global_clist, int num_threads, int *malloc_contig_count)
{
	int i=0, j=0;

        for (i=0; i<num_threads; i++) {
             //printf("Freeing for thread: %d, with count: %d \n", i, global_clist[i].contig_count);
            for (j=0; j<malloc_contig_count[i]; j++) {
                 free(global_clist[i].c_series[j].contig_data);
            }
            free(global_clist[i].c_series);
        } 

	free(global_clist);

}

// check whether valid character
bool is_sigma (char c)
{
	int i = 0;
	int alpha_num = strlen (alphabet_array);

	for (i = 0; i < alpha_num; i++)
		if (alphabet_array[i] == c)
			return true;
	return false;
}

void parse_input_reads_file (FILE *fp, Sequence_r *rd_sequence, int num_reads, size_t size)
{
	int ch;
	int len = 0;
	int read_num=0;

	while ((ch = fgetc (fp)) != EOF)
	{
		len=0, size=ADJUST_SEQ_SIZE; 
		// initialize_read() ?
		while ((ch = fgetc (fp)) != '\n')
		{
			/*if (ch != '>' && flag == 0)
			{
				rd_sequence[read_num].rname[num_chars] = ch;
				num_chars++;
			}
			if (ch == ' ')
				flag = 1;
                        */
		}
		//rd_sequence[read_num].rname[num_chars] = '\0';

		while ((ch = fgetc (fp)) != '\n')
		{
			//if (read_num > 0) {
			//    printf("%c",ch); 
			//} 
			// check if c in alphabet array
			/*if (!is_sigma (ch))
			  {
			  printf ("Character: %c does not match alphabets in sigma \n", ch);
			  exit (0);
			  }
			 */

			if (ch != '\n')
			{
                                if ((ch >= 97) && (ch <= 122))
                                     ch = ch - 32;
				rd_sequence[read_num].read_data[len] = ch;
				len++;
			}

			if (len == size)
			{
				char *ptr = (char *) realloc ((char *) rd_sequence[read_num].read_data,
						sizeof (char) * (size += ADJUST_SEQ_SIZE));

				if (!ptr)
				{
					printf ("Error (re)allocating memory\n");
					exit (-1);
				}
				else
					rd_sequence[read_num].read_data = ptr;


			}
		}
		rd_sequence[read_num].read_data[len] = '\0';
		rd_sequence[read_num].rlen = len;
		//rd_sequence[read_num].read_id = read_num; //+ 1;
		read_num++;


	} // end eof while loop

	if (read_num != num_reads)
		printf("Error in parsing the Reads file \n");
	else
		printf ("All %d reads, stored in memory \n", num_reads);

}

// Hash function from UTHASH library
/*unsigned generate_hash (hash_table_t *hashtab, char *key)
{
	unsigned found_entry;
	unsigned num_buckets = hashtab->size;

	HASH_ID_BUCKET_STR(key, num_buckets, found_entry);

	return found_entry;
}

unsigned hash_to_bucket (char *key, unsigned num_buckets)
{
	unsigned found_entry;

	HASH_ID_BUCKET_STR(key, num_buckets, found_entry);

	return found_entry;
}
*/
/*
unsigned generate_hash (hash_table_t *hashtab, char *key)
{
	unsigned num_buckets = hashtab->size;

        return hash_fcn(key, num_buckets);
}
*/

unsigned hash_to_bucket (char *key, unsigned num_buckets)
{
        return hash_fcn(key, num_buckets);
}


/*htable_entry *Find_succ_exists (hash_table_t *hashtab, char *key, unsigned bucket_id)
{
	htable_entry *list;
 
	for(list = hashtab->buckets[bucket_id].series; list != NULL; list = list->hh_next) {
		if (strcmp(list->kmer_plus_one,key) == 0) 
			return list;
	}

	return NULL;
}
*/

int count_hashtab_entries (CM_type_Heap *cm)
{
	int i=0, tot_count=0;

	for (i=0; i<cm->width; i++)
		tot_count += cm->tab_entries[i].num_coll_kmers;

#ifdef SD_COLL_CHECK
        int *num_kmer_per_buck = (int *) calloc (tot_count, sizeof(int));
        double var=0.00;

        for (i=0; i<cm->width; i++) {
             num_kmer_per_buck[i] = cm->tab_entries[i].num_coll_kmers;
        }
        double coll_rate = (double)tot_count/cm->width;

        for (i=0; i<cm->width; i++)
             var += pow(((double)num_kmer_per_buck[i] - coll_rate), 2); 
        
        printf("Avg collison rate: %f , SD collison: %f \n", coll_rate, sqrt(var/(double)cm->width));

        free(num_kmer_per_buck);
#endif
	return tot_count;
}


/*pred_kmer_series *locate_begin_kmers(CM_type *cm, int *count)
{
   int i=0, j=0, k=0, d=0, loc=0, pred_cntr=0, kmer_count=0;
   cm_cell_list *temp=NULL, *list=NULL, *next=NULL;
   char pred_bp[2];
   char pred_kmer[WIN_SIZE_PLUS];
   int  pred_present[4];
   int def_num_pkmers = DEF_NUM_PRED_KMERS;
   kmer_count = *count;
   pred_kmer_series *pkmers = (pred_kmer_series*) malloc (sizeof(pred_kmer_series) * DEF_NUM_PRED_KMERS);

   for (i=0; i<DEF_NUM_PRED_KMERS; i++)
        memset(pkmers[i].kmer_plus_one, '\0', WIN_SIZE_PLUS);

#pragma omp for schedule(static) private(i) //(, pkmers, kmer_count)
   for (i=0; i<cm->width; i++)
   {
       if (cm->counts[d][i].num_kmers > 0 )
       {
           temp = cm->counts[d][i].kmer_entries;
           while (temp != NULL)
           {
              if (check_kmer_quality(cm, temp) == 0)
              {
#ifdef PRED_DEBUG  
                 printf(" ***** For kmer: %s ********* \n", temp->kmer_plus_one);
#endif

                 pred_cntr = 0;
                 for (j=0; j<4; j++)
                      pred_present[j] = 0;

                 for (j=0; j<4; j++) {
                      loc = 0;
                      list = NULL;
                      change_to_char(pred_bp,j);
                      memcpy(&pred_kmer[0], &pred_bp[0], 1);
                      memcpy(&pred_kmer[1], &temp->kmer_plus_one[0], WIN_SIZE);
                      pred_kmer[WIN_SIZE+1] = '\0'; 
     
                      loc = hash31(cm->hasha[d],cm->hashb[d],hash_str(pred_kmer)) % cm->width;
                      list = find_kmer_exists(cm, pred_kmer, d, loc, &next);
                      if (list == NULL) { 
                         pred_present[j] = 1;
#ifdef PRED_DEBUG  
                         printf(" pred kmer: %s missing \n", pred_kmer);
#endif
                      } else {
#ifdef PRED_DEBUG  
                         printf(" contains pred kmer: %s \n", list->kmer_plus_one);
#endif
                         if (check_kmer_quality(cm, list) == 1)
                             pred_present[j] = 1;  
                      }
                 }
                 for (j=0; j<4; j++)
                      if (pred_present[j] == 1)
                          pred_cntr++;
                 
#ifdef PRED_DEBUG  
                 printf("So for kmer: %s, pred_count: %d \n\n", temp->kmer_plus_one, pred_cntr); 
#endif
                 if (pred_cntr == 4) {
                     strcpy(pkmers[kmer_count].kmer_plus_one, temp->kmer_plus_one);
                     kmer_count++;
              
                     if (kmer_count == def_num_pkmers)
                     {
                         pred_kmer_series *ptr_2 = (pred_kmer_series *) realloc ((pkmers), sizeof(pred_kmer_series) * 
                                                                (def_num_pkmers += DEF_NUM_PRED_KMERS));

                         if (!ptr_2)
                         {
                             printf ("Error (re)allocating memory for Pred kmers\n");
                             exit (-1);
                         }
                         else
                             (pkmers) = ptr_2;
                      
                         for (k=kmer_count; k<def_num_pkmers; k++)
                              memset(pkmers[k].kmer_plus_one, '\0', WIN_SIZE_PLUS);
                             
                     }

                 } 
              } // end of 'if' condition checking is current kmer is viable
              temp = temp->kl_next;
           } // end of while loop
       }
   } // end of for loop

   *count = kmer_count;
   return pkmers;
}
*/

void print_hashtab_entries (CM_type *cm, CM_type_Heap *cm_h, double collison_rate)
{
	int i=0,CM_count=0; //quality=0;
	cm_cell_list list;
        char output_str[WIN_SIZE_PLUS];
        FILE *fp;

        fp = fopen("kmer_count.dat", "w");
        if (fp == NULL)
        {
            fprintf(stdout,"Error opening the htab_count file!\n");
            exit(1);
        }

	for (i=0; i<cm_h->width; i++){
	//	if (cm_h->tab_entries[i].num_coll_kmers > 0 )
	//	{
             
              for(int j=0; j<(int)cm_h->tab_entries[i].kmer_entries.size(); j++)
              {
			list = cm_h->tab_entries[i].kmer_entries[j];
                        //printf("num kmers: %d \n", cm_h->tab_entries[i].num_coll_kmers);

                                memset(output_str, '\0', WIN_SIZE_PLUS);
                                convert_binary_to_char(list.kmer_plus_one, output_str);
                                CM_count = CM_Estimate_str(cm, output_str);
                                /*quality = check_kmer_quality(cm, list, collison_rate);
				printf("Kmer+1 key: %s CM_count: %d Cov_count: %f quality: %d \n", 
                                        list->kmer_plus_one, CM_count, (double)CM_count/collison_rate, quality);
                                */
				printf("Kmer+1 key: %s CM_count: %d pred_flag: %d partial_count: %d \n", 
                                        output_str, CM_count, list.kmer_pred_flag, list.kmer_partial_count);

		}
	}
  
        fclose(fp);

} 

void change_to_num (char *pred_kmer, int *num)
{
     if (strcmp(pred_kmer,"A") == 0)
         *num = 0;
     else if (strcmp(pred_kmer,"C") == 0)
         *num = 1;
     else if (strcmp(pred_kmer,"G") == 0)
         *num = 2;
     else if (strcmp(pred_kmer,"T") == 0)
         *num = 3;
     else {        
            printf (" Error: Predecessor Kmer: %s unrecognizable \n", pred_kmer);
            exit(0);
     }
}

void change_to_char (char *pred_kmer, int num)
{

     switch (num)
     {
           case 0: 
               strcpy(pred_kmer,"A"); 
               pred_kmer[1] = '\0';
               break;

           case 1: 
               strcpy(pred_kmer,"C"); 
               pred_kmer[1] = '\0';
               break;

           case 2: 
               strcpy(pred_kmer,"G"); 
               pred_kmer[1] = '\0';
               break;

           case 3: 
               strcpy(pred_kmer,"T"); 
               pred_kmer[1] = '\0';
               break;
     }
}

int find_max (int array[])
{

  int i=0;
  int max = array[0];   // start with max = first element
  int ind = 0;

  for(i = 1; i<4; i++)
  {
      if(array[i] > max)
      {
        max = array[i];
        ind = i;
      }
  }
  return ind;                    // return location of highest value in array
}

int find_min (int array[], int *min)
{

  int i=0;
  int ind = 0;

  for (i = 0; i<4; i++) {
       if (array[i] >= 0) {
           *min = array[i];
            ind = i;
            break;
       }
  }

  for(i = 0; i<4; i++)
  {
     if (array[i] >= 0) { 
      if(array[i] < *min)
      {
        *min = array[i];
        ind = i;
      }
     }
  }

  return ind;                    // return location of highest value in array
}


int compare (const void * a, const void * b)
{
  //return ( *(int*)a - *(int*)b );
  return ( *(int*)b - *(int*)a );
}

int FindIndex( const int a[], int size, int value )
{
    int index = 0;

    while ( index < size && a[index] != value ) ++index;

    return ( index == size ? -1 : index );
}

/*void calculate_heuristic_cutoff (hash_table_t *hashtable, int num_entries, int num_kmers)
{
    int i=0, j=0, k=0;
    int temp[4];
    int sorted_temp[4];
    char succ_bp[2];
    char succ_kmer[WIN_SIZE_PLUS];
    htable_entry *list, *found_succ;
    unsigned succ_bucket = -1;
    htable_entry *existing_succs[4];

#ifdef _ENABLE_HEURISTIC
    int this_succ=-1;
    char output_file_name[100];
    char num_to_string[20];
    char kwin_to_string[3];
    int read_delta_num = 0;
    double var_sum=0.0, std_dev=0.0;
    int closest_kmer_count = 0;

    sprintf(num_to_string, "%d", num_entries);
    sprintf(kwin_to_string, "%d", WIN_SIZE);
    strcpy(output_file_name,"output_n");
    strcpy(&output_file_name[8],num_to_string);
    strcpy(&output_file_name[strlen(output_file_name)],"_k");
    strcpy(&output_file_name[strlen(output_file_name)],kwin_to_string);
    strcpy(&output_file_name[strlen(output_file_name)],".csv");
    
    FILE *f = fopen(output_file_name, "w+");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
#endif

//#pragma omp parallel for firstprivate(this_succ,succ_bucket) \
                         private(temp,existing_succs,succ_bp,succ_kmer,found_succ,list,i,j,k,closest_kmer_count) \
                         reduction(+:total_delta) 
    for(i=0; i<hashtable->size; i++) {
             
	list = hashtable->buckets[i].series;
        while (list != NULL)
        {
              if (check_kmer_quality(list) == 0)
              {     

                 for (k=0; k<4; k++) {
                     temp[k] = -1;
#ifdef _ENABLE_HEURISTIC
                     this_succ=-1;
#endif
                     existing_succs[k] = NULL;
                 }

#ifdef SUCC_DEBUG
                printf("kmer: %s Count: %d \n", list->kmer_plus_one, list->kmer_count);
#endif
                for (j=0; j<4; j++) {
                    switch (list->flag_acgt[j]) 
                    {
                        succ_bucket = 0;
                        found_succ = NULL;
                        case 0:
                             break;
                        case 1:
                             change_to_char(succ_bp,j);
			     memcpy(&succ_kmer[0],&list->kmer_plus_one[1],WIN_SIZE);
			     memcpy(&succ_kmer[WIN_SIZE],&succ_bp[0],1);
			     succ_kmer[WIN_SIZE+1] = '\0';
			     succ_bucket = generate_hash(hashtable, succ_kmer);
			     found_succ = Find_succ_exists (hashtable, succ_kmer, succ_bucket);

                             if (check_kmer_quality(found_succ) == 0)
                             {
                                 existing_succs[j] = found_succ;
			         temp[j] = abs(list->kmer_count - found_succ->kmer_count);
                             }
#ifdef SUCC_DEBUG
                             printf("In case 1,Succ_kmer: %s Count: %d, j: %d temp[j]: %d \n", found_succ->kmer_plus_one, found_succ->kmer_count, j, temp[j]);
#endif
                             break;
                        case 2:
                             break;
                        case 3:
                             change_to_char(succ_bp,j);
			     memcpy(&succ_kmer[0],&list->kmer_plus_one[1],WIN_SIZE);
			     memcpy(&succ_kmer[WIN_SIZE],&succ_bp[0],1);
			     succ_kmer[WIN_SIZE+1] = '\0';
			     succ_bucket = generate_hash(hashtable, succ_kmer);
			     found_succ = Find_succ_exists (hashtable, succ_kmer, succ_bucket);

                             if (check_kmer_quality(found_succ) == 0)
                             {
                                 existing_succs[j] = found_succ;
			         temp[j] = abs(list->kmer_count - found_succ->kmer_count);
                             }
#ifdef SUCC_DEBUG
                             printf("In case 3,Succ_kmer: %s Count: %d, j: %d temp[j]: %d \n", found_succ->kmer_plus_one, found_succ->kmer_count, j, temp[j]);
#endif
                             break;
                             
                    }
                }

                memcpy(sorted_temp, temp, sizeof(int)*4);
                qsort (sorted_temp, 4, sizeof(int), compare);

                int succ_index = -1;
 
                for (j=0; j<4; j++)
                {
                     if (sorted_temp[j] >= 0)
                     {
                         succ_index = FindIndex(temp, 4, sorted_temp[j]); 
                         if (existing_succs[succ_index]->pred_flag != true)
                            {
                              list->succ_ptr = existing_succs[succ_index];
                              existing_succs[succ_index]->pred_flag = true;
#ifdef SUCC_DEBUG
                              printf("Chosen-- for kmer: %s succ: %s \n", list->kmer_plus_one, list->succ_ptr->kmer_plus_one);
#endif
                              break;
                            }
                      }
                }

#ifdef _ENABLE_HEURISTIC
                //printf("closest kmer succ: %d , kmer: %s, delta: %d \n", this_succ, existing_succs[this_succ]->kmer_plus_one, closest_kmer_count);
                //fprintf(f, "kmer: %s  | closest succ: %s | delta: %d \n", list->kmer_plus_one, existing_succs[this_succ]->kmer_plus_one, closest_kmer_count);
                //fprintf(f, "%s , %d \n", list->kmer_plus_one, closest_kmer_count);
 
                closest_kmer_count = 0;
                this_succ = 0;                
                this_succ = find_min(temp, &closest_kmer_count);

                fprintf(f, "%d\n", closest_kmer_count);
                total_delta += closest_kmer_count;
#endif
             } // end of 'if' condition

             list = list->hh_next;

       } // end of while loop
 
     } // end of for loop 


#ifdef _ENABLE_HEURISTIC
     printf("TOT delta: %ld MEAN delta: %f \n", total_delta, ((double)total_delta/(double)num_kmers));
     double mean_of_delta = (double)total_delta/(double)num_kmers;

     fseek (f, 0, SEEK_SET);
     total_delta = 0;

     while ( fscanf(f,"%d",&read_delta_num) == 1)
     {
         total_delta += read_delta_num;
         var_sum += pow(((double)read_delta_num - mean_of_delta), 2);
     }

     std_dev = sqrt(var_sum/(double)num_kmers);
     printf("Standard Deviation of Delta: %f total delta: %ld \n", std_dev, total_delta);

     fclose(f); 
     int delete_status = remove(output_file_name);
     if (delete_status != 0) printf("File not deleted successfully \n");
#endif

}
*/


void Sliding_window (int idx, CM_type *cm, CM_type_Heap *cm_h)
{
	char kmer_plus_1[WIN_SIZE_PLUS];
	int i=0;
	Sequence_r r_seq;

	r_seq = rd_seq[idx];

	for(i=0; (r_seq.rlen - i) > WIN_SIZE; i++)
	{

		memcpy(kmer_plus_1,&r_seq.read_data[i],WIN_SIZE+1);
		kmer_plus_1[WIN_SIZE+1] = '\0';


                CM_Update_str(cm, cm_h, kmer_plus_1, 1, Heap_Threshold); 
                        //tag->my_kmer_id = __sync_add_and_fetch(&this_kmer_id, 1);


	} 
}

void initialize_contig_array_per_thread (contig_entry **c_seq, size_t def_num_contigs, size_t def_contig_size)
{
     int i=0;

     *c_seq = (contig_entry *) malloc (sizeof(contig_entry) * def_num_contigs);

     for (i=0; i<def_num_contigs; i++) {
          
          (*c_seq)[i].my_contig_id = 0;
          (*c_seq)[i].contig_data = (char *) malloc (sizeof(char) * def_contig_size); 
     }

}

void populate_valid_succ (CM_type *cm_c, CM_type_Heap *cm, char *output_str, int *temp_s, cm_cell_list *existing_succs[4], 
                          double collison_rate)
{

   int k=0, loc=0;
   char succ_bp[2];
   char succ_kmer[WIN_SIZE_PLUS];
   //char output_str[WIN_SIZE_PLUS];
   unsigned char result[K_SIZE];
   cm_cell_list *slist=NULL;

   for (k=0; k<4; k++) {
	temp_s[k] = -1;
	existing_succs[k] = NULL;
   }

   for (k=0; k<4; k++) {
        loc=0;
        slist=NULL;
        succ_kmer[0] = '\0';
	change_to_char(succ_bp,k);

        //memset(output_str, '\0', WIN_SIZE_PLUS);
        //convert_binary_to_char(list->kmer_plus_one, output_str);
	     memcpy(&succ_kmer[0], &output_str[1], WIN_SIZE);
	     //memcpy(&succ_kmer[WIN_SIZE], &succ_bp[0], 1);
             succ_kmer[WIN_SIZE] = succ_bp[0];
	     succ_kmer[WIN_SIZE+1] = '\0';

    
	loc = hash31(cm->hasha,cm->hashb,hash_str(succ_kmer)) % cm->width;
        memset(result, 0x00, sizeof(unsigned char)*K_SIZE);
        convert_char_to_binary(succ_kmer, result);

	slist = find_kmer_exists(cm, result, loc);

	if (slist != NULL) {
	   if (slist->kmer_partial_count > 1) {
	       existing_succs[k] = slist;
	       temp_s[k] = slist->kmer_partial_count;
           }
	     //temp_s[k] = abs(CM_Estimate_str(cm_c, list->kmer_plus_one) - CM_Estimate_str(cm_c, slist->kmer_plus_one));
	}
   }

}

void enumerate_contigs_to_array (CM_type *cm_c, CM_type_Heap *cm, contig_entry **c_seq, int *tot_contigs, int *alloc_contig_num, 
                                 size_t max_contig_length, double collison_rate)
{
     cm_cell_list *entry=NULL;// *temp=NULL, *next=NULL;
     int i, j, k, t, loc=0, len=0, pred_cntr=0, num_contigs = *alloc_contig_num;
     int contig_id = *tot_contigs;
     char pred_bp[2];
     char pred_kmer[WIN_SIZE_PLUS];
     int  pred_present[4];
     int temp_s[4];
     cm_cell_list *existing_succs[4];
     int sorted_temp[4];
     bool found_succ = false; 
     int sum_k;
     char output_str[WIN_SIZE_PLUS];
     unsigned char result[K_SIZE];

#ifdef CM_DEBUG
     FILE *fp_d;
     char thread_id[3];
     char debug_file[15];

     sprintf(thread_id, "%d", omp_get_thread_num()); 
     strcpy(debug_file,"debug_");
     strcpy(&debug_file[strlen(debug_file)],thread_id);
     strcpy(&debug_file[strlen(debug_file)],".log");

     //fp_d = fopen ("debug_succ.log", "w");
     fp_d = fopen (debug_file, "w");
     if (fp_d == NULL) {
         printf ("Error opening the debug dump file \n");
         exit (0);
     }
#endif

#pragma omp for schedule(dynamic) private(i)
     for(i=0; i<cm->width; i++) {
           if (cm->tab_entries[i].num_coll_kmers > 0 ) {

              for (int p=0; p<(int)cm->tab_entries[i].kmer_entries.size(); p++)
              {
                   cm_cell_list &temp = cm->tab_entries[i].kmer_entries[p];
 
 
              //while (temp != NULL) {

                       // checking to see if the kmer qualifies as a 'begin kmer' 
                       pred_cntr = 0;
		       for (j=0; j<4; j++)
			    pred_present[j] = 0;

                       memset(output_str, '\0', WIN_SIZE_PLUS);
                       convert_binary_to_char(temp.kmer_plus_one, output_str);
		       for (j=0; j<4; j++) {
			    loc = 0;
			    entry = NULL;
			    //list = NULL;
			    change_to_char(pred_bp,j);

			    memcpy(&pred_kmer[0], &pred_bp[0], 1);
			    memcpy(&pred_kmer[1], &output_str[0], WIN_SIZE);
			    //memcpy(&pred_kmer[1], &temp->kmer_plus_one[0], WIN_SIZE);
			    pred_kmer[WIN_SIZE+1] = '\0';

			    loc = hash31(cm->hasha,cm->hashb,hash_str(pred_kmer)) % cm->width;
                            memset(result, 0x00, sizeof(unsigned char)*K_SIZE);
                            convert_char_to_binary(pred_kmer, result);

			    entry = find_kmer_exists(cm, result, loc);
			    if (entry == NULL) {
				pred_present[j] = 1;
                            } 
                            /*else {
                                if (check_kmer_quality(cm, list, collison_rate) == 1)
                                    pred_present[j] = 1;
                            }*/
                       }
                       for (j=0; j<4; j++)
                           if (pred_present[j] == 1)
                               pred_cntr++;

                       // Begin contig enumeration process
                       if (pred_cntr == 4) 
                       {
                         //check if the 'begin kmer' has a valid succ
                         populate_valid_succ (cm_c, cm, output_str, temp_s, existing_succs, collison_rate);
                         sum_k=0;

                         for (k=0; k<4; k++)
                              sum_k += temp_s[k];

                         if (sum_k > -4) {
				 //list = NULL;
				 memcpy(&(*c_seq)[contig_id].contig_data[0], output_str, (WIN_SIZE + 1));
				 //memcpy(&(*c_seq)[contig_id].contig_data[0], temp->kmer_plus_one, (WIN_SIZE + 1));
                                 __sync_bool_compare_and_swap(&temp.kmer_pred_flag, false, true);
				 len = WIN_SIZE+1;
				 max_contig_length = DEF_CONTIG_SIZE; 
				 //list = temp;
  
                         do {

                           populate_valid_succ (cm_c, cm, output_str, temp_s, existing_succs, collison_rate);

                           // attempting to pick the best successor, among valid succ's
                           memcpy(sorted_temp, temp_s, sizeof(int)*4);
                           qsort (sorted_temp, 4, sizeof(int), compare);
                           int succ_index = -1;
                           found_succ = false;

#ifdef CM_DEBUG
                           int r=0, sum_r=0;
                           for (r=0; r<4; r++) {
                                if (temp_s[r] > -1)
                                    sum_r++;
                           }
                           
                           if (sum_r > 0) {
                               char temp_str[WIN_SIZE_PLUS];
                               fprintf(fp_d,"Kmer: %s, succ's: ", output_str);
                               for (r=0; r<4; r++) {
                                    if (existing_succs[r] != NULL) {
                                        memset(temp_str, '\0', WIN_SIZE_PLUS);
                                        convert_binary_to_char(existing_succs[r]->kmer_plus_one, temp_str);
                                        fprintf(fp_d,"%s %d |", temp_str, temp_s[r]);
                                    }
                               }
                           }
#endif

                           for (k=0; k<4; k++)
			   {
			        if (sorted_temp[k] >= 0)
			        {
				    succ_index = FindIndex(temp_s, 4, sorted_temp[k]);
                                    
                                    //if ( __sync_bool_compare_and_swap(&existing_succs[succ_index]->kmer_pred_flag, 0, 1) == true) {
                                    if ( __sync_bool_compare_and_swap(&existing_succs[succ_index]->kmer_pred_flag, false, true) == true) {

                                       memset(output_str, '\0', WIN_SIZE_PLUS);
                                       convert_binary_to_char(existing_succs[succ_index]->kmer_plus_one, output_str);

                                       (*c_seq)[contig_id].contig_data[len] = output_str[WIN_SIZE];
                                       //memcpy(&(*c_seq)[contig_id].contig_data[len], &output_str[WIN_SIZE], 1);
                                       //memcpy(&(*c_seq)[contig_id].contig_data[len], &existing_succs[succ_index]->kmer_plus_one[WIN_SIZE], 1);
                                       len += 1;
                                       found_succ = true;
#ifdef CM_DEBUG
                                       fprintf(fp_d, "chosen succ: %s k: %d contig_id: %d\n", output_str, succ_index, contig_id);
#endif 
                                       //list = existing_succs[succ_index]; 
                                       
				       if (len == max_contig_length)
				       {
					   char *ptr_2 = (char *) realloc ((char *) (*c_seq)[contig_id].contig_data, 
                                                                           sizeof(char) * (max_contig_length += DEF_CONTIG_SIZE));

					    if (!ptr_2)
					    {
						 printf ("Error (re)allocating memory for length of contig %d \n", contig_id);
						 exit (-1);
					    }
					    else
						 (*c_seq)[contig_id].contig_data = ptr_2;
				       }

                                       break;
                                    }
                                }
                           }

                        } while (found_succ != false);
                         
                        if (len > (WIN_SIZE+1))  {

			       (*c_seq)[contig_id].contig_data[len] = '\0';
			       (*c_seq)[contig_id].my_contig_id = __sync_add_and_fetch(&this_contig_id, 1);
			       contig_id += 1;

			       if (contig_id == num_contigs)
			       {
				   contig_entry *ptr_2 = (contig_entry *) realloc ((*c_seq), sizeof(contig_entry) *
										(num_contigs += DEF_NUM_CONTIGS));

				    if (!ptr_2)
				    {
					 printf ("Error (re)allocating memory for num contigs\n");
					 exit (-1);
				    }
				    else
					 (*c_seq) = ptr_2;

				    for (t=contig_id; t<num_contigs; t++) {
					 (*c_seq)[t].my_contig_id = 0;
					 (*c_seq)[t].contig_data = (char *) malloc (sizeof(char) * DEF_CONTIG_SIZE);
				    }
			       }
                        } 
                     
                     } // end of 'if condition' checking for contig with more than 1 succ
                    } // end of 'if condition' concluding the enumeration process of a contig

                //temp = temp->kl_next;
            } // end of for loop
         }
       } // end of for loop


     *tot_contigs = contig_id;
     *alloc_contig_num = num_contigs;
#ifdef CM_DEBUG
     fclose(fp_d);
#endif

}

bool check_hashes_file_exists (char *hashes_filename, int depth)
{
    int i=0, j=0;
    char str[100];
    char *num;
    bool hash_file_present = false; 
   
    FILE *fp = fopen(hashes_filename, "r");
    if (fp != NULL) {
       printf("File: %s exists \n", hashes_filename);
 
       for (i=0; i<(depth+1); i++) {
           if (fgets(str, sizeof(str), fp) == NULL)
            {
                fprintf(stderr, "Premature EOF if Hashes file \n");
                exit(1);
            }
            else {
               num = strtok (str,",");
               while (num != NULL) {
                      //printf ("%s\n", num);
                      hashes[j] = atoi(num);
                      j++;
                      num = strtok (NULL, ",");
               }
           }
       }

       fclose(fp);
       hash_file_present = true;

       //for (i=0; i<(depth+1); i++)
         //printf("Depth: %d hashes: %d %d \n", i, hashes[i*2], hashes[i*2+1]);

    }
    else {
       printf("File: %s does not exist \n", hashes_filename);
    }

    return hash_file_present;

}

void print_memory_consumption (CM_type *cm, CM_type_Heap *cm_h, int num_threads, contig_thrd_list *global_clist, int *malloc_contig_count)
{
     long int tot_contig_size = Contigs_size (num_threads, global_clist, malloc_contig_count);
     long int cm_size = CMS_Size (cm);
     long int cm_heap_size = CM_Size_Heap (cm_h);

     printf("Memory consumed by CM sketch: %f MB \n", ((double)cm_size/(double)1000000));
     printf("Memory consumed by the Heap (Hash Table): %f MB \n", ((double)cm_heap_size/(double)1000000));
     printf("Memory consumed by all the contigs: %f MB \n", ((double)tot_contig_size/(double)1000000));

}

int main (int argc, char *argv[])
{

	FILE *ptr_alpha, *ptr_reads;

        //float eps=0.01; /* eps must be in this range: [0.01, 1) */
        //float gamma=0.005; /* gamma must be in this range: (0,1) */


	if (argc != 9)
	{
		printf ("Please provide correct number of input arguments: \n");
		printf ("./exe 1:fasta file with reads  2:alphabet_file  3:Coverage  4:Delta 5:w 6:d 7:Heap_Threshold 8: Hashes file\n");
		exit (0);
	}

	char *reads_file_name = argv[1];
	char *alphabet_file_name = argv[2];
        COVERAGE = atoi(argv[3]);
        delta = atoi(argv[4]);
        int width = atoi(argv[5]);
        int depth = atoi(argv[6]);
        Heap_Threshold = atoi(argv[7]);
        char *hashes_filename = argv[8];

	if (((WIN_SIZE+1)%4) != 0) {
	     printf("Invalid k-mer size!! Please enter a k-mer of size divisible by 4 \n");
	     return 1;
	}

	srand (time (NULL));

        hashes = (unsigned int *)calloc( sizeof(unsigned int), (2*(depth+1)) );
        bool is_hash_present = check_hashes_file_exists (hashes_filename, depth);

	//Counting the number of reads in the reads file
	ptr_reads = fopen (reads_file_name, "r");
	/*check to see if it opened okay */
	if (ptr_reads == NULL)
	{
		printf ("Error opening the reads file \n");
		exit (0);
	}

	int num_lines=0;
	char c;
	while((c = fgetc(ptr_reads)) != EOF)
	{
		if(c == '>')
			num_lines++;
	}
	printf("number of reads in input file:%d\n",num_lines);
	fclose(ptr_reads);

	//initialize the input read sequence data structure
	Initialize_input_read_seq (&rd_seq, num_lines);

	/*open input reads file */
	ptr_reads = fopen (reads_file_name, "r");
	/*check to see if it opened okay */
	if (ptr_reads == NULL)
	{
		printf ("Error opening the reads file \n");
		exit (0);
	}

	/*open alphabet file */
	ptr_alpha = fopen (alphabet_file_name, "r");
	/*check to see if it opened okay */
	if (ptr_alpha == NULL)
	{
		printf ("Error opening alphabet file \n");
		exit (0);
	}

	parse_alphabet_file (ptr_alpha);

	sigma_size = strlen (alphabet_array);

	fclose (ptr_alpha);

	double start_t = omp_get_wtime ();

	// Parsing the fasta while with reads
	parse_input_reads_file(ptr_reads, rd_seq, num_lines, ADJUST_SEQ_SIZE);

	double end_t = omp_get_wtime ();
	printf ("\nTime for loading the reads into memory: %f secs \n", (end_t - start_t));

        /* Initializing the CM_sketch data structure */

	start_t = omp_get_wtime ();

        CM_type *cm = CM_Init(width, depth, time(NULL), is_hash_present, hashes_filename);
        CM_type_Heap *cm_h = CM_Heap_Init(width*depth*100, depth, time(NULL), is_hash_present, hashes_filename);
        
        free(hashes);

#pragma omp parallel shared(rd_seq, cm, cm_h, this_kmer_id) //num_threads(1)
	{
		int i=0;
		if(omp_get_thread_num() == 0) printf("Number of threads: %d Kmer size: %d Coverage: %d Delta: %d CM_w: %d CM_d: %d Heap_width: %d Heap_Threshold: %d \n", 
                   omp_get_num_threads(), WIN_SIZE, COVERAGE, delta, cm->width, cm->depth, cm_h->width, Heap_Threshold);

		//Populating the CM_sketch table
#pragma omp for schedule(static) private(i) firstprivate(num_lines) nowait
		for (i=0; i<num_lines; i++) 
		{
			Sliding_window(i,cm,cm_h);
		}

	} // end of parallel region
	end_t = omp_get_wtime ();
	printf ("\nTime for populating entries to the CM table: %f secs \n", (end_t - start_t));

	int num_entries = count_hashtab_entries(cm_h);
        double collison_rate = (double)num_entries/(double)(width*depth*100);

	printf("There are %d (distinct) entries in the CM table, Collison_rate(k): %f \n", num_entries, collison_rate);

#ifdef CM_DEBUG_PRINT
	//print_hashtab_entries(cm, cm_h, collison_rate);
#endif

	fclose(ptr_reads);
	free_reads(num_lines);

        printf("Total CM entries: %lld \n", cm->count);

        //total_all_kmer_count(my_hash_table);
/* ------------------------------- END OF UPDATING THE CM TABLE ---------------------------------------------------------------------*/

 // ************* begin contig enumeration process ****************************

	start_t = omp_get_wtime ();

        int num_threads=0; 
        // shared array size = # of threads 
        contig_thrd_list *global_clist = NULL;
        int *malloc_contig_count;
        int t=0, tot_num_contigs=0;

#pragma omp parallel shared(cm, cm_h, this_contig_id, global_clist, COVERAGE, delta) //num_threads(1)
{
        int i=0;
        num_threads = omp_get_num_threads();
        contig_entry *contig_list = NULL;

        // master thread initializes the shared array holding pointers for each inidividual thread's contig list
        if (omp_get_thread_num() == 0) {
            global_clist = (contig_thrd_list *) malloc (sizeof(contig_thrd_list) * num_threads);
            malloc_contig_count = (int *) calloc (num_threads, sizeof(int));
  
            for (i=0; i<num_threads; i++) {
                 global_clist[i].contig_count = 0;
                 global_clist[i].c_series = NULL;
            }
        }

        initialize_contig_array_per_thread(&contig_list, DEF_NUM_CONTIGS, DEF_CONTIG_SIZE);
    
        int tot_contigs = 0;
        int alloc_contig_num = DEF_NUM_CONTIGS;

        enumerate_contigs_to_array(cm, cm_h, &contig_list, &tot_contigs, &alloc_contig_num, DEF_CONTIG_SIZE, collison_rate);
       
        global_clist[omp_get_thread_num()].c_series = contig_list;
        global_clist[omp_get_thread_num()].contig_count = tot_contigs;

        malloc_contig_count[omp_get_thread_num()] = alloc_contig_num;

        printf("thread: %d, contig_count: %d, malloc'd count: %d \n", omp_get_thread_num(), tot_contigs, alloc_contig_num);

        //free_contigs(contig_list, tot_contigs);

} // end of second parallel region

	end_t = omp_get_wtime ();

#pragma omp parallel for reduction(+:tot_num_contigs)
        for (t=0; t<num_threads; t++)
            tot_num_contigs += global_clist[t].contig_count;

	printf ("\nTime for populating the Contigs: %f secs, Total number of Contigs: %d \n\n", (end_t - start_t), tot_num_contigs);

#ifdef CM_DEBUG_PRINT
	print_hashtab_entries(cm, cm_h, collison_rate);
#endif

        char output_file_name[250];
        char num_to_string[20];
        char kwin_to_string[3];
        char cov_to_string[4];
        char delta_to_string[3];
        char w_to_string[8];
        char d_to_string[4];
        char thresh_to_string[6];
        char thrd_to_string[4];

	sprintf(num_to_string, "%d", num_lines);
	sprintf(kwin_to_string, "%d", WIN_SIZE);
	sprintf(cov_to_string, "%d", COVERAGE);
	sprintf(delta_to_string, "%d", delta);
	sprintf(w_to_string, "%d", width);
	sprintf(d_to_string, "%d", depth);
	sprintf(thresh_to_string, "%d", Heap_Threshold);
	sprintf(thrd_to_string, "%d", num_threads);

        if (num_threads == 1)
            strcpy(output_file_name,"serial_output_n");
        else
	    strcpy(output_file_name,"contig_output_n");

	strcpy(&output_file_name[15],num_to_string);
	strcpy(&output_file_name[strlen(output_file_name)],"_k");
	strcpy(&output_file_name[strlen(output_file_name)],kwin_to_string);
	strcpy(&output_file_name[strlen(output_file_name)],"_cov");
	strcpy(&output_file_name[strlen(output_file_name)],cov_to_string);
	strcpy(&output_file_name[strlen(output_file_name)],"_del");
	strcpy(&output_file_name[strlen(output_file_name)],delta_to_string);
	strcpy(&output_file_name[strlen(output_file_name)],"_w");
	strcpy(&output_file_name[strlen(output_file_name)],w_to_string);
	strcpy(&output_file_name[strlen(output_file_name)],"_d");
	strcpy(&output_file_name[strlen(output_file_name)],d_to_string);
	strcpy(&output_file_name[strlen(output_file_name)],"_thresh");
	strcpy(&output_file_name[strlen(output_file_name)],thresh_to_string);
	strcpy(&output_file_name[strlen(output_file_name)],"_thrd");
	strcpy(&output_file_name[strlen(output_file_name)],thrd_to_string);
	strcpy(&output_file_name[strlen(output_file_name)],".fasta");

        FILE *f = fopen(output_file_name, "w");
        if (f == NULL)
        {
            printf("Error opening file!\n");
            exit(1);
        }
        
        int i=0, j=0, k=0, l=0;
        char contig_id[15];
        char output_contig_name[100];
        int num_per_line = 0, contig_length = 0, chars_per_line = 60;

        for (i=0; i<num_threads; i++) {
            for (j=0; j<global_clist[i].contig_count; j++) {
                 strcpy(output_contig_name, ">contig_");
                 sprintf(contig_id, "%ld", global_clist[i].c_series[j].my_contig_id);
                 strcpy(&output_contig_name[strlen(output_contig_name)], contig_id);
                 fprintf(f, "%s\n", output_contig_name);

                 contig_length = strlen(global_clist[i].c_series[j].contig_data);
                 num_per_line = contig_length / chars_per_line;
                 for (k=0; k<num_per_line; k++) {
                      for (l=(k*chars_per_line); l<((k+1)*chars_per_line); l++) {
                           fprintf(f, "%c", global_clist[i].c_series[j].contig_data[l]);
                      }
                      fprintf(f, "\n");
                 }
                 for (k=(num_per_line*chars_per_line); k<contig_length; k++) {
                      fprintf(f, "%c", global_clist[i].c_series[j].contig_data[k]);
                 }
                 fprintf(f, "\n");
             }
        } 
      
        print_memory_consumption(cm, cm_h, num_threads, global_clist, malloc_contig_count);

        fclose(f);
        free_contigs(global_clist, num_threads, malloc_contig_count);
        free (malloc_contig_count);

 // end contig enumeration process

        CM_Destroy(cm); 
        CM_Heap_Destroy(cm_h); 

}
