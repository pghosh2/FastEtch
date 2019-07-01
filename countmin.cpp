/********************************************************************
Count-Min Sketches

G. Cormode 2003,2004

Updated: 2004-06 Added a floating point sketch and support for 
                 inner product point estimation
Initial version: 2003-12

This work is licensed under the Creative Commons
Attribution-NonCommercial License. To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc/1.0/ or send a letter
to Creative Commons, 559 Nathan Abbott Way, Stanford, California
94305, USA. 
*********************************************************************/

#include <iostream>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "prng.h"
#include "massdal.h"
#include "countmin.h"
#include <omp.h>

#define min_cm(x,y)	((x) < (y) ? (x) : (y))
//#define max(x,y)	((x) > (y) ? (x) : (y))

double eps;	               /* 1+epsilon = approximation factor */
double Delta;                  /* probability of failure */

//int bits=32;

/************************************************************************/
/* Routines to support Count-Min sketches                               */
/************************************************************************/

CM_type * CM_Init(int width, int depth, int seed, bool is_hash_present, char *hashes_filename)
{     // Initialize the sketch based on user-supplied size
  CM_type * cm;
  int j,k;
  prng_type * prng;

  cm=(CM_type *) malloc(sizeof(CM_type));
  prng=prng_Init(-abs(seed),2); 
  // initialize the generator to pick the hash functions

  if (cm && prng)
    {
      cm->depth=depth;
      cm->width=width;
      cm->count=0;
      cm->counts=(int **)calloc(sizeof(int*),cm->width);
      //cm->counts[0]=(CM_cell_type *)calloc(sizeof(CM_cell_type), cm->depth*cm->width);
      cm->hasha=(unsigned int *)calloc(sizeof(unsigned int),cm->depth);
      cm->hashb=(unsigned int *)calloc(sizeof(unsigned int),cm->depth);
      if (cm->counts && cm->hasha && cm->hashb)  //&& cm->counts[0])
	{
          if (is_hash_present) {
	      for (j=0;j<depth;j++)
	      {
	           cm->hasha[j]= hashes[(2*j)];
	           cm->hashb[j]= hashes[(2*j)+1];
	           // pick the hash functions
	           //cm->counts[j]=(CM_cell_type *)cm->counts[0]+(j*cm->width);
	      }
          } else {
	      for (j=0;j<depth;j++)
	      {
	           cm->hasha[j]=prng_int(prng) & MOD;
	           cm->hashb[j]=prng_int(prng) & MOD;
	           // pick the hash functions
	           //cm->counts[j]=(CM_cell_type *)cm->counts[0]+(j*cm->width);
	      }
          }  
          for (k=0;k<width;k++)
               cm->counts[k] = (int *) calloc(sizeof(int),cm->depth);
        }
      else cm=NULL;
    }

  for (j=0;j<depth;j++)
  {
     printf("For d=%d, hash_a:%d, hash_b:%d \n",j, cm->hasha[j],cm->hashb[j]);
  }

  if (!is_hash_present) {
     FILE *fp_h = fopen (hashes_filename, "w"); 
     if (fp_h == NULL) {
         printf ("Error opening the hashes dump file \n");
         exit (0);
     }

     for (j=0;j<depth;j++)
          fprintf(fp_h,"%d,%d \n", cm->hasha[j],cm->hashb[j]);
    
     fclose (fp_h);
  }

  free(prng);
  return cm;
}

CM_type_Heap * CM_Heap_Init(int width, int depth, int seed, bool is_hash_present, char *hashes_filename)
{     // Initialize the sketch based on user-supplied size
  CM_type_Heap * cm;
  int k;
  prng_type * prng;

  //cm=(CM_type_Heap *) malloc(sizeof(CM_type_Heap));
  cm = new CM_type_Heap;
  prng=prng_Init(-abs(seed),2); 
  // initialize the generator to pick the hash functions

  if (cm && prng)
    {
      cm->width=width;
      cm->tab_entries= new CM_cell_entry[width];
      //cm->tab_entries=(CM_cell_entry *)calloc(sizeof(CM_cell_entry),cm->width);
      //cm->counts[0]=(CM_cell_type *)calloc(sizeof(CM_cell_type), cm->depth*cm->width);
   
      if (!is_hash_present) {
          cm->hasha=prng_int(prng) & MOD;
          cm->hashb=prng_int(prng) & MOD;
      } else {
          cm->hasha=hashes[(2*depth)];
          cm->hashb=hashes[(2*depth)+1];
      }

      if (cm->tab_entries)
	{
          for (k=0;k<width;k++) {
               //cm->tab_entries[k].kmer_entries = NULL;
               cm->tab_entries[k].num_coll_kmers = 0;
               omp_init_lock(&cm->tab_entries[k].writelock);
          }
        }
      else cm=NULL;
    }

    printf("Heap hashes: hash_a:%d, hash_b:%d \n\n", cm->hasha,cm->hashb);
 
    if (!is_hash_present) {
       FILE *fp_h = fopen (hashes_filename, "a"); 
       if (fp_h == NULL) {
           printf ("Error opening the hashes dump file \n");
           exit (0);
       }

       fprintf(fp_h,"%d,%d \n", cm->hasha,cm->hashb);
    
       fclose (fp_h);
    }

    free(prng);
    return cm;
}

/*CM_type * CM_Copy(CM_type * cmold)
{     // create a new sketch with the same parameters as an existing one
  CM_type * cm;
  int j;

  if (!cmold) return(NULL);
  cm=(CM_type *) malloc(sizeof(CM_type));
  if (cm)
    {
      cm->depth=cmold->depth;
      cm->width=cmold->width;
      cm->count=0;
      cm->counts=(int **)calloc(sizeof(int *),cm->depth);
      cm->counts[0]=(int *)calloc(sizeof(int), cm->depth*cm->width);
      cm->hasha=(unsigned int *)calloc(sizeof(unsigned int),cm->depth);
      cm->hashb=(unsigned int *)calloc(sizeof(unsigned int),cm->depth);
      if (cm->counts && cm->hasha && cm->hashb && cm->counts[0])
	{
	  for (j=0;j<cm->depth;j++)
	    {
	      cm->hasha[j]=cmold->hasha[j];
	      cm->hashb[j]=cmold->hashb[j];
	      cm->counts[j]=(int *) cm->counts[0]+(j*cm->width);
	    }
	}
      else cm=NULL;
    }
  return cm;
}
*/

void CM_Destroy(CM_type * cm)
{     // get rid of a sketch and free up the space
  int i;

  if (!cm) return;
  if (cm->counts)
    {
      for (i=0; i<cm->width; i++) {
           free(cm->counts[i]);
      }
                  	
      free(cm->counts);
      cm->counts=NULL;
    }
  if (cm->hasha) free(cm->hasha); cm->hasha=NULL;
  if (cm->hashb) free(cm->hashb); cm->hashb=NULL;
  free(cm);  cm=NULL;
}

void CM_Heap_Destroy(CM_type_Heap * cm)
{     // get rid of a sketch and free up the space
  int i;
  //cm_cell_list *head=NULL, *tmp=NULL;

  if (!cm) return;
  if (cm->tab_entries)
    {
      for (i=0; i<cm->width; i++) {
           /*head = cm->tab_entries[i].kmer_entries;
           while (head != NULL) {
                  tmp = head;
                  head = head->kl_next;
                  free(tmp);
           }*/
           cm->tab_entries[i].kmer_entries.clear();
           omp_destroy_lock(&cm->tab_entries[i].writelock);
      }
        
      delete[] cm->tab_entries;          	
      //free(cm->tab_entries);
      //cm->tab_entries=NULL;
    }
  //if (cm->hasha) free(cm->hasha); cm->hasha=NULL;
  //if (cm->hashb) free(cm->hashb); cm->hashb=NULL;
  //free(cm);  cm=NULL;
  delete[] cm;
}

int CM_Size(CM_type * cm)
{ // return the size of the sketch in bytes
  int counts, hashes, admin;
  if (!cm) return 0;
  admin=sizeof(CM_type);
  counts=cm->width*cm->depth*sizeof(CM_cell_entry);
  hashes=cm->depth*2*sizeof(unsigned int);
  return(admin + hashes + counts);
}

int convert_hash_fcn (const char* word)
{
    unsigned int hash = 0;
    for (int i = 0 ; word[i] != '\0' ; i++)
    {
        hash = 31*hash + word[i];
    }
    return hash;
}

unsigned int hash_str(const char *str)
{
    unsigned long hash = 5381;
    int c;
    while ((c = *str++)) 
    {
           hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
    }

    return hash;
}

/*void CM_Update(CM_type * cm, unsigned int item, int diff)
{
  int j,loc;

  if (!cm) return;
  cm->count+=diff;

  for (j=0;j<cm->depth;j++)
  {
    loc = hash31(cm->hasha[j],cm->hashb[j],item) % cm->width;
    cm->counts[j][loc].aggregate_count+=diff;
  }
  // this can be done more efficiently if the width is a power of two
}
*/

unsigned createMask(unsigned a, unsigned b)
{
   unsigned r = 0;
   for (unsigned i=a; i<=b; i++)
       r |= (1 << i);

   return r;
}

void convert_char_to_binary (char *kmer_name, unsigned char *result)
{

int i=0,j=0,k=0;

for (i=0; i<WIN_SIZE+1; i++) 
 {
     if ((i%4 == 0) && (i>0)) {
         j++;
         k=0;
     }

     switch (kmer_name[i])
     {
           case 'A': 
               result[j] |= 0 << (7 - (2*k));
               result[j] |= 0 << (6 - (2*k));
               //printf("i: %d, j: %d, result: %d  %d  \n", i, j, result[0], result[1]);
               break;

           case 'C': 
               result[j] |= 0 << (7 - (2*k));
               result[j] |= 1 << (6 - (2*k));
               //printf("i: %d, j: %d, result: %d  %d  \n", i, j, result[0], result[1]);
               break;

           case 'G': 
               result[j] |= 1 << (7 - (2*k));
               result[j] |= 0 << (6 - (2*k));
               //printf("i: %d, j: %d, result: %d  %d  \n", i, j, result[0], result[1]);
               break;

           case 'T': 
               result[j] |= 1 << (7 - (2*k));
               result[j] |= 1 << (6 - (2*k));
               //printf("i: %d, j: %d, result: %d  %d  \n", i, j, result[0], result[1]);
               break;
     }
     k++;
 }

}

void convert_binary_to_char (unsigned char *result, char output_str[WIN_SIZE_PLUS])
{
 
int i=0,j=0,len=0;
//len = *out_len;

for (i=0; i<K_SIZE; i++)
{
     for (j=7; j>=0; j-=2) {
          unsigned mask1 = createMask(j-1,j);
          unsigned temp = mask1 & result[i];
          if (j>1)
              temp = temp >> (j-1);
          //printf("i: %d, j: %d, mask1: %d, input: %d, temp: %d \n", i, j, mask1, result[i], temp);
          
	  switch (temp)
	  {
		 case 0: 
                       strcpy(&output_str[len],"A");
		       break;

		 case 1: 
                       strcpy(&output_str[len],"C");
		       break;

		 case 2: 
                       strcpy(&output_str[len],"G");
		       break;

		 case 3: 
                       strcpy(&output_str[len],"T");
		       break;
	  }
          len++;
     }     
}

output_str[len] = '\0';
//*out_len = len;

}

cm_cell_list *find_kmer_exists(CM_type_Heap *cm, unsigned char key[K_SIZE], int loc)
{
     std::vector<cm_cell_list>::iterator it;
     
     it = std::find_if(cm->tab_entries[loc].kmer_entries.begin(), cm->tab_entries[loc].kmer_entries.end(),
                      [key] (const cm_cell_list& d) {
                             return (memcmp(d.kmer_plus_one,key,sizeof(unsigned char)*K_SIZE) == 0);
                             //return d.kmer_plus_one == key;
                      });

     if (it != cm->tab_entries[loc].kmer_entries.end()) {
         return &(cm->tab_entries[loc].kmer_entries[std::distance(cm->tab_entries[loc].kmer_entries.begin(), it)]);
     } 
     
     return 0;

}

/*cm_cell_list *find_kmer_exists(CM_type_Heap *cm, unsigned char key[K_SIZE], int loc, cm_cell_list **next)
{
     cm_cell_list *list , *last_but_one=NULL;

     //for(list = cm->tab_entries[loc].kmer_entries; list != NULL; list = list->kl_next) {
     for(list = cm->tab_entries[loc].kmer_entries; list != NULL; list = list->kl_next) {
                if (list->kl_next == NULL)
                        last_but_one = list;

                if (memcmp(list->kmer_plus_one,key,sizeof(unsigned char)*K_SIZE) == 0)
                    return list;
        }

        *next = last_but_one;
        return NULL;

}*/


void CM_Update_insert (CM_type_Heap *cm, char *key, int loc, int val)
{
    //cm_cell_list *start=NULL, *next=NULL;
    unsigned char result[K_SIZE];

        memset(result, 0x00, sizeof(unsigned char)*K_SIZE);
        convert_char_to_binary(key, result);

#ifdef BIT_CONV
        printf("Insert: key: %s \n", key);
        for (int i=0; i<K_SIZE; i++)
             printf("%d ", result[i]);
        printf("\n");
#endif

        cm_cell_list *start = find_kmer_exists(cm, result, loc);

        if (start == NULL) /* Key does not exist in the Heap, thus add */
        {
          //if (cm_sketch->counts[k_w][k_d].insertion_count < GAMMA_THRESH) 
              cm->tab_entries[loc].num_coll_kmers++;
 
              cm_cell_list first;
              memcpy(first.kmer_plus_one, result, sizeof(unsigned char)*K_SIZE);
              first.kmer_pred_flag = false;
              first.kmer_partial_count = 1;

              cm->tab_entries[loc].kmer_entries.push_back(first);

#ifdef TAB_DEBUG
        printf("Inserted Key: %s \n", key);
#endif
        }
        else {
             start->kmer_partial_count++;
        }

}

void CM_Update_str(CM_type *cm, CM_type_Heap *cm_h, char *key, int val, int Heap_Threshold)
{
    int i=0,loc=0,ans=0,temp=0,h_bucket=0;

    if (!cm) return;

//#pragma omp atomic
//     cm->count+=val;
    __sync_add_and_fetch(&cm->count, val);

    //unsigned int hashval = convert_hash_fcn(key);
    unsigned int hashval = hash_str(key);

    loc = hash31(cm->hasha[0],cm->hashb[0],hashval) % cm->width;
    ans =  __sync_add_and_fetch(&cm->counts[loc][0], val);
#ifdef TAB_DEBUG
    printf("For key: %s, d0: %d,", key, ans);
#endif

    for (i=1; i<cm->depth;i++) {
        loc = hash31(cm->hasha[i],cm->hashb[i],hashval) % cm->width;
        temp = __sync_add_and_fetch(&cm->counts[loc][i], val);
        ans = min(ans,temp);

#ifdef TAB_DEBUG
    printf(" d%d: %d,", i, temp);
#endif
    }

#ifdef TAB_DEBUG
  printf("  Ans: %d \n", ans);
#endif 

    if (ans >= Heap_Threshold) {
        h_bucket = hash31(cm_h->hasha,cm_h->hashb,hashval) % cm_h->width;
        omp_set_lock(&cm_h->tab_entries[h_bucket].writelock);
        CM_Update_insert (cm_h, key, h_bucket, val);
        omp_unset_lock(&cm_h->tab_entries[h_bucket].writelock);
    }

}

int CM_PointEst(CM_type * cm, char *key)
{
  // return an estimate of the count of an item by taking the minimum
  int j, ans, loc=0;
  unsigned int query = hash_str(key);

  if (!cm) return 0;
  ans=cm->counts[hash31(cm->hasha[0],cm->hashb[0],query) % cm->width][0];

#ifdef SUCC_DEBUG
  printf("d:0 ans:%d ", ans);
#endif 

  for (j=1;j<cm->depth;j++) {
       loc = cm->counts[hash31(cm->hasha[j],cm->hashb[j],query)%cm->width][j];
       ans=min(ans,loc);

#ifdef SUCC_DEBUG
       printf("d:%d ans:%d ", j, loc);
#endif 

  }

#ifdef SUCC_DEBUG
  printf("\n");
#endif 

  // this can be done more efficiently if the width is a power of two
  return (ans);
}

int CM_Estimate_str(CM_type *cm, char *key)
{
    //unsigned int hashval = convert_hash_fcn(key);
    unsigned int hashval = hash_str(key);
    //printf("Estimate: hashval: %d \n", hashval);
    //return CM_PointEst(cm, hashval);
    return CM_PointEst(cm, key);
}

