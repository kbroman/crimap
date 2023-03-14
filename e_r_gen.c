#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"

#define MAX_ID_SIZE 80 
read_gen_file(fp, data, g_all)
     FILE *fp;
     struct data *data;
     char g_all;
/* reads a .gen file into a data object. If g_all is 0, ignores genotype info.
 and locusnames (this is for purposes of chrompic)*/ 
{
  SHORT i,j,k, num_fams, num_loci, temp;
  char *our_alloc();
  struct individual *ind_p;      
  char loc_name[MAX_LEN_LOC_NAME];
  SHORT *genos;
  char fam_idvec[MAX_ID_SIZE];

  fscanf(fp, "%hd", &num_fams);
  fscanf(fp, "%hd", &num_loci);
  data->num_loci = num_loci;
  data->num_fams = num_fams;
  if (g_all)
    data->locus_names = (char **)our_alloc((ALLOC)num_loci*sizeof(char *));

  for (k = 0; k < num_loci; k++){
    fscanf(fp,"%s",loc_name);
    if (g_all) {
      data->locus_names[k] = (char *)our_alloc((ALLOC)(1 + strlen(loc_name))
					     * sizeof(char));
      strcpy(data->locus_names[k], loc_name);
    }
  }
  
  data->num_mems = (SHORT *)our_alloc((ALLOC)num_fams * sizeof(SHORT));
  data->fam_id = (char **)our_alloc((ALLOC)num_fams * sizeof(char *));

  data->ind = (struct individual ***)
    our_alloc((ALLOC)num_fams * sizeof(struct individual **));
  for (i = 0; i < num_fams; i++){
    fscanf(fp,"%s", fam_idvec);
    data->fam_id[i] = (char *)our_alloc((ALLOC)(1 + strlen(fam_idvec)) *
					sizeof(char));
    strcpy(data->fam_id[i], fam_idvec);
    fscanf(fp,"%hd", &data->num_mems[i]);
    data->ind[i] = (struct individual **)our_alloc((ALLOC)data->num_mems[i]*sizeof(struct individual *));
    for (j = 0; j < data->num_mems[i]; j++){
      data->ind[i][j] = ind_p = (struct individual *)our_alloc((ALLOC)sizeof(struct individual));
      fscanf(fp, "%ld", &ind_p->id);
      fscanf(fp, "%ld", &ind_p->moth_id);
      fscanf(fp, "%ld", &ind_p->fath_id);
      fscanf(fp, "%hd", &temp);
      ind_p->sex = temp;
      genos = (SHORT *)our_alloc((ALLOC)2 * num_loci * 
						sizeof(SHORT)); 
      if (g_all) {
	ind_p->a = (SHORT **)our_alloc((ALLOC)num_loci*sizeof(SHORT *));
	ind_p->a[0] = genos;
      }
      for (k = 0; k < 2 * num_loci; k++) fscanf(fp,"%hd", genos + k);
      if (g_all) for (k = 1; k < num_loci; k++) ind_p->a[k] = genos += 2;
      else our_free(genos);
    }  /* close j */
    
    for (j = 0; j < data->num_mems[i]; j++){
      ind_p = data->ind[i][j];
      if(ind_p->moth_id){
        for (k = 0; k < data->num_mems[i]; k++){
	  if (ind_p->moth_id == data->ind[i][k]->id){
	    ind_p->moth = k;
	    break;
	  }
	} /* close k */
        for (k = 0; k < data->num_mems[i]; k++)
	  if (ind_p->fath_id == data->ind[i][k]->id){
	    ind_p->fath = k;
	    break;
	  }
      }  /* end if */
    } /* close j */
  } /* close i */
}        


alloc_c_dat(data, chrom_data, pk_chrom_data)
     struct data *data;
     struct chrom_data *chrom_data, *pk_chrom_data;
{
  SHORT i, j, k, num_fams, num_chroms, num_loci;
  char *our_alloc();
  
  chrom_data->num_loci = pk_chrom_data->num_loci = num_loci = data->num_loci;
  chrom_data->num_fams = num_fams = data->num_fams;
  pk_chrom_data->num_fams = 1;
  pk_chrom_data->locus_names = chrom_data->locus_names = data->locus_names;

  /* note that pk_chrom_data & chrom_data do not have their own allocated space for several objects, including locus_names */

  chrom_data->num_chroms = (SHORT *)our_alloc((ALLOC)num_fams*sizeof(SHORT));
  pk_chrom_data->num_chroms = (SHORT *)our_alloc((ALLOC)sizeof(SHORT));
  
  chrom_data->fam_nums_array = data->fam_id;
  chrom_data->phase_choices = (struct phase **)
    our_alloc((ALLOC)num_fams * sizeof(struct phase *));
  chrom_data->chrom_array = (char ***)our_alloc((ALLOC)num_fams * sizeof(char **));
  pk_chrom_data->chrom_array = (char ***)our_alloc((ALLOC)sizeof(char **));
  pk_chrom_data->num_chroms[0] = 0;
  for(i = 0; i < num_fams; i++){
    chrom_data->num_chroms[i] = 0;
    chrom_data->phase_choices[i] = (struct phase *)
      our_alloc((ALLOC)sizeof(struct phase));
    chrom_data->phase_choices[i]->num_switches = 0;
    
    for(j = 0; j < data->num_mems[i]; j++)
      if(data->ind[i][j]->moth_id){
        chrom_data->num_chroms[i] += 2;
        pk_chrom_data->num_chroms[0] += 2;
      }  /* end if */
    num_chroms = chrom_data->num_chroms[i];
    chrom_data->chrom_array[i] = (char **)
      our_alloc((ALLOC)num_chroms * sizeof(char *));
    for (j = 0; j < num_chroms; j++){
      chrom_data->chrom_array[i][j] = (char *)
	our_alloc((ALLOC)num_loci * sizeof(char)); 
      for (k = 0; k < num_loci; k++) chrom_data->chrom_array[i][j][k] = 'X';
    } 
  } /* close i */
  num_chroms = pk_chrom_data->num_chroms[0]; 
  pk_chrom_data->chrom_array[0] = (char **)
    our_alloc((ALLOC)num_chroms * sizeof(char *));
  for (j = 0; j < num_chroms; j++){
    pk_chrom_data->chrom_array[0][j] = (char *)
      our_alloc((ALLOC)num_loci * sizeof(char)); 
    for (k = 0; k < num_loci; k++) pk_chrom_data->chrom_array[0][j][k] = 'X';
  }
}        


