#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"
#include "var1.h"
#include "var2.h"

float **r_ratio;
float *r_ratio_vec;
char **r_r;
char *r_r_vec;
char *rec_list;    /* indicator vector: tells whether corresponding interval in 
		      sorted list has a recombination under the optimal 
		      switch */
extern char gen_file[];

struct cross {
     struct cross *next_cross;
     char *famnum;
     ID sibnum;
     char parent;
};  

double chrompics_mle(ad_num_its, order, num_sub_loci, chrom_data, ad_num_tswitch_elim)
     SHORT *ad_num_its;
     SHORT *order;
     SHORT num_sub_loci;
     SHORT *ad_num_tswitch_elim;
     struct chrom_data *chrom_data;
{
     double log_like, prev_log_like; 
     SHORT num_its, pk_num_its, num_fams, num_tswitch_elim, k, i, j;
     double super_punk(), change_theta();
     SHORT make_fl_sw_int();
     SHORT our_free();

     num_fams = chrom_data->num_fams;
     num_tswitch_elim = make_fl_sw_int(order, num_sub_loci, chrom_data);
     alloc_rec_pos(num_fams);
     num_its = 0;
     log_like = super_punk(num_fams, num_tswitch_elim);

     do {

       for(k = 0; k < theta->num_types; k++){
        	for(i = 0; i < theta->n; i++){
       	        	for(j = 0; j < theta->n-i; j++){
       		        	(nrecs->data)[k][i][j] = num_mei->data[k][i][j] - recs->
                                   data[k][i][j];
     	        	}
               	}
       }

	change_theta(&pk_num_its);
     	prev_log_like = log_like;
        log_like = super_punk(num_fams, num_tswitch_elim);
     	num_its++;

     } while( (log_like - prev_log_like) > TOL );
       
     if(log_like < prev_log_like - MAX_DECREASE)
     	printf("\n\nLOG_LIKE DECREASED IN MLE, DIFF = %f\n\n",
     		log_like - prev_log_like);
     *ad_num_its = num_its;
     *ad_num_tswitch_elim = num_tswitch_elim;
     our_free(r_prod);
    our_free(r_prod_vec);
    chrompics_alloc_rec_pos(num_fams);
    chrompics(order, num_sub_loci, chrom_data);
    free_fl_sw_int(num_fams);
    return(log_like);
}

chrompics_alloc_rec_pos(num_fams)
   SHORT num_fams;
{
   char *our_alloc();
   SHORT i_fam, num_flanks;
   SHORT i_fl, sup_nfl;
   LINDEX n_rprod, sup_rprod;

   sup_nfl = sup_rprod = 0;
   for (i_fam = 0; i_fam < num_fams; i_fam++){
      if ((num_flanks = all_flanks[i_fam].num_flanks) <= 1) continue;
      if (sup_nfl < num_flanks)	 sup_nfl = num_flanks;
      n_rprod = 0;
      for (i_fl = 0; i_fl < num_flanks; i_fl++)
         n_rprod += 1 <<  all_flanks[i_fam].num_tswitchs[i_fl];
      if (sup_rprod < n_rprod) sup_rprod = n_rprod;
   }
   r_ratio = (float **)our_alloc((ALLOC)sup_nfl * sizeof(float *));
   r_ratio_vec = (float *)our_alloc((ALLOC)sup_rprod * sizeof(float));
   r_r = (char **)our_alloc((ALLOC)sup_nfl * sizeof(char *));
   rec_list = (char *)our_alloc((ALLOC)sup_nfl * sizeof(char));
   r_r_vec = (char *)our_alloc((ALLOC)sup_rprod * sizeof(char));
}

chrompics(order, num_loc, chrom_data)
   SHORT *order;        /* order of loci being tested */
   SHORT num_loc;    /* # loci in order */
   struct chrom_data *chrom_data;
{
   SHORT num_fams, num_chroms, i_fam, i_fl, num_flanks, n_loc, i_chrom, i_loc;
   SHORT i, i1, i_cross, i_sw, num_sw, num_cross;
   ID sibnum;
   SHORT *locus_nums;
   LINDEX m, m_offset, pj, j, k, r_ratio_ind;
   LINDEX term_ind[2];
   float *t_ratio;
   float sup_rat;
   char *famnum;
   double *t_prod;
   double term[2];
   double t1t;
   char *t_r;
   char r_term[2];
   char i_rec, sr, rj;
   char e1, e2, pk_flag, prev_pk_flag, pred_val, new_val, prev_val, pk_anchor;
   char parent, prev_cross, sep, prev_name;
   char **sw_array, **array, **chrom_array;  /* array of phase chromosomes */
   char *our_alloc();
   SHORT our_free();
   struct phase *phase;
   struct flank_tswitchs *i_fl_sw; 
   struct intervals *i_int;
   struct cross ***cross_mat;
   struct cross *cross;
   SHORT j_mem, j_count;
   double find_phase_change();
   double plike;
   struct data *data;
   FILE *fp;
   FILE *fopen();
   char pchar[5];
   LINDEX sh_mask, id_mask, rec_mask;
   char sh_dir;

/* globals needed: r_ratio, r_ratio_vec; prod1, rec, prod, pos, offset
   r_r, r_r_vec;
*/
   pchar[0] = '0';
   pchar[1] = '1';
   pchar[2] = 'o';
   pchar[3] = 'i';
   pchar[4] = '-';
/* above are the printing characters for the phase chromosomes:
 'o' and 'i' are the phase unknown versions. */

   printf("\ngen_file = %s",gen_file);
   if (!(fp = fopen(gen_file, "r"))) {
     printf("\nERROR: No .gen file named %s\n",gen_file);
     exit(1);
   }
   data = (struct data *)our_alloc((ALLOC)sizeof(struct data));
   read_gen_file(fp, data, (char)0);
   make_theta_1_t();
   num_fams = chrom_data->num_fams;
   cross_mat = (struct cross ***)our_alloc
                ((ALLOC)(num_loc - 1)*sizeof(struct cross **));
   for (i = 0; i < num_loc - 1; i++){
      cross_mat[i] = (struct cross **)our_alloc
                ((ALLOC)(num_loc - (1 + i))*sizeof(struct cross *));      
      for (j = 0; j < num_loc - (1 +i); j++) cross_mat[i][j] = 0;
   }
    
   for (i_fam = 0; i_fam < num_fams; i_fam++){
      chrom_array = chrom_data->chrom_array[i_fam];
      num_chroms = chrom_data->num_chroms[i_fam];
      phase = chrom_data->phase_choices[i_fam];
      array = (char **)our_alloc((ALLOC)num_chroms * sizeof(char *));

      for (i_chrom = 0; i_chrom < num_chroms; i_chrom++) 
           array[i_chrom] = (char *)our_alloc((ALLOC)num_loc * sizeof(char));

      num_flanks = all_flanks[i_fam].num_flanks;
      if (num_flanks > 1){
          r_ratio[num_flanks - 1] = r_ratio_vec;
          r_ratio_ind = 1;
          r_ratio[num_flanks - 1][0] = 0;
          r_r[num_flanks - 1] = r_r_vec;
          r_r[num_flanks - 1][0] = 0;
          prod1[0] = 1.0;
      }
      for (i_fl = num_flanks - 2; i_fl >= 0; i_fl--){
          m = 1 << all_flanks[i_fam].num_tswitchs[i_fl];
          t_ratio = r_ratio[i_fl] = r_ratio_vec + r_ratio_ind;
          t_r = r_r[i_fl] = r_r_vec + r_ratio_ind;
          r_ratio_ind += m; 
          m_offset = all_flanks[i_fam].m_right_off[i_fl];
          i_int = sort_intervals[i_fam].list[i_fl];
          sr = i_int->r;
          t1t = theta_1_t->data[i_int->i][i_int->j][i_int->k - i_int->j -1];
	  sh_mask = all_flanks[i_fam].l_sh_mask[i_fl + 1];
	  id_mask = all_flanks[i_fam].l_id_mask[i_fl + 1];
	  sh_dir = all_flanks[i_fam].l_sh_dir[i_fl + 1];
	  rec_mask = all_flanks[i_fam].a_r_int[i_fl];
	  if (sh_dir) sh_mask >>= 1;
	  else sh_mask <<= 1;

          if (!m_offset) 
	      for (j = 0; j < m; j++){
		pj = ((j & sh_mask) >> 1) + (j & id_mask);
		prod[j] = (rj = sr != srec[j & rec_mask]) ? t1t * prod1[pj] : prod1[pj];
	      }
          else 
	    for (j = 0; j < m; j++){
	      if (sh_dir) pj = ((j & sh_mask) << 1) + (j & id_mask);
	      else pj = ((j & sh_mask) >> 1) + (j & id_mask);
	      term[0] = prod1[term_ind[0] = pj];
	      term[1] = prod1[term_ind[1] = pj + m_offset];
	      i_rec = sr == srec[j & rec_mask]; /* i.e. i_rec = !rec[j] */
	      term[i_rec] = t1t * term[i_rec];
	      i_rec = term[1] > term[0];   
	      prod[j] = term[i_rec];
	      if (term[i_rec]) t_ratio[j] = term[!i_rec] / term[i_rec]; 
	      else t_ratio[j] = 0;
	      t_r[j] = i_rec; 
	    }
         t_prod = prod1;
         prod1 = prod;
         prod = t_prod;  
      }

      for (i_chrom = 0; i_chrom < num_chroms; i_chrom++) 
	for (i_loc = 0; i_loc < num_loc; i_loc++){
	  if (chrom_array[i_chrom][order[i_loc]] == 'X')
	    array[i_chrom][i_loc] = 4;
	  else array[i_chrom][i_loc] = 
	    chrom_array[i_chrom][order[i_loc]] == '1';
	}

      j = 0;
      sup_rat = 0.0;
      for (i_fl = 0; i_fl < num_flanks - 1; i_fl++){
          m_offset = all_flanks[i_fam].m_right_off[i_fl];
          i_int = sort_intervals[i_fam].list[i_fl];
          sr = i_int->r;
	  sh_mask = all_flanks[i_fam].l_sh_mask[i_fl + 1];
	  id_mask = all_flanks[i_fam].l_id_mask[i_fl + 1];
	  sh_dir = all_flanks[i_fam].l_sh_dir[i_fl + 1];
	  rec_mask = all_flanks[i_fam].a_r_int[i_fl];
	  if (sh_dir) sh_mask >>= 1;
	  else sh_mask <<= 1;

          if (!m_offset) { 
	    rec_list[i_fl] = sr != srec[j & rec_mask];
	    j = ((j & sh_mask) >> 1) + (j & id_mask);
	  }
          else {
	    if (sh_dir) pj = ((j & sh_mask) << 1) + (j & id_mask);
	    else pj = ((j & sh_mask) >> 1) + (j & id_mask);
	    rj = sr != srec[j & rec_mask];
	    if (r_r[i_fl][j]) {
	      pj += m_offset;
	      rj = !rj;
	    }
	    rec_list[i_fl] = rj;
	    if (r_ratio[i_fl][j] > sup_rat) sup_rat = r_ratio[i_fl][j];
	    j = pj;
	  }
      }

      num_sw = phase->num_switches; 
      locus_nums = phase->locus_nums;
      sw_array = phase->array;
      plike = find_phase_change(array, num_chroms, order, 
			num_loc, all_intervals + i_fam, phase);
      if (num_flanks >= 2) 
	if (fam_likes[i_fam]) plike *= prod1[0] / fam_likes[i_fam];
      printf("\n\nFamily %s", famnum = chrom_data->fam_nums_array[i_fam]);

      for (i_loc = 0; i_loc < num_loc; i_loc++)
	for (i = 0; i < num_sw; i++) {
	  if (order[i_loc] < locus_nums[i]) break;
	  else if (order[i_loc] == locus_nums[i]) 
	    for (i_chrom = 0; i_chrom < num_chroms; i_chrom++) 
	      if (sw_array[i][i_chrom] == '1' && array[i_chrom][i_loc] < 2)
		array[i_chrom][i_loc] += 2;
	}

      printf("   phase likelihood = %.3f, 2d best = %.3f\n",plike, sup_rat * plike);
      for (i_chrom = 0; i_chrom < num_chroms; i_chrom++) {
          if (!(i_chrom%2)) {
	    j_count = 0;
	    for (j_mem = 0; j_mem < data->num_mems[i_fam]; j_mem++) 
	      if (data->ind[i_fam][j_mem]->moth_id) {
		j_count++;
		if (j_count == 1 + i_chrom/2) {
		  sibnum = data->ind[i_fam][j_mem]->id;
		  break;
		}
	      }
	      
	    printf ("\n\n%5ld ",sibnum);
	    parent = 'M';
          }
          else {
               printf ("\n      ");    
               parent = 'P';
          }
          for (i = 0; i < num_loc; i += 10) {
              for (j = i; j < (i + 10) && j < num_loc; j++)
                  printf ("%c", pchar[array[i_chrom][j]]);
              printf(" ");
          }

          num_cross = 0;
          e1 = 4;
          for (i = 0; i < num_loc; i++)
	    if ((e2 = array[i_chrom][i]) != 4){
	      if (e1 != 4 && (e1 & 1) != (e2 & 1)) num_cross++;
	      e1 = e2;
	    }

          printf ("   %d",num_cross);
          e1 = 4;
          prev_cross = 1;
          prev_name = 0;
          n_loc = 0;
          for (i = 0; i < num_loc; i++)
              if ((e2 = array[i_chrom][i]) != 4){
                  n_loc++;
                  if (e1 != 4) {
                       if ((e1 & 1) != (e2 & 1)){
                           if (prev_cross) {
                               if (!prev_name) {
                                  printf("\n         ");
                                  prev_name = 1;
                               }
                               printf("%d %s  ",i1 + 1,chrom_data->locus_names[order[i1]]);
                           }
                           prev_cross = 1;
                           cross = (struct cross *)our_alloc
                                 ((ALLOC)sizeof(struct cross)); 
                           cross->sibnum = sibnum;
                           cross->famnum = famnum;
                           cross->parent = parent;
                           cross->next_cross = cross_mat[i1][i - (i1 + 1)];
                           cross_mat[i1][i - (i1 + 1)] = cross;
                       }
                       else prev_cross = 0;
                  }
                  i1 = i;
                  e1 = e2;
              }
           if (prev_cross && n_loc > 1) {
               if (!prev_name) printf("\n         ");
               printf("%d %s  ",i1 + 1,chrom_data->locus_names[order[i1]]);
           }
     } 

     for (i_chrom = 0; i_chrom < num_chroms; i_chrom++) 
          our_free(array[i_chrom]);

     our_free(array);
  }
  printf("\n\n\n\nCROSSOVER CHROMOSOMES FOR EACH INFORMATIVE INTERVAL");
  for (i = 0; i < num_loc - 1; i++)
     for (j = 0; j < num_loc - (i + 1); j++)
         if (cross_mat[i][j]){
              i_cross = 0;
              printf("\n\n%d  %d", i + 1, i + j + 2);
              for (cross = cross_mat[i][j]; cross; 
                           i_cross++, cross = cross->next_cross){
                    if (!(i_cross%10)) printf("\n");
                    printf("%s-%ld-%c  ",cross->famnum, 
                            cross->sibnum, cross->parent);
              }
         }
  printf("\n\n\n\nCONSECUTIVE LOCI UNSEPARATED BY CROSSOVERS: \n");
  for (i = 0; i < num_loc - 1; i++){
     if (cross_mat[i][0]) continue;
     printf("\n%d ",i + 1);
     for (j = i + 1; j < num_loc; j++) {
         sep = 0;

         for (k = i; k < j; k++) 
             if (cross_mat[k][j - (k + 1)]) {
                  sep = 1;
                  break;
             }
         if (sep) break;
         printf("%d ",j + 1);
     } 
  }
}

double find_phase_change(chrom_array, num_chroms, order, num_loc, ival_list, phase)

/* finds the phase change corresponding to a given interval change for a family's chromosomes,
*/  
   char **chrom_array;  /* array of phase chromosomes, rearranged by order of loci */
   SHORT *order;        /* order of loci being tested */
   SHORT num_chroms, num_loc;    /* # chromosomes in family, # loci in order */
   struct interval_list *ival_list;
   struct phase *phase;
{
  SHORT i, num_sw, chrom_sw, i_chrom, i_loc, i1, nnsw, j, k, num_tri, m, i_fl;
  SHORT *locus_nums, *tri_sw, *list_sw, *n_locus_nums;
  SHORT **sub_array;
  char e1, e2, add_int_flag;
  char *use_sw, *t_sw, *app_int;
  char **array, **p_array;
  struct intervals *interval;
  char *our_alloc();
  SHORT our_free();
  double prod;
  
  num_sw = phase->num_switches;  /* # of tswitchs (for all loci) in this 
				    family */
  locus_nums = phase->locus_nums;  /* array of locus #s for tswitchs (assumed 
				      below to be in increasing order) */
  array = phase->array;
  n_locus_nums =  (SHORT *)our_alloc((ALLOC)num_sw * sizeof(SHORT));
  list_sw = (SHORT *)our_alloc((ALLOC)num_sw * sizeof(SHORT));
  nnsw = 0;
  for (i_loc = 0; i_loc < num_loc; i_loc++)
    for (i = 0; i < num_sw; i++){
      if (order[i_loc] < locus_nums[i]) break;
      else if (order[i_loc] == locus_nums[i]) {
	n_locus_nums[nnsw] = i_loc;
	list_sw[nnsw] = i;
	nnsw++;
      }
    }
  
  tri_sw = (SHORT *)our_alloc((ALLOC)nnsw * sizeof(SHORT));
  num_tri = 0;
  sub_array = (SHORT **)our_alloc((ALLOC)nnsw * sizeof(SHORT *));
  p_array = (char **)our_alloc((ALLOC)nnsw * sizeof(char *));
  use_sw = (char *)our_alloc((ALLOC)nnsw * sizeof(char));
  t_sw = (char *)our_alloc((ALLOC)nnsw * sizeof(char));
  for (i = 0; i < nnsw; i++){
    tri_sw[i] = i;
    sub_array[i] = (SHORT *)our_alloc((ALLOC)nnsw * sizeof(SHORT));
    p_array[i] = (char *)our_alloc((ALLOC)nnsw * sizeof(char));
    for (j = 0; j < nnsw; j++) {
      sub_array[i][j] = -1;
      p_array[i][j] = 0;
    }
    p_array[i][i] = 1;
    t_sw[i] = 0;
    use_sw[i] = 0;
  }
  app_int = (char *)our_alloc((ALLOC)nnsw * sizeof(char));
  prod = 1.0;

  for (i_chrom = 0; i_chrom < num_chroms; i_chrom++){
    for (i_loc = 0; i_loc < num_loc &&
	 (e1 = chrom_array[i_chrom] [i_loc]) == 4; i_loc++);
    for (i1 = i_loc++; i_loc < num_loc; i_loc++){
      if ((e2 = chrom_array[i_chrom][i_loc]) != 4){
	add_int_flag = 0;
	for (j = 0; j < nnsw; j++){
	  if ((n_locus_nums[j] == i1 
	       || n_locus_nums[j] == i_loc) 
	      && array[list_sw[j]][i_chrom] == '1'){
	    add_int_flag = 1;
	    app_int[j] = 1;
	  }
	  else app_int[j] = 0; 
	}
	if (add_int_flag){
	  for (i = 0; i < num_tri; i++)
	    if (app_int[tri_sw[i]])
	      for (j = 0; j < nnsw; j++) {
		if ((k = sub_array[i][j]) < 0) break;
		else app_int[k] = !app_int[k];
	      }
	  for (; i < nnsw; i++)
	    if (app_int[tri_sw[i]]){
	      j = tri_sw[num_tri];
	      tri_sw[num_tri] = tri_sw[i];
	      tri_sw[i] = j;
	      /* now need to find interval with the properties below (in ival_list)
		 look it up in flanklist, and see if it needs to be switched. */
	      for (interval = ival_list->first_interval ; interval; 
		   interval = interval->next_interval)
		if (interval->chrom_num == 
		    i_chrom && interval->j == i1) {
		  i_fl = interval->sort_index;
		  use_sw[tri_sw[num_tri]] = 
		    interval->r != rec_list[i_fl];
		  break;
		}
	      if (!interval) {
		use_sw[tri_sw[num_tri]] = e1 != e2;
/*
		prod *= 1.0 - theta->data[!SEX_EQ && i_chrom % 2][i1]
		  [i_loc - i1 - 1];
*/
	      }
	      k = 0;
	      for (j = 0; j < nnsw; j++)
		if (app_int[j] && j != tri_sw[num_tri]){
		  sub_array[num_tri][k++] = j;        
		  app_int[j] = 0;
		  for (m = 0; m < nnsw; m++) p_array[j][m] =
		    p_array[j][m] != p_array[tri_sw[num_tri]][m];
		}
	      num_tri++;
	      break;
	    }
	}
	i1 = i_loc;
	e1 = e2;
      }
    }
  } 

  for (i = 0; i < num_tri; i++) 
    if (use_sw[tri_sw[i]]) 
      for (j = 0; j < nnsw; j++)
	t_sw[j] = t_sw[j] != p_array[tri_sw[i]][j];
  /* in the following loop, perform the locus switching */       
  for (j = 0; j < nnsw; j++)
    if (t_sw[j]) { 
      /* apply corresponding locus switch to chrom_array */
      i_loc = n_locus_nums[j];
      for (i_chrom = 0; i_chrom < num_chroms; i_chrom++)
	if (array[list_sw[j]][i_chrom] == '1') {
	  e2 = chrom_array[i_chrom][i_loc];
	  if (e2 == 1) chrom_array[i_chrom][i_loc] = 0;
	  else if (e2 == 0) chrom_array[i_chrom][i_loc] = 1;
	}
    }

  our_free(tri_sw);
  for (i = 0; i < nnsw; i++) our_free(sub_array[i]);
  our_free(sub_array);
  our_free(app_int);
  for (i = 0; i < nnsw; i++) our_free(p_array[i]);
  our_free(p_array);
  our_free(list_sw);
  our_free(use_sw);
  our_free(t_sw);
  our_free(n_locus_nums);

  return (prod); /* prod gives correction factor for singleton intervals */
}
