#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"
#include "var1.h"
#include <math.h>

#define XRATIO 10
LINDEX *pos;
char *rec, *srec;
double *prod, *prod1, *r_prod_vec;
double **r_prod;

struct int_array *sort_intervals;
struct interval_list *all_intervals;
struct tswitch_list *all_tswitchs;
struct flank_list *all_flanks;


twopoint(ordered_loci, num_ordered, inserted_loci, num_inserted, chrom_data)
     SHORT num_ordered, num_inserted;
     SHORT *ordered_loci, *inserted_loci;
     struct chrom_data *chrom_data;
{
     SHORT pair[2];
     double thet[12];
     SHORT i_loc, j_loc, i_type, i, num_fams, num_tswitch_elim;
     SHORT num_its, num_types;
     double lods, base, log_like;  
     double super_like(), super_punk();
     double theta_copy[2];
     SHORT make_fl_sw_int();
     SHORT num_loci1, num_loci2;
     SHORT *loci_list1, *loci_list2;
     SHORT *insert_hap();
     SHORT n_loci;
     SHORT *new_order;
     double mle();
     SHORT dum_index;

     num_types = 2 - SEX_EQ;
     num_fams = chrom_data->num_fams;
     TOL = TOL/10.0;
     thet[0]=.001;
     thet[1]=.01;
     for (i = 1; i < 11; i++) thet[i + 1] = i * .05; 

     loci_list1 = num_ordered ? ordered_loci : inserted_loci;
     num_loci1 = num_ordered ? num_ordered : num_inserted;
     loci_list2 = num_inserted ? inserted_loci : ordered_loci;
     num_loci2 = num_inserted ? num_inserted : num_ordered;

     for (i_loc = 0; i_loc < num_loci1; i_loc++){
        pair[0] = loci_list1[i_loc];
	if (loci_list1 == loci_list2) num_loci2 = i_loc;
 	for (j_loc = 0; j_loc < num_loci2; j_loc++){
            pair[1] = loci_list2[j_loc];
	    new_order = insert_hap(pair, 2, &n_loci, &dum_index);
            log_like = mle(&num_its, new_order, n_loci, chrom_data, &num_tswitch_elim);
    	    for(i_type = 0; i_type < num_types; i_type++){
                 theta_copy[i_type] = theta->data[i_type][dum_index][0];
                 theta->data[i_type][dum_index][0] = .5;
		 fill_theta();
            }
            base = super_like(num_fams, num_tswitch_elim);   
            lods = log_like - base;
	    if (lods < PUK_LIKE_TOL && -lods < PUK_LIKE_TOL) continue;
	    /*second test is to allow for very small neg. lods when PUK_LIKE_TOL = 0.0*/
            printf("\n%s   %s  rec. fracs.= ",chrom_data->locus_names[pair[0]],
                  chrom_data->locus_names[pair[1]]);
    	    for(i_type = 0; i_type < num_types; i_type++)
                printf("  %5.2f", theta_copy[i_type]);
            printf(",   lods =  %6.2f", lods);

    	    for(i_type = 0; i_type < num_types; i_type++){
                printf("\n");
     		if (num_types > 1)theta->data[!i_type][dum_index][0] = 
                          theta_copy[!i_type];
                for (i = 0; i < 12; i++){
                    theta->data[i_type][dum_index][0] = thet[i];
		    fill_theta();
                    printf(" %6.2f",
                       super_like(num_fams, num_tswitch_elim) - base);  
                } 
            }
            printf("\n");
        }  /*  j_loc */
     } /*i_loc */

            
} 

/* this function computes the mle for phase unknown data */

double mle(ad_num_its, order, num_sub_loci, chrom_data, ad_num_tswitch_elim)
     SHORT *ad_num_its;
     SHORT *order;
     SHORT num_sub_loci;
     SHORT *ad_num_tswitch_elim;
     struct chrom_data *chrom_data;
{
     double log_like, prev_log_like; 
     SHORT num_its, pk_num_its;
     double super_punk(), change_theta();
     SHORT i, j, k;
     SHORT make_fl_sw_int();
     SHORT num_fams;
     SHORT num_tswitch_elim;

     if((num_fams = chrom_data->num_fams) == 1 && 
     		((chrom_data->phase_choices)[0])->num_switches == 0){
     	find_recs_nrecs(chrom_data,order,num_sub_loci);
     	return(change_theta(ad_num_its));
     }
     num_tswitch_elim = make_fl_sw_int(order, num_sub_loci, chrom_data);
     alloc_rec_pos(num_fams);
     num_its = 0;
     log_like = super_punk(num_fams, num_tswitch_elim);
/*
     printf("\nnum_mei\n");
     print_loci_data(num_mei);
     printf("\npk_recs\n");
     print_loci_data(pk_recs);
     printf("\nrecs\n");
     print_loci_data(recs);
     printf("\nlog10_like = %f\n",log_like);
*/  
     do {

       for(k = 0; k < theta->num_types; k++){
	 for(i = 0; i < theta->n; i++){
	   for(j = 0; j < theta->n-i; j++){
	     nrecs->data[k][i][j] = num_mei->data[k][i][j] - recs->
	       data[k][i][j];
     	        	}
               	}
       }

	change_theta(&pk_num_its);
/*
        printf("\n num_its %d , pk_num_its %d",num_its, pk_num_its);
     	printf("\ncurrent mle log10_like = %f\n",log_like);
*/
     	prev_log_like = log_like;
        log_like = super_punk(num_fams, num_tswitch_elim);
     	num_its++;
     } while( (log_like - prev_log_like) > TOL );
       
     if(log_like < prev_log_like - MAX_DECREASE){
     	printf("\n\nLOG_LIKE DECREASED IN MLE, DIFF = %f\n\n",
     		log_like - prev_log_like);
     }
     *ad_num_its = num_its;
     free_rec_pos();
     free_fl_sw_int(num_fams);
/*
    find_frags();
*/
     *ad_num_tswitch_elim = num_tswitch_elim;
     return(log_like);
}

free_rec_pos()
{
    SHORT our_free();

    our_free(pos);
    our_free(rec);
    our_free(srec);
    our_free(prod);
    our_free(prod1);
    our_free(r_prod);
    our_free(r_prod_vec);
}

alloc_rec_pos(num_fams)
   SHORT num_fams;
{
   char *our_alloc();
   SHORT i_fam, num_flanks, big_fam, i_fl, sup_nfl, sup_num_tswitchs;
   LINDEX n_rprod, m, sup_rprod;

   sup_num_tswitchs = sup_nfl = sup_rprod = 0;
   for (i_fam = 0; i_fam < num_fams; i_fam++){
      if ((num_flanks = all_flanks[i_fam].num_flanks) <= 1) continue;
      if (sup_nfl < num_flanks) sup_nfl = num_flanks;
      n_rprod = 0;

      for (i_fl = 0; i_fl < num_flanks; i_fl++){
         n_rprod += 1 <<  all_flanks[i_fam].num_tswitchs[i_fl];
         if (sup_num_tswitchs <  all_flanks[i_fam].num_tswitchs[i_fl])
             sup_num_tswitchs = all_flanks[i_fam].num_tswitchs[i_fl];
      }
      if (sup_rprod < n_rprod) {
             sup_rprod = n_rprod;
             big_fam = i_fam;
      }
   }
   m = 1 << sup_num_tswitchs;
/*
   printf("\n ALLOCATIONS: m = %ld sup_nfl = %d sup_rprod = %ld, i_fam = %d",
       m, sup_nfl, sup_rprod, big_fam);
*/
   pos = (LINDEX *)our_alloc((ALLOC)m * sizeof(LINDEX));
   rec = (char *)our_alloc((ALLOC)m * sizeof(char));
   srec = (char *)our_alloc((ALLOC)m * sizeof(char));
   make_srec(m);
   prod = (double *)our_alloc((ALLOC)m * sizeof(double));
   prod1 = (double *)our_alloc((ALLOC)m * sizeof(double));
/*
   offset = (LINDEX *)our_alloc((ALLOC)sup_m_offset * sizeof(LINDEX));
*/
   r_prod = (double **)our_alloc((ALLOC)sup_nfl * sizeof(double *));
   r_prod_vec = (double *)our_alloc((ALLOC)sup_rprod * sizeof(double));
}



SHORT make_fl_sw_int(order, num_sub_loci, chrom_data)
     SHORT *order;
     SHORT num_sub_loci;
     struct chrom_data *chrom_data;
{
     char *our_alloc();
     SHORT i_fam, num_fams, num_flanks;
     SHORT num_tswitch_elim;
     SHORT revise_tswitch_list();
     char shuffle();
     SHORT make_intervals_tswitchs();
     char flag;
     LINDEX score, impr_score, prev_score;
     LINDEX compute_score();

    num_fams = chrom_data->num_fams;
    clear(num_mei->data, num_sub_loci - 1, num_mei->num_types);
    clear(pk_recs->data, num_sub_loci - 1, num_mei->num_types);
    clear(pk_nrecs->data, num_sub_loci - 1, num_mei->num_types);

    all_tswitchs = (struct tswitch_list *)our_alloc((ALLOC)num_fams*sizeof(struct tswitch_list));
    all_intervals = (struct interval_list *)our_alloc((ALLOC)num_fams*sizeof(struct interval_list));
    all_flanks = (struct flank_list *)our_alloc((ALLOC)num_fams*sizeof(struct flank_list));
    sort_intervals = (struct int_array *)our_alloc((ALLOC)num_fams*sizeof(struct int_array));
    
    reset_int_block();
    reset_tsw_block();
    reset_int_ptr_block();
/*
    reset_fl_tsw_block();
*/
    num_tswitch_elim = 0;

    for (i_fam = 0; i_fam < num_fams; i_fam++){

       all_tswitchs[i_fam].num_tswitchs = 0;
       all_tswitchs[i_fam].first_tswitch = 0;
       all_intervals[i_fam].num_intervals = 0;
       all_intervals[i_fam].first_interval = 0;
       num_tswitch_elim += make_intervals_tswitchs(chrom_data->chrom_array[i_fam],
           chrom_data->num_chroms[i_fam], order, num_sub_loci, all_tswitchs + i_fam, all_intervals + i_fam, chrom_data->phase_choices[i_fam]);
/*
       num_tswitch_elim += revise_tswitch_list(all_tswitchs + i_fam, all_intervals + i_fam);
*/
       sort_intervals[i_fam].n = num_flanks = all_intervals[i_fam].num_intervals;
       sort_intervals[i_fam].list = (struct intervals **)our_alloc((ALLOC)
              num_flanks * sizeof(struct intervals *));
       all_flanks[i_fam].num_flanks = ++num_flanks;
     
       if (num_flanks > 1){
         sort_interval_list(num_sub_loci, sort_intervals[i_fam].list,
                      all_intervals + i_fam);
         make_equiv_classes(all_tswitchs + i_fam, all_intervals + i_fam,
               sort_intervals[i_fam].list);
         alloc_se_list(num_flanks - 1);
         find_basis(all_tswitchs + i_fam, num_flanks - 1);
         compress_basis(all_tswitchs + i_fam, num_flanks - 1, 
               sort_intervals[i_fam].list);
/* 
         printf("\nShuffling: i_fam = %d ", i_fam);
*/
         prev_score = impr_score = score = compute_score(num_flanks);
         while (impr_score/num_flanks > XRATIO) { 
            flag = shuffle(all_tswitchs + i_fam, num_flanks, 
                  sort_intervals[i_fam].list);
            if (flag) 
                compress_basis(all_tswitchs + i_fam, num_flanks - 1, 
                     sort_intervals[i_fam].list);
            score = compute_score(num_flanks);
            impr_score = prev_score - score;
            prev_score = score;
         }
         free_se_list();
         make_all_flanks(all_flanks + i_fam, all_tswitchs + i_fam);
      }
    }
    make_num_mei_split();
    return(num_tswitch_elim);
}

double super_punk(num_fams, num_tswitch_elim)
   SHORT num_fams, num_tswitch_elim;
{
   SHORT i_fam, i_fl, num_flanks;
   LINDEX m, m_offset, pj, j, r_prod_ind, m_add, jm, pjm;
   char sr;
   double like, log_like, tt, recomb;
   double log10(), get_log_like();
   struct intervals *i_int;
   double *t_prod, *t_prod1, *r_i_fl;
   double t1t;

   LINDEX sh_mask, id_mask, rec_mask;
   char sh_dir;

   log_like = num_tswitch_elim * log10(2.0);
   log_like += get_log_like(pk_recs, pk_nrecs);
/*
   printf("\nlog10_like = %f num_tswitch_elim = %d\n",log_like,num_tswitch_elim);
*/
   copy(pk_recs, recs);
   make_theta_1_t();
    
   for (i_fam = 0; i_fam < num_fams; i_fam++){
      if ((num_flanks = all_flanks[i_fam].num_flanks) <= 1) continue;
      r_prod[num_flanks - 1] = r_prod_vec;
      r_prod_ind = 1;
      r_prod[num_flanks - 1][0] = 1.0;
      for (i_fl = num_flanks - 2; i_fl >= 0; i_fl--){
          m = 1 << all_flanks[i_fam].num_tswitchs[i_fl];
          t_prod = r_prod[i_fl] = r_prod_vec + r_prod_ind;
          r_prod_ind += m; 
          t_prod1 = r_prod[i_fl + 1];
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

          if (m_offset) 
	    for (j = 0; j < m; j++){
	      if (sh_dir) pj = ((j & sh_mask) << 1) + (j & id_mask);
	      else pj = ((j & sh_mask) >> 1) + (j & id_mask);
	      t_prod[j] = (sr != srec[j & rec_mask]) ?
		t1t * t_prod1[pj] + t_prod1[pj + m_offset]
		  : t1t * t_prod1[pj + m_offset] + t_prod1[pj];
	    }
          else   /* sh_dir always 0 in this case */
	    for (j = 0; j < m; j++){
	      pj = ((j & sh_mask) >> 1) + (j & id_mask);
	      t_prod[j] = (sr != srec[j & rec_mask]) ?
		t1t * t_prod1[pj] : t_prod1[pj];
	    }
      }
      like = r_prod[0][0];

      prod1[0] = 1.0;
      for (i_fl = 1; i_fl < num_flanks; i_fl++){
          m = 1 << all_flanks[i_fam].num_tswitchs[i_fl];
          m_offset = all_flanks[i_fam].m_left_off[i_fl];
          i_int = sort_intervals[i_fam].list[i_fl - 1];
          sr = i_int->r;
          t1t = theta_1_t->data[i_int->i][i_int->j][i_int->k - i_int->j - 1];
          recomb = 0.0;
          r_i_fl = r_prod[i_fl];

	  sh_mask = all_flanks[i_fam].l_sh_mask[i_fl];
	  id_mask = all_flanks[i_fam].l_id_mask[i_fl];
	  sh_dir = all_flanks[i_fam].l_sh_dir[i_fl];
	  rec_mask =all_flanks[i_fam].a_l_int[i_fl];

          if (m_offset) 
	      for (j = 0; j < m; j++) {
		if (sh_dir) pj = ((j & sh_mask) >> 1) + (j & id_mask);
		else pj = ((j & sh_mask) << 1) + (j & id_mask);
		if (sr != srec[j & rec_mask]) {
		  tt = prod1[pj] * t1t;
		  prod[j] = tt + prod1[pj + m_offset];
		  recomb += tt * r_i_fl[j];
		}
		else {
		  tt = prod1[pj + m_offset] * t1t;
		  prod[j] = tt + prod1[pj];
		  recomb += tt * r_i_fl[j];
		}
	      }

	  else  /* sh_dir always 1 in this case */
	      for (j = 0; j < m; j++) {
		pj = ((j & sh_mask) >> 1) + (j & id_mask);

		if (sr != srec[j & rec_mask]) {
		  prod[j] = tt = t1t * prod1[pj];
		  recomb += tt * r_i_fl[j];
		}
		else prod[j] = prod1[pj];
	      }

          recs->data[i_int->i][i_int->j][i_int->k - i_int->j - 1] += recomb / like;
          t_prod = prod1;
          prod1 = prod;
          prod = t_prod;
      }
      if(like == 0.0) {
	printf("\n\nERROR: 0 likelihood. Check hap_sys0 systems for intralocus recombinants.\n");
	exit(1);
      }
      log_like += log10(like);
      fam_likes[i_fam] = like;
   }
   return(log_like);
}

double super_like(num_fams, num_tswitch_elim)
   SHORT num_fams, num_tswitch_elim;

/* this is the part of super_punk which computes the likelihood */
{
   SHORT i_fam, i_fl, num_flanks;
   LINDEX m, m_offset, pj, j;
   double log_like;
   double log10(), get_log_like();
   struct intervals *i_int;
   double *t_prod;
   double t1t;
   char sr;
   LINDEX sh_mask, id_mask, rec_mask;
   char sh_dir;

   log_like = num_tswitch_elim * log10(2.0);
   log_like += get_log_like(pk_recs, pk_nrecs);
   make_theta_1_t();
    
   for (i_fam = 0; i_fam < num_fams; i_fam++){
      if ((num_flanks = all_flanks[i_fam].num_flanks) <= 1) continue;
      prod1[0] = 1.0;

      for (i_fl = num_flanks - 2; i_fl >= 0; i_fl--){
          m = 1 << all_flanks[i_fam].num_tswitchs[i_fl];
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
		prod[j] = (sr != srec[j & rec_mask]) ? t1t * prod1[pj] : prod1[pj];
	      }
          else 
	      for (j = 0; j < m; j++){
		if (sh_dir) pj = ((j & sh_mask) << 1) + (j & id_mask);
		else pj = ((j & sh_mask) >> 1) + (j & id_mask);
		prod[j] = (sr != srec[j & rec_mask]) ?
		  t1t * prod1[pj] + prod1[pj + m_offset]
		    : t1t * prod1[pj + m_offset] + prod1[pj];
	      }
	  
          t_prod = prod1;
          prod1 = prod;
          prod = t_prod;
      }
      if(prod1[0] == 0.0) {
	printf("\n\nERROR: 0 likelihood. Check hap_sys0 systems for intralocus recombinants.\n");
	exit(1);
      }

      log_like += log10(prod1[0]);
   }
   return(log_like);
}


SHORT make_intervals_tswitchs(chrom_array,
      num_chroms, order, num_loc, sw_list, ival_list, phase)

/* makes the (unsorted) list of phase unknown intervals for a family's chromosomes,
  the tswitch list, and increments pk_recs, pk_nrecs, and num_mei for the phase known intervals. 
*/  
   char **chrom_array;  /* array of phase chromosomes */
   SHORT *order;        /* order of loci being tested */
   SHORT num_chroms, num_loc;    /* # chromosomes in family, # loci in order */
   struct tswitch_list *sw_list;
   struct interval_list *ival_list;
   struct phase *phase;
/* also uses global variables num_mei, pk_recs, pk_nrecs, SEX_EQ */

{
   SHORT i, num_sw, chrom_sw;
   SHORT i_chrom, sex, i_loc, i1, index, i_chrom1;
   char e1, e2, add_int_flag;
   struct tswitchs *i_sw;
   struct tswitchs *append_tswitch();
   struct intervals *interval;
   struct intervals *append_interval();
   struct interval_ptrs *int_ptr;
   struct interval_ptrs *append_interval_ptr();
   SHORT *locus_nums;
   char **array;
   char rank;
   SHORT nnsw, j, k, num_tri;
   char *app_int;
   SHORT *tri_sw;
   SHORT **sub_array;
   char *our_alloc();
   SHORT num_tswitch_elim;
   SHORT our_free();

/* initialize tswitch list */
   num_sw = phase->num_switches;  /* # of tswitchs (for all loci) in this family */
   locus_nums = phase->locus_nums;  /* array of locus #s for tswitchs (assumed below to be in increasing order */

   for (i_loc = 0; i_loc < num_loc; i_loc++){
      for (i = 0; i < num_sw; i++){
         if (order[i_loc] < locus_nums[i]) break;
         else if (order[i_loc] == locus_nums[i]){
            i_sw = append_tswitch(sw_list);
            i_sw->chrom_tswitch_index = i;
         }
      }
   }

   nnsw = sw_list->num_tswitchs;
   tri_sw = (SHORT *)our_alloc((ALLOC)nnsw * sizeof(SHORT));
   num_tri = 0;
   sub_array = (SHORT **)our_alloc((ALLOC)nnsw * sizeof(SHORT *));
   for (i = 0; i < nnsw; i++){
        tri_sw[i] = i;
        sub_array[i] = (SHORT *)our_alloc((ALLOC)nnsw * sizeof(SHORT));
        for (j = 0; j < nnsw; j++) sub_array[i][j] = -1;
   }
   app_int = (char *)our_alloc((ALLOC)nnsw * sizeof(char));
   array = phase->array;
   index = 0;
   for (i_chrom = 0, sex = 0; i_chrom < num_chroms; i_chrom++, sex = !(sex || SEX_EQ)){
       for (i_loc = 0; i_loc < num_loc && (e1 = chrom_array[i_chrom]
           [order[i_loc]]) == 'X'; i_loc++);
       for (i1 = i_loc++; i_loc < num_loc; i_loc++){
           if ((e2 = chrom_array[i_chrom][order[i_loc]]) != 'X'){
              add_int_flag = 0;
              rank = 5;
              for (j = 0, i_sw = sw_list->first_tswitch; i_sw; 
                       j++, i_sw = i_sw->next_tswitch){
                  chrom_sw = i_sw->chrom_tswitch_index;
                  if ((locus_nums[chrom_sw] == order[i1] || locus_nums[chrom_sw] == order[i_loc]) && array[chrom_sw][i_chrom] == '1'){
                    if (!add_int_flag){
                       interval = append_interval(ival_list);
                       interval->chrom_num = i_chrom;
                       interval->i = sex;
                       interval->j = i1;
                       interval->k = i_loc;
                       interval->r = e1 != e2;
                       interval->index = index++;
                       interval->sort_index = -1;
                       add_int_flag = 1;
                       if (i_chrom % 2) i_chrom1 = i_chrom - 1;
                       else i_chrom1 = i_chrom + 1;
                    }
/*
                    int_ptr = append_interval_ptr(i_sw);
                    int_ptr->interval = interval;
*/
                    app_int[j] = 1;
                    if (array[chrom_sw][i_chrom1] == '1'){
                       if (locus_nums[chrom_sw] == order[i_loc]) rank |= 8;
                       else rank &= 11;
                    } 
                    else{
                       if (locus_nums[chrom_sw] == order[i_loc]) rank |= 2;
                       else rank &= 14;
                    } 
                  }
                  else app_int[j] = 0; 
                
               }
               if (!add_int_flag){
                   if (e1 == e2) pk_nrecs->data[sex][i1][i_loc-i1-1] += 1;
                   else pk_recs->data[sex][i1][i_loc-i1-1] += 1;
                   num_mei->data[sex][i1][i_loc-i1-1] += 1;
               }
               else {
                  if ((rank & 12) != 4) rank = rank & 12;
                  interval->rank = rank;
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
                      k = 0;
                      for (j = 0; j < nnsw; j++)
                         if (app_int[j] && (j != tri_sw[num_tri])){
                            sub_array[num_tri][k++] = j;        
                            app_int[j] = 0;
                         }
                      num_tri++;
                      break;
                    }
                 for (j = 0, i_sw = sw_list->first_tswitch; i_sw; 
                       j++, i_sw = i_sw->next_tswitch)
                    if (app_int[j]){
                       int_ptr = append_interval_ptr(i_sw);
                       int_ptr->interval = interval;
                    }
               }
               i1 = i_loc;
               e1 = e2;
           }
       }
   } 
   num_tswitch_elim = 0;

   for (i_sw = sw_list->first_tswitch; i_sw; i_sw = i_sw->next_tswitch){

       if (!i_sw->num_intervals){     /* eliminate empty tswitch */
          elim_tswitch(i_sw, sw_list);
          num_tswitch_elim++;  /* use this to adjust likelihood later */
       }
       else
         if (i_sw->num_intervals == 1){
            elim_interval(i_sw->first_ptr->interval, ival_list);
            elim_tswitch(i_sw, sw_list);
         }  
   }   /* i_sw */

   our_free(tri_sw);
   for (i = 0; i < nnsw; i++) our_free(sub_array[i]);
   our_free(sub_array);
   our_free(app_int);
   return(num_tswitch_elim);
}

make_all_flanks(fl_list, sw_list)
    struct flank_list *fl_list;
    struct tswitch_list *sw_list;
{
    SHORT i, j, index, min_index, max_index, num_flanks;
    LINDEX m, m1;  
    struct tswitchs *i_sw;
    struct interval_ptrs *i_ptr;
    struct flank_tswitchs *i_fl_sw, *i1_fl_sw;
    char *our_alloc();
/*
    struct flank_tswitchs *append_flank_tswitch();
*/
    if ((num_flanks = fl_list->num_flanks) <= 1) return;

    fl_list->num_tswitchs = (SHORT *)our_alloc((ALLOC)num_flanks * sizeof(SHORT));
    fl_list->m_left_off = (LINDEX *)our_alloc((ALLOC)num_flanks * sizeof(LINDEX));
    fl_list->m_right_off = (LINDEX *)our_alloc((ALLOC)num_flanks * sizeof(LINDEX));
    fl_list->a_l_int = (LINDEX *)our_alloc((ALLOC)num_flanks * sizeof(LINDEX));
    fl_list->a_r_int = (LINDEX *)our_alloc((ALLOC)num_flanks * sizeof(LINDEX));
    fl_list->n_in_l_list = (SHORT *)our_alloc((ALLOC)num_flanks * sizeof(SHORT));
    fl_list->n_in_r_list = (SHORT *)our_alloc((ALLOC)num_flanks * sizeof(SHORT));
    fl_list->l_id_mask = (LINDEX *)our_alloc((ALLOC)num_flanks * sizeof(LINDEX));
    fl_list->l_sh_mask = (LINDEX *)our_alloc((ALLOC)num_flanks * sizeof(LINDEX));
    fl_list->l_sh_dir = (char *)our_alloc((ALLOC)num_flanks * sizeof(char));

    for (i = 0; i < num_flanks; i++) {
      fl_list->num_tswitchs[i] = fl_list->m_left_off[i] = fl_list->m_right_off[i] = 0; 
      fl_list->a_l_int[i] = fl_list->a_r_int[i] = 0; 
      fl_list->n_in_l_list[i] = fl_list->n_in_r_list[i] = 0; 
      fl_list->l_id_mask[i] = fl_list->l_sh_mask[i] = fl_list->l_sh_dir[i] = 0; 
    }

    for (i_sw = sw_list->first_tswitch; i_sw; i_sw = i_sw->next_tswitch){
        i_ptr = i_sw->first_ptr;
        min_index = max_index = i_ptr->interval->sort_index;
        for (i_ptr = i_ptr->next_ptr; i_ptr; i_ptr = i_ptr->next_ptr){
            index = i_ptr->interval->sort_index;
            if (index < min_index){
               for (j = min_index; j > index; j--) {
		 fl_list->num_tswitchs[j] += 1;
		 fl_list->a_l_int[j] <<= 1;
		 fl_list->a_r_int[j] <<= 1;
	       }
               fl_list->a_l_int[index+1] += 1;
               fl_list->a_r_int[min_index] += 1;
               min_index = index;
            }
            else if (index > max_index){
               for (j = index; j > max_index; j--) {
		 fl_list->num_tswitchs[j] += 1;
		 fl_list->a_l_int[j] <<= 1;
		 fl_list->a_r_int[j] <<= 1;
	       }
               fl_list->a_r_int[index] += 1;
               fl_list->a_l_int[max_index + 1] += 1;
               max_index = index;
            }
            else{
               fl_list->a_r_int[index] += 1;
               fl_list->a_l_int[index + 1] += 1;
            }
        }
        fl_list->n_in_l_list[min_index+1] = fl_list->num_tswitchs[min_index + 1];
        fl_list->n_in_r_list[max_index] = fl_list->num_tswitchs[max_index];
   }

/* note: we're assuming now that there is only one switch with a given start or  
   end */

   for (i = 1; i < num_flanks; i++){
     fl_list->n_in_l_list[i] = fl_list->num_tswitchs[i] - fl_list->n_in_l_list[i];
     if (fl_list->n_in_l_list[i] < fl_list->num_tswitchs[i]) 
       fl_list->m_right_off[i - 1] = 1 << fl_list->n_in_l_list[i];
     fl_list->n_in_r_list[i - 1] = fl_list->num_tswitchs[i - 1] - fl_list->n_in_r_list[i - 1];
     if (fl_list->n_in_r_list[i - 1]  < fl_list->num_tswitchs[i - 1]) 
       fl_list->m_left_off[i] = 1 << fl_list->n_in_r_list[i - 1];

     if (fl_list->n_in_l_list[i] < fl_list->n_in_r_list[i - 1]) {
       fl_list->l_id_mask[i] = ~(0L) << fl_list->n_in_r_list[i - 1] + 1 |
	 ~(~(0L) << fl_list->n_in_l_list[i]);
       fl_list->l_sh_dir[i] = 1;
     }
     else {
       fl_list->l_id_mask[i] = ~(0L) << fl_list->n_in_l_list[i] + 1|
	 ~(~(0L) << fl_list->n_in_r_list[i - 1]);
       fl_list->l_sh_dir[i] = 0;
     }

     fl_list->l_sh_mask[i] = ~fl_list->l_id_mask[i] -
       (1 << fl_list->n_in_l_list[i]);
   }
} 


sort_interval_list(num_loc, sort_list, ival_list)
/* following sorting loop may be modified later */
    SHORT num_loc;
    struct intervals **sort_list;
    struct interval_list *ival_list;
/* uses global variable num_mei */
{
   SHORT sort_index, i_loc;
   struct intervals *i_int;
   char rank;

   sort_index = 0;
   for (i_loc = 0; i_loc < num_loc - 1; i_loc++)
     for (rank = 0; rank < 16; rank++)
      for (i_int = ival_list->first_interval; i_int; i_int = i_int->next_interval)
          if ((i_int->j == i_loc) && (i_int->rank == rank)){
             i_int->sort_index = sort_index;
             sort_list[sort_index] = i_int;
             num_mei->data[i_int->i][i_int->j][i_int->k - i_int->j - 1] += 1;
             pk_nrecs->data[i_int->i][i_int->j][i_int->k - i_int->j - 1] += 1;
             sort_index++;
          }
}

/*
make_l_rec_pos(i_fl_sw)
       struct flank_tswitchs *i_fl_sw;
{
       LINDEX m, j, m_add, jm;
       char r;

       pos[0] = 0;	
       for (m = 1; i_fl_sw; i_fl_sw = i_fl_sw->next_tswitch, m <<= 1){
             if (i_fl_sw->in_l_list){
                 r = i_fl_sw->affects_l_interval;
                 m_add = i_fl_sw->m_in_l_list;
                 for (j = 0; j < m; j++){
                    
                     rec[jm = j + m] = r != rec[j];
                     pos[jm] = pos[j] + m_add;
                 }
             }
             else
                 for (j = 0; j < m; j++){
                     rec[jm = j + m] = !rec[j];
                     pos[jm] = pos[j];
                 }
       }
}

make_r_rec_pos(i_fl_sw)
       struct flank_tswitchs *i_fl_sw;
{
       LINDEX m, j, m_add, jm;
       char r;

       pos[0] = 0;	
       for (m = 1; i_fl_sw; i_fl_sw = i_fl_sw->next_tswitch, m <<= 1){
             if (i_fl_sw->in_r_list){
                 r = i_fl_sw->affects_r_interval;
                 m_add = i_fl_sw->m_in_r_list;
                 for (j = 0; j < m; j++){
                     rec[jm = j + m] = r != rec[j];
                     pos[jm] = pos[j] + m_add;
                 }
             }
             else
                 for (j = 0; j < m; j++){
                     rec[jm = j + m] = !rec[j];
                     pos[jm] = pos[j];
                 }

       }
}
*/

make_srec(tm)
     LINDEX tm;
{
  LINDEX m, j;
  
  srec[0] = 0;
  for (m = 1; m < tm; m <<= 1) 
    for (j = 0; j < m; j++)
      srec[j + m] = !srec[j];
}

free_fl_sw_int(num_fams)
    SHORT num_fams;
{
    SHORT i_fam;
    SHORT our_free();

    for (i_fam = num_fams - 1; i_fam >= 0; i_fam--){
       our_free (sort_intervals[i_fam].list);
       if (all_flanks[i_fam].num_flanks <= 1) continue;
       our_free(all_flanks[i_fam].num_tswitchs);
       our_free(all_flanks[i_fam].m_left_off);
       our_free(all_flanks[i_fam].m_right_off);
       our_free(all_flanks[i_fam].a_l_int); 	
       our_free(all_flanks[i_fam].n_in_l_list); 	
       our_free(all_flanks[i_fam].l_id_mask); 	
       our_free(all_flanks[i_fam].a_r_int); 	
       our_free(all_flanks[i_fam].n_in_r_list); 	
       our_free(all_flanks[i_fam].l_sh_mask); 	
       our_free(all_flanks[i_fam].l_sh_dir); 	
    }
    our_free(all_flanks);
    our_free(all_tswitchs);
    our_free(all_intervals);
    our_free(sort_intervals);
}

