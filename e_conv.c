#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"

#define  MAX_NUM_SWITCHS 1000
#define  MAX_NUM_CHROMS 400
SHORT num_alleles = 0;

derive_gen(data, chrom_data, pk_chrom_data)
    struct data *data;
    struct chrom_data *chrom_data, *pk_chrom_data;
{
    SHORT all_alleles, t;
    SHORT c[2], d[2];
    SHORT *s, *m, *p, *inda;
    SHORT *par[2];
    SHORT mm, pp, flag, i_chr, index;
    SHORT i, j, k, i_loc, i_fam, i_mem, num_mems, i_sw, i_chrom, isw;
    SHORT r, rs, rm0, rm1, num_chroms;
    char **chrom_array, **pk_chrom_array, **ch_vec;
    struct individual *ind;
    struct phase *phase;
    char phase_code();
    SHORT i_pkchroms;
    char *sw_vec;
    char *our_alloc();
    SHORT our_free();
    SHORT *p_locus_nums;
    char **p_array;
    SHORT aa, bb;
    SHORT compare();
    SHORT find_common();
    SHORT convert();

    i_pkchroms = 0;
    pk_chrom_array = pk_chrom_data->chrom_array[0];
    p_locus_nums = (SHORT *)our_alloc((ALLOC)MAX_NUM_SWITCHS * sizeof(SHORT));
    p_array = (char **)our_alloc((ALLOC)MAX_NUM_SWITCHS * sizeof(char *));
    for (i_sw = 0; i_sw < MAX_NUM_SWITCHS; i_sw++)
        p_array[i_sw] = (char *)our_alloc((ALLOC)MAX_NUM_CHROMS * sizeof(char));  
    for (i_fam = 0; i_fam < data->num_fams; i_fam++){
        printf("\nfamily id %s",data->fam_id[i_fam]);
        chrom_array = chrom_data->chrom_array[i_fam];

        num_mems = data->num_mems[i_fam];
        num_chroms = chrom_data->num_chroms[i_fam];
        if (num_chroms > MAX_NUM_CHROMS) 
             printf("\n\nERROR: num_chroms = %d, exceeds max.\n",num_chroms);
        phase = chrom_data->phase_choices[i_fam]; 

        ch_vec = (char **)our_alloc((ALLOC)num_mems * sizeof(char *));
        for (i_mem = 0; i_mem < num_mems; i_mem++) {
             ch_vec[i_mem] = (char *)our_alloc((ALLOC)num_chroms * sizeof(char));
             for (j = 0; j < num_chroms; j++) ch_vec[i_mem][j] = '0';
        }

        i_chrom = 0; 
        for (i_mem = 0; i_mem < num_mems; i_mem++){
                  ind = data->ind[i_fam][i_mem];
                  if (!(ind->moth_id && ind->fath_id)) continue;  
	          ch_vec[ind->moth][i_chrom++] = '1' ;
	          ch_vec[ind->fath][i_chrom++] = '1' ;
/*
printf("\n%ld %ld %d %d ",ind->moth, ind->fath,i_chrom, i_mem);
*/
        }

        for (i_loc = 0; i_loc < data->num_loci; i_loc++){
/* convert alleles in two passes: first non-zero, then zero */    
           num_alleles = 0;
           for (i_mem = 0; i_mem < num_mems; i_mem++){
                inda = data->ind[i_fam][i_mem]->a[i_loc];
                for (k = 0; k < 2; k++) inda[k] = convert(inda[k]);
           }
           all_alleles = 0;
           for (j = 0; j <= num_alleles; j++) all_alleles += 1 << j;        
           for (i_mem = 0; i_mem < num_mems; i_mem++){
                inda = data->ind[i_fam][i_mem]->a[i_loc];
                for (k = 0; k < 2; k++) if (!inda[k]) inda[k] = all_alleles;
           }
           do {
	       flag = 0;
               for (i_mem = 0; i_mem < num_mems; i_mem++){
                  ind = data->ind[i_fam][i_mem];
                  if (!(ind->moth_id && ind->fath_id)) continue;  

	          s = ind->a[i_loc];
	          m = data->ind[i_fam][ind->moth]->a[i_loc];
	          p = data->ind[i_fam][ind->fath]->a[i_loc];

	          mm = m[0] | m[1];
	          pp = p[0] | p[1];
	          for (j = 0; j < 2; j++) d[j] = s[j] & (mm | pp);
	          flag += compare(s, d);
 
	          r = find_common(s, c, mm, pp);
                  flag += compare(s, c);

                  if ((!flag) && r == 2){
                     t = s[0];
                     s[0] = s[1];
                     s[1] = t;
                  }   
	          r = find_common(m, d, c[0], all_alleles);
                  if (!r) {
                       printf("\n NONINHERITANCE: family %s,  individ %ld, locus %d",
                       data->fam_id[i_fam], data->ind[i_fam][i_mem]->id, i_loc);
                  }
                  else flag += compare(m, d); 
	          r = find_common(p, d, c[1], all_alleles);
                  if (!r) {
                       printf("\n NONINHERITANCE: family %s,  individ %ld, locus %d",
                       data->fam_id[i_fam], data->ind[i_fam][i_mem]->id, i_loc);
                  }
                  else flag += compare(p, d);
	       }
	    } while (flag);

            i_chrom = 0;
            for (i_mem = 0; i_mem < num_mems; i_mem++){
                  ind = data->ind[i_fam][i_mem];
                  if (!(ind->moth_id && ind->fath_id)) continue;

	          s = ind->a[i_loc];
	          par[0] = m = data->ind[i_fam][ind->moth]->a[i_loc];
	          par[1] = p = data->ind[i_fam][ind->fath]->a[i_loc];
	          mm = m[0] | m[1];
	          pp = p[0] | p[1];
 
	          rs = find_common(s, c, mm, pp);

                  if (!rs) {
                       printf("\n NONINHERITANCE: family %s,  individ %ld, locus %d",
                       data->fam_id[i_fam], data->ind[i_fam][i_mem]->id, i_loc);
                       printf("\n rs = %d, s = %d, c = %d, mm = %d, pp = %d\n",
                         rs,s,c,mm,pp);
                  }
                  if (rs == 2) printf("\n PROGRAM ERROR: family %s, individ %ld, locus %d",
                       data->fam_id[i_fam], data->ind[i_fam][i_mem]->id, i_loc);
                  if ((rs == 3) && (s[0] != s[1])) isw = 1;
                  else isw = 0; 

                  for (k = 0; k < 2; k++){
  	             rm0 = find_common(par[k], d, s[k], all_alleles);
                     chrom_array[i_chrom + k][i_loc] = phase_code(rm0);
                     if ((rm0 < 3) && isw) {
        	          rm1 = find_common(par[k], d, s[!k], all_alleles);
                          if (rm1 == 3) chrom_array[i_chrom + k][i_loc] = 'X';
                     }
                     mm = par[!k][0] | par[!k][1];
/*
     printf("\n test: %d %d %d %d",par[0][0],par[0][1],par[1][0],par[1][1]);   

*/
                     if ((aa = !(par[k][0] & mm)) || (bb = !(par[k][1] & mm)))
                        pk_chrom_array[i_pkchroms + i_chrom + k][i_loc]
                           = chrom_array[i_chrom + k][i_loc] ;

                  }
                 i_chrom += 2;
	    }
            i_chrom = 0;
            for (i_mem = 0; i_mem < num_mems; i_mem++){
               ind = data->ind[i_fam][i_mem];
               if (phase->num_switches >= MAX_NUM_SWITCHS) 
                 printf("\n\nERROR: num_switches = %d, exceeds max.\n",
                       phase->num_switches);
               sw_vec = p_array[phase->num_switches];

               if (!(ind->moth_id && ind->fath_id)) {
                       copy_n (ch_vec[i_mem], sw_vec, num_chroms);
                       p_locus_nums[phase->num_switches] = i_loc; 
                       phase->num_switches += 1;
               }

               else {
	          s = ind->a[i_loc];
	          par[0] = m = data->ind[i_fam][ind->moth]->a[i_loc];
	          par[1] = p = data->ind[i_fam][ind->fath]->a[i_loc];
	          mm = m[0] | m[1];
	          pp = p[0] | p[1];
 
	          rs = find_common(s, c, mm, pp);

                  if ((rs == 3) && (s[0] != s[1])) {
                       copy_n (ch_vec[i_mem], sw_vec, num_chroms);
                       p_locus_nums[phase->num_switches] = i_loc; 
                       phase->num_switches += 1;
                       for (k = 0; k < 2; k++){
     	                  rm0 = find_common(par[k], d, s[k], all_alleles);
                          if (rm0 < 3) {
        	               rm1 = find_common(par[k], d, s[!k], all_alleles);
                               if (rm1 < 3 && rm0 != rm1)
                                       sw_vec[i_chrom + k] = '1';
                          }
                       }
                       i_chrom += 2;
                  }
                  else {
                       i_chrom += 2;
                       continue;
                  }
               }
               flag = 0;
               for (i_chr = 0; i_chr < num_chroms; i_chr++) 
                     if (sw_vec[i_chr] != '0' && 
                             chrom_array[i_chr][i_loc] != 'X'){
                          flag++;
                          index = i_chr;
                     }   
               if (flag == 1){
                     chrom_array[index][i_loc] = 'X';
                     pk_chrom_array[i_pkchroms + index][i_loc] = 'X';
                     phase->num_switches -= 1;
               } 
               else if (flag == 0) phase->num_switches -= 1;
/*
      	printf("\n%d \n",flag);
        for(i_chr=0;i_chr < num_chroms; i_chr++) printf("%c",sw_vec[i_chr]);
*/
              }
      } /* i_loc */        
      for (i = 0; i < phase->num_switches; i++) 
          for (j = 0; j < num_chroms; j++)
              if (p_array[i][j] == '1')
                   pk_chrom_array[i_pkchroms + j][p_locus_nums[i]] = 'X';
      i_pkchroms += num_chroms;
      for (i_mem = 0; i_mem < num_mems; i_mem++) our_free(ch_vec[i_mem]);
      our_free(ch_vec);
      phase->array = (char **)our_alloc((ALLOC)phase->num_switches * sizeof(char *));
      for (i = 0; i < phase->num_switches; i++){
           phase->array[i] = (char *)our_alloc((ALLOC)num_chroms * sizeof(char));
           for (j = 0; j < num_chroms; j++)
               phase->array[i][j] = p_array[i][j];
      } 
      phase->locus_nums = (SHORT *)our_alloc((ALLOC)phase->num_switches * sizeof(SHORT));
      for (i = 0; i < phase->num_switches; i++) 
           phase->locus_nums[i] = p_locus_nums[i];  
   } /* i_fam */
}

copy_n (s, t, n)
    char *s, *t;
    SHORT n;
{
    for (n--; n >= 0; n--) t[n] = s[n]; 
}

char phase_code(r)
    SHORT r;
{

    if (r == 1) return ('0'); 
    if (r == 2) return ('1');
    if (r == 3) return ('X');
    else return('X');
}
 
SHORT find_common(a, c, b0, b1)
    SHORT *a, *c;
    SHORT b0, b1;
{
    SHORT j, t0, t1;
    SHORT r;

    c[0] = c[1] = r = 0;
    for (j = 0; j < 2; j++) {
        if ((t0 = a[j] & b0) && (t1 = a[!j] & b1)){
            c[0] |= t0;
            c[1] |= t1;
            r += 1 <<j;
        }
    }
    return (r);
}

SHORT compare(a, c)
    SHORT *a, *c;
{
    SHORT anb[2], cnb[2];

    if ((a[0] == c[0]) && (a[1] == c[1])) return (0);
    num_bits(a, anb);
    num_bits(c, cnb);
    if (!cnb[0]) return (0); 
    if ((cnb[0] < anb[0]) || ((cnb[0] == anb[0]) && cnb[1] < anb[1])){
            a[0] = c[0];
            a[1] = c[1];
            return (1);
    }
    else return (0);
}

num_bits(a, nb)
/* nb (a 2-element array) on return gives number of bits in each element of a 
   with smaller number first */
    SHORT *a, *nb;
{
    SHORT i, j, t;
    for (i = 0; i < 2; i++){
         nb[i] = 0;
         for (j = 0; j < 15; j++) nb[i] += (1 << j) & a[i];
    }
    if (nb[0] > nb[1]) {
         t = nb[0];
         nb[0] = nb[1];
         nb[1] = t;
    }
}

SHORT convert(a)
    SHORT a;
{
    SHORT j;
    static SHORT alleles[16];

    if (!a) return (0);
    for (j = 0; j < num_alleles; j++) if (a == alleles[j]) return (1 << j);
    alleles[num_alleles++] = a;
    return (1 << (num_alleles - 1));
}



