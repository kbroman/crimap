#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"
#include "var1.h" /*need FIXED_INTERVALS, theta, SEX_EQ */

struct haplos{
    SHORT num_in_hap;
    SHORT *loc_list;
    char fix_zero;
    struct haplos *next;
};

typedef struct haplos GLOBE_HAPLOS;
static GLOBE_HAPLOS *haplos;

struct fixed_dists{
    SHORT locus1;
    SHORT locus2;
    SHORT sex;
    double dist;
    struct fixed_dists *next;
};

typedef struct fixed_dists GLOBE_FIXED;
static GLOBE_FIXED *fixed_dists;

SHORT *insert_hap(order, num_loci, ad_num_loci, ad_index)
    SHORT *order;
    SHORT num_loci;
    SHORT *ad_num_loci, *ad_index;
{
    SHORT new_nloci, i, j, k, i_type, j_orig, jz;
    char *our_orders_alloc();
    SHORT *new_order;
    struct haplos *i_hap;
    struct fixed_dists *i_fixed;

    new_nloci = num_loci;
    if (haplos)
      for (i = 0; i < num_loci; i++)
    for (i_hap = haplos; i_hap; i_hap = i_hap->next)
      if (order[i] == i_hap->loc_list[0]){
        new_nloci += i_hap->num_in_hap - 1;
        break;
      }

    new_order = (SHORT *)our_orders_alloc((ALLOC)new_nloci * sizeof(SHORT));

    if (new_nloci != theta->n + 1 || 2 - SEX_EQ != theta->num_types) {
      free_global(theta->num_types, theta->n);
      malloc_global(2-SEX_EQ, new_nloci - 1);
    }

    for(i_type = 0; i_type < theta->num_types; i_type++)
      for(j = 0; j < theta->n; j++) FIXED_INTERVALS[i_type][j] = 0;

    for (i = j = 0; i < num_loci; i++){
      j_orig = j;
      new_order[j++] = order[i];
      for (i_hap = haplos; i_hap; i_hap = i_hap->next)
    if (order[i] == i_hap->loc_list[0]){
      for (k = 1; k < i_hap->num_in_hap; k++)
        new_order[j++] = i_hap->loc_list[k];
      for(i_type = 0; i_type < theta->num_types; i_type++)
        for(jz = j_orig; jz < j - 1; jz++) {
          if (i_hap->fix_zero) {
        theta->data[i_type][jz][0] = 0.0;
        FIXED_INTERVALS[i_type][jz] = 1;
          }
          else theta->data[i_type][jz][0] = 0.1;
        }
      break;
    }
      if (i < num_loci - 1) {
    for(i_type = 0; i_type < theta->num_types; i_type++)
      theta->data[i_type][j - 1][0] = 0.1;
    *ad_index = j - 1; /*only relevant for twopoint */
    for (i_fixed = fixed_dists; i_fixed; i_fixed = i_fixed->next)
      if ((order[i] == i_fixed->locus1 && order[i+1] == i_fixed->locus2
           || order[i] == i_fixed->locus2 && order[i+1] == i_fixed->locus1)
          && i_fixed->sex < theta->num_types) {
        theta->data[i_fixed->sex][j - 1][0] = i_fixed->dist;
        FIXED_INTERVALS[i_fixed->sex][j - 1] = 1;
      }
      }
    }
    fill_theta();

    *ad_num_loci = new_nloci;
    return (new_order);
}

SHORT delete_hap(order, num_loci)
    SHORT *order;
    SHORT num_loci;
{
    SHORT i, j, k;
    struct haplos *i_hap;

    if (!haplos) return (num_loci);
    for (i = j = 0; i < num_loci; i++) {
      for (i_hap = haplos; i_hap; i_hap = i_hap->next)
    for (k = 1; k < i_hap->num_in_hap; k++)
      if (order[i] == i_hap->loc_list[k]) goto nexti;
      order[j++] = order[i];
    nexti:;
    }
    return (j);
}

read_haps(hap_list, num_in_hap, fix_zero)
     SHORT *hap_list;
     SHORT num_in_hap;
     char fix_zero;
{
    SHORT j;
    char *our_orders_alloc();
    struct haplos *i_hap;

    i_hap = (struct haplos *)our_orders_alloc((ALLOC)sizeof(struct haplos));
    i_hap->next = haplos;
    haplos = i_hap;
    haplos->num_in_hap = num_in_hap;
    haplos->loc_list = (SHORT *)our_orders_alloc((ALLOC)num_in_hap * sizeof(SHORT));
    for (j = 0; j < num_in_hap; j++) haplos->loc_list[j] = hap_list[j];
    haplos->fix_zero = fix_zero;
}

read_fixed(fixed_list, fixed_dist, num_fixed)
     SHORT *fixed_list;
     SHORT num_fixed;
     double fixed_dist;
{
    struct fixed_dists *i_fixed;
    char flag;
    char *our_orders_alloc();

    flag = 0;
    do {
      i_fixed = (struct fixed_dists *)our_orders_alloc((ALLOC)sizeof(struct fixed_dists));
      i_fixed->next = fixed_dists;
      fixed_dists = i_fixed;

      fixed_dists->dist = fixed_dist;
      fixed_dists->locus1 = fixed_list[0];
      fixed_dists->locus2 = fixed_list[1];
      if (num_fixed >= 4) fixed_dists->sex = fixed_list[2];
      else fixed_dists->sex = flag = !flag;
    } while (flag);
}

set_null_hap()
{
    haplos = 0;
}

set_null_fixed()
{
    fixed_dists = 0;
}

print_haps(locus_names)
     char **locus_names;
{
  struct haplos *i_hap;

    if (!haplos) return;
    printf("\n\n");
    for (i_hap = haplos; i_hap; i_hap = i_hap->next) {
      printf("Haplotyped system (distances%sforced to 0.0):",
         i_hap->fix_zero ? " " : " not ");
      print_names(locus_names, i_hap->loc_list, i_hap->num_in_hap, 3, 1);
    }
    printf("N.B. Only the first locus in each set is retained in the orders");
    printf("\nobjects, but the remaining loci are used in all likelihood calculations.\n\n");

}

print_fixed()
{
  struct fixed_dists *i_fixed;

  if(!fixed_dists) return;
  printf("\n\nFixed recombination fractions (rec. frac., locus numbers, sex):");
  for (i_fixed = fixed_dists; i_fixed; i_fixed = i_fixed->next)
    printf("\n%.3f  %3d %3d  %3d", i_fixed->dist, i_fixed->locus1,
       i_fixed->locus2, i_fixed->sex);
  printf("\n\nN.B. These are fixed only when the loci in question are adjacent.\n");
}

write_haps(fp)
     FILE *fp;
{
  struct haplos *i_hap;
  SHORT i;

  if (!haplos) return;
  fprintf(fp,"\n\n");
  for (i_hap = haplos; i_hap; i_hap = i_hap->next) {
    fprintf(fp,"\nhap_sys%c ",i_hap->fix_zero ? '0' : ' ');
    for (i = 0; i< i_hap->num_in_hap; i++)
      fprintf(fp,"%d ", i_hap->loc_list[i]);
    fprintf(fp,"*");
  }
}

write_fixed(fp)
     FILE *fp;
{
  struct fixed_dists *i_fixed;

  if(!fixed_dists) return;
  for (i_fixed = fixed_dists; i_fixed; i_fixed = i_fixed->next)
    fprintf(fp, "\nfixed_dist %.3f  %3d %3d  %3d *",
        i_fixed->dist, i_fixed->locus1, i_fixed->locus2, i_fixed->sex);
}
