#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"
#include "var1.h"

find_recs_nrecs(chrom_data,loci_order, num_sub_loci)
     struct chrom_data *chrom_data;
     SHORT *loci_order;
     SHORT num_sub_loci;

{
     char e1, e2;
     char **chrom_ptr;
     SHORT i, j, k, gender, n, num_types;

     num_types = theta->num_types;
     n = num_sub_loci - 1;

     /* first initialize recs and nrecs */

     for (k = 0; k < num_types; k++)
       for (i = 0; i < n; i++)
     for (j = 0; j < n - i; j++){
       recs->data[k][i][j] = 0;
       nrecs->data[k][i][j] = 0;
     }

     chrom_ptr = chrom_data->chrom_array[0];

     for (k = 0, gender = 0; k < chrom_data->num_chroms[0]; k++){
       for (i = 0; i < num_sub_loci; i++){
     if ((e1 = chrom_ptr[k][loci_order[i]]) != 'X'){
       for (j = i + 1; j < num_sub_loci; j++){
         if ((e2 = chrom_ptr[k][loci_order[j]]) != 'X'){
           if (e1 != e2) recs->data[gender][i][j - i - 1] += 1;
           else nrecs->data[gender][i][j - i - 1] += 1;
           i = j;
           e1 = e2;
         }
       }              /* close j loop */
       break;            /* break out of i loop */
     }                 /* close if e1 = 'x' */
       }                    /* close i loop */
       gender = !(gender || SEX_EQ);
     }                   /* close k loop */

     /* set num_mei = recs + nrecs */

     for (k = 0; k < num_types; k++)
       for (i = 0; i < n; i++)
     for (j = 0; j < n - i; j++)
       num_mei->data[k][i][j] = recs->data[k][i][j] + nrecs->data[k][i][j];
     make_num_mei_split();
     return;
}
