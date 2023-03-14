#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"
#include "var1.h"

/* this routine takes a structure of type loci_orders, a chrom_array,  */
/* num_choms, and num_loci                             */

double *get_likelihoods(mle_func, orders,chrom_data)
     double (*mle_func)();
     struct loci_orders *orders;
     struct chrom_data *chrom_data;
{
     LINDEX i;
     double *like_array;
     SHORT num_its, num_loci;
     char *our_alloc();
     SHORT *new_order;
     SHORT *insert_hap();
     SHORT our_orders_free();
     SHORT dum_index;
     SHORT num_tswitch_elim;

     like_array = (double *)our_alloc((ALLOC)orders->num_orders * sizeof(double));

     for(i = 0; i < orders->num_orders; i++){
        new_order = insert_hap(orders->orders[i], orders->num_loci, &num_loci, &dum_index);
        like_array[i] = (*mle_func)(&num_its,new_order,num_loci,chrom_data, &num_tswitch_elim);
        our_orders_free(new_order);
    }

     return(like_array);
}
