#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"

/* this routines takes a loci_orders structure, an array of log_likes */
/* corresponding to each order, and a tolerence. It returns a loci_orders */
/* structure with the orders which have log_likes < (max_log_like-log_tol) */
/* deleted.								*/


screen_orders(orders,likelihoods,like_tol)
     struct loci_orders *orders;
     double *likelihoods;
     double like_tol;
{
     LINDEX i, j;
     double max_like;
     SHORT our_orders_free();

     max_like = likelihoods[0];
     for(i = 1; i < orders->num_orders; i++)
        if (likelihoods[i] > max_like) max_like = likelihoods[i];

     /* now select orders */
     for(i = 0, j = 0; i < orders->num_orders; i++){
        if( likelihoods[i] >= (max_like - like_tol)) {
            orders->orders[j] = orders->orders[i];
            j++;
        }
    else our_orders_free(orders->orders[i]);

     }

     orders->num_orders = j;

}
