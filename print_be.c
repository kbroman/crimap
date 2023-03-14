#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"

/* this routine takes a loci_orders structure, an array of log_likes */
/* corresponding to each order, and a tolerance. It prints  */
/* the orders which have log_likes >= (max_log_like-log_tol) */
/* in decreasing order of likelihood */
print_best_orders(orders,likelihoods,like_tol)
     struct loci_orders *orders;
     double *likelihoods;
     double like_tol;
{
     LINDEX i, n;
     SHORT j;
     double max_like;
     char *our_orders_alloc();
     LINDEX *v;

     n = orders->num_orders;
     v = (LINDEX *)our_orders_alloc((ALLOC)n * sizeof(LINDEX));
     short_sort_x(v, n, likelihoods);
     printf("\n");
     max_like = likelihoods[v[0]];
     for(i = 0; i < n; i++)
        if(likelihoods[v[i]] >= max_like - like_tol){
           for(j = 0; j < orders->num_loci; j++)
            printf("%3d ", orders->orders[v[i]][j] );
           printf("  %8.3f\n", likelihoods[v[i]]);
        }
     our_orders_free(v);
}

short_sort_x(v, n, x)
     LINDEX *v;
     LINDEX n;
     double *x;
/* Modified from Shell sort on p.108 of Kernighan-Ritchie */
{
  LINDEX gap,i,j;
  LINDEX temp;

  for (i = 0; i < n; i++) v[i] = i;
  for(gap = n/2; gap > 0; gap /= 2)
    for(i = gap; i < n; i++)
      for(j = i-gap; j >= 0; j -= gap){
    if (x[v[j]] >= x[v[j + gap]]) break;
    temp = v[j];
    v[j] = v[j + gap];
    v[j + gap] = temp;
      }
}
