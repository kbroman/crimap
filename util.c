#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"


clear(mat,n,num_types)
   double ***mat;
   SHORT n;
   SHORT num_types;
{

   SHORT i,j,k;

   for(k = 0; k < num_types; k++)
      for(i = 0; i < n; i++)
     for(j = 0; j < n-i; j++)
        mat[k][i][j] = 0;
}

copy(object1, object2)
     struct loci_data *object1, *object2;
{
     SHORT i,j,k;

     for(k = 0; k < object1->num_types; k++)
       for(i = 0; i < object1->n; i++)
         for(j = 0; j < object1->n - i; j++)
            object2->data[k][i][j] = object1->data[k][i][j];

     object2->n = object1->n;
     object2->num_types = object1->num_types;
}

copy_orders(orders1, orders2)
     struct loci_orders *orders1;
     struct loci_orders *orders2;
{
     LINDEX i,j;
     char *our_orders_alloc();

     orders2->num_loci = orders1->num_loci;
     orders2->num_orders = orders1->num_orders;
     orders2->orders = (SHORT **)our_orders_alloc((ALLOC)orders2->num_orders *
          sizeof(SHORT *));

     for(i = 0; i < orders2->num_orders; i++){
        orders2->orders[i] = (SHORT *)our_orders_alloc
                 ((ALLOC)orders2->num_loci * sizeof(SHORT));
        for(j = 0; j < orders2->num_loci; j++)
            orders2->orders[i][j] = orders1->orders[i][j];
     }
     orders2->next_orders = orders1->next_orders;
}
