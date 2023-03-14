#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"

print_loci_orders(orders)
     struct loci_orders *orders;
{
     LINDEX i,j;

     for(i = 0; i< orders->num_orders; i++){
        for(j = 0; j<orders->num_loci; j++){
            printf("%d", orders->orders[i][j] );
        if(orders->orders[i][j] > 9 ) printf("  ");
        else  printf("   ");
        }
        printf("\n");
     }
}


/* print_array is used in flips */

print_array(array,n)
     SHORT *array;
     SHORT n;
{
     SHORT i;

     for(i = 0; i<n; i++){
    printf("%d ", *(array+i));
     }
     printf("\n");

}

print_loci_data(object)
     struct loci_data *object;
{
     SHORT i,j,k;

     for(k = 0; k<object->num_types; k++){
       for(i = 0; i<object->n; i++){
         for(j = 0; j<(object->n-i); j++){
            printf("%.4f ", (object->data)[k][i][j]);
          }
          printf("\n");
       }
       printf("\n\n");
     }

}
