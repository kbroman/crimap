#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"

struct loci_orders *make_loci_orders(num_loci,num_orders,orders)
     SHORT num_loci;
     LINDEX num_orders;
     SHORT **orders;
{
     struct loci_orders *result_orders;
     char *our_orders_alloc();

     result_orders = (struct loci_orders *)our_orders_alloc((ALLOC)sizeof(struct loci_orders));

     result_orders->num_loci = num_loci;
     result_orders->num_orders = num_orders;
     result_orders->orders = orders;
     result_orders->next_orders = 0;
     return(result_orders);
}

/* orderings of that list */
struct loci_orders *get_all_orders(orders,array,num_loci)
     struct loci_orders *orders;
     SHORT *array;
     SHORT num_loci;
{
     SHORT i;
     SHORT free_orders();

     struct loci_orders *insert();
     struct loci_orders *orders_temp;

     /* repeatedly call insert until a list of orders is built */
     for(i = 0; i<num_loci; i++){
        orders_temp = orders;
        orders = insert(orders_temp, array[i]);

        free_orders(orders_temp);
     }

     return(orders);
}

struct loci_orders *insert(orders,element)
     struct loci_orders *orders;
     SHORT element;
{
     SHORT **new_orders;
     LINDEX i, j;
     char *our_orders_alloc();
     SHORT **new_orders_ptr;

     new_orders = (SHORT **)our_orders_alloc
        ((ALLOC)(orders->num_orders)*(orders->num_loci +1)*sizeof(SHORT *));


     new_orders_ptr = new_orders;
     for(j = 0; j<orders->num_loci+1; j++){
         for(i = 0; i<orders->num_orders; i++){
        *new_orders_ptr =
            (SHORT *)our_orders_alloc((ALLOC) (orders->num_loci+1)
                            *sizeof(SHORT));

            put_index( orders->orders[i], *new_orders_ptr,
                            j,element,orders->num_loci);
            new_orders_ptr ++;
        }
     }

     return(make_loci_orders( (orders->num_loci+1),(orders->num_loci+1)*
            (orders->num_orders), new_orders) );

}



put_index(array1,array2,slot,element,n)
     SHORT *array1, *array2;
     SHORT slot;
     SHORT  element;
     SHORT n;
{
     SHORT i;

     for(i = 0; i < slot; i++) array2[i] = array1[i];

     array2[slot] = element;

     for(i = slot + 1; i < n + 1; i++) array2[i] = array1[i-1];
}
