#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"


struct loci_orders *malloc_orders(num_orders, num_loci)
     LINDEX num_orders;
     SHORT num_loci;
{
     LINDEX i;
     char *our_orders_alloc();
     SHORT flag;
     struct loci_orders *ptr;

     flag = 0;

     ptr = (struct loci_orders *)our_orders_alloc((ALLOC)sizeof(struct loci_orders));
     if(!ptr) flag = 1;


     ptr->orders = (SHORT **)our_orders_alloc((ALLOC)num_orders*sizeof(SHORT *));
     if(!ptr->orders) flag = 1;

     for(i = 0; i < num_orders; i++){
        (ptr->orders)[i] = (SHORT *)our_orders_alloc((ALLOC)num_loci*sizeof(SHORT));
        if( !( (ptr->orders)[i] ) ) flag =1;
     }

     ptr->num_orders = num_orders;
     ptr->num_loci = num_loci;
     ptr->next_orders = 0;

     if(flag){
        printf("\n\n\nERROR IN MALLOC_ORDERS\n\n");
        return(0);
     }

     return(ptr);

}

SHORT free_orders(orders)
     struct loci_orders *orders;
{
     LINDEX i;
     SHORT flag;
     SHORT our_orders_free();

     flag = 0;
     for(i = 0; i<orders->num_orders; i++){
        if(orders->orders != NULL)
          flag += our_orders_free(orders->orders[i]);
    if(!flag) orders->orders[i] = NULL;
     }
     flag += our_orders_free(orders->orders);
     if(!flag) orders->orders = NULL;
     if(orders != NULL)
       flag += our_orders_free(orders);
       if(!flag) orders = NULL;

     if(flag){
        printf("\n\n\nFREE_ORDERS FAILED\n\n");
        return(-1);
     }

     return(0);
}
