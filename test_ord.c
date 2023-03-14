#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"

SHORT test_order(order,n_loci1,orderstruc)
/* tests an order against a linked series of loci_order structures */
    SHORT *order;
    SHORT n_loci1;
    struct loci_orders *orderstruc;

{

    SHORT i_1, i_2, i_loc1, i_loc2, n_loci2, order_found;
    LINDEX i_order, n_orders;
    SHORT **orders;

    for (;orderstruc;orderstruc=orderstruc->next_orders){
         order_found=0;
         n_orders=orderstruc->num_orders;
         n_loci2=orderstruc->num_loci;
         orders = orderstruc->orders;

     for (i_order=0; i_order<n_orders; i_order++){
            i_1 = i_2 = -1;
            for (i_loc1=0; i_loc1<n_loci1; i_loc1++){
                for (i_loc2=0; i_loc2<n_loci2; i_loc2++)
                    if (order[i_loc1]==orders[i_order][i_loc2]){
                        if((i_1>=0)&&(i_2>=0))
                            if((i_1<i_2)!=(i_2<i_loc2)){
                            goto nextorder;
                        }
                    i_1 = i_2;
                        i_2 = i_loc2;

                }
            }
            order_found = 1;  /* a compatible order has been found */
            break;

            nextorder: ;  /* current order is not compatible */
     }
     if(!order_found) return (0);  /* no order is compatible */
    }
    return (1); /* every structure in list has a compatible order */
}

compress_orderstruc(orderstruc,i)
/* deletes the i-th order from orderstruc */
    struct loci_orders *orderstruc;
    LINDEX i;
{
    LINDEX j;
    SHORT **orders;
    SHORT our_orders_free();

    orders= orderstruc->orders;
    our_orders_free (orders[i]);
    orderstruc->num_orders -= 1;
    for (j = i;j < orderstruc->num_orders;j++)
        orders[j] = orders[j+1];

}

append_orderstruc(new_orders,orderstruc)
/* adds an orderstruc (linked to 0) to a linked list of orderstrucs */
    struct loci_orders *new_orders;
    struct loci_orders **orderstruc;

{
    struct loci_orders *old_orders, *prev_orders, *temp_orders;
    SHORT included, in_order, i_loc, j_loc;
    SHORT free_orders();

    new_orders->next_orders = 0;
    prev_orders = new_orders;

    for(old_orders = *orderstruc;old_orders;old_orders=old_orders->next_orders){
                test_and_compress(old_orders,new_orders);

    }

    for(old_orders = *orderstruc;old_orders;){
/*    test inclusion;  */
          included = 1;
          for (i_loc=0; i_loc<old_orders->num_loci; i_loc++){
                in_order = 0;
                for(j_loc=0; j_loc<new_orders->num_loci; j_loc++)
                   if(old_orders->orders[0][i_loc]==
            new_orders->orders[0][j_loc]){
                       in_order=1;
                       break;
                   }
                if(!in_order){
                   included = 0;
                   break;
                }
          }
          if(!included){
                prev_orders->next_orders = old_orders;
                prev_orders = old_orders;
        old_orders = old_orders->next_orders;
          }

          else {
        temp_orders = old_orders;
        old_orders = old_orders->next_orders;
            free_orders(temp_orders);
          }
     }
     prev_orders->next_orders = 0;
     *orderstruc = new_orders;
}



test_and_compress (orders, orderstruc)
    struct loci_orders *orders, *orderstruc;
{
    SHORT **ordermat;
    LINDEX i_order;
    SHORT test_order();

    ordermat = orders -> orders;
    for (i_order=0;i_order<orders->num_orders;i_order++)
        if (!test_order(ordermat[i_order],orders->num_loci,orderstruc)){
            compress_orderstruc(orders,i_order);
            i_order--;
        }
}
