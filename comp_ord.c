#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"
#include "var1.h"

/* this function takes a list of loci, a chrom_array, num_loci,    */
/* num_chroms.    It returns an object of type loci_orders         */
/* in building the map it tries to optimize by always keeping the */
/* orders object with the smallest number of orders                */


struct loci_orders *comp_ords_build_map(ord_file,orders,pk_chrom_data, 
     		puk_chrom_data,	array_orig, num_sub_loci, ad_orders_list)

     char *ord_file;
     struct loci_orders *orders;
     struct chrom_data *pk_chrom_data, *puk_chrom_data;
     SHORT *array_orig;
     SHORT num_sub_loci;
     struct loci_orders **ad_orders_list;
{
     SHORT i, j;
     SHORT *array;
     char *our_alloc(); 
     SHORT our_free();
     char *our_orders_alloc();	
     SHORT free_orders();
     struct loci_orders *insert();
     struct loci_orders  *orders_orig,*orders_temp, *best_orders, *orders_list;
     SHORT current_num_loci, del_ind;
     LINDEX current_num_orders;
     double *get_likelihoods();
     double *likelihoods;
     double mle();

     array = (SHORT *)our_alloc((ALLOC)num_sub_loci * sizeof(SHORT));
     for(i = 0;i<num_sub_loci; i++) array[i] = array_orig[i];

     orders_orig= (struct loci_orders *)our_orders_alloc((ALLOC)sizeof(struct loci_orders));
     copy_orders(orders,orders_orig);

     orders_list = *ad_orders_list;
      best_orders = orders_orig;
     
     for (current_num_loci = num_sub_loci; current_num_loci && PK_NUM_ORDERS_TOL
                  ; current_num_loci--){

     	printf("\n\ncurrent orders\n");
     	print_loci_orders(orders);

     	for (i=0;  i<current_num_loci; i++){
     		orders_temp = insert(orders,array[i]);
                test_and_compress(orders_temp,orders_list);
		if(orders_temp->num_orders > 1){
         	   likelihoods = get_likelihoods(mle,orders_temp,pk_chrom_data);
        	   screen_orders(orders_temp,likelihoods,PK_LIKE_TOL);
                   our_free(likelihoods);
                   
		}
     		append_orderstruc(orders_temp,&orders_list);
                test_and_compress(orders,orders_list);
     		if( i==0 || orders_temp->num_orders < best_orders->num_orders){
     			best_orders = orders_temp;
     			del_ind = i;
                        if (best_orders->num_orders <= orders->num_orders)
                               break;
     		}
     	
     	}

	write_orders_file(ord_file, orders_list);     
        if (best_orders->num_orders > PK_NUM_ORDERS_TOL)
                break; 
     	for(j = del_ind; j < current_num_loci-1; j++)
     		array[j] = array[j+1];
     	
	free_orders(orders);

        orders= (struct loci_orders *)our_orders_alloc((ALLOC)sizeof(struct loci_orders));
	copy_orders(best_orders, orders);
/*
        find_frags();
        find_orders_frags();
*/

     }

     for(i = 0;i < num_sub_loci; i++) array[i] = array_orig[i];

     free_orders(orders);
     best_orders = orders = orders_orig;
     for (current_num_loci = num_sub_loci; current_num_loci
                  ; current_num_loci--){
     	printf("\n\ncurrent orders\n");
     	print_loci_orders(orders);

     	for (i=0;  i<current_num_loci; i++){
     		orders_temp = insert(orders,array[i]);
                test_and_compress(orders_temp,orders_list);
		
	     	printf("\n\norders_temp\n");
     		print_loci_orders(orders_temp);

     		if( i==0 || orders_temp->num_orders < best_orders->num_orders){
     			best_orders = orders_temp;
     			del_ind = i;
                        if (best_orders->num_orders <= orders->num_orders)
                               break;
     		}
     	
     	}
        if (best_orders->num_orders > 1) break; 
     	for(j = del_ind; j < current_num_loci - 1; j++) array[j] = array[j+1];
     	
     	orders = best_orders;
     
     }    

     best_orders= (struct loci_orders *)our_orders_alloc((ALLOC)sizeof(struct loci_orders));
     copy_orders(orders,best_orders);
     orders = best_orders;

     for (; current_num_loci; current_num_loci--){

     	printf("\n\ncurrent orders\n");
     	print_loci_orders(orders);

     	for (i=0;  i<current_num_loci; i++){    	
     		orders_temp = insert(orders,array[i]);
                test_and_compress(orders_temp,orders_list);
                current_num_orders = orders_temp->num_orders;

		if(orders_temp->num_orders > 1){
         	  likelihoods = get_likelihoods(mle,orders_temp,puk_chrom_data);
        	  screen_orders(orders_temp,likelihoods,PUK_LIKE_TOL);
                  our_free(likelihoods);
		}
     		append_orderstruc(orders_temp,&orders_list);
                test_and_compress(orders,orders_list);
     	     	printf("\n\norders_temp\n");
     		print_loci_orders(orders_temp);

     		if( i==0 || orders_temp->num_orders < best_orders->num_orders){
     			best_orders = orders_temp;
     			del_ind = i;
                        if (best_orders->num_orders <= orders->num_orders)
                               break;
          		}
/* commented out for run pk&puk tol=1.0 to speed things up
                if(current_num_orders != orders->num_orders) write_orders_file(ord_file,orders_list);     
*/
     	}

	write_orders_file(ord_file,orders_list);     

        if (best_orders->num_orders > PUK_NUM_ORDERS_TOL) break; 
     	for(j = del_ind; j<current_num_loci-1; j++) array[j] = array[j+1];
     	
	free_orders(orders);
        orders= (struct loci_orders *)our_orders_alloc((ALLOC)sizeof(struct loci_orders));
	copy_orders(best_orders, orders);
/*
        find_frags();
        find_orders_frags();     	  	
*/     
     }    
	*ad_orders_list = orders_list; 
     	return;
}
