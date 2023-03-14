#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"
#include "var1.h"

/* this routine prints the final results in the form of maps */

print_results1(orders,puk_chrom_data,loci_list,num_sub_loci,orders_list,like,sorted_list,sex)
     struct loci_orders *orders;
     struct loci_orders *orders_list;
     SHORT num_sub_loci;
     SHORT *loci_list;
     struct chrom_data *puk_chrom_data;
     SHORT like;
     SHORT *sorted_list;
     char sex;
{
     SHORT i, j, k;
     LINDEX l;
     struct loci_orders *orders_temp, *best_orders;
     struct loci_orders *insert();
     double kosambi();
     SHORT flag;
     SHORT current_num_loci;
     SHORT del_ind;	
     double *likelihoods;
     double *get_likelihoods();  
     double mle();
     SHORT *insert_hap();
     SHORT *new_order;
     SHORT n_loci;
     SHORT our_free();
     SHORT free_orders();
     SHORT dum_index;

     print_names(puk_chrom_data->locus_names,sorted_list,
        orders->num_loci + num_sub_loci,1,1);

     for (current_num_loci = num_sub_loci; current_num_loci
                  ; current_num_loci--){

     	for (i=0;  i < current_num_loci; i++){    	
     		orders_temp = insert(orders,loci_list[i]);
                test_and_compress(orders_temp,orders_list);
     		if( i==0 || orders_temp->num_orders < best_orders->num_orders){
     			best_orders = orders_temp;
     			del_ind = i;
                        if (best_orders->num_orders <= orders->num_orders)
                               break;
     		}
     	}
        if (best_orders->num_orders > 1)
                break; 
     	for(j = del_ind; j<current_num_loci-1; j++)
     		loci_list[j] = loci_list[j+1];
     	orders = best_orders;
     
     }    

     printf("\n\n\n\n\n");
     new_order = insert_hap(orders->orders[0],orders->num_loci,&n_loci, &dum_index);

     /* sex averaged */
     SEX_EQ = 1;
     likelihoods = get_likelihoods(mle,orders,puk_chrom_data);
     print_map(puk_chrom_data->locus_names,new_order,0);
     printf("\n\nlog10_like = %.2f\n\n",likelihoods[0]);
     our_free(likelihoods);
     printf("\n\n\n");

     /* sex specific */
     SEX_EQ = 0;
     likelihoods = get_likelihoods(mle,orders,puk_chrom_data);
     print_map(puk_chrom_data->locus_names,new_order,0);
     printf("\nlog10_like = %.2f\n\n",likelihoods[0]);
     our_free(likelihoods);
     printf("\n\n\n");

     	for(j = 0; current_num_loci-- ; j++){
     		printf("\n\n%s\n",(puk_chrom_data->locus_names)[loci_list[j]] );

		orders_temp = insert(orders,loci_list[j]);
     		test_and_compress(orders_temp,orders_list);

     		printf("    ");
     		print_loci_orders(orders);
     		for(k = 0; k<orders_temp->num_loci; k++){
     		   for(l = 0, flag=0; l<orders_temp->num_orders;l++){
     	              if( orders_temp->orders[l][k] == loci_list[j]){
     			 flag = 1;
     			 break;
     		      }
     		   }
 
    		   if(flag){
     		        printf("  X ");
     		   }
     		   else{
       		      	printf("    ");
     		   }
     		   
     		}

                if(like){
                  SEX_EQ = sex;
     		  likelihoods = get_likelihoods(mle,orders_temp,puk_chrom_data);
     	
     		  printf("\n\n");
     		  for(l = 0; l<orders_temp->num_orders; l++){
     			printf("%.2f\n",likelihoods[l]);
     		  }
                  our_free(likelihoods);
                }
       
     		free_orders(orders_temp);
     	}

     return;

}
