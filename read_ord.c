#if vms
#include stdio
#include ctype
#else
#include <stdio.h>
#include <ctype.h>
#endif

#include "defs.h"

/* this routine reads in a previously computed orders list as a single */
/* orders object */

read_orders_file(fp,ad_orders_list)
     FILE *fp;
     struct loci_orders **ad_orders_list;
{
     SHORT j,k;
     SHORT num_objects;
     SHORT num_loci;
     LINDEX num_orders, i;
     struct loci_orders *orders_list, *orders, *last_orders;
     struct loci_orders *malloc_orders();

   fscanf(fp,"%hd",&num_objects);
   fscanf(fp,"%ld",&num_orders);
   fscanf(fp,"%hd",&num_loci);

   orders_list = malloc_orders(num_orders,num_loci);

   for(i= 0; i < num_orders; i++){
     	for(j = 0; j < num_loci; j++){
     		fscanf(fp,"%hd",&(orders_list->orders[i][j]));
     	}
   }

   last_orders = orders_list;
   for(k = 1; k < num_objects; k++){
     fscanf(fp,"%ld",&num_orders);
     fscanf(fp,"%hd",&num_loci);
     orders = malloc_orders(num_orders,num_loci);

     for(i= 0; i < num_orders; i++){
     	for(j = 0; j < num_loci; j++){
     		fscanf(fp,"%hd",&(orders->orders[i][j]));
     	}
     }

     last_orders->next_orders = orders;
     last_orders = orders;

   }
     
     *ad_orders_list = orders_list;

}

/* this routine writes an orders_list to a file */

write_orders_file(ord_file, orders_list)
     char *ord_file;
     struct loci_orders *orders_list;
{
     FILE *fp;
     FILE *fopen();
     LINDEX i;
     SHORT j;
     SHORT num_objects;
     struct loci_orders *orders;

     if (!strcmp(ord_file, "")) return;
     if(!(fp = fopen(ord_file, "w"))){
       printf("\nERROR: can't write to the .ord file %s\n",ord_file);
       exit(1);
     }
     orders = orders_list;
     num_objects = 1;

     while (orders = orders->next_orders) num_objects++;
     fprintf(fp,"%d\n\n\n",num_objects); 

     orders = orders_list;
     while(orders){
     	fprintf(fp, "%ld  %d\n",orders->num_orders, orders->num_loci);
     	for(i = 0; i < orders->num_orders; i++){
     	    for(j = 0; j < orders->num_loci; j++)
     		fprintf(fp,"%d  ",orders->orders[i][j]);
     	    fprintf(fp,"\n");
     	}
     	fprintf(fp,"\n\n\n");
     	orders = orders->next_orders;
     }
     fclose(fp);
}
