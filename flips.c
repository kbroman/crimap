#if vms
#include stdio
#else
#include <stdio.h>
#endif
  
#include "defs.h"
#include "var1.h"
  
flipsn(loci_list,num_loci,chrom_data,n_to_flip,use_ord_file,orders_list)
     SHORT *loci_list;
     SHORT num_loci,n_to_flip,use_ord_file;
     struct chrom_data *chrom_data;
     struct loci_orders *orders_list;
{
  SHORT **orders_array;
  SHORT start, end;
  LINDEX i, j, n_fac, k, m;
  struct loci_orders *orders, *orders_full, *orders_temp;
  struct loci_orders *make_loci_orders(),*get_all_orders();
  double mle(),kosambi();
  double *likelihoods;
  double like_orig;
  double *get_likelihoods();
  char *our_orders_alloc(),*our_alloc();
  SHORT our_free(), free_orders();
  
  for (i = 2, n_fac = 1; i <= n_to_flip; i++) n_fac *= i;
  orders_full = (struct loci_orders *)our_orders_alloc((ALLOC)sizeof
	       (struct loci_orders));
  orders_full->orders = (SHORT **)our_orders_alloc((ALLOC)
  	       (n_fac) * sizeof(SHORT *));
  for(i = 0; i < n_fac; i++){
    orders_full->orders[i] = (SHORT *)our_orders_alloc((ALLOC)
			       (num_loci) * sizeof(SHORT));
    for (j = 0; j < num_loci; j++) orders_full->orders[i][j] = loci_list[j];
  }
  orders_full->num_loci = num_loci;
  orders_full->next_orders = 0;
  
  for(start = 0, end = n_to_flip; end <= num_loci; start++, end++){
    orders_array = (SHORT **)our_orders_alloc((ALLOC)sizeof(SHORT *)); 
    *orders_array = (SHORT *)our_orders_alloc((ALLOC)sizeof(SHORT)); 
    orders_array[0][0] = loci_list[start];
    orders = make_loci_orders(1,(LINDEX)1,orders_array);
    orders = get_all_orders(orders,loci_list + start + 1, n_to_flip - 1);
    
    if (start) orders_full->num_orders = n_fac - n_fac/n_to_flip;
    else orders_full->num_orders = n_fac;
    
    for (i = 0; i < orders_full->num_orders; i++){
      for (j = 0; j < start; j++) orders_full->orders[i][j] = loci_list[j];
      for(k = 0; k < orders->num_loci; k++,j++)
	orders_full->orders[i][j] = orders->orders[i][k];
    }
    orders_temp = (struct loci_orders *)our_orders_alloc((ALLOC)
                    sizeof(struct loci_orders));
    copy_orders(orders_full, orders_temp);
    if (use_ord_file) test_and_compress(orders_temp, orders_list);
    likelihoods = get_likelihoods(mle, orders_temp, chrom_data);
    if(!start){
      printf("\n\nOriginal order, & its log10_likelihood,  followed by");
      printf("\nflipped orders, with their relative log10_likelihoods \n");
      printf("(= log10_like[orig] - log10_like[curr])\n\n");
      for (i = 0; i < num_loci; i++)  {
	m = orders_temp->orders[orders_temp->num_orders - 1][i];
	printf("%3d",m);
	if (m != loci_list[i]){
	  printf("\n\nERROR: input ordered_loci are incompatible with the orders database\n");
	  exit(1);
	}
      }
      like_orig = likelihoods[orders_temp->num_orders - 1];
      printf("%10.2f\n",like_orig);
      orders_temp->num_orders -= 1;
    }
    
    for(i = 0; i < orders_temp->num_orders; i++){
      if(n_to_flip == 2 || likelihoods[i] >= like_orig - PUK_LIKE_TOL){
	printf("\n");
	for (j = 0; j < num_loci; j++){
	  if (orders_temp->orders[i][j] == loci_list[j]) printf("  -");
	  else printf("%3d", orders_temp->orders[i][j]);
	}
	printf("%10.2f",like_orig - likelihoods[i]);
      }
    }
    our_free(likelihoods);
    free_orders(orders);
    free_orders(orders_temp);
  }
  printf("\n\n");
}
/*
  struct loci_orders *prepend(orders,start,loci_list)
  struct loci_orders *orders;
  SHORT *loci_list;
  SHORT start;
  {
  SHORT i,j,k;
  char *our_orders_alloc();
  struct loci_orders *new_list;
  
  new_list = (struct loci_orders *)our_orders_alloc((ALLOC)sizeof
  (struct loci_orders));
  new_list->orders = (SHORT **)our_orders_alloc((ALLOC)
  (orders->num_orders) * sizeof(SHORT *));
  for(i = 0; i < orders->num_orders; i++){
  new_list->orders[i] = (SHORT *)our_orders_alloc((ALLOC)
  (orders->num_loci + start) * sizeof(SHORT));
  for (j = 0; j < start; j++) new_list->orders[i][j] = loci_list[j];
  for(k = 0; k < orders->num_loci; k++,j++)
  new_list->orders[i][j] = orders->orders[i][k];
  }
  new_list->num_orders = orders->num_orders;
  new_list->num_loci = orders->num_loci + start;
  new_list->next_orders = 0;
  return (new_list);
  }
  
  struct loci_orders *append(orders,end,loci_list,n_loci)
  struct loci_orders *orders;
  SHORT *loci_list;
  SHORT end, n_loci;
  {
  SHORT i,j;
  char *our_orders_alloc();
  struct loci_orders *new_list;
  
  new_list = (struct loci_orders *)our_orders_alloc((ALLOC)sizeof 
  (struct loci_orders));
  new_list->orders = (SHORT **)our_orders_alloc((ALLOC)(orders->num_orders)*
  sizeof(SHORT *));
  for(i=0; i< orders->num_orders;i++){
  (new_list->orders)[i] = (SHORT *)our_orders_alloc((ALLOC)n_loci*sizeof(SHORT));
  for(j=0; j< orders->num_loci;j++)
  (new_list->orders)[i][j] = (orders->orders)[i][j];
  for(j=end;j< n_loci;j++)
  (new_list->orders)[i][j] = loci_list[j];
  }
  new_list->num_orders = orders->num_orders;
  new_list->num_loci = n_loci;
  new_list->next_orders = 0;
  return(new_list);
  }
  
  
  
  screen_and_print(orders,like_orig,likelihoods,like_tol)
  struct loci_orders *orders;
  double like_orig,like_tol;
  double *likelihoods;
  {
  SHORT i;
  for(i=0;i < orders->num_orders;i++){
  if((orders->num_orders == 1) || 
  likelihoods[i] >= (like_orig - like_tol)){
  printf("\n");
  print_array(orders->orders[i],orders->num_loci);
  printf("log10 of L[orig]/L[curr] = %.2f\n",like_orig - 
  likelihoods[i]);
  }
  }
  }
  */
