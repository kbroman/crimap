#if vms
#include stdio
#else
#include <stdio.h>
#include <math.h>
#endif

#include "defs.h"
#include "var1.h"

double get_log_like(recs_t, nrecs_t)
     struct loci_data *recs_t, *nrecs_t;
{
     SHORT i, j,k, n, num_types;
     double log10();
     double num_recs, num_nrecs, theta_val;
     double log_like;

     n = theta->n;
     num_types = theta->num_types;
     
     log_like = 0.0;
     for(k = 0; k < num_types; k++){
       for(i = 0; i < n; i++){
	 for(j = 0; j < n-i; j++){
	   num_recs = recs_t->data[k][i][j];
	   num_nrecs = nrecs_t->data[k][i][j];
	   theta_val = theta->data[k][i][j];
	   if(num_recs) {
	     if(theta_val)
	       log_like += num_recs*log10(theta_val);
/*
	     else {
	       printf("\n\nERROR: 0 likelihood. Check hap_sys0 systems for intralocus recombinants.\n");
	       exit();
	     }
*/
	   }
	   if(num_nrecs)
	     log_like += num_nrecs*log10(1-theta_val);
	   
	 }     	     		
       }
     }
     return(log_like);
}
