#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"
#include "var1.h"

/* these routines allocate memory and free memory */

malloc_global(num_types,n)
     SHORT num_types;
     SHORT n;
{
     SHORT i;

     struct loci_data *make_space_loci_data();
     char *our_alloc();

     recs = make_space_loci_data(num_types,n);
     nrecs = make_space_loci_data(num_types,n);
     pk_recs = make_space_loci_data(num_types,n);
     pk_nrecs = make_space_loci_data(num_types,n);
     theta = make_space_loci_data(num_types,n);
     theta_t = make_space_loci_data(num_types,n);
     theta_1_t = make_space_loci_data(num_types,n);
     recs_temp = make_space_loci_data(num_types,n);
     nrecs_temp = make_space_loci_data(num_types,n);
     num_mei = make_space_loci_data(num_types,n);
     num_mei_split = make_space_loci_data(num_types,n);
     FIXED_INTERVALS = (SHORT **)our_alloc((ALLOC)num_types
                       *sizeof(SHORT *));
     for (i = 0; i < num_types; i++)
       FIXED_INTERVALS[i] = (SHORT *)our_alloc((ALLOC)n*sizeof(SHORT));
}

/* alloc space for 4-D array */
/*
double ****make_array(d1,d2)
     SHORT d1;
     SHORT d2;
{
     SHORT i1,i2,i3;
     double ****ptr;
     char *our_alloc();

     ptr = (double ****)our_alloc((ALLOC)d1*sizeof(double ***));

     for(i1 = 0; i1<d1; i1++){
        ptr[i1] = (double ***)our_alloc((ALLOC)d2*sizeof(double **));
        for(i2 = 0; i2<d2; i2++){
            ptr[i1][i2]= (double **)our_alloc((ALLOC)(d2-i2)*sizeof(double *));
            for(i3 = 0; i3<(d2-i2); i3++){
                ptr[i1][i2][i3] = (double *)our_alloc((ALLOC)NUM_PRE_COMP*
                            sizeof(double));
            }
        }
     }

     return(ptr);
}
*/
/* our_alloc space for 3-D array */
double ***make_mat(d1,d2)
     SHORT d1;
     SHORT d2;
{
     SHORT i1,i2;
     double ***ptr;
     char *our_alloc();

     ptr = (double ***)our_alloc((ALLOC)d1*sizeof(double ***));

     for(i1 = 0; i1<d1; i1++){
        ptr[i1] = (double **)our_alloc((ALLOC)d2*sizeof(double *));
        for(i2 = 0; i2<d2; i2++){
            ptr[i1][i2]= (double *)our_alloc((ALLOC)(d2-i2)*sizeof(double));
        }
     }

     return(ptr);
}


struct loci_data *make_space_loci_data(num_types,n)
     SHORT num_types;
     SHORT n;
{
     char *our_alloc();
     struct loci_data *result;
     double ***make_mat();

     result = (struct loci_data *)our_alloc((ALLOC)sizeof(struct loci_data));

     result->num_types = num_types;
     result->n = n;
     result->data = make_mat(num_types,n);

     return(result);
}



free_global(num_types,n)
     SHORT num_types;
     SHORT n;
{
     SHORT flag, i;
     SHORT free_space_loci_data();
     SHORT our_free();

     flag = 0;

     flag += free_space_loci_data(recs,num_types,n);
     flag += free_space_loci_data(nrecs,num_types,n);
     flag += free_space_loci_data(pk_recs,num_types,n);
     flag += free_space_loci_data(pk_nrecs,num_types,n);
     flag += free_space_loci_data(theta, num_types,n);
     flag += free_space_loci_data(theta_t,num_types,n);
     flag += free_space_loci_data(recs_temp,num_types,n);
     flag += free_space_loci_data(nrecs_temp,num_types,n);
     flag += free_space_loci_data(num_mei,num_types,n);
     flag += free_space_loci_data(num_mei_split,num_types,n);
     for (i = 0; i < num_types; i++)
       flag += our_free(FIXED_INTERVALS[i]);
     flag += our_free(FIXED_INTERVALS);
     if(flag){
        printf("\n\n\n\nFAILED IN FREE_GLOBAL\n\n\n");
      }

}

SHORT free_space_loci_data(ptr, num_types,n)
     struct loci_data *ptr;
     SHORT num_types;
     SHORT n;
{
     SHORT flag;
     SHORT free_mat();
     SHORT our_free();

     flag = free_mat(ptr->data,num_types,n);
     flag += our_free(ptr);

     if(flag){
        printf("\n\n\n\nFREE_SPACE_LOCI_DATA FAILED\n\n");
        return(-1);
     }
     return(0);

}

/* free space for 3-D array */
SHORT free_mat(ptr,d1,d2)
     double ***ptr;
     SHORT d1;
     SHORT d2;
{
     SHORT i1,i2;
     SHORT flag;
     SHORT our_free();

     flag = 0;

     for(i1 = 0; i1<d1; i1++){
        for(i2 = 0; i2<d2; i2++){
            flag += our_free(ptr[i1][i2]);
        }
        flag += our_free(ptr[i1]);
     }
     flag += our_free(ptr);

     if(flag){
        printf("\n\n\n\nFREE FAILED IN FREE_MAT\n");
        return(-1);
     }

     return(0);
}
