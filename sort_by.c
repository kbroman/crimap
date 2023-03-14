#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"

/* this routine sorts the loci in chrom_data by number of informative meioses */
/* the sorting algorithm is very inefficient */

SHORT *sort_by_info(chrom_data)
     struct chrom_data *chrom_data;
{
     SHORT num_loci;
     SHORT i,j;
     SHORT *info;
     char *our_alloc();
     SHORT info_max, index;
     SHORT *new_array;
     SHORT our_free();
     SHORT *count_meioses();

     num_loci = chrom_data->num_loci;
     new_array = (SHORT *)our_alloc((ALLOC)num_loci*sizeof(SHORT));
     info = count_meioses(chrom_data);

     for (i = 0; i < num_loci; i++){
     	index = 0;
     	info_max = -1;
     	for (j = 0; j < num_loci; j++){
     		if (info_max < info[j]){
     			info_max = info[j];
     			index = j;
     		}
     	}
     	info[index] = -1;
     	new_array[i] = index;
     }
     our_free(info);
     return(new_array);
}

SHORT *count_meioses(chrom_data)
     struct chrom_data *chrom_data;
{
     char **chrom_array;
     SHORT i, j, i_fam, num_loci;
     SHORT *info;
     char *our_alloc();
     SHORT our_free();

     num_loci = chrom_data->num_loci;
     info = (SHORT *)our_alloc((ALLOC)num_loci*sizeof(SHORT));
     for (i = 0; i < num_loci; i++)  info[i] = 0;

     for (i_fam = 0; i_fam < chrom_data->num_fams; i_fam++) {
       chrom_array = chrom_data->chrom_array[i_fam];
       for (i = 0; i < chrom_data->num_chroms[i_fam]; i++)
	 for (j = 0; j < num_loci; j++)
	   if (chrom_array[i][j] != 'X') info[j] += 1;
     }
     return(info);
}

