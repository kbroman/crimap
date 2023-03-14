#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"

write_file(fp, data, chrom_data, pk_chrom_data)
     FILE *fp;
     struct data *data;
     struct chrom_data *chrom_data, *pk_chrom_data;
{
     SHORT i,j,k;
     SHORT m,n;
     struct phase *phase;

     fprintf(fp, "%d  ", data->num_loci);
     fprintf(fp, "1\n\n");
     fprintf(fp,"1\n\n");

     for(k = 0; k<data->num_loci; k++) fprintf(fp, "%s\n", data->locus_names[k]);
     fprintf(fp,"\n\n\n%d", pk_chrom_data->num_chroms[0]);
     for (j = 0; j < pk_chrom_data->num_chroms[0]; j++){
         fprintf(fp, "\n");
         for (k = 0; k < pk_chrom_data->num_loci; k++) 
               fprintf(fp, "%c",pk_chrom_data->chrom_array[0][j][k]);
     }
     fprintf(fp,"\n\n0\n\n\n");

     fprintf(fp, "%d\n", data->num_loci);
     fprintf(fp, "%d\n\n", data->num_fams);

     for(i = 0; i<data->num_fams; i++) fprintf(fp,"%s\n",data->fam_id[i]);
     fprintf(fp,"\n\n");
     for(k = 0; k<data->num_loci; k++) fprintf(fp, "%s\n", data->locus_names[k]);

     for(i=0; i<data->num_fams; i++){
      fprintf(fp,"\n\n\n%d", chrom_data->num_chroms[i]);
      for(j = 0; j < chrom_data->num_chroms[i]; j++){
        fprintf(fp,"\n");
        for(k = 0; k<chrom_data->num_loci; k++)
        	fprintf(fp, "%c",chrom_data->chrom_array[i][j][k]);
      }
      phase = chrom_data->phase_choices[i];
      fprintf(fp,"\n\n%d\n",phase->num_switches);
      for(m = 0; m < phase->num_switches; m++){
	fprintf(fp, "\n%d ",phase->locus_nums[m] + 1);
	for(n = 0; n < chrom_data->num_chroms[i]; n++)
  	   fprintf(fp, "%c",phase->array[m][n]);
      }

     } /* close i */

}
/*
write_data(fp, data)
     FILE *fp;
     struct data *data;
{
     SHORT i,j,k;

     fprintf(fp, "%d \n %d ", data->num_fams, data->num_loci);
     fprintf(fp, "1\n\n");
     fprintf(fp,"1\n\n");

     for(k = 0; k<data->num_loci; k++)
     	fprintf(fp, "%s\n", data->locus_names[k]);

     for (i = 0; i < data->num_fams; i++){
         fprintf(fp, "\n%s\n%d", data->fam_id[i], data->num_mems[i]);
         for (j = 0; j < data->num_mems[i]; j++){
             fprintf(fp, "\n%ld %ld %ld %d\n", data->ind[i][j]->id,
                data->ind[i][j]->moth_id, data->ind[i][j]->fath_id,
                      data->ind[i][j]->sex);
             for (k = 0; k < data->num_loci; k++) 
                fprintf(fp,"%d %d ", data->ind[i][j]->a[k][0],
                    data->ind[i][j]->a[k][1]);
         }
     }

}
*/
