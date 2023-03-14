#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"
#define MAX_ID_LEN 80
/* this routine reads the the data file given by fp */


read_file(fp,chrom_data)
     FILE *fp;
     struct chrom_data *chrom_data;
{
     SHORT i,j, num_loci, num_fams;
     SHORT *num_chroms;
     struct phase **phase_choices;
     char ***chrom_array;
     char **fam_nums_array;
     char **locus_names;
     char *our_alloc();
     char fam_ids[MAX_ID_LEN];
     char loc_name[MAX_LEN_LOC_NAME];

     fscanf(fp,"%hd",&num_loci);
     fscanf(fp,"%hd",&num_fams);
     fam_nums_array = (char **)our_alloc((ALLOC)num_fams*sizeof(char *));

     for (i = 0; i < num_fams; i++) {
       fscanf(fp,"%s",fam_ids);
       fam_nums_array[i] = (char *)our_alloc((ALLOC)(1 + strlen(fam_ids))
                             * sizeof(char));
       strcpy(fam_nums_array[i],fam_ids);
     }
     locus_names = (char **)our_alloc((ALLOC)num_loci*sizeof(char *));

     for(i = 0; i < num_loci; i++){
        fscanf(fp,"%s", loc_name);
    locus_names[i] = (char *)our_alloc((ALLOC)(1 + strlen(loc_name))
                         * sizeof(char));
    strcpy(locus_names[i], loc_name);
     }

     num_chroms = (SHORT *)our_alloc((ALLOC)num_fams * sizeof(SHORT));
     chrom_array = (char ***)our_alloc((ALLOC)num_fams*sizeof(char **));
     phase_choices = (struct phase **)our_alloc((ALLOC)num_fams*sizeof(struct phase *));

     for(i = 0; i<num_fams; i++){
        fscanf(fp, "%hd", num_chroms + i);
        chrom_array[i] = (char **)our_alloc((ALLOC)num_chroms[i]*sizeof(char *));

        for(j = 0; j<num_chroms[i]; j++){
        chrom_array[i][j] = (char *)our_alloc((ALLOC)(num_loci+1)*sizeof(char));
            fscanf(fp,"%s",chrom_array[i][j]);
        }
    phase_choices[i]=(struct phase *)our_alloc((ALLOC)num_fams*sizeof(struct phase));
        read_phase(fp,phase_choices[i], num_chroms[i]);
     }

     chrom_data->num_loci = num_loci;
     chrom_data->num_fams = num_fams;
     chrom_data->num_chroms = num_chroms;
     chrom_data->phase_choices = phase_choices;
     chrom_data->chrom_array = chrom_array;
     chrom_data->fam_nums_array = fam_nums_array;
     chrom_data->locus_names = locus_names;
}

read_phase(fp, phase, num_chroms)
     FILE *fp;
     struct phase *phase;
     SHORT num_chroms;
{
     SHORT i;
     char *our_alloc();
     SHORT num_switches;
     SHORT *locus_nums;
     char **array;

     fscanf(fp,"%hd",&num_switches);
     array = (char **)our_alloc((ALLOC)num_switches*sizeof(char *));
     locus_nums = (SHORT *)our_alloc((ALLOC)num_switches*sizeof(SHORT));

     for(i = 0; i<num_switches; i++){
        array[i] = (char *)our_alloc((ALLOC)(num_chroms+1)*sizeof(char));
        fscanf(fp,"%hd", &locus_nums[i]);
        locus_nums[i]--;
        fscanf(fp,"%s", array[i]);
     }

     phase->array = array;
     phase->num_switches = num_switches;
     phase->locus_nums = locus_nums;
}
