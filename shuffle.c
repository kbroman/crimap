#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"
#include "var1.h"
#include "var2.h"

#define NMAX 9

static struct tswitchs **s_list, **e_list;
static struct tswitchs *swit_list[100];
static SHORT *s_inds, *e_inds, *cvec;

find_score(iorder, num_flanks, sw_list, s_vec)
   SHORT *iorder, *s_vec;
   SHORT num_flanks;
   struct tswitch_list *sw_list;
{
   struct tswitchs *i_sw;
   SHORT i, min_index, max_index, index, jind;
   struct interval_ptrs *i_ptr;
   SHORT *e_vec;
   char *our_alloc();
   SHORT our_free();

   e_vec = (SHORT *)our_alloc((ALLOC)num_flanks * sizeof(SHORT));

   for (i = 0; i < num_flanks; i++) s_vec[i] = e_vec[i] = 0;

   for (i_sw = sw_list->first_tswitch; i_sw; i_sw = i_sw->next_tswitch){
             i_ptr = i_sw->first_ptr;
             min_index = max_index = iorder[i_ptr->interval->sort_index];
             for (i_ptr = i_ptr->next_ptr; i_ptr; i_ptr = i_ptr->next_ptr){
               index = iorder[i_ptr->interval->sort_index];
               if (index < min_index) min_index = index;
               else if (index > max_index) max_index = index;
             }
/*
             for (jind = max_index; jind > min_index; jind--) s_vec[jind] += 1;
*/
/*
             s_vec[max_index + 1] -= 1;
             s_vec[min_index + 1] += 1;
*/
             for (i = min_index + 1; s_vec[i]; i++);
             s_vec[i] = 1;
             for (i = max_index + 1; e_vec[i]; i--);
             e_vec[i] = 1;

   }
   jind = 0;
/*
   for (i = 0; i < num_flanks; i++) s_vec[i] = jind += s_vec[i];
*/
   for (i = 0; i < num_flanks; i++) s_vec[i] = jind += s_vec[i] - e_vec[i];
   our_free(e_vec);
}

quick_find_score(num_flanks)
   SHORT num_flanks;
{
   SHORT i, jind;

   jind = 0;
   cvec[0] = 0;
   for (i = 1; i < num_flanks; i++) cvec[i] = jind +=
              (s_list[i - 1] != 0) - (e_list[i - 1] != 0);
}

LINDEX compute_score(num_flanks)
   SHORT num_flanks;
{
   SHORT i;
   LINDEX x;

   quick_find_score(num_flanks);

   x = 0;
   for (i = 0; i < num_flanks; i++) {
             if (cvec[i] > 20){
                x = 2147483647;
                break;
             }
             x += 1 << cvec[i];
   }
   return (x);

}

change_se(start, end)
   SHORT start, end;
{
   SHORT i, nsw;
   char *our_alloc();
   SHORT our_free();
   nsw = 0;
   for (i = start; i < end; i++){
       if (s_list[i]){
             swit_list[nsw++] = s_list[i];
             s_list[i] = 0;
             e_list[e_inds[i]] = 0;
             s_inds[e_inds[i]] = 0;
             e_inds[i] = 0;
       }
       if (e_list[i]){
             swit_list[nsw++] = e_list[i];
             e_list[i] = 0;
             s_list[s_inds[i]] = 0;
             e_inds[s_inds[i]] = 0;
             s_inds[i] = 0;
       }

   }
   for (i = 0; i < nsw; i++) align_se(swit_list[i]);
}

char shuffle(sw_list, num_flanks, list)
   struct tswitch_list *sw_list;
   SHORT num_flanks;
   struct intervals **list;
{
   SHORT i, k, m, n, start, offset,  n_ints;
   SHORT perm[NMAX];
   struct intervals *i_int;
   struct intervals *temp[NMAX];
   SHORT *itemp, *ctemp;
   char *our_alloc();
   SHORT our_free();
   SHORT nmax, num_its;
   char flag, loc_flag;
   LINDEX cdiff;

   num_its = 0;
   n_ints = num_flanks - 1;
   nmax = (n_ints > NMAX) ? NMAX : n_ints - 1;
   itemp = (SHORT *)our_alloc((ALLOC)n_ints * sizeof(SHORT));
   ctemp = (SHORT *)our_alloc((ALLOC)num_flanks * sizeof(SHORT));
   quick_find_score(num_flanks);
   flag = 0;
restart:
   for (n = nmax; n >= 2; n--) {
           for (k = 0; k < n; k++) perm[k] = n - 1 - k;
           for (start = 0; start < n; start++){
               loc_flag = 0;
               for (i = 0; i < start; i++) itemp[i] = i;
               for (i = start; i <= n_ints - n; i += n)
                   for (k = 0; k < n; k++)
                       itemp[k + i] = i + perm[k];
               for (; i < n_ints; i++) itemp[i] = i;

               find_score(itemp, num_flanks, sw_list, ctemp);
               num_its++;
               for (i = start; i <=n_ints - n; i += n){
                   cdiff = 0;
                   for (k = i + 1; k < i + n; k++){
                       if (ctemp[k] > 25){
                           if (ctemp[k] < cvec[k]){
                               cdiff -= 1 << 25;
                           }
                           else if (ctemp[k] > cvec[k]) {
                               cdiff = 0;
                               break;
                           }
                       }
                       else cdiff += (1 << ctemp[k]) - (1 << cvec[k]);
                   }
                   if (cdiff < 0) {
                       loc_flag = 1;
/*
                       for (k = i + 1; k < i + n; k++) cvec[k] = ctemp[k];
*/
                       for (k = 0; k < n; k++){
                           temp[k] = i_int = list[k + i];
                           i_int->sort_index = i + perm[k];
                       }
                       for (k = 0; k < n; k++) list[i + perm[k]] = temp[k];
                       change_se(i, i + n);
                   }
               }
               if (loc_flag) {
                   flag = 1;
                   quick_find_score(num_flanks);
               }
/*
               if (flag) {
                   printf("\n %ld %d %d %d", ctot, n, start, offset);
                   goto restart;
               }
*/
          }
      }

 finish:

 freeup:
      our_free((char *)itemp);
      our_free((char *)ctemp);
      return (flag);
}


make_equiv_classes(sw_list, ival_list, sort_list)
     struct tswitch_list *sw_list;
     struct interval_list *ival_list;
     struct intervals **sort_list;
{
    SHORT num_sw, num_int, i, j, k, index, max_ind, temp;
    SHORT *ind_sw, *ind_int, *count;
    struct tswitchs *i_sw;
    struct interval_ptrs *i_int;
    struct intervals *i_ival;
    SHORT our_free();
    char *our_alloc();

    num_sw = sw_list->num_tswitchs;
    num_int = ival_list->num_intervals;
    ind_sw = (SHORT *)our_alloc((ALLOC)num_sw * sizeof(SHORT));
    ind_int = (SHORT *)our_alloc((ALLOC)num_int * sizeof(SHORT));
    count = (SHORT *)our_alloc((ALLOC)num_sw * sizeof(SHORT));

    for (i = 0; i < num_sw; i++) ind_sw[i] = i;
    for (i = 0; i < num_int; i++) ind_int[i] = -1;
    for (j = 0, i_sw = sw_list->first_tswitch; i_sw;
                i_sw = i_sw->next_tswitch, j++)
        for (i_int = i_sw->first_ptr; i_int; i_int = i_int->next_ptr){
            index = i_int->interval->sort_index;
            if ( (k = ind_int[index]) >= 0) {
                i = k;
                while (i != ind_sw[i]){
                   temp = i;
                   i = ind_sw[i];
                   ind_sw[temp] = j;
                }
                ind_sw[i] = j;
            }
            ind_int[index] = j;
        }
    for (j = 0; j < num_sw; j++){
            count[j] = 0;
            while (ind_sw[j] != ind_sw[ind_sw[j]])
                   ind_sw[j] = ind_sw[ind_sw[j]];
    }
    for (i = 0; i < num_int; i++) {
            ind_int[i] = ind_sw[ind_int[i]];
            count[ind_int[i]] += 1;
    }
    max_ind = 0;
    for (j = 0; j < num_sw; j++)
        if (count[j]) {
            temp = count[j];
            count[j] = max_ind;
            max_ind += temp;
        }
    for (i = 0; i < num_int; i++) {
        temp = ind_int[i];
        ind_int[i] = count[temp];
        count[temp] += 1;
    }
    for (i_ival = ival_list->first_interval; i_ival;
               i_ival = i_ival->next_interval) {
        i_ival->sort_index = ind_int[i_ival->sort_index];
        sort_list[i_ival->sort_index] = i_ival;
    }

    our_free(ind_int);
    our_free(count);
    our_free(ind_sw);
}

add_tswitchs(i_sw, j_sw)
/* "adds" switch i_sw to j_sw */
    struct tswitchs *i_sw, *j_sw;
{
    struct interval_ptrs *i_int, *j_int, *int_ptr;
    char inc_flag;
    struct interval_ptrs *append_interval_ptr();

    for (i_int = i_sw->first_ptr; i_int; i_int = i_int->next_ptr) {
       inc_flag = 1;
       for (j_int = j_sw->first_ptr; j_int; j_int = j_int->next_ptr)
           if (i_int->interval == j_int->interval){
              elim_interval_ptr(j_int, j_sw);
              inc_flag = 0;
              break;
           }
       if (inc_flag) {
           int_ptr = append_interval_ptr(j_sw);
           int_ptr->interval = i_int->interval;
       }
   }
}

find_start_end(i_sw, start, end)
     struct tswitchs *i_sw;
     SHORT *start, *end;
{
     SHORT max, min, index;
     struct interval_ptrs *i_int;

     i_int = i_sw->first_ptr;
     max = min = i_int->interval->sort_index;

     for (i_int = i_int->next_ptr; i_int; i_int = i_int->next_ptr) {
         index = i_int->interval->sort_index;
         if (max < index) max = index;
         else if (min > index) min = index;
     }
     *start = min;
     *end = max;
}

alloc_se_list(num_ints)
     SHORT num_ints;
{
     char *our_alloc();

     s_list = (struct tswitchs **)our_alloc((ALLOC)num_ints *
             sizeof(struct tswitchs *));
     e_list = (struct tswitchs **)our_alloc((ALLOC)num_ints *
             sizeof(struct tswitchs *));
     s_inds = (SHORT *)our_alloc((ALLOC)num_ints * sizeof(SHORT));
     e_inds = (SHORT *)our_alloc((ALLOC)num_ints * sizeof(SHORT));
     cvec = (SHORT *)our_alloc((ALLOC)(num_ints + 1) * sizeof(SHORT));

}

free_se_list()
{
     SHORT our_free();

     our_free(s_list);
     our_free(e_list);
     our_free(e_inds);
     our_free(s_inds);
     our_free(cvec);

}

find_basis(sw_list, num_ints)
     SHORT num_ints;
     struct tswitch_list *sw_list;
{
     SHORT i;
     struct tswitchs *i_sw;

     for (i = 0; i < num_ints; i++) s_list[i] = e_list[i] = 0;
     for (i_sw = sw_list->first_tswitch; i_sw; i_sw = i_sw->next_tswitch)
           align_se(i_sw);

}

align_se(i_sw)
     struct tswitchs *i_sw;
{
     SHORT start, end;
     struct tswitchs *j_sw;
     char reset;

           find_start_end(i_sw, &start, &end);
restart:
           reset = 0;
           if (j_sw = s_list[start]) {
               if (e_inds[start] > end) {
                   add_tswitchs(i_sw, j_sw);
                   s_list[start] = 0;
                   e_list[e_inds[start]] = 0;
                   s_inds[e_inds[start]] = 0;
                   e_inds[start] = 0;
                   align_se(j_sw);
               }
               else {
                   add_tswitchs(s_list[start], i_sw);
                   find_start_end(i_sw, &start, &end);
                   reset = 1;
               }
           }

           if (j_sw = e_list[end]) {
               if (s_inds[end] < start) {
                   add_tswitchs(i_sw, j_sw);
                   e_list[end] = 0;
                   s_list[s_inds[end]] = 0;
                   e_inds[s_inds[end]] = 0;
                   s_inds[end] = 0;
                   align_se(j_sw);
               }
               else {
                   add_tswitchs(e_list[end], i_sw);
                   find_start_end(i_sw, &start, &end);
                   reset = 1;
               }
           }

           if (!reset) {
               s_list[start] = e_list[end] = i_sw;
               e_inds[start] = end;
               s_inds[end] = start;
           }
           else goto restart;
}

compress_basis(sw_list, num_ints, inverse)
     struct tswitch_list *sw_list;
     SHORT num_ints;
     struct intervals **inverse;
{
     SHORT i, j, k, m, last_e, new_start, new_end, index;
     SHORT e_count, s_count, max, maxk, temp;
     SHORT count[4];
     struct intervals **temp_ival;
     char *our_alloc();
     char *ilist;
     SHORT our_free();
     struct intervals *ival;
     struct interval_ptrs *i_int;
     struct tswitchs *i_sw, *j_sw;

     ilist = (char *)our_alloc((ALLOC)num_ints * sizeof(char));
     temp_ival = (struct intervals **)our_alloc((ALLOC)num_ints *
            sizeof(struct intervals *));
     for (i = 0; i < num_ints; i++) ilist[i] = 0;

     for (i = num_ints - 1; i >= 0; i--)
         if (i_sw = s_list[i]) {
              e_count = s_count = 0;
              for (j = i; j < num_ints; j++) {
                  if (s_list[j]) {
                      if (e_count || s_count == 2) break;
                      if (s_count && e_list[j]) break;
                      if (s_count) j_sw = s_list[j];
                      ++s_count;
                  }
                  if (e_list[j]) {
                      if (s_count == 2) break;
                      if (2 <= ++e_count) {
                          j++;
                          break;
                      }
                  }
              }
              if (j <= i + 1) continue;
              for (i_int = i_sw->first_ptr; i_int; i_int = i_int->next_ptr) {
                  index = i_int->interval->sort_index;
                  if ((index >= i) && index < j) ilist[index] = 1;
              }
              if (s_count == 2) {
                 for (i_int = j_sw->first_ptr; i_int; i_int = i_int->next_ptr) {
                     index = i_int->interval->sort_index;
                     if ((index >= i) && index < j) ilist[index] += 2;
                 }
                 for (k = 0; k < 4; k++) count[k] = 0;
                 for (k = i; k < j; k++)  count[ilist[k]] += 1;
                 max = count[1];
                 maxk = 1;
                 for (k = 2; k < 4; k++)
                     if (count[k] > max) {
                          max = count[k];
                          maxk = k;
                     }
                 max += count[maxk] = i + count[0];
                 count[0] = i;
                 for (k = 1; k < 4; k++)
                     if (k != maxk) {
                          temp = count[k];
                          count[k] = max;
                          max += temp;
                     }
                 for (k = i; k < j; k++) {
                     temp = count[ilist[k]];
                     temp_ival[temp] = inverse[k];
                     inverse[k]->sort_index = temp;
                     count[ilist[k]] += 1;
                 }
                 for (k = i; k < j; k++) inverse[k] = temp_ival[k];
                 change_se(i, j);
              }
              else {
                 new_start = i;

                 for (k = i + 1; k < j; k++)
                     if (!ilist[k]) {
                         ival = inverse[k];
                         for (m = k - 1; m >= new_start; m--) {
                               inverse[m]->sort_index = m + 1;
                               inverse[m + 1] = inverse[m];
                         }
                         ival->sort_index = new_start;
                         inverse[new_start] = ival;
                         new_start++;
                     }
                 if (new_start > i) {
                     if (e_count) change_se(i, j);
                     else {
                         s_list[new_start] = s_list[i];
                         e_inds[new_start] = e_inds[i];
                         s_inds[e_inds[i]] = new_start;
                         s_list[i] = 0;
                         e_inds[i] = 0;
                     }
                 }
              }
              for (k = i; k < j; k++) ilist[k] = 0;
         }

     last_e = -1;
     for (i = 0; i < num_ints; i++)
         if (i_sw = e_list[i]) {
              e_count = s_count = 0;
              for (j = i; j >= 0; j--) {
                  if (e_list[j]) {
                      if (s_count || e_count == 2) break;
                      if (e_count && s_list[j]) break;
                      if (e_count) j_sw = e_list[j];
                      ++e_count;
                  }
                  if (s_list[j]) {
                      if (e_count == 2) break;
                      if (2 <= ++s_count) {
                          j--;
                          break;
                      }
                  }
              }
              if (j >= i - 1) continue;
              for (i_int = i_sw->first_ptr; i_int; i_int = i_int->next_ptr) {
                  index = i_int->interval->sort_index;
                  if ((index <= i) && index > j) ilist[index] = 1;
              }
              if (e_count == 2) {
                 for (i_int = j_sw->first_ptr; i_int; i_int = i_int->next_ptr) {
                     index = i_int->interval->sort_index;
                     if ((index <= i) && index > j) ilist[index] += 2;
                 }
                 for (k = 0; k < 4; k++) count[k] = 0;
                 for (k = i; k > j; k--)  count[ilist[k]] += 1;
                 max = count[1];
                 maxk = 1;
                 for (k = 2; k < 4; k++)
                     if (count[k] > max) {
                          max = count[k];
                          maxk = k;
                     }
                 count[maxk] = i - count[0];
                 max = count[maxk] - max;
                 count[0] = i;
                 for (k = 1; k < 4; k++)
                     if (k != maxk) {
                          temp = count[k];
                          count[k] = max;
                          max -= temp;
                     }
                 for (k = i; k > j; k--) {
                     temp = count[ilist[k]];
                     temp_ival[temp] = inverse[k];
                     inverse[k]->sort_index = temp;
                     count[ilist[k]] -= 1;
                 }
                 for (k = i; k > j; k--) inverse[k] = temp_ival[k];
                 change_se(j + 1, i + 1);
              }
              else {
                 new_end = i;

                 for (k = i - 1; k > j; k--)
                     if (!ilist[k]) {
                         ival = inverse[k];
                         for (m = k + 1; m <= new_end; m++) {
                               inverse[m]->sort_index = m - 1;
                               inverse[m - 1] = inverse[m];
                         }
                         ival->sort_index = new_end;
                         inverse[new_end] = ival;
                         new_end--;
                     }
                 last_e = new_end;
                 if (new_end < i) {
                     if (s_count) change_se(j + 1, i + 1);
                     else {
                         e_list[new_end] = e_list[i];
                         s_inds[new_end] = s_inds[i];
                         e_inds[s_inds[i]] = new_end;
                         e_list[i] = 0;
                         s_inds[i] = 0;
                     }
                 }
              }
              for (k = i; k > j; k--) ilist[k] = 0;
         }
     our_free(ilist);
     our_free(temp_ival);
}
