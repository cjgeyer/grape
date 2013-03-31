
/*
*   input variables:
*       nnode  : scalar int, number of nodes in graph
*       nspac  : scalar int, spacing of iterates
*       nbatch : scalar int, number of batches (for batch means)
*       blen   : scalar int, length of batches (in iterations)
*       state  : int nnode by nnode matrix, 0-or-1-valued,
*                    1 in (i, j) cell indicates edge between nodes i and j
*                    must be symmetric (undirected graph) with 0's on diagonal
*                    only strict lower triangle is actually used
*       theta  : double vector of length nnode, the canonical parameter vector
*                    the corresponding statistics are the number of nodes
*                    of each degree, that is the log unnormalized density
*                    is sum(t * theta) where t[i] is the number of nodes
*                    of degree i (in C where we use 0-origin indexing) or
*                    of degree i - 1 (in R where we use 1-origin indexing)
*       batch  : double nbatch by nnode matrix, each row being the batch
*                    means of the canonical statistic vector for one batch 
*/

#include <R.h>

void grape(int *nnode_in, int *nspac_in, int *nbatch_in, int *blen_in,
    int *state, double *theta, double *batch)
{
    int nnode = nnode_in[0];
    int nspac = nspac_in[0];
    int nbatch = nbatch_in[0];
    int blen = blen_in[0];

    int *d = (int *) R_alloc(nnode, sizeof(int));
    int *dpro = (int *) R_alloc(nnode, sizeof(int));
    int *t = (int *) R_alloc(nnode, sizeof(int));
    int *tpro = (int *) R_alloc(nnode, sizeof(int));
    int *tbat = (int *) R_alloc(nnode, sizeof(int));

    for (int i = 0; i < nnode; ++i)
        d[i] = 0;

    for (int icol = 0; icol < nnode; ++icol)
        for (int irow = icol + 1; irow < nnode; ++irow) {
            int tmp = state[irow + nnode * icol];
            if (! (tmp == 0 || tmp == 1))
                error("bad initial value of state: upper tri not 0 or 1\n");
            if (tmp == 1) {
                d[irow]++;
                d[icol]++;
            }
        }

#ifdef CHECK_IT
    for (int i = 0; i < nnode; ++i)
        if (d[i] < 0 || d[i] >= nnode)
            error("initial check of d failed\n");
#endif /* CHECKIT */

    for (int i = 0; i < nnode; ++i)
        t[i] = 0;
    for (int i = 0; i < nnode; ++i)
        t[d[i]]++;

#ifdef BLATHER
    printf("initial :\n");
    printf("    d :");
    for (int i = 0; i < nnode; ++i)
        printf(" %d", d[i]);
    printf("\n");
    printf("    t :");
    for (int i = 0; i < nnode; ++i)
        printf(" %d", t[i]);
    printf("\n");
#endif /* BLATHER */

#ifdef CHECK_IT
    {
        for (int i = 0; i < nnode; ++i)
            if (t[i] < 0)
                error("initial check of sign(t) failed\n");
        int sum = 0;
        for (int i = 0; i < nnode; ++i)
            sum += t[i];
        if (sum != nnode)
            error("initial check of sum(t) failed\n");
    }
#endif /* CHECKIT */

    for (int i = 0; i < nnode; ++i) {
        dpro[i] = d[i];
        tpro[i] = t[i];
    }

    GetRNGstate();

    for (int ibatch = 0, kbatch = 0; ibatch < nbatch; ++ibatch) {

        for (int i = 0; i < nnode; ++i)
            tbat[i] = 0;

        for (int jbatch = 0; jbatch < blen; ++jbatch) {

        for (int kbatch = 0; kbatch < nspac; ++kbatch) {

            int i, j;
            do {
                i = nnode * unif_rand();
                j = nnode * unif_rand();
            } while (i == j);
            if (i < j) {
                int tmp = i;
                i = j;
                j = tmp;
            }
            /* propose to flip (i,j)-th edge */
            int state_ij = state[i + nnode * j];

            /* loop invariant: tpro[i] = t[i] & dpro[i] = d[i] for all i */

#ifdef BLATHER
    printf("iteration %d of batch %d, attempt to flip (%d, %d) edge\n",
        jbatch + 1, ibatch + 1, i + 1, j + 1);
    printf("    d :");
    for (int i = 0; i < nnode; ++i)
        printf(" %d", d[i]);
    printf("\n");
    printf("    t :");
    for (int i = 0; i < nnode; ++i)
        printf(" %d", t[i]);
    printf("\n");
#endif /* BLATHER */

            if (state_ij == 1) {
                dpro[i]--;
                dpro[j]--;
            } else {
                dpro[i]++;
                dpro[j]++;
            }

#ifdef BLATHER
    printf("    dpro :");
    for (int i = 0; i < nnode; ++i)
        printf(" %d", dpro[i]);
    printf("\n");
#endif /* BLATHER */
#ifdef CHECK_IT
    for (int i = 0; i < nnode; ++i)
        if (dpro[i] < 0 || dpro[i] >= nnode)
            error("check of dpro in %d-th update of %d-th batch failed\n", jbatch, ibatch);
#endif /* CHECKIT */

            tpro[d[i]]--;
            tpro[d[j]]--;
            tpro[dpro[i]]++;
            tpro[dpro[j]]++;

#ifdef BLATHER
    printf("    tpro :");
    for (int i = 0; i < nnode; ++i)
        printf(" %d", tpro[i]);
    printf("\n");
#endif /* BLATHER */

#ifdef CHECK_IT
    {
        for (int i = 0; i < nnode; ++i)
            if (tpro[i] < 0)
                error("check of sign(tpro) in %d-th update of %d-th batch failed\n", jbatch, ibatch);
        int sum = 0;
        for (int i = 0; i < nnode; ++i)
            sum += tpro[i];
        if (sum != nnode)
            error("check of sum(tpro) in %d-th update of %d-th batch failed\n", jbatch, ibatch);
    }
#endif /* CHECKIT */

            double log_odds = theta[i] * (tpro[i] - t[i])
                + theta[j] * (tpro[j] - t[j]);

            int accept = FALSE;
            if (log_odds >= 0.0 || unif_rand() < exp(log_odds))
                accept = TRUE;

#ifdef BLATHER
    if (accept)
        printf("    accept : TRUE\n");
    else
        printf("    accept : FALSE\n");
#endif /* BLATHER */

            if (accept) {
                state[i + nnode * j] = 1 - state_ij;
                d[i] = dpro[i];
                d[j] = dpro[j];
                for (int k = 0; k < nnode; ++k)
                    t[k] = tpro[k];
            } else {
                /* restore loop invariants */
                dpro[i] = d[i];
                dpro[j] = d[j];
                for (int k = 0; k < nnode; ++k)
                    tpro[k] = t[k];
            }

        }

            for (int i = 0; i < nnode; ++i)
                tbat[i] += t[i];
        }

        for (int i = 0; i < nnode; ++i)
                batch[kbatch++] = ((double) tbat[i]) / ((double) blen);
    }

    PutRNGstate();
}

