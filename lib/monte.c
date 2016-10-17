#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mcce.h"

/* normal pw is garanteed to be smaller than 2000. When it is is bigger than 5000, it is 
 * represents a value refenced to a clashed pair of conformers.
 */

/* self energy in prot, prot_free, prot_fixed and copied to conflist
 * before sampling
 * pairwise in pairwise
 * counter in conflist
 */
/* This Monte Carlo Program adds entropy correction to each conformer type different by
 * proton and electron, so that the ionized and neutral will get equal chance regardless of
 * the population of conformers.
 * It has been found that the drift of calculated pKa with rotamer number is caused by this
 * entropy effect
 */
typedef struct {
    float vdw0;
    float vdw1;
    float tors;
    float ebkb;
    float dsol;
    float offset;
    float pHpK0;
    float EhEm0;
    float TS;
    float residues;
    float total;
    float *mfePair;
    float *crg;
    float *vdw;
    float *ele;
} MFE;  // all the mfe energy terms

typedef struct  {
    int n;      /* Number of free confs on this residue */
    int on;
    int *conf;  /* idices of free confs */
    MFE mfe;   // let residue include mfe
} RESIDUE;

typedef struct  {
    int n;
    int *res;
} BIGLIST;

typedef struct {
    char  resName[4];
    char  chainID;
    char  iCode;
    int   resSeq;
    int   e;
    int   H;
    char  DM;
    float E_TS;
} CONFTYPE;

typedef struct {
   int n;
   CONFTYPE *conftype;
}  TYPES;

/* public variables */
PROT     prot;
RES      conflist, conflist_bak;
float    **pairwise, **pairwise_vdw, **pairwise_ele;
float    E_base, E_state;
float    ph, eh;
int      *state, *state_bak;

RESIDUE  *free_res, *fixed_res, *all_res, *mfe_res; // mfe_res is all the mfe residue
int      n_free, n_fixed, n_all, n_mfe;
BIGLIST  *biglist;   /* same length as free residues */
TYPES    Sconverge, SconvergeBak;
float    E_minimum;
FILE     *fp;
time_t   timerA, timerB, timerC;

float    **MC_occ;    /* occ of 3 parallel MC */
float    **occ_table;   /* occ of conformers at various pH/Eh */
float    **entropy_table;


int   verify_flag();
int   load_pairwise();
int   load_pairwise_vdw();
void  group_confs();
float get_base();
float get_E();
void  mk_neighbors();
void MC(int n);
int reduce_conflist();
PROT monte_load_conflist(char *fname);
int enumerate(int i_ph_eh);
int load_conflist();
int group_conftype();
int cmp_conftype(CONFTYPE t1, CONFTYPE t2);
CONFTYPE get_conftype(CONF conf);
void update_Sconvergence();
float s_stat();
float get_entropy(int iconf, int start, int end);
char **shead;

#define mev2Kcal 0.0235  // to keep consistent with mfe.py, use this constant instead of PH2KCAL/58, just for mfe part

int monte()
{   int i, j, k, counter, k_run;
    float *sigma, sigma_max, t;
    int N_smp;
    float S_max;
    PROT prot2;
    timerA = time(NULL);

    /* Load conformer list from FN_CONFLIST3 */
    printf("   Load conformer list from file \"%s\" ...\n", FN_CONFLIST3);
    fflush(stdout);
    prot = new_prot();
    prot2 = monte_load_conflist(FN_CONFLIST3);
    load_conflist();
    if (conflist.n_conf == 0) {
        printf("   FATAL: error in reading conformer list %s", FN_CONFLIST3);
        return USERERR;
    }
    printf("   Scaling factors:\n");
    printf("   VDW0  = %.3f:\n", env.scale_vdw0);
    printf("   VDW1  = %.3f:\n", env.scale_vdw1);
    printf("   VDW   = %.3f:\n", env.scale_vdw);
    printf("   TORS  = %.3f:\n", env.scale_tor);
    printf("   ELE   = %.3f:\n", env.scale_ele);
    printf("   DSOLV = %.3f:\n", env.scale_dsolv);
    printf("   Done\n\n"); fflush(stdout);


    /* load pairwise */
    printf("   Load pairwise interactions ...\n"); fflush(stdout);
    if (load_pairwise()) {
        printf("   FATAL: pairwise interaction not loaded\n");
        return USERERR;
    }
    printf("   Done\n\n");
    fflush(stdout);

    /* Verify self-consistancy between occupancy and flag */
    printf("   Verify self consistancy of conformer flag ...\n");
    if (verify_flag()) {
        printf("   Conformer flags updated due to self consistancy.\n\n");
    }
    else {
        printf("   Done\n\n");
    }
    fflush(stdout);

    /* backup */
    conflist_bak.n_conf = 0; conflist_bak.conf = NULL;
    cpy_res(&conflist_bak, &conflist);

    all_res = NULL;
    free_res = NULL;
    fixed_res = NULL;
    biglist  = NULL;
    state = NULL;


    /* titration */
    printf("   Do titration at %d points...\n", env.titr_steps);
    printf("   Detailed progress is in file \"%s\"\n", MC_OUT);
    fflush(stdout);
    if (env.monte_seed < 0) srand(time(NULL));
    else srand(env.monte_seed);
    
    MC_occ = (float **) malloc(env.monte_runs * sizeof(float *));
    for (i=0; i<env.monte_runs; i++) {
        if (!(MC_occ[i] = (float *) malloc(conflist.n_conf * sizeof(float)))) {
            printf("   FATAL: memory error.\n");
            return USERERR;
        }
        
    }
    
    occ_table = (float **) malloc(conflist.n_conf*sizeof(float *));
    for (i=0; i<conflist.n_conf; i++) occ_table[i] = (float *) malloc(env.titr_steps * sizeof(float));

    entropy_table = (float **) malloc(conflist.n_conf*sizeof(float *));
    for (i=0; i<conflist.n_conf; i++) entropy_table[i] = (float *) malloc(env.titr_steps * sizeof(float));

    sigma = (float *) malloc(conflist.n_conf*sizeof(float *));
    
    if (!env.minimize_size) {
       if (!(fp = fopen(MC_OUT, "w"))) {
          printf("   FATAL: Can not write to file \"%s\"\n", MC_OUT);
          return USERERR;
       }
    }
    else {
       fp = tmpfile();
    }
    
    timerB = time(NULL);
    printf("   Monte Carlo set up time: %ld seconds.\n", timerB-timerA); fflush(stdout);
    
    
    /* entropy sampling and correction */
    /* group conformers */
    group_conftype();  /* This ininitializes Sconverge and SconvergeBak */

    for (i=0; i<env.titr_steps; i++) {
        timerB = time(NULL);
        
        ph = env.titr_ph0;
        eh = env.titr_eh0;
        if (env.titr_type == 'p') ph = env.titr_ph0 + i*env.titr_phd;
        else eh = env.titr_eh0 + i*env.titr_ehd;
        fprintf(fp, "Titration (%dth) at pH = %6.2f Eh = %6.2f:\n", i+1, ph, eh); fflush(fp);

        
        /* reduce prot to free residues, conformers fixed as 0 occ are excluded from this free residues */
        free(conflist.conf);
        conflist.n_conf = 0; conflist.conf = NULL;
        cpy_res(&conflist, &conflist_bak);
        
        group_confs();
        fprintf(fp, "Identified %d free residues from %d residues.\n", n_free, n_all);
        fflush(fp);
        
        /* get big list */
        mk_neighbors();
        /* Compute base energy, self and mfe */
        E_base = get_base();
        
        /* get a microstate */
        state = realloc(state, n_free*sizeof(int));
        for (j=0; j<n_free; j++)
            state[j] = free_res[j].conf[rand() / (RAND_MAX/free_res[j].n + 1)];
        
        /* DEBUG
        for (j=0; j<n_free; j++) printf("%03d ", state[j]);
        printf("\n");
        */

        /* Do annealing */
        /*
        fprintf(fp, "Conformer list before reduction: E_base = %10.3f\n", E_base);
        for (j=0; j<conflist.n_conf; j++) {
            fprintf(fp, "%s %c %4.2f self = %8.3f mfe = %8.3f\n", conflist.conf[j].uniqID, 
            conflist.conf[j].on,
            conflist.conf[j].occ,
            conflist.conf[j].E_self0,
            conflist.conf[j].E_mfe);
        }
        fprintf(fp, "\n"); fflush(fp);
        */

        counter = 0;
        for (j=0; j<n_free; j++) counter+=free_res[j].n;
        N_smp = env.monte_nstart * counter;
        fprintf(fp, "Doing annealing... \n"); fflush(fp);
        if (N_smp) MC(N_smp);
        fprintf(fp, "Done\n\n");
        
        /* memcpy(state_check, state, sizeof(int)*n_free);  DEBUG */
        
        /* Do equalibriation */
        counter = 0;
        for (j=0; j<n_free; j++) counter+=free_res[j].n;
        N_smp = env.monte_neq * counter;
        fprintf(fp, "Doing equalibration ... \n"); fflush(fp);
        if (N_smp) MC(N_smp);
        fprintf(fp, "Done\n\n");
        
        /* do reduction */
        fprintf(fp, "Reducing conflist ... %d conformers are marked as fixed.\n\n", reduce_conflist());
        group_confs();
        fprintf(fp, "%d free residues from %d residues, after reduction.\n\n", n_free, n_all);
        fflush(fp);
        /* reset */
        mk_neighbors();
        j = 0; S_max = 999.9;
        for (k=0; k<Sconverge.n; k++) SconvergeBak.conftype[k].E_TS = 0.0;
        if (env.monte_tsx) {
           while (S_max > 0.5) {
             E_base = get_base();
             if (enumerate(i) == -1) { /* Non-ANALYTICAL SOLOTION */
                counter = 0;
                for (k=0; k<n_free; k++) counter+=free_res[k].n;
                N_smp = env.monte_niter * counter;

                fprintf(fp, "Doing Entropy sampling cycle %d...\n", j+1); fflush(fp);
                if (env.monte_nstart * counter) MC(env.monte_nstart * counter);
                if (N_smp) MC(N_smp);
             }
             update_Sconvergence(); /* calculate entropy from occupancy */
             S_max = s_stat();
             fprintf(fp, "Max delta = %8.3f\n", S_max);
           
             /* backup */
             for (k=0; k<Sconverge.n; k++) SconvergeBak.conftype[k] = Sconverge.conftype[k];
          
             j++;
           }
           fprintf(fp, "Done, exit at max entropy convergence %.3f \n\n", S_max);
        }
        
        
        E_base = get_base();
        /* ANALYTICAL SOLOTION */
        if (enumerate(i) == -1) {
            /*
            fprintf(fp, "Conformer list after reduction: E_base = %10.3f\n", E_base);
            for (j=0; j<conflist.n_conf; j++) {
                fprintf(fp, "%s %c %4.2f self = %8.3f mfe = %8.3f\n", conflist.conf[j].uniqID,
                conflist.conf[j].on,
                conflist.conf[j].occ,
                conflist.conf[j].E_self0,
                conflist.conf[j].E_mfe);
            }
            fprintf(fp, "\n"); fflush(fp);
            */
            counter = 0;
            for (j=0; j<n_free; j++) counter+=free_res[j].n;
            N_smp = env.monte_niter * counter;
            
            for (j=0; j<env.monte_runs; j++) {
                /* a new state */
                for (k=0; k<n_free; k++)
                    state[k] = free_res[k].conf[rand() / (RAND_MAX/free_res[k].n + 1)];
                
                fprintf(fp, "Doing annealing of MC %2d ...\n", j+1); fflush(fp);
                if (env.monte_nstart * counter) MC(env.monte_nstart * counter);
                
                fprintf(fp, "Doing MC %2d ... \n", j+1); fflush(fp);
                if (N_smp) MC(N_smp);
                for (k=0; k<conflist.n_conf; k++) {
                    MC_occ[j][k] = conflist.conf[k].occ;
                }
                fprintf(fp, "Done\n\n");
                 /* average -Yifan */
                for (k=0; k<conflist.n_conf; k++) conflist.conf[k].occ = 0.0;
                for (k_run=0; k_run<j+1; k_run++) {
                    for (k=0; k<conflist.n_conf; k++) {
                        conflist.conf[k].occ +=  MC_occ[k_run][k]/(j+1);
                    }
                }
                if (env.monte_tsx) {
                    update_Sconvergence(); /* calculate entropy from occupancy */
                    E_base = get_base();
                }
            }
            
            /* average */
            for (k=0; k<conflist.n_conf; k++) conflist.conf[k].occ = 0.0;
            for (j=0; j<env.monte_runs; j++) {
                for (k=0; k<conflist.n_conf; k++) {
                    conflist.conf[k].occ +=  MC_occ[j][k]/env.monte_runs;
                }
            }
            
            /* standard deviation */
            sigma_max = 0.0;
            for (j=0; j<conflist.n_conf; j++) {
                t = 0.0;
                for (k=0; k<env.monte_runs; k++) {
                    t += (MC_occ[k][j] - conflist.conf[j].occ) * (MC_occ[k][j] - conflist.conf[j].occ) ;
                }
                if (env.monte_runs > 1) sigma[j] = sqrt(t/(env.monte_runs-1));
                else sigma[j] = 999.00;
                
                if (sigma_max < sigma[j]) sigma_max = sigma[j];
            }
        }    
        /* write statstics */
        fprintf(fp, "Conformer     flag   E_self");
        for (j=0; j<env.monte_runs; j++) fprintf(fp, "   mc%02d", j+1);
        fprintf(fp, "      occ Sgm(n-1)\n");
        for (j=0; j<conflist.n_conf; j++) {
            occ_table[j][i] = conflist.conf[j].occ;
            entropy_table[j][i] = conflist.conf[j].E_TS;
            fprintf(fp, "%s   %c %8.3f", conflist.conf[j].uniqID,
            conflist.conf[j].on,
            conflist.conf[j].E_self0);
            for (k=0; k<env.monte_runs; k++) fprintf(fp, " %6.3f", MC_occ[k][j]);
            fprintf(fp, " Av=%5.3f Sg=%5.3f\n", conflist.conf[j].occ, sigma[j]);
        }
        fprintf(fp, "\n"); fflush(fp);
        
        timerC = time(NULL);
        printf("   Titration %2d: %5ld seconds, biggest stdev of conformer occ = %5.3f\n",
        i+1, timerC-timerB, sigma_max);
        fflush(stdout);
  
    }

    fclose(fp);
    printf("   Done\n\n"); fflush(stdout);


    /* writing occ table */
    if (!(fp = fopen(OCC_TABLE, "w"))) {
        printf("   FATAL: Can not write occupancy to file \"%s\"\n", OCC_TABLE);
        return USERERR;
    }
    if (env.titr_type == 'p') {
        fprintf(fp, " ph           ");
        for (i=0; i<env.titr_steps; i++) fprintf(fp, " %5.1f", env.titr_ph0+i*env.titr_phd);
    }
    else {
        fprintf(fp, " eh           ");
        for (i=0; i<env.titr_steps; i++) fprintf(fp, " %5.0f", env.titr_eh0+i*env.titr_ehd);
    }
    fprintf(fp, "\n");
    for (i=0; i<conflist.n_conf; i++) {
        fprintf(fp, "%s", conflist.conf[i].uniqID);
        for (j=0; j<env.titr_steps; j++) {
            fprintf(fp, " %5.3f", occ_table[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    /* writing entropy table */
    if (!(fp = fopen("entropy.out", "w"))) {
        printf("   FATAL: Can not write occupancy to file \"%s\"\n", "entropy.out");
        return USERERR;
    }
    if (env.titr_type == 'p') {
        fprintf(fp, " ph           ");
        for (i=0; i<env.titr_steps; i++) fprintf(fp, " %5.1f", env.titr_ph0+i*env.titr_phd);
    }
    else {
        fprintf(fp, " eh           ");
        for (i=0; i<env.titr_steps; i++) fprintf(fp, " %5.0f", env.titr_eh0+i*env.titr_ehd);
    }
    fprintf(fp, "\n");
    for (i=0; i<conflist.n_conf; i++) {
        fprintf(fp, "%s", conflist.conf[i].uniqID);
        for (j=0; j<env.titr_steps; j++) {
            fprintf(fp, " %5.3f", entropy_table[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);


    
    
    /* write self energy terms in right unit */
    
    
    
    
    timerC=time(NULL);
    printf("   Total time on MC: %ld seconds\n\n", timerC-timerA);
    
    printf("   Output files:\n");
    printf("      %-16s: MC progress and convergence.\n", MC_OUT);
    printf("      %-16s: Occupancy table.\n", OCC_TABLE);
    
    
    /* free MC stat */
    for (i=0; i<env.monte_runs; i++) {
        free(MC_occ[i]);
    }
    free(MC_occ);
    free(sigma);
    
    /* free conflist */
    free(conflist.conf);
    free(conflist_bak.conf);
    
    /* free pairwise */
    for (i=0; i<conflist.n_conf; i++)
        free(pairwise[i]);
    free(pairwise);
    
    /* free residues */
    for (i=0; i<n_free; i++)
        free(free_res[i].conf);
    free(free_res);
    
    for (i=0; i<n_fixed; i++)
        free(fixed_res[i].conf);
    free(fixed_res);
    
    for (i=0; i<n_all; i++)
        free(all_res[i].conf);
    free(all_res);

    free(Sconverge.conftype);
    //free(SconvergeBak.conftype);

        
    return 0;
}

int load_conflist()
{  FILE *fp;
    char sbuff[MAXCHAR_LINE];
    char stemp[MAXCHAR_LINE];
    CONF conf_temp;
    char notfound;
    int iconf, ires;
    int kr;
    int counter;
    
    conflist.n_conf = 0;
    conflist.conf   = NULL;
       
    if (!(fp=fopen(FN_CONFLIST3, "r"))) {
        printf("   FATAL: Can't open file %s\n", FN_CONFLIST3);
        return USERERR;
    }
    fgets(sbuff, sizeof(sbuff), fp); /* skip the first line */
    counter = 0;
    while(fgets(sbuff, sizeof(sbuff), fp)) {
        /* load this line to a conf template */
        if (strlen(sbuff) < 20) continue;
        sscanf(sbuff, "%d %s %c %f %f %f %f %d %d %f %f %f %f %f %f %s", &conf_temp.iConf,
        conf_temp.uniqID,
        &conf_temp.on,
        &conf_temp.occ,
        &conf_temp.netcrg,
        &conf_temp.Em,
        &conf_temp.pKa,
        &conf_temp.e,
        &conf_temp.H,
        &conf_temp.E_vdw0,
        &conf_temp.E_vdw1,
        &conf_temp.E_tors,
        &conf_temp.E_epol,
        &conf_temp.E_dsolv,
        &conf_temp.E_extra,
        conf_temp.history);

        conf_temp.E_TS = 0.0; /* initialize entropy effect at the time of loading conflist */

        /* rescale */
        conf_temp.E_vdw0 *= env.scale_vdw0;
        conf_temp.E_vdw1 *= env.scale_vdw1;
        conf_temp.E_epol *= env.scale_ele;
        conf_temp.E_tors *= env.scale_tor;
        conf_temp.E_dsolv*= env.scale_dsolv;

        strncpy(conf_temp.resName, conf_temp.uniqID, 3); conf_temp.resName[3] = '\0';
        strncpy(conf_temp.confName, conf_temp.uniqID, 5); conf_temp.confName[5] = '\0';
        conf_temp.chainID = conf_temp.uniqID[5];
        strncpy(stemp, conf_temp.uniqID+6, 4); stemp[4] = '\0';
        conf_temp.resSeq = atoi(stemp);
        conf_temp.iCode = conf_temp.uniqID[10];
        conf_temp.n_atom = 0;
        if (conf_temp.on == 't' || conf_temp.on == 'T') conf_temp.on = 't';
        else conf_temp.on = 'f';
        conf_temp.iConf = counter;
        /* creating conflist */
        iconf = ins_conf(&conflist, conflist.n_conf, 0);
        cpy_conf(&conflist.conf[iconf], &conf_temp);
        counter++;
        
        notfound = 1;
        for (kr=0; kr<prot.n_res; kr++) {
            if (conf_temp.chainID != prot.res[kr].chainID ||
                strncmp(conf_temp.resName, prot.res[kr].resName, 3) ||
            conf_temp.resSeq != prot.res[kr].resSeq ||
            conf_temp.iCode != prot.res[kr].iCode)
            continue;
            /* residue found */
            iconf = ins_conf(&prot.res[kr], prot.res[kr].n_conf, 0);
            cpy_conf(&prot.res[kr].conf[iconf], &conf_temp);
            
            notfound = 0;
            break;
        }
        if (notfound) { /* belongs to new residue */
            ires  = ins_res(&prot, prot.n_res);
            iconf = ins_conf(&prot.res[kr], prot.res[kr].n_conf, 0);

            cpy_conf(&prot.res[kr].conf[iconf], &conf_temp);

            strncpy(prot.res[ires].resName, conf_temp.resName,3);
            prot.res[ires].chainID = conf_temp.chainID;
            prot.res[ires].resSeq = conf_temp.resSeq;
            prot.res[ires].iCode = conf_temp.iCode;
        }
    }
    return 0;
}

int verify_flag()
/* this program will update conflist, the backup is in conflist_bak */
{  PROT old_prot;
    int kr, kc;
    float socc;
    int changed = 0;
    int n_free;

    
    memset(&old_prot, 0, sizeof(PROT));
    
    cpy_prot(&old_prot, &prot);
    
    for (kr=0; kr<prot.n_res; kr++) {
        /* set occ of free conformer to be 0 */
        for (kc=0; kc<prot.res[kr].n_conf; kc++) {
            if (prot.res[kr].conf[kc].on == 'f') {
                prot.res[kr].conf[kc].occ = 0.0;
                if (fabs(old_prot.res[kr].conf[kc].occ) > 0.001) {
                    printf("   %s %c %4.2f -> %c %4.2f (free conformer has 0.0 occ)\n",
                    old_prot.res[kr].conf[kc].uniqID,
                    old_prot.res[kr].conf[kc].on,
                    old_prot.res[kr].conf[kc].occ,
                    prot.res[kr].conf[kc].on,
                    prot.res[kr].conf[kc].occ);
                    changed++;
                }
            }
        }

        /* set t for single free conformer */
        socc = 0.0;
        n_free = prot.res[kr].n_conf;
        for (kc=0; kc<prot.res[kr].n_conf; kc++) {
            if (prot.res[kr].conf[kc].on == 't') {
                n_free --;
                socc += prot.res[kr].conf[kc].occ;
            }
        }
        if (n_free == 1) {         /* single free conf, convert it to t */
            for (kc=0; kc<prot.res[kr].n_conf; kc++) {
                if (prot.res[kr].conf[kc].on == 'f') {
                    prot.res[kr].conf[kc].on = 't';
                    prot.res[kr].conf[kc].occ = 1.0-socc;
                    if (prot.res[kr].conf[kc].occ < 0.0) prot.res[kr].conf[kc].occ = 0.0;
                    if (prot.res[kr].conf[kc].occ > 1.0) prot.res[kr].conf[kc].occ = 1.0;
                    printf("   %s %c %4.2f -> %c %4.2f (single conformer treated as fixed)\n",
                    old_prot.res[kr].conf[kc].uniqID,
                    old_prot.res[kr].conf[kc].on,
                    old_prot.res[kr].conf[kc].occ,
                    prot.res[kr].conf[kc].on,
                    prot.res[kr].conf[kc].occ);
                    changed++;
                }
            }
        }
        else if (n_free>1) {      /* multiple free conf, free confs take 1.00 occ */
            for (kc=0; kc<prot.res[kr].n_conf; kc++) {
                if (prot.res[kr].conf[kc].on == 't' && prot.res[kr].conf[kc].occ > 0.001) {
                    prot.res[kr].conf[kc].occ = 0.0;
                    printf("   %s %c %4.2f -> %c %4.2f (with multiple free conformers)\n",
                    old_prot.res[kr].conf[kc].uniqID,
                    old_prot.res[kr].conf[kc].on,
                    old_prot.res[kr].conf[kc].occ,
                    prot.res[kr].conf[kc].on,
                    prot.res[kr].conf[kc].occ);
                    changed++;
                }
            }
        }
    }

    del_prot(&old_prot);

    /* update conflist */
    for (kr=0; kr<prot.n_res; kr++) {
        for (kc=0; kc<prot.res[kr].n_conf; kc++) {
            conflist.conf[prot.res[kr].conf[kc].iConf].occ = prot.res[kr].conf[kc].occ;
            conflist.conf[prot.res[kr].conf[kc].iConf].on  = prot.res[kr].conf[kc].on;
        }
    }

    return changed;
}

int load_pairwise()
{   int i, j, kc;
    EMATRIX ematrix;

    ematrix.n = 0;  
    if (load_energies(&ematrix, ".", 0)<0) {
        printf("   File %s not found\n", ENERGY_TABLE);
        return USERERR;
    }

    /* scan conformers to see if all pw were calculated, DISABLED TEMPERORY 
    for (i=0; i<ematrix.n; i++) {
        if (!ematrix.conf[i].on && (ematrix.conf[i].uniqID[3] != 'D' || ematrix.conf[i].uniqID[4] != 'M')) {
           printf("      Incompleted delphi run, the first place detected at %s\n ", ematrix.conf[i].uniqID);
           return USERERR;
        }
    }
    */
    
    
    if (!(pairwise = (float **) malloc(ematrix.n * sizeof(float *)))) {
        printf("   FATAL: memory error in make_matrices()\n");
        return USERERR;
    }
    

    for (kc=0; kc<ematrix.n; kc++) {
        if (!(pairwise[kc] = (float *) malloc(ematrix.n * sizeof(float)))) {
            printf("   FATAL: memory error in make_matrices()\n");
            return USERERR;
        }
    }
    
    for (i=0; i<ematrix.n; i++) {
       for (j=0; j<ematrix.n; j++) {
           /*
           if (ematrix.pw[i][j].vdw > 990. && ematrix.pw[j][i].vdw > 990.) {
               pairwise[i][j] = pairwise[j][i] = 999.;
           }
           else {
               pairwise[i][j] = pairwise[j][i] = ((ematrix.pw[i][j].ori + ematrix.pw[j][i].ori)*env.scale_ele \
                   +(ematrix.pw[i][j].vdw + ematrix.pw[j][i].vdw)*env.scale_vdw)/2.0 ;
           }
           */
           pairwise[i][j] = pairwise[j][i] = ((ematrix.pw[i][j].ele + ematrix.pw[j][i].ele)*env.scale_ele
                                            +(ematrix.pw[i][j].vdw + ematrix.pw[j][i].vdw)*env.scale_vdw)/2.0 ;
           //pairwise[i][j] = fround3(ematrix.pw[i][j].ele) *env.scale_ele + fround3(ematrix.pw[i][j].vdw)*env.scale_vdw;
           /* proprocessing 
           if ((ematrix.pw[i][j].vdw + ematrix.pw[j][i].vdw) > 999.0) {
              pairwise[i][j] = pairwise[j][i] = 999.0;
           }
           else {
              pairwise[i][j] = pairwise[j][i] = ((ematrix.pw[i][j].ele + ematrix.pw[j][i].ele)*env.scale_ele \
                                                +(ematrix.pw[i][j].vdw + ematrix.pw[j][i].vdw)*env.scale_vdw)/2.0 ;
           }
           */
        }
    }

    /* DEBUG print pairwise table 
    int kr;
    printf("                ");
    for (kc=0; kc<conflist.n_conf; kc++) {
        printf("%16s", conflist.conf[kc].uniqID);
    }
    printf("\n");
    for (kr=0; kr<conflist.n_conf; kr++) {
        printf("%s ", conflist.conf[kr].uniqID);
        for (kc=0; kc<conflist.n_conf; kc++) {
            printf("%16.3f", pairwise[kr][kc]);
        }
        printf("\n");
    }
    */

    
    /* free memory */
    free_ematrix(&ematrix);

    return 0;
}

int load_pairwise_vdw()
{
    int i, j, kc;
    EMATRIX ematrix;

    ematrix.n = 0;  
    if (load_energies(&ematrix, ".", 0)<0) {
        printf("   File %s not found\n", ENERGY_TABLE);
        return USERERR;
    }

    if (!(pairwise_vdw = (float **) malloc(ematrix.n * sizeof(float *)))) {
        printf("   FATAL: memory error in make_matrices()\n");
        return USERERR;
    }
    for (kc=0; kc<ematrix.n; kc++) {
        if (!(pairwise_vdw[kc] = (float *) malloc(ematrix.n * sizeof(float)))) {
            printf("   FATAL: memory error in make_matrices()\n");
            return USERERR;
        }
    }
    
    for (i=0; i<ematrix.n; i++) {
       for (j=0; j<ematrix.n; j++) {
           pairwise_vdw[i][j] = pairwise_vdw[j][i] = (ematrix.pw[i][j].vdw + ematrix.pw[j][i].vdw)*env.scale_vdw/2.0 ;
        }
    }
    
    /* free memory */
    free_ematrix(&ematrix);

    return 0;
}

void group_confs()
{  int  kr, kc;
    int free_conf;
    int ir, ic;
    
    /* delete existing big list when regrouping */
    if (biglist != NULL) {
        for (ir=0; ir<n_free; ir++) {
            if (biglist[ir].res != NULL)
                free(biglist[ir].res);
        }
        free(biglist);
    }

    
    /* create all_res */
    if (all_res != NULL) {
        for (ir=0; ir<n_all; ir++) {
            if (all_res[ir].conf != NULL)
                free(all_res[ir].conf);
        }
        free(all_res);
    }
    all_res = (RESIDUE *) malloc(prot.n_res * sizeof(RESIDUE));
    n_all = prot.n_res;
    for (ir=0; ir<n_all; ir++) {
        all_res[ir].n    = prot.res[ir].n_conf;
        all_res[ir].conf = (int *) malloc(prot.res[ir].n_conf * sizeof(int));
        all_res[ir].mfe.crg = (float *) calloc(env.titr_steps, sizeof(float));
    }
    for (ir=0; ir<n_all; ir++) {
        for (ic=0; ic<all_res[ir].n; ic++) {
            all_res[ir].conf[ic] = prot.res[ir].conf[ic].iConf;
        }
    }
    
    
    /* create free residues */
    if (free_res != NULL) {
        for (ir=0; ir<n_free; ir++) {
            free(free_res[ir].conf);
        }
        free(free_res);
        free_res = NULL;
    }
    
    if (fixed_res != NULL) {
        for (ir=0; ir<n_fixed; ir++) {
            free(fixed_res[ir].conf);
        }
        free(fixed_res);
        fixed_res = NULL;
    }
    
    n_free = 0; n_fixed = 0;
    for (kr=0; kr<n_all; kr++) {
        free_conf = 0;
        for (kc=0; kc<all_res[kr].n; kc++) {
            if (conflist.conf[all_res[kr].conf[kc]].on == 'f') free_conf ++;
            /*
            printf("%s, %c, %f\n", conflist.conf[all_res[kr].conf[kc]].uniqID,
            conflist.conf[all_res[kr].conf[kc]].on,
            conflist.conf[all_res[kr].conf[kc]].occ);
            */
        }
        if (free_conf) {       /* free_conf stores number of free conformers in this residue */
            n_free++;
            free_res = realloc(free_res, n_free * sizeof(RESIDUE));
            free_res[n_free-1].n = free_conf;
            free_res[n_free-1].conf = (int *) malloc(free_conf*sizeof(int));
            
            ic = 0;
            for (kc=0; kc<all_res[kr].n; kc++) {
                if (conflist.conf[all_res[kr].conf[kc]].on == 'f') {
                    free_res[n_free-1].conf[ic] = conflist.conf[all_res[kr].conf[kc]].iConf;
                    ic++;
                }
                else { /* 't' must be with occ = 0.0, ignored */
                    if (fabs(conflist.conf[all_res[kr].conf[kc]].occ) > 0.000000001)  {
                        printf("   WARNING: %s has a fixed occ %.3f in a free residue\n",
                        conflist.conf[all_res[kr].conf[kc]].uniqID, conflist.conf[all_res[kr].conf[kc]].occ);
                        fflush(stdout);
                    }
                }
            }
        }
        else {
            n_fixed++;
            fixed_res = realloc(fixed_res, n_fixed * sizeof(RESIDUE));
            fixed_res[n_fixed-1].n = all_res[kr].n;
            fixed_res[n_fixed-1].conf = (int *) malloc(all_res[kr].n*sizeof(int));
            
            ic = 0;
            for (kc=0; kc<all_res[kr].n; kc++) {
                fixed_res[n_fixed-1].conf[ic] = conflist.conf[all_res[kr].conf[kc]].iConf;
                ic++;
            }
        }
    }
    
    return;
}

float get_base()
{  float E;
    int kr, kc, ir, ic;
    float mfe;
    
    E = 0.0;
    
    /* compute pH and Eh effect for each conformer */
    for (kc=0; kc<conflist.n_conf; kc++) {
        conflist.conf[kc].E_ph =
        env.monte_temp/ROOMT * conflist.conf[kc].H * (ph-conflist.conf[kc].pKa) * PH2KCAL;
        conflist.conf[kc].E_eh =
        env.monte_temp/ROOMT * conflist.conf[kc].e * (eh-conflist.conf[kc].Em) * PH2KCAL/58.0;
    }
    
    /* self without mfe */
    for (kc=0; kc<conflist.n_conf; kc++) {
        conflist.conf[kc].E_self0 = conflist.conf[kc].E_vdw0
        + conflist.conf[kc].E_vdw1
        + conflist.conf[kc].E_epol
        + conflist.conf[kc].E_tors
        + conflist.conf[kc].E_dsolv
        + conflist.conf[kc].E_extra
        + conflist.conf[kc].E_ph
        + conflist.conf[kc].E_eh
        + conflist.conf[kc].E_TS; /* Monte Carlo "knows" entropy, we want to "correct" it, so the signe is "+" */
/*        if (conflist.conf[kc].E_self0 > 900.0)
           conflist.conf[kc].E_self0 = 999.0;
*/
    }

    /* mfe on all conformers from fixed conformers */
    for (kr=0; kr<n_all; kr++) {
        for (kc=0; kc<all_res[kr].n; kc++) {
            mfe = 0.0;
            for (ir=0; ir<n_fixed; ir++) {
                for (ic=0; ic<fixed_res[ir].n; ic++) {
                    mfe += pairwise[all_res[kr].conf[kc]][fixed_res[ir].conf[ic]]
                    *conflist.conf[fixed_res[ir].conf[ic]].occ; /* self - self was already set to 0 */
                }
            }
            conflist.conf[all_res[kr].conf[kc]].E_mfe = mfe;
            conflist.conf[all_res[kr].conf[kc]].E_self = mfe + conflist.conf[all_res[kr].conf[kc]].E_self0;
        }
    }
    
    /* base energy of fixed conformers */
    for (ir=0; ir<n_fixed; ir++) {
        for (ic=0; ic<fixed_res[ir].n; ic++) {
            E+=conflist.conf[fixed_res[ir].conf[ic]].E_self * conflist.conf[fixed_res[ir].conf[ic]].occ;
        }
    }
    /* substract double counted mfe on fixed conformers */
    for (ir=0; ir<n_fixed; ir++) {
        for (kr=ir+1; kr<n_fixed; kr++) {
            for (ic=0; ic<fixed_res[ir].n; ic++) {
                for (kc=0; kc<fixed_res[kr].n; kc++) {
                    E-= pairwise[fixed_res[ir].conf[ic]][fixed_res[kr].conf[kc]]
                    *conflist.conf[fixed_res[ir].conf[ic]].occ
                    *conflist.conf[fixed_res[kr].conf[kc]].occ;
                }
            }
        }
    }
    
    return E;
}


float get_E()
{  float E = 0.0;
    int kr, ir;
    
    /* Triangl calculation, no bias but expensive */
    for (kr=0; kr<n_free; kr++) E += conflist.conf[state[kr]].E_self;
    for (kr=0; kr<n_free; kr++) {
        for (ir=0; ir<kr; ir++)
            E += pairwise[state[kr]][state[ir]];
    }
    
    return E;
}

void mk_neighbors()
{  int kr, ir, kc, ic;
    char big;
    
    /* here */
    biglist = (BIGLIST *) malloc(n_free * sizeof(BIGLIST));
    
    for (kr=0; kr<n_free; kr++) {
        biglist[kr].n = 0;
        biglist[kr].res = NULL;
        for (ir=0; ir<n_free; ir++) {
            if (kr==ir) continue;
            big = 0;
            for (kc=0; kc<free_res[kr].n; kc++) {
                if (big) break;
                for (ic=0; ic<free_res[ir].n; ic++) {
                    if (big) break;
                    if (fabs(pairwise[free_res[kr].conf[kc]][free_res[ir].conf[ic]])>env.big_pairwise)
                        big = 1;
                }
            }
            if (big) {
                biglist[kr].n++;
                biglist[kr].res = realloc(biglist[kr].res, biglist[kr].n*sizeof(int));
                biglist[kr].res[biglist[kr].n-1] = ir;
            }
        }
    }

    return;
}


void MC(int n)
{  int cycles, n_total, n_cycle;
    int i, j, k;
    register int iters;
    int mem;
    int *old_state;
    float  old_E;
    float dE;
    float b;
    int nflips;
    double H_average;

    int iflip, ires, iconf; /* iconf is 0 to n of the conf in a res */
    int old_conf, new_conf; /* old_conf and new_conf are from 0 to n_conf in conflist */


    b = -KCAL2KT/(env.monte_temp/ROOMT);

    mem =n_free * sizeof(int);
    old_state = (int *) malloc(mem);
    E_minimum = E_state = get_E();

    /* number of cycles and iters in each cycle */
    if (env.monte_trace > 0) {
        cycles = (n-1)/env.monte_trace + 1; /* round up */
        n_total = cycles*env.monte_trace;
        n_cycle = env.monte_trace;
    }
    else {
        cycles = 1;
        n_total = n_cycle = n;
    }

    /* clear counters */
    for (i=0; i<conflist.n_conf; i++) conflist.conf[i].counter = 0;
    H_average = 0.0;

    for (i=0; i<cycles; i++) {
        /*
        fprintf(fp, "Step %10d, E_minimum = %10.2f, E_running = %10.2f, E_reset = %10.2f\n",
        i*n_cycle, E_minimum+E_base, E_state+E_base, get_E()+E_base);
        */
        fprintf(fp, "Step %10d, E_minimum = %10.2f, E_running = %10.2f\n",
        i*n_cycle, E_minimum+E_base, E_state+E_base);
        fflush(fp);
        iters = n_cycle;
        while (iters) {
            /*  save state */
            old_E = E_state;
            memcpy(old_state, state, mem);

            /* 1st flip */
            ires  = rand()/(RAND_MAX/n_free + 1);
            while (1) {
                iconf = rand()/(RAND_MAX/free_res[ires].n + 1);
                old_conf = state[ires];
                new_conf = free_res[ires].conf[iconf];
                if (old_conf != new_conf) break;
            }
            state[ires] = new_conf;
            E_state += conflist.conf[new_conf].E_self - conflist.conf[old_conf].E_self;
            for (j=0; j<n_free; j++) {
                E_state += pairwise[new_conf][state[j]] - pairwise[old_conf][state[j]];
            }

            /* now multiple flip */
            /*   1st flip -> No (50% probablity)
             *     | Yes (50% probablity)
             *   2nd flip (any res in big list)
             *     |
             *   3rd flip (any res in big list)
             *     |
             *   4th flip (any res in big list)
             */
            if (rand() & 1) {   /* do multiple flip if odd number */
                if (biglist[ires].n) {
                    nflips = env.monte_flips > (biglist[ires].n+1) ? biglist[ires].n+1: env.monte_flips;
                    for (k=1; k<nflips; k++) {
                        
                        iflip = biglist[ires].res[rand()/(RAND_MAX/biglist[ires].n + 1)];
                        iconf = rand()/(RAND_MAX/free_res[iflip].n + 1);
                        old_conf = state[iflip];
                        new_conf = free_res[iflip].conf[iconf];

                        state[iflip] = new_conf;
                        E_state += conflist.conf[new_conf].E_self - conflist.conf[old_conf].E_self;
                        for (j=0; j<n_free; j++) {
                            E_state += pairwise[new_conf][state[j]] - pairwise[old_conf][state[j]];
                        }
                    }
                }
            }


            /* DEBUG
            for (j=0; j<n_free; j++) printf("%03d ", state[j]);
            printf("\n");
            */

            if (E_minimum > E_state) E_minimum = E_state;

            dE = E_state - old_E;
            if (dE < 0.0) {                                 /* go to new low */
            }
            /*<<< Boltman distribution >>>*/
            else if ((float) rand()/RAND_MAX < exp(b*dE)) { /* Go */
            }
            else {                                                    /* stay, restore the state */
                memcpy(state, old_state, mem);
                E_state = old_E;
            }

            /*
            if (!memcmp(state, state_check, sizeof(int)*n)) {
                printf("E=%12.5f\n", E_state);
            }  DEBUG */

            /* count this state energy */
            H_average += E_state/n_total;

            /*<<< Do statistics >>>*/
            for (j=0; j<n_free; j++) {
                conflist.conf[state[j]].counter++;
                //+= pairwise_vdw[j]
            }
            
            iters --;
        }
    }

    fprintf(fp, "Exit %10d, E_minimum = %10.2f, E_running = %10.2f\n", n_total, E_minimum+E_base, E_state+E_base);
    fprintf(fp, "The average running energy, corresponding to H, is %8.3f kCal/mol\n", H_average+E_base);
    fflush(fp);

    /* compute occ */
    for (i=0; i<conflist.n_conf; i++) {
        if (conflist.conf[i].on == 't') continue;
        conflist.conf[i].occ = (float) conflist.conf[i].counter / n_total;
    }

    free(old_state);

    return;
}

int reduce_conflist()
{  int i, j, t;
    int counter = 0;

    for (i=0; i<n_free; i++) {
        t = 0;
        for (j=0; j<free_res[i].n; j++) {
            if (conflist.conf[free_res[i].conf[j]].on == 'f'
            && conflist.conf[free_res[i].conf[j]].occ < env.monte_reduce) {
                conflist.conf[free_res[i].conf[j]].on = 't';
                conflist.conf[free_res[i].conf[j]].occ = 0.0;
                counter++;
                t++;
            }
        }
        if ((free_res[i].n-t) == 1) {  /* leave only one free conf, got to be 't' 1.00 */
            for (j=0; j<free_res[i].n; j++) {
                if (conflist.conf[free_res[i].conf[j]].on == 'f') {
                    conflist.conf[free_res[i].conf[j]].on = 't';
                    conflist.conf[free_res[i].conf[j]].occ = 1.0;
                    counter++;
                }
            }
        }
    }
    
    return counter;
}


int Nx;       /* number of titration points */
float *xp, *yp;   /* titration points */

struct STAT {
    float a;
    float b;
    float chi2;
};




PROT monte_load_conflist(char *fname) {
    PROT prot;
    FILE *fp;
    char sbuff[MAXCHAR_LINE];
    char stemp[MAXCHAR_LINE];
    CONF conf;
    int  k_res, k_conf;
    
    memset(&prot,0,sizeof(PROT));
    if (!(fp=fopen(fname, "r"))) {
        printf("   FATAL: Can't open file %s\n", fname);fflush(stdout);
        return prot;
    }
    
    fgets(sbuff, sizeof(sbuff), fp); /* skip the first line */
    
    while(fgets(sbuff, sizeof(sbuff), fp)) {
        if (strlen(sbuff) < 20) continue;
        
        /*
        iConf CONFORMER     FL  occ    crg   Em0  pKa0 ne nH    vdw0    vdw1    tors    epol   dsolv   extra    self
        01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
        00001 NTR01_0001_001 f 0.00  0.000     0  0.00  0  0  -0.002   1.335   0.000   0.005   0.000   0.000 01O000M000t
        */
        
        /* load this line to a conf */
        memset(&conf,0,sizeof(CONF));
        strncpy(stemp, sbuff,    5); stemp[5] = '\0'; conf.i_conf_prot = atoi(stemp);
        strncpy(conf.uniqID, sbuff+6, 14); conf.uniqID[14]  = '\0';
        conf.toggle = sbuff[21];
        strncpy(stemp, sbuff+23, 4); stemp[4] = '\0'; conf.occ     = atof(stemp);
        strncpy(stemp, sbuff+28, 6); stemp[6] = '\0'; conf.netcrg  = atof(stemp);
        strncpy(stemp, sbuff+35, 5); stemp[5] = '\0'; conf.Em      = atof(stemp);
        strncpy(stemp, sbuff+41, 5); stemp[5] = '\0'; conf.pKa     = atof(stemp);
        strncpy(stemp, sbuff+47, 2); stemp[2] = '\0'; conf.e       = atof(stemp);
        strncpy(stemp, sbuff+50, 2); stemp[2] = '\0'; conf.H       = atof(stemp);
        strncpy(stemp, sbuff+53, 7); stemp[7] = '\0'; conf.E_vdw0  = atof(stemp) * env.scale_vdw0;
        strncpy(stemp, sbuff+61, 7); stemp[7] = '\0'; conf.E_vdw1  = atof(stemp) * env.scale_vdw1;
        strncpy(stemp, sbuff+69, 7); stemp[7] = '\0'; conf.E_tors  = atof(stemp) * env.scale_tor;
        strncpy(stemp, sbuff+77, 7); stemp[7] = '\0'; conf.E_epol  = atof(stemp) * env.scale_ele;
        strncpy(stemp, sbuff+85, 7); stemp[7] = '\0'; conf.E_dsolv  = atof(stemp)* env.scale_dsolv;
        strncpy(stemp, sbuff+93, 7); stemp[7] = '\0'; conf.E_extra  = atof(stemp);
        strncpy(conf.history, sbuff+101, 11); conf.history[11]  = '\0';

        strncpy(conf.resName,  conf.uniqID, 3); conf.resName[3] = '\0';
        strncpy(conf.confName, conf.uniqID, 5); conf.confName[5] = '\0';
        conf.chainID = conf.uniqID[5];
        strncpy(stemp, conf.uniqID+6, 4); stemp[4] = '\0'; conf.resSeq = atoi(stemp);
        conf.iCode = conf.uniqID[10];

        /* search for the residue: using same procedure as in load_pdb() */
        for (k_res = prot.n_res - 1; k_res >= 0; k_res--) {
            if (!strcmp(conf.resName, prot.res[k_res].resName) &&
                conf.chainID == prot.res[k_res].chainID &&
            conf.resSeq  == prot.res[k_res].resSeq  &&
            conf.iCode   == prot.res[k_res].iCode  )
            {
                break;
            }
        }
        /* If couldn't find the residue, add a new one */
        if (k_res == -1) {
            k_res = ins_res(&prot, prot.n_res);
            strcpy(prot.res[k_res].resName, conf.resName);
            prot.res[k_res].chainID = conf.chainID;
            prot.res[k_res].resSeq  = conf.resSeq;
            prot.res[k_res].iCode   = conf.iCode;
            
            /* Insert a backbone conformer */
            ins_conf(&prot.res[k_res], 0, 0);
        }
        
        k_conf = ins_conf(&prot.res[k_res], prot.res[k_res].n_conf, 0);
        prot.res[k_res].conf[k_conf] = conf;
        if (conf.toggle == 't') prot.res[k_res].conf[k_conf].toggle = 'f';
        else if (conf.toggle == 'f') prot.res[k_res].conf[k_conf].toggle = 't';
        else if (conf.toggle == 'a') prot.res[k_res].conf[k_conf].toggle = 'a';
        else prot.res[k_res].conf[k_conf].toggle = 't';
    }
    fclose(fp);
    
    for (k_res=0;k_res<prot.n_res;k_res++) {
        for (k_conf=1;k_conf<prot.res[k_res].n_conf;k_conf++) {
            prot.nc++;
            prot.conf = realloc(prot.conf, prot.nc*sizeof(void *));
            prot.res[k_res].conf[k_conf].i_conf_prot = prot.nc-1;
            prot.conf[prot.res[k_res].conf[k_conf].i_conf_prot] = &prot.res[k_res].conf[k_conf];
        }
    }
    return prot;
}


int enumerate(int i_ph_eh)
{
    long    nstate,istate;
    float   b;
    int     ires,jres,iconf,old_conf,new_conf;
    float   *E_states,*occ_states;
    float   E_min, tot_occ;
    b = -KCAL2KT/(env.monte_temp/ROOMT);
    
    if (!n_free) {
        for (iconf=0; iconf<conflist.n_conf; iconf++) {
            occ_table[iconf][i_ph_eh] = conflist.conf[iconf].occ;
        }
        return 0;
    }
    
    /* each residue state with first conformer */
    nstate = 1;
    for (ires=0;ires<n_free;ires++) {
        free_res[ires].on=0;
        state[ires]=free_res[ires].conf[0];
        nstate=nstate*free_res[ires].n;
        if (nstate>env.nstate_max) return -1;
    }
    if (nstate>env.nstate_max) return -1;
    
    
    fprintf(fp, "%ld <= %d microstates, will use analytical solution\n\n", nstate, env.nstate_max);
    fflush(fp);
    
    E_states = malloc(nstate*sizeof(float));
    occ_states = malloc(nstate*sizeof(float));
    
    istate = 0;
    E_min = E_state = get_E();  /* get microstate energy */
    E_states[istate]=E_state;
    
    while (1) {
        istate++;
        ires = 0;
        old_conf = state[ires];
        free_res[ires].on++;
        while (free_res[ires].on>=free_res[ires].n) {
            free_res[ires].on=0;
            new_conf = state[ires] = free_res[ires].conf[0];
            
            E_state += conflist.conf[new_conf].E_self - conflist.conf[old_conf].E_self;
            for (jres=0; jres<n_free; jres++) {
                E_state += pairwise[new_conf][state[jres]] - pairwise[old_conf][state[jres]];
            }
            
            ires++;
            if (ires==n_free) break;
            old_conf = state[ires];
            free_res[ires].on++;
        }
        if (ires==n_free) break;
        
        new_conf = state[ires] = free_res[ires].conf[free_res[ires].on];
        E_state += conflist.conf[new_conf].E_self - conflist.conf[old_conf].E_self;
        for (jres=0; jres<n_free; jres++) {
            E_state += pairwise[new_conf][state[jres]] - pairwise[old_conf][state[jres]];
        }
        
        E_states[istate]=E_state;
        if (E_state<E_min) E_min = E_state;
    }
    
    tot_occ=0.;
    for (istate=0;istate<nstate;istate++) {
        occ_states[istate] = exp(b*(E_states[istate]-E_min));
        tot_occ += occ_states[istate];
    }
    for (istate=0;istate<nstate;istate++) {
        occ_states[istate] = occ_states[istate]/tot_occ;
    }
    
    /* get conformer occ */
    /* reset occ to 0 */
    for (iconf=0;iconf<conflist.n_conf;iconf++) {
        if (conflist.conf[iconf].on == 't') continue;
        conflist.conf[iconf].occ = 0.;
    }
    
    /* first microstate */
    istate = 0;
    for (ires=0;ires<n_free;ires++) {
        free_res[ires].on=0;
        state[ires]=free_res[ires].conf[0];
        
        iconf = free_res[ires].conf[free_res[ires].on];
        conflist.conf[iconf].occ += occ_states[istate];
    }
    
    while (1) {
        istate++;
        /* flipping (looping over) microstate */
        ires = 0;
        free_res[ires].on++;
        while (free_res[ires].on>=free_res[ires].n) {
            free_res[ires].on=0;
            new_conf = free_res[ires].conf[0];
            
            ires++;
            if (ires==n_free) break;
            free_res[ires].on++;
        }
        if (ires==n_free) break;
        
        /* collect occ for each working conformer */
        for (ires=0;ires<n_free;ires++) {
            iconf = free_res[ires].conf[free_res[ires].on];
            conflist.conf[iconf].occ += occ_states[istate];
        }
    }
    
    for (iconf=0; iconf<conflist.n_conf; iconf++) {
        occ_table[iconf][i_ph_eh] = conflist.conf[iconf].occ;
    }
    
    free(E_states);
    free(occ_states);
    return 0;
}

int group_conftype()
{  int i;
   int n = 0;
   CONFTYPE oldID, newID;
   
   oldID = get_conftype(conflist.conf[0]);
   for (i=1; i<conflist.n_conf; i++) {
      newID = get_conftype(conflist.conf[i]);
      if (cmp_conftype(oldID, newID) == 0); /* same group of conformers */
      else {
         n++;
         oldID = newID;
      }
   }
   n++;
   
   /* now we have number of conformer groups */
   Sconverge.n = SconvergeBak.n = n;
   if (!(Sconverge.conftype = (CONFTYPE *) malloc(n*sizeof(CONFTYPE)))) {
      printf("   Memory allocation error in \"group_conftype\"\n");
      return USERERR;
   }
   if (!(SconvergeBak.conftype = (CONFTYPE *) malloc(n*sizeof(CONFTYPE)))) {
      printf("   Memory allocation error in \"group_conftype\"\n");
      return USERERR;
   }

   /* real grouping */
   oldID = get_conftype(conflist.conf[0]);
   n = 0;
   Sconverge.conftype[0] = oldID;
   for (i=1; i<conflist.n_conf; i++) {
      newID = get_conftype(conflist.conf[i]);
      if (cmp_conftype(oldID, newID) != 0) {
         n++;
         oldID = newID;
         Sconverge.conftype[n] = oldID;
      }
   }

   /* backup */
   for (i=0; i<Sconverge.n; i++) SconvergeBak.conftype[i] = Sconverge.conftype[i];

   return 0;
}

int cmp_conftype(CONFTYPE t1, CONFTYPE t2)
{  if (!strcmp(t1.resName, t2.resName) && \
        t1.chainID==t2.chainID && \
        t1.iCode==t2.iCode && \
        t1.resSeq==t2.resSeq && \
        t1.e == t2.e && \
        t1.H == t2.H &&\
        t1.DM == t2.DM)
      return 0;
   return -1;
}

CONFTYPE get_conftype(CONF conf)
{  CONFTYPE typeid;
   strcpy(typeid.resName, conf.resName);
   typeid.chainID = conf.chainID;
   typeid.iCode = conf.iCode;
   typeid.resSeq = conf.resSeq;
   typeid.e = conf.e;
   typeid.H = conf.H;
   if (conf.uniqID[3] == 'D' && conf.uniqID[4] == 'M') typeid.DM = 't';
   else typeid.DM = '\0';
   return typeid;
}

void update_Sconvergence()
/* This function actually computes the entropy of each comformer type */
{  int i, pstart, pend;
   CONFTYPE newID;
   pstart = 0; /* pointer to the line of conflist of the first conformer */
   pend = 1;   /* pointer to the line of the next new conf type */

   for (i=0; i<Sconverge.n; i++) {
      while (pend < conflist.n_conf) {
         newID = get_conftype(conflist.conf[pend]);
         if (cmp_conftype(Sconverge.conftype[i], newID) == 0) {
            pend ++;
            continue; /* same group of conformers */
         }
         else break;
      }
      Sconverge.conftype[i].E_TS = get_entropy(i, pstart, pend);
      /* Have done that in get_entropy */
      /* for (j=pstart; j<pend; j++) conflist.conf[j].E_TS = Sconverge.conftype[i].E_TS; */

      /* done with is type, go to the next */
      pstart = pend;
      pend++;
   }

   return;
}


#define TINY 1.0E-10
#define NMAX 5000
#define GET_PSUM for (j=0; j<ndim; j++) {\
                        for (sum=0.0, i=0; i<mpts; i++) sum += p[i][j];\
                         psum[j] = sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void dhill(float **p, float *y, int ndim, float ftol, float (*funk)(float []), int *nfunk);


void dhill(float **p, float *y, int ndim, float ftol, float (*funk)(float []), int *nfunk)
{
    float dhtry(float **p, float *y, float *psum, int ndim, float (*funk)(float []), int ihi, float fac);
    int i, ihi, ilo, inhi,j, mpts = ndim+1;
    float rtol,sum,swap,ysave,ytry,*psum;
    
    
    psum = (float *)malloc(ndim * sizeof(float));
    *nfunk = 0;
    GET_PSUM
    
    for (;;) {
        ilo = 0;
        ihi = y[0] > y[1] ? (inhi = 1,0):(inhi = 0,1);
        for(i=0; i<mpts;i++) {
            if(y[i] <= y[ilo]) ilo=i;
            if(y[i] > y[ihi]) {
                inhi = ihi;
                ihi  = i;
            }
            else if (y[i] > y[inhi] && i != ihi) inhi = i;
        }
        
        /* DEBUG
        for (i=0; i<mpts; i++) {
            printf("%8.3f at (%8.3f, %8.3f)\n", y[i], p[i][0], p[i][1]);
        }
        printf("\n");
        */
        
        rtol = 2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
        
        if(rtol < ftol || *nfunk >= NMAX) {
            SWAP(y[0],y[ilo])
            for (i=0; i<ndim; i++) SWAP(p[0][i],p[ilo][i])
                break;
        }
        
        *nfunk += 2;
        
        ytry = dhtry(p,y,psum,ndim,funk,ihi,-1.0);
        
        if (ytry <= y[ilo]) ytry=dhtry(p,y,psum,ndim,funk,ihi,2.0);
        else if (ytry >= y[inhi]) {
            ysave = y[ihi];
            ytry=dhtry(p,y,psum,ndim,funk,ihi,0.5);
            if (ytry >= ysave) {
                for (i=0; i<mpts; i++) {
                    if (i !=ilo) {
                        for (j=0; j<ndim; j++)
                            p[i][j] = psum[j] = 0.5*(p[i][j]+p[ilo][j]);
                        y[i] = (*funk)(psum);
                    }
                }
                *nfunk += ndim;
                GET_PSUM
            }
        }
        else --(*nfunk);
    }
    free(psum);
}

float dhtry(float **p, float *y, float *psum, int ndim, float (*funk)(float []), int ihi, float fac)
{
    int j;
    float fac1, fac2, ytry, *ptry;
    
    ptry = (float *)malloc(ndim*sizeof(float));
    fac1 = (1.0-fac)/ndim;
    fac2 = fac1 - fac;
    for (j=0; j<ndim; j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
    ytry = (*funk)(ptry);   /* evaluate the function at the trial point */
    if (ytry < y[ihi]) {    /* if it's better than the highest, then replace the highest */
        y[ihi] = ytry;
        for (j=0; j<ndim; j++) {
            psum[j] += ptry[j] - p[ihi][j];
            p[ihi][j] = ptry[j];
        }
    }
    free(ptry);
    return ytry;
}

float s_stat()
/* This function returns the maximum difference of Sconvergence and the backup */
{  int i;
   float max = 0.0;
   float diff;
   
   for (i=0; i<Sconverge.n; i++) {
      diff = fabs(Sconverge.conftype[i].E_TS - SconvergeBak.conftype[i].E_TS);
      if (diff > max) max = diff;
   }
   return max;
}

float get_entropy(int iconf, int start, int end)
{  int i;
   float TS = 0.0;
   float sum = 0.0;
   
   for (i=start; i<end; i++) {
      sum += conflist.conf[i].occ;
   }

   if (sum<0.0001) {
      for (i=start; i<end; i++) {
         conflist.conf[i].E_TS = 0.0;
      }
      return 0.0;
   }

   for (i=start; i<end; i++) {
      if (conflist.conf[i].occ/sum > 0.00001)
         TS -= conflist.conf[i].occ/sum*log(conflist.conf[i].occ/sum)/1.688;
   }

   /* assign the calculated entropy to conformers */
   for (i=start; i<end; i++) {
      conflist.conf[i].E_TS = TS;
   }


   return TS;
}
