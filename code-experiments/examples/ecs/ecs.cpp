/* --------------------------------------------------------- */
/* --- File: ecs.c  -------- Author: ACMO                --- */
/* --------------------------------------------------------- */
/*   
    ECS for non-linear function minimization. 
    Copyright 2020 ACMO
    e-mail: alexandre.cesar .AT. ufma.br
    SOURCE: 
        https://github.com/cavalcantigor/ecs
*/

#define GARE "ECS15"
#define NOXLS
#define NODUMP
#define CONSO
#define NOCONVER
#define PLOT 0

/* PLOTS E CONVERGENCIA */
#define PLOTGERA 10
#define CONVAVAL 10

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* Problema */
#define MAXGER 180000
#define NUMCRU 200
#define TAXERR 0.0001

/*Clusters */
#define AAPALFA 0.05
//define AAPALFA 0.25F

//define INICLUS 3
#define MORTUS 0
#define GELADO 1
#define QUENTE 2
#define ASSIMPLES 0
#define NASRECOMBI 1
#define NASCAMINHO 2

/* Populacao*/
#define PPIOR 0.1

/* Operadores Evolutivos */
#define NUMSELS 2
#define ROLGIRO 6
#define BLXALFA 0.35
#define NI 2
#define MUTNUNI 1
#define PREBATI (rand() % 101 / 100.F)

/* Parametros comuns Busca Local */
#define ESCALA 0.1F
#define PSGRI 0.05F

/* Constantes */
#define TRUE 1
#define FALSE 0
#define INFINITO 9999999
#define PLOTODOS -1

#define max(x, y) (x > y ? x : y)
#define randi(x, y) (x + ((rand() % 1001 / 1000.) * (y - x)))

/*---------------------------------------------------------------------------#
Author : ACMO

ALGORITMO GENETICO COM POPULACAO FIXA, BLX E GEOMETRICO, MUTACAO NAO UNIFORME
        clock_t clock(void);

*/

/*************** G L O B A I S *******************/

struct
{
    int numAval;
} FuncaoTeste;

typedef struct
{
    double *var;
    double fit;
    int sel;
} Cromossomo;

typedef struct
{
    Cromossomo *indiv;
    Cromossomo centr;
    double sumFit;
    double media;
    double dvpad;
    int tamPop;
    int tamInd;
    int melhor;
    int pior;
    int numMuta;
    int iguais;
    int gerMelhor;
    int pai[NUMSELS];
} Populacao;

Populacao P;

FILE *saida = NULL;

#include "./funtesv5i.h"

int funcao;
double SOLUCAO;
int MAXVAR;
int ROLETA = 1;
int MAXAVA;
int PASSOS = 10;
double MUTPROB = 1;
int MAXPOP = 10;
double PPROMIS = 2.6F;
int NUMCLUS = 20;

typedef struct
{
    Cromossomo ponto;
    int conta;
    double alert;
    int stats;
} Centro;

// pos = ultima posica; num = total = pos - MORTUS; max = limite
typedef struct
{
    Centro *grupos;
    int posGrp;
    int numGrp;
    int maxGrp;
    double limiar;
    int densid;
} Prototipos;

void PlotPop(Populacao *Prox, double *p, Prototipos C, int numGeracoes, int flag)
{
    int i, j, ini, fim, index, cor;

    fprintf(saida, "\n figure ");
    fprintf(saida, "\n [c,h] = contour(X,Y,Z);");
    fprintf(saida, "\n hold on;");

    // POPULACAO INICIAL A CADA N CRUZAMENTOS/ A CADA GERACAO
    if (p == NULL)
    {
        fprintf(saida, "\nP%d=[", numGeracoes);
        for (i = 0; i < Prox->tamPop; i++)
        {
            fprintf(saida, "\n");
            for (j = 0; j < 2; j++)
                fprintf(saida, "%.10f ", Prox->indiv[i].var[j]);
        }
        fprintf(saida, "\n]; \n plot (P%d(:,1),P%d(:,2),'.k');", numGeracoes, numGeracoes);
    }
    else
    {
        fprintf(saida, "\np%d=[", numGeracoes);
        for (j = 0; j < 2; j++)
            fprintf(saida, "%.10f ", p[j]);
        fprintf(saida, "\n]; \nplot (p%d(:,1),p%d(:,2),'sb');", numGeracoes, numGeracoes);
    }
    if (flag == PLOTODOS)
    {
        ini = 0;
        fim = C.posGrp - 1;
    }
    else
    {
        ini = fim = flag;
    }
    for (index = ini; index <= fim; index++)
    {
        if (C.grupos[index].stats > MORTUS)
        {
            fprintf(saida, "\nB%d=[", index);
            for (j = 0; j < 2; j++)
                fprintf(saida, "%.10f ", C.grupos[index].ponto.var[j]);
            fprintf(saida, "\n]; \nplot (B%d(:,1),B%d(:,2),'ob');", index, index);
            fprintf(saida, "\n circle(%.2f,%.2f,%.2f,'r');", C.limiar, C.grupos[index].ponto.var[0], C.grupos[index].ponto.var[1]);
        }
    }
}

int CorrigeInviavel(double *xr, double linf, double lsup)
{
    if ((*xr) > lsup)
        *xr = lsup - PREBATI * ((*xr) - lsup) / ((*xr) - linf);
    else if ((*xr) < linf)
        *xr = linf + PREBATI * (linf - (*xr)) / (lsup - (*xr));
    return (1);
}

int MutaNaoUni(double *indiv, int tamind, int tampop, int ger, int expo, float pmut)
{
    int mutou = FALSE;
    float fRandVal, fFactor;
    float fNewt, fNewT;
    int iExponent, iIndex;

    for (iIndex = 0; iIndex < tamind; iIndex++)
    {
        if (rand() % 100 < pmut)
        {
            mutou = TRUE;
            fRandVal = (rand() % 101 / 100.);
            /* pick either the max or min. limit */
            if (fRandVal < 0.5) /* lower */
            {
                fNewt = ger;
                fNewT = MAXGER;
                fRandVal = (rand() % 101 / 100.);
                fFactor = pow((1.0F - (fNewt / fNewT)), expo) * fRandVal;
                if (fFactor < TAXERR / 10.0F)
                    fFactor = TAXERR / 10.0F;
                fFactor = fFactor * (indiv[iIndex] - FuncoesTeste[funcao].inf);
                indiv[iIndex] = indiv[iIndex] - fFactor;
            }
            else
            {
                fNewt = ger;
                fNewT = MAXGER;
                fRandVal = (rand() % 101 / 100.);
                fFactor = pow((1.0F - (fNewt / fNewT)), expo) * fRandVal;
                if (fFactor < TAXERR / 10.0F)
                    fFactor = TAXERR / 10.0F;
                fFactor = fFactor * (FuncoesTeste[funcao].sup - indiv[iIndex]);
                indiv[iIndex] = indiv[iIndex] + fFactor;
            }
        } // if prob
    }     // for

    return mutou;
}

double Simplex(double (*func)(double[], int n), double start[], double fstart, int n,
               double EPSILON, double SCALE, int MAX_IT, double ALPHA, double BETA, double GAMMA)
{
    // by ACMO
    double ajuste, nureal, interv;

    //
    int vs; /* vertex with smallest value */
    int vh; /* vertex with next smallest value */
    int vg; /* vertex with largest value */

    int i, j, m, row;
    int k;   /* track the number of function evaluations */
    int itr; /* track the number of iterations */

    double **v;    /* holds vertices of simplex */
    double pn, qn; /* values used to create initial simplex */
    double *f;     /* value of function at each vertex */
    double fr;     /* value of function at reflection point */
    double fe;     /* value of function at expansion point */
    double fc;     /* value of function at contraction point */
    double *vr;    /* reflection - coordinates */
    double *ve;    /* expansion - coordinates */
    double *vc;    /* contraction - coordinates */
    double *vm;    /* centroid - coordinates */
    double min;

    double fsum, favg, s, cent;

    /* dynamically allocate arrays */

    /* allocate the rows of the arrays */
    v = (double **)malloc((n + 1) * sizeof(double *));
    f = (double *)malloc((n + 1) * sizeof(double));
    vr = (double *)malloc(n * sizeof(double));
    ve = (double *)malloc(n * sizeof(double));
    vc = (double *)malloc(n * sizeof(double));
    vm = (double *)malloc(n * sizeof(double));

    /* allocate the columns of the arrays */
    for (i = 0; i <= n; i++)
    {
        v[i] = (double *)malloc(n * sizeof(double));
    }

    /* create the initial simplex */
    /* assume one of the vertices is 0,0 */

    pn = SCALE * (sqrt(n + 1) - 1 + n) / (n * sqrt(2));
    qn = SCALE * (sqrt(n + 1) - 1) / (n * sqrt(2));

    for (i = 0; i < n; i++)
    {
        v[0][i] = start[i];
    }

    for (i = 1; i <= n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i - 1 == j)
            {
                v[i][j] = pn + start[j];
            }
            else
            {
                v[i][j] = qn + start[j];
            }
        }
    }
    /* find the initial function values */
    f[0] = fstart;
    for (j = 1; j <= n; j++)
    {
        f[j] = func(v[j], n);
    }

    k = n + 1;

    /* begin the main loop of the minimization */
    for (itr = 1; itr <= MAX_IT; itr++)
    {
        /* find the index of the largest value */
        vg = 0;
        for (j = 0; j <= n; j++)
        {
            if (f[j] > f[vg])
            {
                vg = j;
            }
        }

        /* find the index of the smallest value */
        vs = 0;
        for (j = 0; j <= n; j++)
        {
            if (f[j] < f[vs])
            {
                vs = j;
            }
        }

        /* find the index of the second largest value */
        vh = vs;
        for (j = 0; j <= n; j++)
        {
            if (f[j] > f[vh] && f[j] < f[vg])
            {
                vh = j;
            }
        }

        /* calculate the centroid */
        for (j = 0; j <= n - 1; j++)
        {
            cent = 0.0;
            for (m = 0; m <= n; m++)
            {
                if (m != vg)
                {
                    cent += v[m][j];
                }
            }
            vm[j] = cent / n;
        }

        /* reflect vg to new vertex vr */
        for (j = 0; j <= n - 1; j++)
        {
            /*vr[j] = (1+ALPHA)*vm[j] - ALPHA*v[vg][j];*/
            vr[j] = vm[j] + ALPHA * (vm[j] - v[vg][j]);
        }
        fr = func(vr, n);
        k++;

        if (fr < f[vh] && fr >= f[vs])
        {
            for (j = 0; j <= n - 1; j++)
            {
                v[vg][j] = vr[j];
            }
            f[vg] = fr;
        }

        /* investigate a step further in this direction */
        if (fr < f[vs])
        {
            for (j = 0; j <= n - 1; j++)
            {
                /*ve[j] = GAMMA*vr[j] + (1-GAMMA)*vm[j];*/
                ve[j] = vm[j] + GAMMA * (vr[j] - vm[j]);
            }
            fe = func(ve, n);
            k++;

            /* by making fe < fr as opposed to fe < f[vs],
			   Rosenbrocks function takes 63 iterations as opposed
			   to 64 when using double variables. */

            if (fe < fr)
            {
                for (j = 0; j <= n - 1; j++)
                {
                    v[vg][j] = ve[j];
                }
                f[vg] = fe;
            }
            else
            {
                for (j = 0; j <= n - 1; j++)
                {
                    v[vg][j] = vr[j];
                }
                f[vg] = fr;
            }
        }

        /* check to see if a contraction is necessary */
        if (fr >= f[vh])
        {
            if (fr < f[vg] && fr >= f[vh])
            {
                /* perform outside contraction */
                for (j = 0; j <= n - 1; j++)
                {
                    /*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
                    vc[j] = vm[j] + BETA * (vr[j] - vm[j]);
                }
                fc = func(vc, n);
                k++;
            }
            else
            {
                /* perform inside contraction */
                for (j = 0; j <= n - 1; j++)
                {
                    /*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
                    vc[j] = vm[j] - BETA * (vm[j] - v[vg][j]);
                }
                fc = func(vc, n);
                k++;
            }

            if (fc < f[vg])
            {
                for (j = 0; j <= n - 1; j++)
                {
                    v[vg][j] = vc[j];
                }
                f[vg] = fc;
            }
            /* at this point the contraction is not successful,
			   we must halve the distance from vs to all the
			   vertices of the simplex and then continue.
			   10/31/97 - modified to account for ALL vertices.
			*/
            else
            {
                for (row = 0; row <= n; row++)
                {
                    if (row != vs)
                    {
                        for (j = 0; j <= n - 1; j++)
                        {
                            v[row][j] = v[vs][j] + (v[row][j] - v[vs][j]) / 2.0;
                        }
                    }
                }
                f[vg] = func(v[vg], n);
                k++;
                f[vh] = func(v[vh], n);
                k++;
            }
        }

        /* test for convergence */
        fsum = 0.0;
        for (j = 0; j <= n; j++)
        {
            fsum += f[j];
        }
        favg = fsum / (n + 1);
        s = 0.0;
        for (j = 0; j <= n; j++)
        {
            s += pow((f[j] - favg), 2.0) / (n);
        }
        s = sqrt(s);
        if (s < EPSILON)
            break;
    }
    /* end main loop of the minimization */

    /* find the index of the smallest value */
    vs = 0;
    for (j = 0; j <= n; j++)
    {
        if (f[j] < f[vs])
        {
            vs = j;
        }
    }
    for (j = 0; j < n; j++)
    {
        start[j] = v[vs][j];
        // rebate ponto invi�vel de volta ao espa�o de busca com concentra��o
        // perto da fronteira tanto maior qto for MENOR o valor de REBATIMENTO
        CorrigeInviavel(&start[j], FuncoesTeste[funcao].inf, FuncoesTeste[funcao].sup);
    }
    min = func(v[vs], n);
    k++;

    for (i = 0; i <= n; i++)
    {
        free(v[i]);
    }

    free(f);
    free(vr);
    free(ve);
    free(vc);
    free(vm);
    free(v);
    return min;
}

float randgen(float fLlim, float fUlim)
{
    float fRandomVal;

    fRandomVal = rand() % 101 / 100.; // rand entre 0 e 1

    return (fLlim + (float)(fRandomVal * (fUlim - fLlim)));
}

/* ************************  GENETICS **************************** */

void GeraIndividuos(Populacao *p, int mp, int mv, int melhor, int nfun)
{
    int i, j, pior;
    double soma, fit;

    // inicializa o centroide
    for (j = 0; j < mv; j++)
        p->centr.var[j] = 0.0F;
    for (i = 0, soma = 0, pior = 0; i < mp; i++)
    {
        // gera individuo e acumula-o no centroide (nao ocorre a divisao - centroide = soma)
        for (j = 0; j < mv; j++)
        {
            p->indiv[i].var[j] = (double)randgen(FuncoesTeste[nfun].inf, FuncoesTeste[nfun].sup);
            p->centr.var[j] += p->indiv[i].var[j];
        }
        fit = funccod[nfun](p->indiv[i].var, mv);
        p->indiv[i].fit = fit;
        p->indiv[i].sel = 0;
        if (fit > p->indiv[pior].fit)
            pior = i;
        if (fit < p->indiv[melhor].fit)
            melhor = i;
        soma += (fit);
    }
    // retira do centroide a parte do pior individuo, antecipando, pois ele sera substituido
    // na primeira atualiza��o
    for (j = 0; j < mv; j++)
        p->centr.var[j] -= p->indiv[pior].var[j];

    p->sumFit = soma;
    p->melhor = melhor;
    p->pior = pior;
    p->media = p->sumFit / p->tamPop;
}

void IniciaCls(Prototipos *c, int mc, int mv)
{
    int i;
    c->grupos = (Centro *)malloc(sizeof(Centro) * NUMCLUS);

    for (i = 0; i < mc; i++)
    {
        c->grupos[i].conta = 0;
        c->grupos[i].stats = MORTUS;
        c->grupos[i].alert = 0.0F;
        c->grupos[i].ponto.var = (double *)malloc(sizeof(double) * mv);
        if (c->grupos[i].ponto.var == NULL)
        {
            fprintf(saida, "ERRO/IniciaCls(#1): Memoria!!!");
            exit(-1);
        }
    }
    c->maxGrp = mc;
    c->posGrp = 0;
    c->numGrp = 0;
    c->limiar = 0.0F;
}

void IniciaPop(Populacao *p, int mp, int mv)
{
    int i;

    p->centr.var = (double *)malloc(sizeof(double) * mv);
    if (p->centr.var == NULL)
    {
        fprintf(saida, "ERRO(1): Problemas de memoria!!!");
        exit(-1);
    }

    p->indiv = (Cromossomo *)malloc(sizeof(Cromossomo) * mp);
    if (p->indiv == NULL)
    {
        fprintf(saida, "ERRO(2a): Problemas de memoria!!!");
        exit(-1);
    }

    for (i = 0; i < mp; i++)
    {
        p->indiv[i].var = (double *)malloc(sizeof(double) * mv);
        if (p->indiv[i].var == NULL)
        {
            fprintf(saida, "ERRO(2b): Problemas de memoria!!!");
            exit(-1);
        }
    }

    p->tamPop = mp;
    p->tamInd = mv;
    p->numMuta = 0;
    p->iguais = 0;
}

void RoletaPressaoSeletiva(Populacao *p, double melhorfit)
{
    int i, pos, sel, fator;
    double z, gatilho, acum;

    sel = 0;
    while (sel < NUMSELS)
    {
        pos = rand() % p->tamPop;
        fator = (p->indiv[pos].sel > 3 ? 3 : p->indiv[pos].sel);
        //        fator = 2;
        z = (1.0F / pow((p->indiv[pos].fit - melhorfit + 1), fator));
        gatilho = z * (ROLGIRO - (ROLGIRO - 1) * z);
        acum = 0;
        i = 0;
        while (acum < gatilho && i <= ROLGIRO)
        {
            pos = (pos < p->tamPop - 1 ? pos + 1 : 0);
            fator = (p->indiv[pos].sel > 3 ? 3 : p->indiv[pos].sel);
            //                fator = 2;
            z = 1.0F / pow((p->indiv[pos].fit - melhorfit + 1), fator);
            acum += z;
            i++;
        }
        p->indiv[pos].sel++;
        p->pai[sel] = pos;
        sel++;
    }
}

void CruzaBlend(Populacao *p, int pai, int mae, int filho, float alfa)
{
    double a, b, r;
    int i;

    a = -alfa;
    b = 1 + alfa;

    for (i = 0; i < p->tamInd; i++)
    {
        r = a + (rand() % 101 / 100.) * (b - a);
        // gera filho
        p->indiv[filho].var[i] = p->indiv[pai].var[i] + r * (p->indiv[mae].var[i] - p->indiv[pai].var[i]);
        // rebate se invi�vel
        CorrigeInviavel(&(p->indiv[filho].var[i]), FuncoesTeste[funcao].inf, FuncoesTeste[funcao].sup);
    }
}

int AtualizaGrp(Prototipos *c, Populacao *p)
{
    int i, j, k;
    int indice, pertice, dispice;
    double dist, menorDist;
    int maiorCont = 0;
    int numsels = NUMSELS;
    char pertence;

    for (i = 0; i < numsels; i++)
    {
        pertence = FALSE;
        menorDist = INFINITO;
        // insercao a priori na ultima posicao
        dispice = c->posGrp;
        for (j = 0; j < c->posGrp; j++)
        {
            if (c->grupos[j].stats != MORTUS)
            {
                dist = sqrt(DistEucl(p->indiv[p->pai[i]].var, c->grupos[j].ponto.var, p->tamInd));
                if (dist < c->limiar && !pertence)
                {
                    pertence = TRUE;
                    pertice = j; // relativo a pertincia
                }
                if (dist < menorDist)
                {
                    menorDist = dist;
                    indice = j; // relativo a assimila��o
                }
            }
            else
                dispice = j;
        }
        // se n�o pertencer a ninguem o ultimo mortus sera usado
        // se n�o entrar nenhuma vez no if, vale a inicializacao (ultimo centro)
        if (!pertence && !(dispice >= c->maxGrp))
        {
            // salva o primeiro ponto e seu fitness
            memcpy(c->grupos[dispice].ponto.var, p->indiv[p->pai[i]].var, p->tamInd * sizeof(double));
            c->grupos[dispice].ponto.fit = p->indiv[p->pai[i]].fit;
            // Cluster come�a gelado
            c->grupos[dispice].conta = 1;
            c->grupos[dispice].stats = GELADO;
            // aumenta posGrp se inseriu na ultima posicao
            c->posGrp += (dispice == c->posGrp);
            // independente disso, aumenta o num
            c->numGrp++;
        }
        else
        { // ou pertence ou nao cabe mais clusters, entao assimila o cluster mais proximo
            if (menorDist > TAXERR)
            {
//--------------------------------------------------
// protege de assimila��o pontos muito pr�ximos ao centro
// TRES TIPOS DE ASSIMILA��O
//--------------------------------------------------
#ifdef ASSIMPLES
                // ASSIMPLES usa um AAPALFA para gerar um novo centros na reta entre o antigo e o ponto assimilado
                // N�o precisa avaliar o novo centro, isso � feito antes de Hooke
                // Se for executar sem BLOCAL, mas n�o deixar de avalia-lo em IIIntenso

                {
                    int flag = 0;
                    if (rand() % 1 == 0 && PLOT)
                    {
                        PlotPop(p, p->indiv[p->pai[i]].var, *c, indice, indice);
                        flag = 1;
                    }
                    for (k = 0; k < p->tamInd; k++)
                    {
                        c->grupos[indice].ponto.var[k] += AAPALFA * (p->indiv[p->pai[i]].var[k] - c->grupos[indice].ponto.var[k]);
                    }
                    if (flag && PLOT)
                        PlotPop(p, p->indiv[p->pai[i]].var, *c, indice, indice);

                } // auto
#endif
#ifdef ASRECOMBI
                // ASRECOMBI usa n - AAPALFA's para gerar um novo centros no hiperplano entre o antigo e o ponto assimilado
                // N�o precisa avaliar o novo centro, isso � feito antes de Hooke
                // Se for executar sem BLOCAL, mas n�o deixar de avalia-lo em IIIntenso

                {
                    double a, b, r;
                    int flag = 0;
                    a = -AAPALFA;
                    b = 1 + AAPALFA;

                    if (rand() % 1 == 0 && PLOT)
                    {
                        PlotPop(p, p->indiv[p->pai[i]].var, *c, indice, indice);
                        flag = 1;
                    }

                    for (k = 0; k < p->tamInd; k++)
                    {
                        r = a + (rand() % 101 / 100.) * (b - a);
                        c->grupos[indice].ponto.var[k] += r * (p->indiv[p->pai[i]].var[k] - c->grupos[indice].ponto.var[k]);
                        CorrigeInviavel(&(c->grupos[indice].ponto.var[k]), FuncoesTeste[funcao].inf, FuncoesTeste[funcao].sup);
                    }
                    if (flag && PLOT)
                        PlotPop(p, p->indiv[p->pai[i]].var, *c, indice, indice);

                } // auto
#endif

#ifdef ASCAMINHO
                // ASCAMINHO usa v�rios AAPALFA para gerar v�rios novos centros entre o antigo e o ponto assimilado
                // Ele j� avalia o novo centro. Em Hooke n�o precisa avaliar de novo

                {
                    double *aux, fitaux, *sau, fitsau;
                    int namost, flag = 0;

                    aux = (double *)malloc(p->tamInd * sizeof(double));
                    sau = (double *)malloc(p->tamInd * sizeof(double));
                    memcpy(aux, c->grupos[indice].ponto.var, p->tamInd * sizeof(double));
                    fitsau = c->grupos[indice].ponto.fit;
                    namost = (int)floor(1.0F / AAPALFA);
                    while (--namost)
                    {
                        if (rand() % 1 == 0 && PLOT)
                        {
                            PlotPop(p, p->indiv[p->pai[i]].var, *c, indice, indice);
                            flag = 1;
                        }
                        for (k = 0; k < p->tamInd; k++)
                        {
                            aux[k] += AAPALFA * (p->indiv[p->pai[i]].var[k] - c->grupos[indice].ponto.var[k]);
                        }
                        fitaux = funccod[funcao](aux, p->tamInd);
                        if (fitaux < fitsau)
                        {
                            memcpy(sau, aux, p->tamInd * sizeof(double));
                            fitsau = fitaux;
                            if (flag && PLOT)
                            {
                                PlotPop(p, p->indiv[p->pai[i]].var, *c, indice, indice);
                                flag = 0;
                            }
                        }
                    }
                    memcpy(aux, p->indiv[p->pai[i]].var, p->tamInd * sizeof(double));
                    do
                    {
                        if (rand() % 1 == 0 && PLOT)
                        {
                            PlotPop(p, p->indiv[p->pai[i]].var, *c, indice, indice);
                            flag = 1;
                        }

                        for (k = 0; k < p->tamInd; k++)
                        {
                            aux[k] += AAPALFA * (p->indiv[p->pai[i]].var[k] - c->grupos[indice].ponto.var[k]);
                            CorrigeInviavel(&aux[k], FuncoesTeste[funcao].inf, FuncoesTeste[funcao].sup);
                        }
                        fitaux = funccod[funcao](aux, p->tamInd);
                        if (fitaux < fitsau)
                        {
                            memcpy(sau, aux, p->tamInd * sizeof(double));
                            fitsau = fitaux;
                            if (flag && PLOT)
                            {
                                PlotPop(p, p->indiv[p->pai[i]].var, *c, indice, indice);
                                flag = 0;
                            }
                        }
                        else
                            break;
                    } while (1);
                    if (fitsau < c->grupos[indice].ponto.fit)
                    {
                        memcpy(c->grupos[indice].ponto.var, sau, p->tamInd * sizeof(double));
                        c->grupos[indice].ponto.fit = fitsau;
                    }
                    free(aux);
                    free(sau);
                } // auto
#endif

            } //if
            if (pertence)
            { // pertinente � o cluster proximo (limiar) mais antigo (primeiro a ser encontrado)
                c->grupos[pertice].conta++;
                // guarda a maior contagem
                if (c->grupos[pertice].conta > maiorCont)
                    maiorCont = c->grupos[pertice].conta;

                // independente de ser o mesmo ou n�o, o cluster � esquentado
                c->grupos[pertice].stats = QUENTE;

                // se houver dispice do lado esquerdo de pos, pos � decrementado
                c->posGrp -= (c->posGrp == dispice + 1);
            } //if
        }     // else (pertence)
    }
    return (maiorCont);
}

void AtualizaPop(Populacao *p, int pos, double fit, int ger)
{
    int i, j, salto, maxIt;

    /*
        Tira de soma o fit de quem vai ser substituido
        Poe em soma o fit de quem vai substituir
*/
    p->sumFit -= (p->indiv[pos].fit);
    p->sumFit += (fit);
    p->indiv[pos].fit = fit;
    p->indiv[pos].sel = 0;

    p->media = p->sumFit / p->tamPop;

    if (fit < p->indiv[p->melhor].fit)
    {
        p->melhor = pos;
        p->gerMelhor = ger;
    }
    /* ***** procura um outro pior ******** */
    maxIt = (int)ceil(PPIOR * p->tamPop);
    salto = (int)ceil(maxIt / 3.0F);
    //        p->pior=p->melhor == p->tamPop ? p->melhor - 1: p->melhor +1;
    for (i = 0; i < maxIt; i++)
    {
        j = rand() % p->tamPop;
        if (p->indiv[j].fit > p->indiv[p->pior].fit)
        {
            p->pior = j;
            i += salto;
        }
    }
    /*
        Antecipando que o pior vai sair, retira-o do centroide e poe o novo(POS)
        Na hora de considerar o criterio, calcula o centroide efetivo para tampop-1
        individuos
*/

    for (j = 0; j < p->tamInd; j++)
        p->centr.var[j] = p->centr.var[j] + p->indiv[pos].var[j] - p->indiv[p->pior].var[j];
}

/*
        H O O K E - J E E V E S  ROUTINES
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*/

double HookExplore(double (*fobj)(double[], int n), double *xr, double fp, double dx, int n)
{
    int i, j;
    double fr;
    double linf, lsup, salvo;

    linf = FuncoesTeste[funcao].inf;
    lsup = FuncoesTeste[funcao].sup;

    for (i = 0; i < n; i++)
    {
        // first direction
        salvo = xr[i];
        xr[i] = salvo + dx;
        // viability
        CorrigeInviavel(&xr[i], linf, lsup);
        // evaluate
        fr = fobj(xr, n);
        if (fr < fp)
        {
            // success
            fp = fr;
        }
        else
        {
            // failure: other direction
            dx = -dx;
            xr[i] = salvo + dx;
            // viability
            CorrigeInviavel(&xr[i], linf, lsup);
            // evaluate
            fr = fobj(xr, n);
            if (fr < fp)
            {
                // success
                fp = fr;
            }
            else
            {
                // reset direction: ACMO bichado por que houve corre��o
                xr[i] = salvo;
            }
        }
    }
    return (fp);
}

double HookeJeeves(double (*fobj)(double[], int n), double xc[], double fc, int n,
                   double epsilon, int passos, double scala, double step)
{
    double linf, lsup, dx, err, fp, inif;
    static double mdx = 1.0F, melf = 0.0F;
    static int cont = 100;
    int i, m;
    char reduz;

    double *xr = (double *)NULL;
    double *xp = (double *)NULL;

    linf = FuncoesTeste[funcao].inf;
    lsup = FuncoesTeste[funcao].sup;

    xp = (double *)malloc(n * sizeof(double));
    xr = (double *)malloc(n * sizeof(double));
    if (xr == NULL || xp == NULL)
    {
        fprintf(saida, "ERRO(5): Problemas de memoria!!!");
        fprintf(saida, "ABEND 102");
    }

    inif = fc;
    if (cont > 0)
    {
        dx = step;
        cont--;
    }
    else
    {
        dx = step = mdx;
    }

    m = 0;

    while (m <= passos && FuncaoTeste.numAval < MAXAVA)
    {
        // Assign base point
        fp = fc;
        memcpy(xr, xc, n * sizeof(double));
        fp = HookExplore(fobj, xr, fp, dx, n);
        // if it doesnt get into; it must be reduced
        reduz = TRUE;
        while (fp < fc && fabs(fp - fc) > epsilon && fabs(fp - SOLUCAO) > epsilon && FuncaoTeste.numAval < MAXAVA)
        {
            reduz = FALSE;
            // set base point
            fc = fp;
            memcpy(xp, xc, n * sizeof(double));
            memcpy(xc, xr, n * sizeof(double));
            for (i = 0; i < n; i++)
            {
                xr[i] = xr[i] + (xr[i] - xp[i]);
                CorrigeInviavel(&xr[i], linf, lsup);
            }
            fp = fobj(xr, n);
            fp = HookExplore(fobj, xr, fp, dx, n);
        }
        if (reduz && fabs(fp - SOLUCAO) > epsilon)
        {
            /*     for (i=0;i<n;i++) {
                                dx[i] = scala * dx[i];
                        }
                        m++; s� incrementa m se houver redu��es*/
            dx = scala * dx;
        }
        // difere do original -- sempre incrementa m
        m++;
    }
    if (inif - fc > melf)
    {
        mdx = step;
        melf = inif - fc;
    }
    free(xr);
    free(xp);
    return (fc);
}

/*******************************************************************/

void EstratIIIntenso(Prototipos &C, Populacao &P, int &achou, int &bLocOk, int &bLocTot, double step, int funcao, double taxerr, int passos, double escala, double solucao)
{
    int index;
    double erro, fitant, fitdep;

    int j;

    for (index = 0; (!achou) && (index < C.posGrp); index++)
    {

        if (C.grupos[index].conta >= C.densid)
        {
#ifdef ASCAMINHO
            // N�o precisa avaliar o centro
            fitant = C.grupos[index].ponto.fit;
#else
            // precisa avaliar o centro
            fitant = funccod[funcao](C.grupos[index].ponto.var, P.tamInd);
#endif
            fitdep = HookeJeeves(funccod[funcao], C.grupos[index].ponto.var, fitant, P.tamInd, taxerr, passos, escala, step);
            if (fitdep < P.indiv[P.melhor].fit)
            {
                memcpy(P.indiv[P.melhor].var, C.grupos[index].ponto.var, P.tamInd * sizeof(double));
                P.indiv[P.melhor].fit = fitdep;
                erro = (double)max(P.indiv[P.melhor].fit - solucao, 0.0F);
                achou = (erro <= taxerr ? TRUE : FALSE);
            }
            bLocTot++;
            if (fitdep < fitant)
                bLocOk++;
            //reinicia so se foi verificado
            C.grupos[index].conta = 1;
        }
    }
    return;
}

void EsfriaGrupos(Prototipos &C)
{
    int index, j;

    for (index = 0; index < C.posGrp; index++)
    {
        if (C.grupos[index].stats > MORTUS)
        {
            C.grupos[index].stats--;
            C.numGrp -= (C.grupos[index].stats == MORTUS);
        }
    }
    return;
}

int main(int argc, char *argv[])
{

    Prototipos C;

    clock_t start, end;

    char nomarq[100];

    int numGeracoes, numCruza, bLocTot, bLocOk, nGrupos = 0;

    double erro, step, fit, dist;

    int i, j, achou, mutou;

    int semente;

    if (argc < 12)
    {
        printf("\n ERRO(3): Falta argumentos! \nTente %s <#prob> <#var> < .ext > <#semente>", argv[0]);
        puts("\nOnde ...");
        puts("#prob = {dj1=0,  dj3, dj4, dj5, uno, bmp=5,  gol(2), fea(2), sph, ros}");
        puts("#prob = {sch=10, ras, gri, mic(5/10), lan(5/10), she=15, ack, rot}");
        exit(0);
    }
    funcao = atoi(argv[1]);
    strcpy(nomarq, GARE);
    strcat(nomarq, FuncoesTeste[funcao].nom);
    strcat(&nomarq[strlen(FuncoesTeste[funcao].nom)], argv[3]);
    if (!(saida = fopen(nomarq, "w")))
    {
        perror("");
        exit(-1);
    }
    MAXVAR = atoi(argv[2]);
    if (MAXVAR < 2)
    {
        puts("ERRO(4): Problemas com o numero de variaveis!!");
        exit(-1);
    }

    MAXAVA = atoi(argv[4]);   // max calls to objective function
    MAXPOP = atoi(argv[5]);   // max population size (10 - 20 - 60)
    MUTPROB = atof(argv[6]);  // mutation likelihood (0.01 - 0.05)
    NUMCLUS = atoi(argv[7]);  // max clusters (20)
    PASSOS = atoi(argv[8]);   // max steps to perform local search (10 * DIMENSION)
    PPROMIS = atof(argv[9]);  // promiscuous likelihood (1.6 - 2.6)
    ROLETA = atoi(argv[10]);  // uses roulette wheel selection (0 - false, >= 1 true)
    semente = atoi(argv[11]); // execution seed (for repeatability)

    /* randomico ou nao */
    srand((unsigned)time(0) + semente);

    /***************** SOLUCAO das funcoes *************/

    if ((funcao == fea || funcao == gol || funcao == dj5) && MAXVAR != 2)
        MAXVAR = 2;
    if ((funcao == she) && MAXVAR != 4)
        MAXVAR = 4;
    if ((funcao == har) && MAXVAR != 6)
        MAXVAR = 6;
    if ((funcao == lan0 || funcao == lan1 || funcao == mic) && MAXVAR != 5 && MAXVAR != 10)
        MAXVAR = 2;

    SOLUCAO = FuncoesTeste[funcao].opt[0];
#ifdef COCOCIPLOT
    SOLUCAO = -INFINITO;
#endif

    /* *********     Schwefel ***************************/
    if (funcao == sch)
    {
        SOLUCAO = MAXVAR * FuncoesTeste[funcao].opt[0];
    }
    /* ***********  Michalewicz 5 ou 10 *************/
    if (funcao == mic && MAXVAR == 10)
        SOLUCAO = FuncoesTeste[funcao].opt[1];

    IniciaPop(&P, MAXPOP, MAXVAR);

    IniciaCls(&C, NUMCLUS, MAXVAR);
    C.limiar = (FuncoesTeste[funcao].sup - FuncoesTeste[funcao].inf) / (2.0F * pow(NUMCLUS, (1.0F / P.tamInd)));
    C.densid = (int)PPROMIS * NUMSELS * NUMCRU / NUMCLUS;

    GeraIndividuos(&P, MAXPOP, MAXVAR, 0, funcao);

#ifdef CONSO
    printf("\n***   (%s) por ACMO/CAP/LAC/INPE    ***", GARE);
    printf("\n      Evoluindo %s com %d vars", FuncoesTeste[funcao].nom, MAXVAR);
#endif

    erro = (double)max(P.indiv[P.melhor].fit - SOLUCAO - TAXERR, 0.0F);
    achou = (erro <= TAXERR ? TRUE : FALSE);
    bLocTot = bLocOk = 0;
    numGeracoes = 1;
    start = clock();

    while (FuncaoTeste.numAval < MAXAVA && !achou && numGeracoes < MAXGER)
    {

        numCruza = NUMCRU;

        while ((!achou) && (numCruza--) && FuncaoTeste.numAval < MAXAVA)
        {
            if (ROLETA)
            {
                RoletaPressaoSeletiva(&P, P.indiv[P.melhor].fit);
            }
            else
            {
                P.pai[0] = rand() % P.tamPop;
                P.pai[1] = (rand() % (P.pai[0] + 1)) + rand() % (P.tamPop - P.pai[0]);
            }

            if (AtualizaGrp(&C, &P) >= C.densid)
            {
                step = PSGRI * C.limiar;
                //                        step = randi(0.51, 1.59);
                PASSOS = MAXVAR;
                EstratIIIntenso(C, P, achou, bLocOk, bLocTot, step, funcao, TAXERR / 10.0F, PASSOS, ESCALA, SOLUCAO);
            }

            if (P.pai[0] != P.pai[1])
                CruzaBlend(&P, P.pai[0], P.pai[1], P.pior, BLXALFA);
            else
                P.iguais++;
            fit = funccod[funcao](P.indiv[P.pior].var, P.tamInd);
            if (fit >= P.indiv[P.melhor].fit)
            {
                mutou = MutaNaoUni(P.indiv[P.pior].var, P.tamInd, P.tamPop, numGeracoes, MUTNUNI, MUTPROB);
                if (mutou)
                {
                    P.numMuta++;
                    fit = funccod[funcao](P.indiv[P.pior].var, P.tamInd);
                }
            }
            AtualizaPop(&P, P.pior, fit, numGeracoes);
            erro = (double)max(P.indiv[P.melhor].fit - SOLUCAO - TAXERR, 0.0F);
            achou = (erro <= TAXERR ? TRUE : FALSE);
        } // while NUM CRUZA

        if (!achou)
        {
            EsfriaGrupos(C);
            C.numGrp += !C.numGrp;
            C.limiar = (FuncoesTeste[funcao].sup - FuncoesTeste[funcao].inf) / (2.0F * pow((C.numGrp), (1.0F / P.tamInd)));
            C.densid = (int)PPROMIS * NUMSELS * NUMCRU / (C.numGrp);
            nGrupos = C.numGrp;
        }
        //        if (!((numGeracoes-1)%PLOTGERA)){
        //                PlotPop (&P, NULL, C, numGeracoes, PLOTODOS);
        //        }

        if (numGeracoes > PLOTGERA && PLOT)
        {
            PlotPop(&P, NULL, C, numGeracoes, PLOTODOS);
        }

        numGeracoes++;
    } // while Geracoes

    end = clock();
    /* *******************RESULTADOS ******************** */

    P.media = P.sumFit / P.tamPop;
    P.dvpad = 0;
    for (i = 0; i < P.tamPop; i++)
        P.dvpad += pow(P.indiv[i].fit - P.media, 2);
    P.dvpad = sqrt(P.dvpad / (P.tamPop - 1));

#ifdef XLS
    // fprintf(saida,"  prob;  n;  esp;  enc;  aval; tempo; media;  dvp; ger; muta;\n");
    fprintf(saida, "    %s; %d;  %.3f;  %.3f;    %d;  %.3f;  %d;   %d;  %d; %d; %.2f; %d; %d; %.2f; %s\n",
            FuncoesTeste[funcao].nom, MAXVAR, SOLUCAO + TAXERR, P.indiv[P.melhor].fit, FuncaoTeste.numAval,
            (float)(end - start) / CLOCKS_PER_SEC, numGeracoes, P.numMuta,
            MAXAVA, MAXPOP, MUTPROB, NUMCLUS, PASSOS, PPROMIS, (ROLETA ? "ROLETA" : "ALEATA"));

#endif

#ifdef CONSO
    printf("\n\t Variaveis ...");
    for (i = 0; i < P.tamInd; i++)
    {
        printf("\n\t\t%.6f", P.indiv[P.melhor].var[i]);
    }
    // se nao houve Blocal, n�o d� erro
    bLocTot += !bLocTot;

    printf("\n\t Min = %.6f ( %.6f ??? %s!!!); \n\t#Ger = %d; \n\t#Mut?= %d; \n\t#Aval = %d;\n\t(%.4f, %.4f)\n; \n\t#BLT = %d;\n\t#BLO = %d;\n\t#GRP = %d; \n\tdFIT = %.4f;",
           P.indiv[P.melhor].fit, SOLUCAO + TAXERR, (achou ? "SUCESSO" : "FRACASSO"), numGeracoes, P.numMuta, FuncaoTeste.numAval, P.media, P.dvpad, bLocTot, bLocOk, nGrupos, (double)C.limiar);
    printf("\n\tTempo(s) = %.4f, ", (float)(end - start) / CLOCKS_PER_SEC);
#endif

#ifdef DUMP

    fprintf(saida, "\nD U M P **************");
    fprintf(saida, "\n\tPopulacao{Tpop=%d, Tind=%d, Som=%.4f, X*=%d, X =%d, #m=%d, #i=%d, (%.4f, %.4f)",
            P.tamPop, P.tamInd, P.sumFit, P.melhor, P.pior, P.numMuta, P.iguais, P.media, P.dvpad);

    for (i = 0; i < P.tamPop; i++)
    {
        fprintf(saida, "\n\tIndividuo (%d) = %.10f", i, P.indiv[i].fit);
        for (j = 0; j < P.tamInd; j++)
            fprintf(saida, "\n\t\t%.4f", P.indiv[i].var[j]);
    }
    fprintf(saida, "\n\tMelhor=%.10f,Media=%.10f, Desvio=%.10f", P.indiv[P.melhor].fit, P.media, P.dvpad);
#endif

    fclose(saida);
    strcat(nomarq, ".doc");
    if (!(saida = fopen(nomarq, "w")))
    {
        perror("");
        exit(-1);
    }

    fprintf(saida, "\nSETTINGS  **************");
    fprintf(saida, "\n\t    MAXAVA=%d,MAXPOP=%d,MUTPROB=%.2f, \n NUMCLUS=%d, PASSOS=%d, PPROMIS=%.2f ,  %d , %s",
            MAXAVA, MAXPOP, MUTPROB, NUMCLUS, PASSOS, PPROMIS, ROLETA, (ROLETA ? "ROLETA" : "ALEATA"));
    fprintf(saida, "\n************** SETTINGS  ");
}
