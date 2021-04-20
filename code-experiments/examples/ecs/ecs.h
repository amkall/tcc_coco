#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "coco.h"

#define GARE "ECS15"
#define NOXLS
#define NODUMP
#define CONSO
#define NOCONVER
#define PLOT 0

/* PLOTS E CONVERGENCIA */
#define PLOTGERA 10
#define CONVAVAL 10

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Problema */
#define MAXGER 180000
#define NUMCRU 200
#define TAXERR 0.0001

/*Clusters */
#define AAPALFA 0.05

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

#define max(a, b) ((a) > (b) ? (a) : (b))
#define randi(x, y) (x + ((rand() % 1001 / 1000.) * (y - x)))

/**
 * A function type for evaluation functions, where the first argument is the
 * vector to be evaluated and the second argument the vector to which the
 * evaluation result is stored.
 */
typedef double (*single_evaluate_function_t)(const double *x, double *y);

struct {
  int numAval;
} FuncaoTeste;

typedef struct {
  double *var;
  double fit;
  int sel;
} Cromossomo;

typedef struct {
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

typedef struct {
  Cromossomo ponto;
  int conta;
  double alert;
  int stats;
} Centro;

typedef struct {
  Centro *grupos;
  int posGrp;
  int numGrp;
  int maxGrp;
  double limiar;
  int densid;
} Prototipos;

void initPopulation(Populacao *p, uint max_population, uint problem_dimension) {
  int i;

  p->centr.var = coco_allocate_vector(problem_dimension);

  p->indiv = (Cromossomo *)malloc(sizeof(Cromossomo) * max_population);
  if (p->indiv == NULL) {
    fprintf(stderr, "Error (memory allocation): p->indiv");
    exit(-1);
  }

  for (i = 0; i < max_population; i++) {
    p->indiv[i].var = coco_allocate_vector(problem_dimension);
  }

  p->tamPop = max_population;
  p->tamInd = problem_dimension;
  p->numMuta = 0;
  p->iguais = 0;
}

void initClusters(Prototipos *c, uint n_clusters, uint problem_dimension) {
  int i;
  c->grupos = (Centro *)malloc(sizeof(Centro) * n_clusters);

  for (i = 0; i < n_clusters; i++) {
    c->grupos[i].conta = 0;
    c->grupos[i].stats = MORTUS;
    c->grupos[i].alert = 0.0F;
    c->grupos[i].ponto.var = coco_allocate_vector(problem_dimension);
  }
  c->maxGrp = n_clusters;
  c->posGrp = 0;
  c->numGrp = 0;
  c->limiar = 0.0F;
}

float randgen(float fLlim, float fUlim) {
  float fRandomVal;

  fRandomVal = rand() % 101 / 100.; // rand entre 0 e 1

  return (fLlim + (float)(fRandomVal * (fUlim - fLlim)));
}

double DistEucl(double x1[], double x2[], int n) {
  double dist = 0.0, d;
  int i;

  dist = 0;

  for (i = 0; i < n; i++) {
    d = x1[i] - x2[i];
    dist += d * d;
  }

  return (dist);
}

int CorrigeInviavel(double *xr, double linf, double lsup) {
  if ((*xr) > lsup)
    *xr = lsup - PREBATI * ((*xr) - lsup) / ((*xr) - linf);
  else if ((*xr) < linf)
    *xr = linf + PREBATI * (linf - (*xr)) / (lsup - (*xr));
  return (1);
}

double HookExplore(single_evaluate_function_t single_evaluation_function,
                   double *y, double *xr, double fp, double dx, int n,
                   const double lower_bound, const double upper_bound) {
  int i, j;
  double fr;
  double salvo;

  for (i = 0; i < n; i++) {
    // first direction
    salvo = xr[i];
    xr[i] = salvo + dx;
    // viability
    CorrigeInviavel(&xr[i], lower_bound, upper_bound);
    // evaluate
    fr = single_evaluation_function(xr, y);
    if (fr < fp) {
      // success
      fp = fr;
    } else {
      // failure: other direction
      dx = -dx;
      xr[i] = salvo + dx;
      // viability
      CorrigeInviavel(&xr[i], lower_bound, upper_bound);
      // evaluate
      fr = single_evaluation_function(xr, y);
      if (fr < fp) {
        // success
        fp = fr;
      } else {
        // reset direction: ACMO bichado por que houve corre��o
        xr[i] = salvo;
      }
    }
  }
  return (fp);
}

void generateIndividualsPopulation(
    Populacao *p, uint max_population, uint problem_dimension, int best,
    single_evaluate_function_t single_evaluation_function, double *y,
    const double lower_bound, const double upper_bound) {
  uint i, j, pior;
  double soma, fit;

  // inicializa o centroide
  for (j = 0; j < problem_dimension; j++)
    p->centr.var[j] = 0.0F;

  for (i = 0, soma = 0, pior = 0; i < max_population; i++) {
    // gera individuo e acumula-o no centroide (nao ocorre a divisao - centroide
    // = soma)
    for (j = 0; j < problem_dimension; j++) {
      p->indiv[i].var[j] = (double)randgen(lower_bound, upper_bound);
      p->centr.var[j] += p->indiv[i].var[j];
    }
    fit = single_evaluation_function(p->indiv[i].var, y);
    p->indiv[i].fit = fit;
    p->indiv[i].sel = 0;
    if (fit > p->indiv[pior].fit)
      pior = i;
    if (fit < p->indiv[best].fit)
      best = i;
    soma += (fit);
  }
  // retira do centroide a parte do pior individuo, antecipando, pois ele sera
  // substituido na primeira atualizacao
  for (j = 0; j < problem_dimension; j++)
    p->centr.var[j] -= p->indiv[pior].var[j];

  p->sumFit = soma;
  p->melhor = best;
  p->pior = pior;
  p->media = p->sumFit / p->tamPop;
}

void rouletteSelectionPressure(Populacao *p, double melhorfit) {
  int i, pos, sel, fator;
  double z, gatilho, acum;

  sel = 0;
  while (sel < NUMSELS) {
    pos = rand() % p->tamPop;
    fator = (p->indiv[pos].sel > 3 ? 3 : p->indiv[pos].sel);
    z = (1.0F / pow((p->indiv[pos].fit - melhorfit + 1), fator));
    gatilho = z * (ROLGIRO - (ROLGIRO - 1) * z);
    acum = 0;
    i = 0;
    while (acum < gatilho && i <= ROLGIRO) {
      pos = (pos < p->tamPop - 1 ? pos + 1 : 0);
      fator = (p->indiv[pos].sel > 3 ? 3 : p->indiv[pos].sel);
      z = 1.0F / pow((p->indiv[pos].fit - melhorfit + 1), fator);
      acum += z;
      i++;
    }
    p->indiv[pos].sel++;
    p->pai[sel] = pos;
    sel++;
  }
}

int updateGrp(Prototipos *c, Populacao *p) {
  int i, j, k;
  int indice, pertice, dispice;
  double dist, menorDist;
  int maiorCont = 0;
  int numsels = NUMSELS;
  char pertence;

  for (i = 0; i < numsels; i++) {
    pertence = FALSE;
    menorDist = INFINITO;
    // insercao a priori na ultima posicao
    dispice = c->posGrp;
    for (j = 0; j < c->posGrp; j++) {
      if (c->grupos[j].stats != MORTUS) {
        dist = sqrt(DistEucl(p->indiv[p->pai[i]].var, c->grupos[j].ponto.var,
                             p->tamInd));
        if (dist < c->limiar && !pertence) {
          pertence = TRUE;
          pertice = j; // relativo a pertincia
        }
        if (dist < menorDist) {
          menorDist = dist;
          indice = j; // relativo a assimila��o
        }
      } else
        dispice = j;
    }
    // se n�o pertencer a ninguem o ultimo mortus sera usado
    // se n�o entrar nenhuma vez no if, vale a inicializacao (ultimo centro)
    if (!pertence && !(dispice >= c->maxGrp)) {
      // salva o primeiro ponto e seu fitness
      memcpy(c->grupos[dispice].ponto.var, p->indiv[p->pai[i]].var,
             p->tamInd * sizeof(double));
      c->grupos[dispice].ponto.fit = p->indiv[p->pai[i]].fit;
      // Cluster come�a gelado
      c->grupos[dispice].conta = 1;
      c->grupos[dispice].stats = GELADO;
      // aumenta posGrp se inseriu na ultima posicao
      c->posGrp += (dispice == c->posGrp);
      // independente disso, aumenta o num
      c->numGrp++;
    } else { // ou pertence ou nao cabe mais clusters, entao assimila o cluster
             // mais proximo
      if (menorDist > TAXERR) {
//--------------------------------------------------
// protege de assimila��o pontos muito pr�ximos ao centro
// TRES TIPOS DE ASSIMILA��O
//--------------------------------------------------
#ifdef ASSIMPLES
        // ASSIMPLES usa um AAPALFA para gerar um novo centros na reta entre o
        // antigo e o ponto assimilado N�o precisa avaliar o novo centro, isso �
        // feito antes de Hooke Se for executar sem BLOCAL, mas n�o deixar de
        // avalia-lo em IIIntenso

        {
          int flag = 0;
          if (rand() % 1 == 0 && PLOT) {
            flag = 1;
          }
          for (k = 0; k < p->tamInd; k++) {
            c->grupos[indice].ponto.var[k] +=
                AAPALFA *
                (p->indiv[p->pai[i]].var[k] - c->grupos[indice].ponto.var[k]);
          }

        } // auto
#endif
#ifdef ASRECOMBI
        // ASRECOMBI usa n - AAPALFA's para gerar um novo centros no hiperplano
        // entre o antigo e o ponto assimilado N�o precisa avaliar o novo
        // centro, isso � feito antes de Hooke Se for executar sem BLOCAL, mas
        // n�o deixar de avalia-lo em IIIntenso

        {
          double a, b, r;
          int flag = 0;
          a = -AAPALFA;
          b = 1 + AAPALFA;

          if (rand() % 1 == 0 && PLOT) {
            PlotPop(p, p->indiv[p->pai[i]].var, *c, indice, indice);
            flag = 1;
          }

          for (k = 0; k < p->tamInd; k++) {
            r = a + (rand() % 101 / 100.) * (b - a);
            c->grupos[indice].ponto.var[k] +=
                r *
                (p->indiv[p->pai[i]].var[k] - c->grupos[indice].ponto.var[k]);
            CorrigeInviavel(&(c->grupos[indice].ponto.var[k]),
                            FuncoesTeste[funcao].inf, FuncoesTeste[funcao].sup);
          }
          if (flag && PLOT)
            PlotPop(p, p->indiv[p->pai[i]].var, *c, indice, indice);

        } // auto
#endif

#ifdef ASCAMINHO
        // ASCAMINHO usa v�rios AAPALFA para gerar v�rios novos centros entre o
        // antigo e o ponto assimilado Ele j� avalia o novo centro. Em Hooke n�o
        // precisa avaliar de novo

        {
          double *aux, fitaux, *sau, fitsau;
          int namost, flag = 0;

          aux = (double *)malloc(p->tamInd * sizeof(double));
          sau = (double *)malloc(p->tamInd * sizeof(double));
          memcpy(aux, c->grupos[indice].ponto.var, p->tamInd * sizeof(double));
          fitsau = c->grupos[indice].ponto.fit;
          namost = (int)floor(1.0F / AAPALFA);
          while (--namost) {
            if (rand() % 1 == 0 && PLOT) {
              PlotPop(p, p->indiv[p->pai[i]].var, *c, indice, indice);
              flag = 1;
            }
            for (k = 0; k < p->tamInd; k++) {
              aux[k] += AAPALFA * (p->indiv[p->pai[i]].var[k] -
                                   c->grupos[indice].ponto.var[k]);
            }
            fitaux = funccod[funcao](aux, p->tamInd);
            if (fitaux < fitsau) {
              memcpy(sau, aux, p->tamInd * sizeof(double));
              fitsau = fitaux;
              if (flag && PLOT) {
                PlotPop(p, p->indiv[p->pai[i]].var, *c, indice, indice);
                flag = 0;
              }
            }
          }
          memcpy(aux, p->indiv[p->pai[i]].var, p->tamInd * sizeof(double));
          do {
            if (rand() % 1 == 0 && PLOT) {
              PlotPop(p, p->indiv[p->pai[i]].var, *c, indice, indice);
              flag = 1;
            }

            for (k = 0; k < p->tamInd; k++) {
              aux[k] += AAPALFA * (p->indiv[p->pai[i]].var[k] -
                                   c->grupos[indice].ponto.var[k]);
              CorrigeInviavel(&aux[k], FuncoesTeste[funcao].inf,
                              FuncoesTeste[funcao].sup);
            }
            fitaux = funccod[funcao](aux, p->tamInd);
            if (fitaux < fitsau) {
              memcpy(sau, aux, p->tamInd * sizeof(double));
              fitsau = fitaux;
              if (flag && PLOT) {
                PlotPop(p, p->indiv[p->pai[i]].var, *c, indice, indice);
                flag = 0;
              }
            } else
              break;
          } while (1);
          if (fitsau < c->grupos[indice].ponto.fit) {
            memcpy(c->grupos[indice].ponto.var, sau,
                   p->tamInd * sizeof(double));
            c->grupos[indice].ponto.fit = fitsau;
          }
          free(aux);
          free(sau);
        } // auto
#endif

      }               // if
      if (pertence) { // pertinente � o cluster proximo (limiar) mais antigo
                      // (primeiro a ser encontrado)
        c->grupos[pertice].conta++;
        // guarda a maior contagem
        if (c->grupos[pertice].conta > maiorCont)
          maiorCont = c->grupos[pertice].conta;

        // independente de ser o mesmo ou n�o, o cluster � esquentado
        c->grupos[pertice].stats = QUENTE;

        // se houver dispice do lado esquerdo de pos, pos � decrementado
        c->posGrp -= (c->posGrp == dispice + 1);
      } // if
    }   // else (pertence)
  }
  return (maiorCont);
}

double HookeJeeves(coco_problem_t *problem,
                   single_evaluate_function_t single_evaluation_function,
                   double *y, double xc[], double fc, int n, double epsilon,
                   int passos, double scala, double step,
                   uint problem_dimension, const double lower_bound,
                   const double upper_bound) {
  double dx, err, fp, inif;
  static double mdx = 1.0F, melf = 0.0F;
  static int cont = 100;
  int i, m;
  char reduz;

  double *xr = (double *)NULL;
  double *xp = (double *)NULL;

  xp = coco_allocate_vector(n);
  xr = coco_allocate_vector(n);

  inif = fc;
  if (cont > 0) {
    dx = step;
    cont--;
  } else {
    dx = step = mdx;
  }

  m = 0;

  while (m <= passos && FuncaoTeste.numAval < problem_dimension) {
    // Assign base point
    fp = fc;
    memcpy(xr, xc, n * sizeof(double));
    fp = HookExplore(single_evaluation_function, y, xr, fp, dx, n, lower_bound,
                     upper_bound);
    // if it doesnt get into; it must be reduced
    reduz = TRUE;
    while (fp < fc && fabs(fp - fc) > epsilon && test_solution_found(problem) &&
           FuncaoTeste.numAval < problem_dimension) {
      reduz = FALSE;
      // set base point
      fc = fp;
      memcpy(xp, xc, n * sizeof(double));
      memcpy(xc, xr, n * sizeof(double));
      for (i = 0; i < n; i++) {
        xr[i] = xr[i] + (xr[i] - xp[i]);
        CorrigeInviavel(&xr[i], lower_bound, upper_bound);
      }
      // TODO: pq sobrescreve fp?
      fp = single_evaluation_function(xr, y);
      fp = HookExplore(single_evaluation_function, y, xr, fp, dx, n,
                       lower_bound, upper_bound);
    }
    if (reduz && test_solution_found(problem)) {
      dx = scala * dx;
    }
    // difere do original -- sempre incrementa m
    m++;
  }
  if (inif - fc > melf) {
    mdx = step;
    melf = inif - fc;
  }
  free(xr);
  free(xp);
  return (fc);
}

void EstratIIIntenso(Prototipos *C, Populacao *P, int *achou, int *bLocOk,
                     int *bLocTot, double step,
                     single_evaluate_function_t single_evaluation_function,
                     coco_problem_t *problem, double *y, double taxerr,
                     int passos, double escala, uint problem_dimension,
                     const double lower_bound, const double upper_bound) {
  int index;
  double erro, fitant, fitdep;

  int j;

  for (index = 0; (!achou) && (index < C->posGrp); index++) {

    if (C->grupos[index].conta >= C->densid) {
#ifdef ASCAMINHO
      // N�o precisa avaliar o centro
      fitant = C.grupos[index].ponto.fit;
#else
      // precisa avaliar o centro
      fitant = single_evaluation_function(C->grupos[index].ponto.var, y);
#endif
      fitdep = HookeJeeves(problem, single_evaluation_function, y,
                           C->grupos[index].ponto.var, fitant, P->tamInd,
                           taxerr, passos, escala, step, problem_dimension,
                           lower_bound, upper_bound);
      if (fitdep < P->indiv[P->melhor].fit) {
        memcpy(P->indiv[P->melhor].var, C->grupos[index].ponto.var,
               P->tamInd * sizeof(double));
        P->indiv[P->melhor].fit = fitdep;
        *achou = test_solution_found(problem);
      }
      bLocTot++;
      if (fitdep < fitant)
        bLocOk++;
      // reinicia so se foi verificado
      C->grupos[index].conta = 1;
    }
  }
  return;
}

void CruzaBlend(Populacao *p, int pai, int mae, int filho, float alfa,
                const double lower_bound, const double upper_bound) {
  double a, b, r;
  int i;

  a = -alfa;
  b = 1 + alfa;

  for (i = 0; i < p->tamInd; i++) {
    r = a + (rand() % 101 / 100.) * (b - a);
    // gera filho
    p->indiv[filho].var[i] = p->indiv[pai].var[i] +
                             r * (p->indiv[mae].var[i] - p->indiv[pai].var[i]);
    // rebate se invi�vel
    CorrigeInviavel(&(p->indiv[filho].var[i]), lower_bound, upper_bound);
  }
}

int MutaNaoUni(double *indiv, int tamind, int tampop, int ger, int expo,
               float pmut, const double lower_bound, const double upper_bound) {
  int mutou = FALSE;
  float fRandVal, fFactor;
  float fNewt, fNewT;
  int iExponent, iIndex;

  for (iIndex = 0; iIndex < tamind; iIndex++) {
    if (rand() % 100 < pmut) {
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
        fFactor = fFactor * (indiv[iIndex] - lower_bound);
        indiv[iIndex] = indiv[iIndex] - fFactor;
      } else {
        fNewt = ger;
        fNewT = MAXGER;
        fRandVal = (rand() % 101 / 100.);
        fFactor = pow((1.0F - (fNewt / fNewT)), expo) * fRandVal;
        if (fFactor < TAXERR / 10.0F)
          fFactor = TAXERR / 10.0F;
        fFactor = fFactor * (upper_bound - indiv[iIndex]);
        indiv[iIndex] = indiv[iIndex] + fFactor;
      }
    }
  }

  return mutou;
}

void updatePopulation(Populacao *p, int pos, double fit, int ger) {
  int i, j, salto, maxIt;

  p->sumFit -= (p->indiv[pos].fit);
  p->sumFit += (fit);
  p->indiv[pos].fit = fit;
  p->indiv[pos].sel = 0;

  p->media = p->sumFit / p->tamPop;

  if (fit < p->indiv[p->melhor].fit) {
    p->melhor = pos;
    p->gerMelhor = ger;
  }

  maxIt = (int)ceil(PPIOR * p->tamPop);
  salto = (int)ceil(maxIt / 3.0F);

  for (i = 0; i < maxIt; i++) {
    j = rand() % p->tamPop;
    if (p->indiv[j].fit > p->indiv[p->pior].fit) {
      p->pior = j;
      i += salto;
    }
  }

  for (j = 0; j < p->tamInd; j++)
    p->centr.var[j] =
        p->centr.var[j] + p->indiv[pos].var[j] - p->indiv[p->pior].var[j];
}

void groupsCoolDown(Prototipos *C) {
  int index, j;

  for (index = 0; index < C->posGrp; index++) {
    if (C->grupos[index].stats > MORTUS) {
      C->grupos[index].stats--;
      C->numGrp -= (C->grupos[index].stats == MORTUS);
    }
  }
  return;
}

int test_solution_found(coco_problem_t *problem) {
  if ((coco_problem_final_target_hit(problem) &&
       coco_problem_get_number_of_constraints(problem) == 0)) {
    return TRUE;
  } else {
    return FALSE;
  }
}
