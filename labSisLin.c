#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"
#include "Metodos.h"
#define QTDE_MATRIZES 10

#define T_EGP_POS   1
#define NR_EGP_POS  2
#define T_GS_POS    3
#define IT_GS_POS   4
#define NR_GS_POS   5
#define T_REF_POS   6
#define IT_REF_POS  7
#define NR_REF_POS  8

int main ()
{
  // inicializa gerador de números aleatóreos
  srand(202202);
  //declarando sistemas de tamanhos  10, 30, 50, 128, 256, 512, 1000, 2000, 3000
  double tempo; real_t norma_egp, norma_ref, norma_gs;
  int tamanhoMatrizes [QTDE_MATRIZES] = {10, 30, 50, 128, 256, 512, 1000, 2000, 3000};
  real_t Resultados[QTDE_MATRIZES][9];
  //declarando vetor solução, sistema linear e os vetores para guardar os dois
  real_t *vetorSolucao;
  real_t **vetorDeVetorSolucao = malloc (sizeof(vetorSolucao) * QTDE_MATRIZES);
  SistLinear_t *SL;
  SistLinear_t **vetorSL = malloc (sizeof(*SL) *  QTDE_MATRIZES );
  SistLinear_t **vetorSLcopias = malloc (sizeof(*SL) *  QTDE_MATRIZES );

  real_t testevalor = 0.0;
  
  alocaVetoresSolucao(vetorDeVetorSolucao, QTDE_MATRIZES, tamanhoMatrizes);
  iniVetoresSolucao(vetorDeVetorSolucao, QTDE_MATRIZES, tamanhoMatrizes);
  alocaVetSisLin(vetorSL, QTDE_MATRIZES, tamanhoMatrizes);
  alocaVetSisLin(vetorSLcopias, QTDE_MATRIZES, tamanhoMatrizes);
  for (int tipo = 0; tipo< 3; tipo++)
  {
    iniVetSisLin(vetorSL, QTDE_MATRIZES, tipo);
    copiaVetSisLin(vetorSL, vetorSLcopias, QTDE_MATRIZES);
    for (int i = 0; i < QTDE_MATRIZES; i++)
    {
      printf("Computando matriz %d...\n", i );
      eliminacaoGauss(vetorSL[i], vetorDeVetorSolucao[i], &tempo );
      Resultados[i][T_EGP_POS] = tempo;
      copiaSisLin(vetorSLcopias[i], vetorSL[i]);
      Resultados[i][IT_GS_POS] = gaussSeidel(vetorSL[i], vetorDeVetorSolucao[i], ERRO, &tempo, &norma_gs);
      Resultados[i][T_GS_POS] = tempo;
      Resultados[i][NR_GS_POS] = norma_gs;
      Resultados[i][IT_REF_POS] = refinamento(vetorSL[i], vetorDeVetorSolucao[i], ERRO, &tempo, &norma_egp, &norma_ref);
      Resultados[i][NR_EGP_POS] = norma_egp;
      Resultados[i][T_REF_POS] = tempo;
      Resultados[i][NR_REF_POS] = norma_ref;
    
    }
    switch (tipo)
    {
    case 0:
      printf("|||GENÉRICO|||\n");
      tamanhoDasMatrizes(Resultados, QTDE_MATRIZES, tamanhoMatrizes);
      imprimeCabecario();
      imprimeResultado(Resultados, QTDE_MATRIZES);
      printf("\n\n\n\n\n");
      break;
    case 1:
      printf("|||Hilbert|||\n");
      tamanhoDasMatrizes(Resultados, QTDE_MATRIZES, tamanhoMatrizes);
      imprimeCabecario();
      imprimeResultado(Resultados, QTDE_MATRIZES);
      printf("\n\n\n\n\n");
      break;
    case 2:
      printf("|||DIAG DOMINANTE|||\n");
      tamanhoDasMatrizes(Resultados, QTDE_MATRIZES, tamanhoMatrizes);
      imprimeCabecario();
      imprimeResultado(Resultados, QTDE_MATRIZES);
      break;
    
    default:
      break;
    }
  }
  


  tamanhoDasMatrizes(Resultados, QTDE_MATRIZES, tamanhoMatrizes);
  imprimeCabecario();
  imprimeResultado(Resultados, QTDE_MATRIZES);


  

  liberaVetSisLin(vetorSL, QTDE_MATRIZES);  
  liberaVetoresSolucao(vetorDeVetorSolucao, QTDE_MATRIZES, tamanhoMatrizes);


}

