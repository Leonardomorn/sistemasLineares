#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"
#include "Metodos.h"
#define QTDE_MATRIZES 10

int main ()
{
  // inicializa gerador de números aleatóreos
  srand(202202);
  //declarando sistemas de tamanhos  10, 30, 50, 128, 256, 512, 1000, 2000, 3000
  double testetempo;
  int tamanhoMatrizes [QTDE_MATRIZES] = {10, 30, 50, 128, 256, 512, 1000, 2000, 3000};

  //declarando vetor solução, sistema linear e os vetores para guardar os dois
  real_t *vetorSolucao;
  real_t **vetorDeVetorSolucao = malloc (sizeof(vetorSolucao) * QTDE_MATRIZES);
  SistLinear_t *SL;
  SistLinear_t **vetorSL = malloc (sizeof(*SL) *  QTDE_MATRIZES );
  real_t testevalor = 0.0;
  
  alocaVetoresSolucao(vetorDeVetorSolucao, QTDE_MATRIZES, tamanhoMatrizes);
  iniVetoresSolucao(vetorDeVetorSolucao, QTDE_MATRIZES, tamanhoMatrizes);
  alocaVetSisLin(vetorSL, QTDE_MATRIZES, tamanhoMatrizes);
  iniVetSisLin(vetorSL, QTDE_MATRIZES);
  
  prnSisLin(vetorSL[0]);

  eliminacaoGauss(vetorSL[0], vetorDeVetorSolucao[0],&testetempo);
  prnSisLin(vetorSL[0]);

  prnVetor(vetorDeVetorSolucao[0], vetorSL[0]->n);

  


  liberaVetSisLin(vetorSL, QTDE_MATRIZES);  
  liberaVetoresSolucao(vetorDeVetorSolucao, QTDE_MATRIZES, tamanhoMatrizes);

  // código do programa aqui
  
}

