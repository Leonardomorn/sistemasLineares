#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "sislin.h"
#include "Metodos.h"

/*!
  \brief Método da Eliminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param tTotal tempo total em milisegundos gastos pelo método

  \return código de erro. 0 em caso de sucesso.
*/
int eliminacaoGauss (SistLinear_t *SL, real_t *vetorSolucao, double *tTotal)
{
  *tTotal = timestamp();
  triangularizaSistema(SL);
  retrosubs(SL, vetorSolucao);
  *tTotal = timestamp() - *tTotal;

}

void triangularizaSistema(SistLinear_t *SL)
{
  int tam = SL->n;
  for (int i = 0; i < tam; i++)
  {
    pivoteamentoParcial(SL, i);
    for (int k = i+1; k < tam; k++)
    {
      real_t m = SL->A[k][i] / SL->A[i][i];
      SL->A[k][i] = 0.0;
      for(int j=i+1; j<tam; j++)
      {
          SL->A[k][j] -= SL->A[i][j] * m;
          SL->b[k] = SL->b[k] - SL->b[i] * m;
      }
    } 
  }  
}

void pivoteamentoParcial(SistLinear_t *SL, int linhaAtual)
{
  int tam = SL->n;
  int linhaMax = linhaAtual;
  //descobre a linha que possui o maior número na coluna da linha atual
  for (int i = linhaAtual+1; i < tam; i++)
  {
    if (SL->A[i][linhaAtual] > SL->A[linhaMax][linhaAtual])
      linhaMax = i;
  }

  if (linhaMax == linhaAtual)
    return;
  real_t *vetAux = malloc (sizeof(real_t) * SL->n);

  for (int i = linhaAtual; i < tam; i++)
  {
    vetAux[i] = SL->A[linhaAtual][i];
  }
  for (int i = linhaAtual; i < tam; i++)
  {
    SL->A[linhaAtual][i] = SL->A[linhaMax][i];
    SL->A[linhaMax][i] = vetAux[i];
  }
  
  free(vetAux);
}

void retrosubs(SistLinear_t *SL, real_t *vetorSolucao)
{
  int tam = SL->n;

  for (int i = tam-1; i >= 0; i--)
  {
    vetorSolucao[i] = SL->b[i];
    for (int j = tam - 1 ; j > i; j--)
    {
      vetorSolucao[i] = vetorSolucao[i] - SL->A[i][j] * vetorSolucao[j];
  
    }
    vetorSolucao[i] = vetorSolucao[i] / SL->A[i][i];
    
  }
  
  
}

/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
*/
real_t normaL2Residuo(SistLinear_t *SL, real_t *x)
{
  
}


/*!
  \brief Método de Gauss-Seidel

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param erro menor erro aproximado para encerrar as iterações
  \param tTotal tempo total em milisegundos gastos pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro.
  */
int gaussSeidel (SistLinear_t *SL, real_t *x, real_t erro, double *tTotal)
{
  
}


/*!
  \brief Método de Refinamento

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param erro menor erro aproximado para encerrar as iterações
  \param tTotal tempo total em milisegundos gastos pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro.
  */
int refinamento (SistLinear_t *SL, real_t *x, real_t erro, double *tTotal)
{

}


