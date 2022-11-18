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
real_t normaL2Residuo(SistLinear_t *SL, real_t *residuo)
{
 int tam = SL->n;
 real_t *auxResiduo = malloc (SL->n * sizeof(real_t));

 real_t normaL2 = 0.0;
 for (int i = 0; i < tam; i++)
 {
  auxResiduo[i] = residuo[i] * residuo[i];
 }
 
 normaL2 = ABS(somaKahan(auxResiduo, tam));

 free(auxResiduo);
 return normaL2; 
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
  int tam = SL->n;
  int iteracoes = 0;
  real_t *xAux;
  zeraVetorSolucao(SL, x);
  xAux = (real_t*) malloc(tam * sizeof(real_t));


  for (iteracoes = 0; iteracoes < MAXIT; iteracoes++)
  {
    atribuiAuxSolucao(SL, xAux, x);
  
    for (int i = 0; i<tam; i++)
    {
      x[i] = SL->b[i];
      for (int j = i+1; j < tam; j++)
      {
        x[i] = x[i] - SL->A[i][j] * x[j];
      }
      for (int j = i-1; j>= 0; j--)
      {
        x[i] = x[i] - SL->A[i][j] * x[j];
      }
      x[i] = x[i] / SL->A[i][i];
    }


    if (erroMaximo(SL, xAux, x) <= erro)
    {
      return iteracoes;
      free(xAux);
    }
  }
  free(xAux);
  return iteracoes;
}

real_t erroMaximo(SistLinear_t* SL,real_t *xAux, real_t *x)
{
  real_t maxErro = 0;
  for (int i = 0; i < SL->n; i++)
  {
    if (maxErro < (abs(xAux[i] - x[i])))
      maxErro = abs(xAux[i] - x[i]);
  }
  return maxErro;
}

int atribuiAuxSolucao(SistLinear_t *SL,real_t *xAux, real_t *x)
{
  for (int i = 0; i < SL->n; i++)
  {
    xAux[i] = x[i];
  }

  return 1;
  
}
int zeraVetorSolucao(SistLinear_t *SL, real_t *x)
{
  int tam = SL->n;
  for (int i = 0; i<tam ; i++)
  {
    x[i] = 0;
  }  
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
int refinamento (SistLinear_t *SL, real_t *vetorSolucao, real_t erro, double *tTotal)
{
*tTotal = timestamp();
real_t *residuo = malloc (SL->n * sizeof (real_t));
real_t variavelDescartavel; // apenas para poder chamar a eliminacao de gauss
real_t normaL2;
int tam = SL->n;
int iteracoes = 0;
//passo 1 
//obter solucao inicial
eliminacaoGauss(SL, vetorSolucao, &variavelDescartavel);

for (int iteracoes = 0; iteracoes < tam; iteracoes++)
{
  //passo 2
  // Calcular o resíduo r = b - Ax(i) e testar critério de parada (a)
  calculaEAtribuiResiduo(SL, vetorSolucao, residuo);
  normaL2 = normaL2Residuo(SL, residuo); 
  if (normaL2 < ERRO)
  {
    free(residuo);
    return iteracoes;
  }
  
  
}

//passo 3
//


//passo 4
//

//passo 5
//

free(residuo);
}

//funcao que calcula o residuo e atribui na variável residuo
void calculaEAtribuiResiduo(SistLinear_t* SL, real_t *vetorSolucao, real_t *residuo)
{
  int tam = SL->n; 
  real_t *resAux = malloc(tam * sizeof(real_t)); 

  for (int i = 0; i < tam; i++)
  {
    residuo[i] = 0.0;
    for (int j = i; j < tam; j++)
    {
      resAux[j] = SL->A[i][j] * vetorSolucao[j];
    }
    residuo[i] = somaKahan(resAux, tam); //soma do lado esquerdo do Sislin (Ax)
    residuo[i] = SL->b[i] - residuo[i]; // r = b - Ax

  }
  
  free (resAux);

}

//algoritmo de soma para evitar muitos cancelamentos subtrativos
real_t somaKahan( real_t *dados, int tam )
{
    real_t soma = 0.0; // Prepara o acumulador
    real_t compensador = 0.0;   // compensador para a perda de bits de baixa ordem

    for(int i = 0; i < tam; i++)
    {
        real_t y = dados[i] - compensador;
        real_t t = soma + y;
        compensador = (t - soma) - y;
        soma = t;
    }
    return soma;
}