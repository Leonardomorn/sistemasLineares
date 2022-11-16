#ifndef __METODOS_H__
#define __METODOS_H__

// Parâmetros para teste de convergência
#define MAXIT 100     // número máximo de iterações para métodos iterativos
#define ERRO 1.0e-10  // Tolerância para critérios de parada em métodos iterativos

// Calcula a normaL2 do resíduo
real_t normaL2Residuo(SistLinear_t *SL, real_t *x);

// Método da Eliminação de Gauss
int eliminacaoGauss (SistLinear_t *SL, real_t *vetorSolucao, double *tTotal);
void triangularizaSistema(SistLinear_t *SL);
void pivoteamentoParcial(SistLinear_t *SL, int linhaAtual);
void retrosubs(SistLinear_t *SL, real_t *vetorSolucao);

// Método de Refinamento
int refinamento (SistLinear_t *SL, real_t *x, real_t erro, double *tTotal);

// Método de Gauss-Seidel
int gaussSeidel (SistLinear_t *SL, real_t *x, real_t erro, double *tTotal);

#endif // __METODOS_H__

