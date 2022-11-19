#ifndef __SISLIN_H__
#define __SISLIN_H__


#define COEF_MAX 32.0 // Valor máximo usado para gerar valores aleatórios de
		      // coeficientes nos sistemas lineares.

// Estrutura para definiçao de um sistema linear qualquer
typedef struct {
  real_t **A; // coeficientes
  real_t *b; // termos independentes
  unsigned int n; // tamanho do SL
} SistLinear_t;

// Tipos de matrizes de coeficientes usados pela função 'inicializaSistLinear()'
typedef enum {
    generico = 0,
    hilbert,
    diagDominante
} tipoSistLinear_t;


// Alocaçao e desalocação de matrizes
SistLinear_t* alocaSisLin (unsigned int n);
void liberaSisLin (SistLinear_t *SL);
void iniSisLin (SistLinear_t *SL, tipoSistLinear_t tipo, real_t coef_max);

//Aloca e desaloca vetor de sistemas lineares e vetor de solucoes
void alocaVetSisLin(SistLinear_t** vetSL, int tam, int *vetTam);
void liberaVetSisLin(SistLinear_t **vetSL, int tam);
void iniVetSisLin(SistLinear_t **vetSL, int tam, int tipo);
void alocaVetoresSolucao(real_t **vetorDeVetorSolucao, int QTDE_MATRIZES, int *tamanhoMatrizes);
void iniVetoresSolucao(real_t **vetorDeVetorSolucao, int numeroMatrizes, int *tamanhoMatrizes);
void liberaVetoresSolucao(real_t **vetorDeVetorSolucao, int numeroMatrizes, int *tamanhoMatrizes);
void copiaSisLin (SistLinear_t* SL, SistLinear_t* SLcopia);
void copiaVetSisLin(SistLinear_t** vetSL, SistLinear_t** vetSLcopia, int quantidadeMatrizes);


// Leitura e impressão de sistemas lineares
SistLinear_t *lerSisLin ();
void prnSisLin (SistLinear_t *SL);
void prnVetor (real_t *vet, unsigned int n);
void imprimeResultado(real_t R[10][9], int quantidadeLinha);
void imprimeCabecario();
void tamanhoDasMatrizes(real_t R[10][9], int quantidadeLinha, int tamanhoMatrizes[]);


#endif // __SISLIN_H__

