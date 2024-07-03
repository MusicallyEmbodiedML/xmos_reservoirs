// #define THIS_RUNS_ON_XMOS 
// #define EIGEN_DONT_ALIGN
#include <platform.h>
#include <stdio.h>
// #include <xcore/port.h>
// #include <xcore/hwtimer.h>
//#include <xcore/select.h>
// #include <xcore/clock.h>
#include <vector>
#include <random>
// using namespace std;
// #include "eigen.h"

// extern "C" {
// #include "EmbeddedLapack/src/LinearAlgebra/declareFunctions.h"
// #include "CControl/src/CControl/Sources/LinearAlgebra/linearalgebra.h"

// }

#include "jacobi_pd.hpp"
#include "matrix_alloc_jpd.hpp"

using RESTYPE = float;

struct reservoir {
  std::vector<RESTYPE> activations;
  std::vector<RESTYPE> inWeights;
  std::vector<RESTYPE> resWeights;
  std::vector<RESTYPE> inputs;
  RESTYPE noise;
  RESTYPE scale;
  size_t N;
  size_t N_INS;
};

struct matrixOps {
  static void randomiseMatrix(std::vector<RESTYPE> &mat, float connectionProb, RESTYPE low, RESTYPE high) {
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(low, high);
    std::uniform_real_distribution<> cxdis(0.f,1.f);
    for(size_t i=0; i < mat.size(); i++) {
      mat[i] = cxdis(gen) < connectionProb  ? dis(gen) : 0;
    }
  }

  template<typename T=RESTYPE>
  static void printMatrix(std::vector<T> &mat, size_t rows, size_t cols) {
    for(size_t i=0; i < rows; i++) {
      for(size_t j=0; j < cols; j++) {
        printf("%f\t", mat[(i*cols) + j]);
      }
      printf("\n");
    }
  }

  static std::vector<double> floatToDouble(std::vector<float> &mat) {
    std::vector<double> out(mat.size());
    for(size_t i=0; i < mat.size(); i++) {
      out[i] = static_cast<double>(mat[i]);
    }
    return std::move(out);
  }
};

struct reserviorOps {
  static void initReservoir(reservoir &res, size_t nInputs, size_t nRes,
  float connProb, RESTYPE resLow, RESTYPE resHigh) {
    res.N = nRes;
    res.N_INS = nInputs;
    res.noise=0.f;
    res.scale = 1.0f;
    res.activations.resize(nRes);
    res.inWeights.resize(nRes * nInputs);
    res.resWeights.resize(nRes * nRes);
    res.inputs.resize(nInputs);
    //init reservoir
    matrixOps::randomiseMatrix(res.resWeights, connProb, resLow, resHigh);
    //scaling
    // std::vector<float> Ereal(res.N); // Eigenvalues real
    // std::vector<float> Eimag(res.N); // Eigenvalues imag part
    // std::vector<float> Vreal(res.N * res.N);
    // std::vector<float> Vimag(res.N * res.N);

    // int n = 3;       // Matrix size
    float **M;      // A symmetric n x n matrix you want to diagonalize
    float *evals;   // Store the eigenvalues here.
    float **evecs;  // Store the eigenvectors here.

    matrix_alloc_jpd::Alloc2D(res.N,res.N, &M);
    size_t idx=0;
    for(size_t i=0; i < res.N; i++) {
      for(size_t j=0; j < res.N; j++) {
        M[i][j] = res.resWeights[idx];
        idx++;
      }
    }
    matrix_alloc_jpd::Alloc2D(res.N,res.N, &evecs);
    evals = new float[res.N];


    // M[0][0] = 2.0; M[0][1] = 1.0; M[0][2] = 1.0;
    // M[1][0] = 1.0; M[1][1] = 2.0; M[1][2] =-1.0;  //Note: The matrix
    // M[2][0] = 1.0; M[2][1] =-1.0; M[2][2] = 2.0;  //must be symmetric.

    jacobi_pd::Jacobi<float, float*, float**> eigen_calc(res.N);

    eigen_calc.Diagonalize(M, evals, evecs);  //(successful if return value is > 0)

// If you have many matrices to diagonalize, you can re-use "eigen_calc". (This
// is more efficient than creating a new "Jacobi" class instance for each use.)

    printf("eigenvalues:  ");
    for (int i=0; i < res.N; i++)
      printf("%f, ", evals[i]);
    printf("\n");

    delete M;
    delete evecs; 
    delete evals;
    // cout << endl;
    // for (int i=0; i < n; i++) {
    //   cout << "eigenvector" <<i+1<< ": ";
    //   for (int j=0; j < n; j++)
    //     cout << evecs[i][j] << " ";
    //   cout << endl;
    // }
    // std::vector<double> tmp = matrixOps::floatToDouble(res.resWeights);
    
    // Solve
    // eig(res.resWeights.data(),Ereal.data(),Eimag.data(),Vreal.data(),
    // Vimag.data(),res.N);
    
    // matrixOps::printMatrix<float>(Ereal, 1, res.N);
  };


  static void dump(reservoir &res) {
    printf("Reservoir:\n");
    printf("%d inputs\n", res.N_INS);
    printf("%d nodes\n", res.N);
    printf("Res weights: \n");
    matrixOps::printMatrix(res.resWeights, res.N, res.N);
  };
};

int main() {
  printf("XMOS Reservoir Test\n");
  reservoir res;
  reserviorOps::initReservoir(res, 2, 10, 0.2, -1.0, 0.4);
  reserviorOps::dump(res);
}
