#include <platform.h>
#include <stdio.h>
// #include <xcore/port.h>
// #include <xcore/hwtimer.h>
//#include <xcore/select.h>
// #include <xcore/clock.h>
#include <vector>
#include <random>

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
  RESTYPE alpha;
  RESTYPE eig;
};

struct readout {
  std::vector<RESTYPE> weights;
  size_t nOutputs;
};

struct matrixOps {

  void mul(std::vector<RESTYPE> &A, std::vector<RESTYPE> &B, std::vector<RESTYPE> &C, int row_a, int column_a, int column_b) {

    // Data matrix
    RESTYPE* data_a = A.data();
    RESTYPE* data_b = B.data();
    RESTYPE* data_c = C.data();

    for (int i = 0; i < row_a; i++) {

      // Then we go through every column of b
      for (size_t j = 0; j < column_b; j++) {
        data_a = &A[i * column_a];
        data_b = &B[j];

        *data_c = 0; // Reset
        // And we multiply rows from a with columns of b
        for (size_t k = 0; k < column_a; k++) {
          *data_c += *data_a * *data_b;
          data_a++;
          data_b += column_b;
        }
        data_c++;
      }
    }
  }

  static void randomiseMatrix(std::vector<RESTYPE> &mat, float prob, RESTYPE low, RESTYPE high) {
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(low, high);
    std::uniform_real_distribution<> cxdis(0.f,1.f);
    for(size_t i=0; i < mat.size(); i++) {
      mat[i] = cxdis(gen) < prob  ? dis(gen) : 0;
    }
  }

  static void scalarMul(std::vector<RESTYPE> &mat, RESTYPE x) {
    for(size_t i=0; i < mat.size(); i++) {
      mat[i] *= x;
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
  float connProb, RESTYPE resLow, RESTYPE resHigh, RESTYPE inConnProb, RESTYPE inLow, RESTYPE inHigh, RESTYPE alpha) {
    res.N = nRes;
    res.N_INS = nInputs;
    res.noise=0.f;
    res.scale = 1.0f;
    res.activations.resize(nRes);
    res.inWeights.resize(nRes * nInputs);
    res.resWeights.resize(nRes * nRes);
    res.inputs.resize(nInputs);
    res.alpha = alpha;
    //init Win
    matrixOps::randomiseMatrix(res.inWeights, inConnProb, inLow, inHigh);
    //init reservoir
    matrixOps::randomiseMatrix(res.resWeights, connProb, resLow, resHigh);
    //scaling

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

    jacobi_pd::Jacobi<float, float*, float**> eigen_calc(res.N);

    eigen_calc.Diagonalize(M, evals, evecs);  //(successful if return value is > 0)

    printf("eigenvalues:  ");
    for (int i=0; i < res.N; i++)
      printf("%f, ", evals[i]);
    printf("\n");

    //scale by largest eigenvalue
    res.eig = evals[0];
    RESTYPE scaling = res.alpha / fabs(res.eig); 
    matrixOps::scalarMul(res.resWeights, scaling);

    //clear up
    delete M;
    delete evecs; 
    delete evals;
  };


  static void dump(reservoir &res) {
    printf("Reservoir:\n");
    printf("%d inputs\n", res.N_INS);
    printf("%d nodes\n", res.N);
    printf("Res weights: \n");
    matrixOps::printMatrix(res.resWeights, res.N, res.N);
  };

  static void iterateReservoir(reservoir &res, std::vector<RESTYPE> &inputs) {
    for(size_t i=0; i < res.N_INS; i++) {
      res.inputs[i] = inputs[i];
    }

  }
};

struct readoutOps {
  static void initReadout(readout &ro, reservoir &res, size_t nOutputs) {
    ro.nOutputs = nOutputs;
    ro.weights.resize(res.N * nOutputs);
  };
};

struct simulatorOps {
};

int main() {
  printf("XMOS Reservoir Test\n");
  reservoir res;
  reserviorOps::initReservoir(res, 2, 10, 0.2, -1.0, 0.4, 0.5, -1, 1, 1.1);
  // reserviorOps::dump(res);
  readout ro;
  readoutOps::initReadout(ro, res, 1);
}
