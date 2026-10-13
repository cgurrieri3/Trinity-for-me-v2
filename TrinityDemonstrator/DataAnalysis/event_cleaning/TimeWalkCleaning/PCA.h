#ifndef PCA_H
#define PCA_H

#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>
#include <vector>

class PCA { 
    public:
        PCA(int features);
        void AddRow(const std::vector<double>& row);
        void ComputeCovarianceMatrix();
        void ComputeEigenDecomposition();
        TVectorD GetEigenValues() const { return eigenValues; }
        TMatrixD GetEigenVectors() const { return eigenVectors; }
        std::vector<double> GetSigmas() const;
        // Matrix dumps are per-event and noisy, so they follow the caller's diagnostics setting
        // rather than always printing.
        static void SetVerbose(bool on) { verbose = on; }

    private:
        static bool verbose;
        int numFeatures;
        std::vector<std::vector<double>> data;
        TMatrixD matrix;
        TMatrixDSym covarianceMatrix;
        TVectorD eigenValues;
        TMatrixD eigenVectors;
};
#endif


