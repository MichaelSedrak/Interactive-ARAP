#pragma once

#include "Eigen/Dense"
#include "Eigen/Sparse"
// TODO Eigen Cholesky

class Solver {
    public:

        // Solver constructor
        Solver(const Eigen::MatrixXd& v, const Eigen::MatrixXd& f, int iter = 5);

        // Precomputes the weights used in the paper
        void PrecomputeCotangentWeights();

        // Compute Cholesky factorization of L
        bool PrecomputeLaplaceBeltrami();

        // Solve ARAP
        void Solve();

        // Set a new constraint
        void SetConstraint(int idx, bool fixed, const Eigen::Vector3d& pos);

        Eigen::MatrixXd GetTransformedVertices();

   private:

        // Cotangent for every vertex in triangle
        Eigen::Vector3d ComputeFaceCotangent(int face);

        // Energy function - as scalar
        double ComputeEnergyFunction();

        // max number of optimization iterations
        int maxIter;

        // constraints, i.e. free = 0, fixed = 1
        Eigen::MatrixXd constraints;

        // vertices of the mesh
        Eigen::MatrixXd vertices;

        // vertices of the mesh after ARAP is applied
        Eigen::MatrixXd vertTransformed;

        // faces of the mesh
        Eigen::MatrixXd faces;

        // rotation matrices, stored row-wise
        Eigen::MatrixXd rotations;

        // contangent weight matrix 
        Eigen::SparseMatrix<double> weights;


        // Triangle index map
        Eigen::MatrixXi map(3, 2);
}


