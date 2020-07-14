#pragma once

#include "Eigen/Dense"
#include "Eigen/Sparse"
// TODO Eigen Cholesky
#include "../include/mesh.h"

class Solver {
    public:

        // Solver constructor
        Solver();
        
        // Loads the Mesh object
        // Needed for: vertices, faces, constaints
        // TODO Add constraint bool to every vertex in Mesh class
        // 0 - free
        // 1 - fixed
        bool LoadData(const Mesh& baseMesh);
        
        // Precompute 1-Ring-Neigborhood
        bool PrecomputeNeighbors();

        // Precomputes the weights used in the paper
        void PrecomputeCotangentWeights();

        // Compute Cholesky factorization of L
        bool PrecomputeLaplace-Beltrami();

        // Solve ARAP
        void Solve(/* TODO */);

        // TODO
        //
        // GetFaces()
        // GetFreeVertices()
        // GetFixedVertices()
        // GetFreeIndices()
        // GetFixedIndices()
    private:

        // Cotangent for every vertex in triangle
        Eigen::Vector3d ComputeFaceCotangent(int face);

        // Energy function - as scalar
        double ComputeEnergyFunction();

        // max number of optimization iterations
        int maxIter;

        // TODO
        // free
        // fixed
        // cholesky solver
        // rotations
        // L

        Eigen::MatrixXd vertices;
        Eigen::MatrixXd faces;
        Eigen::SparseMatrix<double> weights;
}


