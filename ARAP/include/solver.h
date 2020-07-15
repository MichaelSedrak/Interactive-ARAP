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
        map << 1, 2, 0, 2, 0, 1;

}


