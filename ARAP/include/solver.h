#pragma once

#include "Eigen/Dense"
#include "Eigen/Sparse"

class Solver {
public:

    Solver() {}

    // Solver constructor
    Solver(const Eigen::MatrixXd& v, const Eigen::MatrixXi& f, int iter = 15);

    // Precomputes the weights used in the paper
    void PrecomputeCotangentWeights();

    // Energy function - as scalar
    double ComputeEnergyFunction();

    // Cotangent for every vertex in triangle
    Eigen::Vector3d ComputeFaceCotangent(int face);

    // Solve ARAP
    void Solve();

    // Set a new constraint
    void SetConstraint(int idx, bool fixed);

    // Set updated position
    void SetPosition(int idx, const Eigen::Vector3d& pos);

    // Get the solvers solution
    Eigen::MatrixXd GetTransformedVertices();

private:
    // max number of optimization iterations
    int maxIter;

    // constraints, i.e. free = 0, fixed = 1
    Eigen::MatrixXi constraints;

    // vertices of the mesh
    Eigen::MatrixXd vertices;

    // vertices of the mesh after ARAP is applied
    Eigen::MatrixXd vertTransformed;

    // faces of the mesh
    Eigen::MatrixXi faces;

    // rotation matrices, stored row-wise
    Eigen::MatrixXd rotations;

    // contangent weight matrix 
    Eigen::SparseMatrix<double> weights;

    // weight matrix with only free vertices
    Eigen::SparseMatrix<double> freeWeights;

    // Covariance matrices
    Eigen::MatrixXd covarianceMatrices;

    // Triangle index map
    Eigen::MatrixXi map;
    
    // Solver used for rhs computation
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> testSolver;

    // check if constrains are updated
    bool updated;

    // free/ fixed vertices count
    int free;
    int fixed;

    //indices of the free vertices inside the vertices array
    Eigen::VectorXi free_indices;
};


