#include <iostream>
#include "../include/solver.h"
//#include <igl/slice.h>

// Constructor - sets vertices, faces, max number of iterations, constraints
Solver::Solver(const Eigen::MatrixXd& v, const Eigen::MatrixXi& f, int iter){
        vertices = v;
        faces = f;
        maxIter = iter;
        map.resize(3, 2);
        map << 1, 2, 2, 0, 0, 1;

        // init updated to true
        updated = true;

        // All vertices are "free" with initialization
        constraints.resize(vertices.rows(), 1);
        constraints.setZero();

        // vertTransformed initialization
        vertTransformed.resize(vertices.rows(), vertices.cols());
        vertTransformed.setZero();

        // initialize vertTransformed with vertices
        vertTransformed = vertices;

        // resize covariance matrix to match number of vertices
        covarianceMatrices.resize(vertices.rows() * 3, 3);
        covarianceMatrices.setZero();
        for(int i = 0; i < vertices.rows(); i++){
             covarianceMatrices.block<3, 3>(i * 3, 0) = Eigen::Matrix3d::Zero(); 
        }

        // resize rotation matrix to match number of vertices
        rotations.resize(vertices.rows() * 3, 3);
        rotations.setZero();
        for(int i = 0; i < vertices.rows(); i++){
             rotations.block<3, 3>(i * 3, 0) = Eigen::Matrix3d::Identity(); 
        }

        std::cout << "The matrix vertices is of size "
        << vertices.rows() << "x" << vertices.cols() << std::endl;

        std::cout << "The matrix faces is of size "
        << faces.rows() << "x" << faces.cols() << std::endl;

        std::cout << "The matrix rotations is of size "
        << rotations.rows() << "x" << rotations.cols() << std::endl;

        // Precompute function calls
        PrecomputeCotangentWeights();

        // Set free/ fixed
        free = vertices.rows();
        fixed = 0;
        /*
        SetPosition(855, Eigen::Vector3d(-0.3443436, -0.3025245, -0.4035977));
        SetConstraint(855, true);
        SetConstraint(109, true);
        SetConstraint(328, true);
        SetConstraint(518, true);
        SetConstraint(316, true);
        SetConstraint(740, true);
        SetConstraint(216, true);
        SetConstraint(206, true);
        */
        ComputeEnergyFunction();
}

// Solve ARAP
void Solver::Solve() {
    //solving iterations
    for (int k = 0; k < maxIter; k++) {
        std::cout << "Solving iteration" << std::endl;
        // set up covariance matrix S
        // This is equation (5) from the paper
        Eigen::MatrixXd covarianceMatrices;
        covarianceMatrices.resize(vertices.rows() * 3, 3);
        covarianceMatrices.setZero();

        for (int face = 0; face < faces.rows(); face++) {
            for (int edge = 0; edge < 3; edge++) {
                int i = faces(face, map(edge, 0));
                int j = faces(face, map(edge, 1));

                Eigen::Vector3d v_i = vertices.row(i);
                Eigen::Vector3d v_j = vertices.row(j);

                Eigen::Vector3d v_transformed_i = vertTransformed.row(i);
                Eigen::Vector3d v_transformed_j = vertTransformed.row(j);

                Eigen::Vector3d e_ij = v_i - v_j;
                Eigen::Vector3d e_ij_prime = v_transformed_i - v_transformed_j;

                covarianceMatrices.block<3, 3>(i * 3, 0) += weights.coeff(i, j) * e_ij * e_ij_prime.transpose();
            }
        }

        //solve for rotation
        for (int i = 0; i < vertices.rows(); i++) {
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(covarianceMatrices.block<3, 3>(i * 3, 0), Eigen::ComputeThinU | Eigen::ComputeThinV);
            Eigen::Matrix3d u = svd.matrixU();
            Eigen::Matrix3d v = svd.matrixV();

            // This is equation (6) from the paper
            rotations.block<3, 3>(i * 3, 0) = (v * u.transpose()).transpose();
        }

        // Construct weight matrix that only contains free vertices
        // only if free/ fixed has changed
        if(updated){
            // constraints have changed we need to update Laplace-Beltrami operator
            updated = false;

            //get the indices of the free vertices
            free_indices.resize(free);
            int counter = 0;
            for (int i = 0; i < vertices.rows(); i++) {
                if (constraints(i, 0) == 0) {
                    free_indices[counter++] = i;
                }
            }

            //set up the Laplace-Beltrami operator
            //igl::slice(weights, free_indices, free_indices, freeWeights);
            freeWeights = Eigen::MatrixXd(weights)(free_indices, free_indices).sparseView();
            testSolver.compute(freeWeights * -1.0);
            if (testSolver.info() != Eigen::Success) {
                std::cout << "Fail to do Cholesky factorization." << std::endl;
                return;
            }
        }

        Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(vertices.rows(), 3);
        for (int face = 0; face < faces.rows(); face++) {
            for (int edge = 0; edge < 3; edge++) {
                int i = faces(face, map(edge, 0));
                int j = faces(face, map(edge, 1));

                if (constraints(i, 0) == 1)
                    continue;
                rhs.row(i) += weights.coeff(i, j) / 2.0 * ((rotations.block<3, 3>(i * 3, 0)
                    + rotations.block<3, 3>(j * 3, 0))
                    * (vertices.row(i) - vertices.row(j)).transpose()).transpose();
                if (constraints(j, 0) == 1){
                    rhs.row(i) += weights.coeff(i, j) * vertTransformed.row(j); 
                }
            }
        }

        /*
        for (int i = 0; i < constraints.rows(); i++) {
            // drop if vertex constraint set to fixed
            if (constraints(i, 0) == 1) {
                removeRow(rhs, i);
            }
        }
        */
        Eigen::MatrixXd b(free, 3);
        b = rhs(free_indices, Eigen::all);
        Eigen::MatrixXd x;
        Eigen::VectorXd solution;
        // outer loop for x, y, z solving
        for (int k = 0; k < 3; k++) {
            solution = testSolver.solve(b.col(k));
            if (testSolver.info() != Eigen::Success) {
                std::cout << "Fail to solve the sparse linear system." << std::endl;
                return;
            }
                
            int j = 0;
            // if vertex is fixed --> continue
            for (int i = 0; i < vertices.rows(); i++) {
                if(constraints(i, 0) == 1){
                    continue;
                }
                vertTransformed(i, k) = solution(j);
                j++;
            }
        }

        ComputeEnergyFunction();
    }
}

// Precomputes the weights used in the paper
// http://rodolphe-vaillant.fr/?e=69
void Solver::PrecomputeCotangentWeights(){
    std::cout << "Precomputing cotangent weights ..." << std::endl;
    
    // weights is a symetric matrix for vertex tuples
    weights.resize(vertices.rows(), vertices.rows());
    for(int i = 0; i < faces.rows(); i++){
        Eigen::Vector3d cot = ComputeFaceCotangent(i);
        for(int j = 0; j < 3; j++){

            // https://wikimedia.org/api/rest_v1/media/math/render/svg/8231849c9a676c7dc50c5ce348de162a19e411b2
            // Edge obviously exists
            weights.coeffRef(faces(i, map(j, 0)), faces(i, map(j, 1))) += (cot(j) / 2.0); 
            weights.coeffRef(faces(i, map(j, 1)), faces(i, map(j, 0))) += (cot(j) / 2.0); 

            // "i = j"
            weights.coeffRef(faces(i, map(j, 0)), faces(i, map(j, 0))) -= (cot(j) / 2.0); 
            weights.coeffRef(faces(i, map(j, 1)), faces(i, map(j, 1))) -= (cot(j) / 2.0); 

            // else 0 -> already satisfied with sparse matrix
        }
    }
    std::cout << "Precomputing cotangent weights done" << std::endl;
}

// Cotangent for every vertex in triangle
Eigen::Vector3d Solver::ComputeFaceCotangent(int face){
    // If you consider the angle between two vectors (v and w), 
    // you can also obtain the cotangent as follow (using Eigen::Vector3d):
    // cot(theta) = cos(theta) / sin(theta) = (v . w) / |v x w|

    Eigen::Vector3d cot;
    Eigen::Vector3d v1 = vertices.row(faces(face, 0));
    Eigen::Vector3d v2 = vertices.row(faces(face, 1));
    Eigen::Vector3d v3 = vertices.row(faces(face, 2));

    // TODO / epsilon + rest
    // cot at v1
    cot(0) = (v2 - v1).dot(v3 - v1) / ((v2 - v1).cross(v3 - v1).norm());
    // cot at v2
    cot(1) = (v1 - v2).dot(v3 - v2) / ((v1 - v2).cross(v3 - v2).norm());
    // cot at v3
    cot(2) = (v1 - v3).dot(v2 - v3) / ((v1 - v3).cross(v2 - v3).norm());
    
    return cot;
}

// Energy function - as scalar
// This is the implementation of equation (3) of the paper
double Solver::ComputeEnergyFunction(){
    double energy = 0.0;
    for(int face = 0; face < faces.rows(); face++){
        for(int edge = 0; edge < 3; edge++){
            int i = faces(face, map(edge, 0));
            int j = faces(face, map(edge, 1));
            // Maybe some more debug needed
            Eigen::Vector3d deltaPLine = (vertTransformed.row(i) - vertTransformed.row(j)).transpose(); 
            Eigen::Vector3d deltaP = (vertices.row(i) - vertices.row(j)).transpose();
            Eigen::Matrix3d rotI = rotations.block<3, 3>(i * 3, 0);
            //Eigen::Vector3d leastSquares = (deltaPLine - (rotI * deltaP)).squaredNorm();
            double leastSquares = (deltaPLine - (rotI * deltaP)).squaredNorm();
            energy += weights.coeff(i, j) * leastSquares;
        }
    } 
    std::cout << "Energy: " << energy << std::endl;
    return energy;
}

Eigen::MatrixXd Solver::GetTransformedVertices() {
    return vertTransformed;
}

// Set a new constraint
void Solver::SetConstraint(int idx, bool fixed){
    // Change only if the constraint is changed
    // else fixed/ free will get messed up
    if (constraints(idx, 0) != int(fixed)){
        if(fixed){
            constraints(idx, 0) = 1; 
            fixed += 1;
            free -= 1;
        } else {
            constraints(idx, 0) = 0;
            fixed -= 1;
            free +=1;
        }
    }
    updated = true;
    std::cout << "constraint for vertex at " << idx << " set to: " << (int) fixed << std::endl;
}


// Set updated position
void Solver::SetPosition(int idx, const Eigen::Vector3d& pos){
    vertTransformed.row(idx) = pos.transpose();
    std::cout << "position for vertex at " << idx << " set to: " << 
    pos.x() << " "<< pos.y() << " " << pos.z() << std::endl;
}

/*
// https://stackoverflow.com/questions/13290395/how-to-remove-a-certain-row-or-column-while-using-eigen-library-c
void Solver::removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove){
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.bottomRows(numRows-rowToRemove);

    matrix.conservativeResize(numRows,numCols);
}

// https://stackoverflow.com/questions/13290395/how-to-remove-a-certain-row-or-column-while-using-eigen-library-c
void Solver::removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove){
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.rightCols(numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

*/
