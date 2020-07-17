#include "../include/solver.h"

// Solve ARAP
void Solver::Solve(/* TODO */){
    //update constraints
    //TODO

    //solving iterations
    for (int k = 0; k < maxIter; k++) {
        // set up covariance matrix S
        // This is equation (5) from the paper 
        
        for (int face = 0; face < faces.rows(); face++) {
            for (int edge = 0; edge < 3; edge++) {
                int i = faces(face, map(edge, 0));
                int j = faces(face, map(edge, 1));

                v_i = vertices.row(i);
                v_j = vertices.row(j);

                v_transformed_i = vertTransformed.row(i);
                v_transformed_j = vertTransformed.row(j);

                Eigen::Vector3d e_ij = v_i - v_j;
                Eigen::Vector3d e_ij_prime = v_transformed_i - v_transformed_j;

                // this is no longer a vector
                // covarianceMatrices.block<3, 3>(i * 3, 0) += weights.coeff(i, j) * e_ij * e_ij_prime.transpose();
                covarianceMatrices[i] += weights.coeff(i, j) * e_ij * e_ij_prime.transpose();
            }
        }

        //solve for rotation

        //TODO Change code to use the row wise rotation matrix storage instead of vector of matrices
        // still .clear() ?
        rotations.clear();
        for (int i = 0; i < vertices.rows(); i++) {
            JacobiSVD<MatrixXf> svd(covarianceMatrices[i], ComputeThinU | ComputeThinV);
            Matrix3f u = svd.matrixU();
            Matrix3f v = svd.matrixV();

            // This is equation (6) from the paper
            // rotations.block<3, 3>(i * 3, 0) = v * u.transpose();
            rotations[i] = v * u.transpose();
        }

        //compute right hand side
        Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(vertices.rows(), 3);
        for (int i = 0; i < constraints.rows()) {
            if (constraints.row(i) == 1) {
                rhs.row(i) = (L * vertices.row(i).transpose()).transpose();
            }
        }

        for (int face = 0; face < faces.rows(); face++) {
            for (int edge = 0; edge < 3; edge++) {
                int i = faces(face, map(edge, 0));
                int j = faces(face, map(edge, 1));

                if (constraints.row(i) == 1)
                    continue;
                // rhs.row(i) += (weights.coeff(i, j) / 2.0 * (rotations.block<3, 3>(i * 3, 0)  
                //               + rotations.block<3, 3>(j * 3, 0))
                //               * (vertice.row(i) - vertices.row(j)).transpose()).transpose();
                rhs.row(i) += (weights.coeff(i, j) / 2.0 * (rotations[i] + rotations[j]) * (vertice.row(i) - vertices.row(j)).transpose()).transpose();
            }
        }

        //solve for updated positions
        //TODO
    }
}

// Constructor - sets vertices, faces, max number of iterations, constraints
Solver::Solver(const Eigen::MatrixXd& v, const Eigen::MatrixXd& f, int iter = 5){
        vertices = v;
        faces = f;
        maxIter = iter;
        map << 1, 2, 0, 2, 0, 1;

        // All vertices are "free" with initialization
        constraints.resize(vertices.rows(), 1);
        constraints.setZero();

        // resize covariance matrix to match number of vertices
        // TODO this is not a vector
        // covarianceMatrices.resize(vertices.rows() * 3, 3);
        // for(int i = 0; i < vertices.rows(); i++){
        //      covarianceMatrices.block<3, 3>(i * 3, 0) = Eigen::Matrix3d::Zero(); 
        // }

        covarianceMatrices.resize(vertices.rows(), Eigen::Matrix3d::Zero());

        // resize rotation matrix to match number of vertices
        // TODO this is not a vector
        // rotations.resize(vertices.rows() * 3, 3);
        // for(int i = 0; i < vertices.rows(); i++){
        //      rotaions.block<3, 3>(i * 3, 0) = Eigen::Matrix3d::Identity(); 
        // }
        rotations.resize(vertices.rows(), Eigen::Matrix3d::Identity());

        std::cout << "The matrix vertices is of size "
        << vertices.rows() << "x" << vertices.cols() << std::endl;

        std::cout << "The matrix faces is of size "
        << faces.rows() << "x" << faces.cols() << std::endl;

        std::cout << "The matrix constraints is of size "
        << constraints.rows() << "x" << constraints.cols() << std::endl;

        std::cout << "The matrix covarianceMatrices is of size "
        << covarianceMatrices.rows() << "x" << covarianceMatrices.cols() << std::endl;

        std::cout << "The matrix rotations is of size "
        << rotations.rows() << "x" << rotations.cols() << std::endl;

}

bool Solver::PrecomputeLaplaceBeltrami() {
    //TODO check wether this is the correct interpretation of the paper
    // L = weights * -1.0
    L = weights *= -1.0;
}

// Precomputes the weights used in the paper
// http://rodolphe-vaillant.fr/?e=69
void Solver::PrecomputeCotangentWeights(){
    // weights is a symetric matrix for vertex tuples
    weights.resize(vertices.rows(), vertices.rows());
    for(int i = 0; i < faces.rows(); i++){
        Eigen::Vector3d cot = ComputeFaceCotangent(i);
        for(int j = 0; j < 3; j++){

            // https://wikimedia.org/api/rest_v1/media/math/render/svg/8231849c9a676c7dc50c5ce348de162a19e411b2
            // Edge obviously exists
            weights.coeffRef(faces(i, map(j, 0)), face(i, map(j, 1))) += (cot(j) / 2.0); 
            weights.coeffRef(faces(i, map(j, 1)), face(i, map(j, 0))) += (cot(j) / 2.0); 

            // "i = j"
            weights.coeffRef(faces(i, map(j, 0)), face(i, map(j, 0))) -= (cot(j) / 2.0); 
            weights.coeffRef(faces(i, map(j, 1)), face(i, map(j, 1))) -= (cot(j) / 2.0); 

            // else 0 -> already satisfied with sparse matrix
        }
    }
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
            Eigen::Vector3d leastSquares = (deltaPLine - (rotI * deltaP)).squaredNorm();
            energy += weights.coeff(i, j) * leastSquares;
        }
    } 
    return energy;
}
