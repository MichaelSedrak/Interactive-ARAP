#include "../include/solver.h"

// Precomputes the weights used in the paper
// http://rodolphe-vaillant.fr/?e=69
void Solver::PrecomputeCotangentWeights(){
    // weights is a symetric matrix for vertex tuples
    weights.resize(vertices.rows(), vertices.rows());
    // Not expensive but equivalent to if statements
    Eigen::MatrixXi map(3, 2);
    map << 1, 2, 0, 2, 0, 1;
    for(int i = 0; i < faces.rows(); i++){
        Eigen::Vector3d cot = ComputeFaceCotangent(i);
        for(int j = 0; j < 3; j++){
            /*
             * Too expensive
            int oppositeEdgeVertex1 = -1;
            int oppositeEdgeVertex2 = -1;
            if(j == 0){
                oppositeEdgeVertex1 = 1;
                oppositeEdgeVertex2 = 2;
            }
            if(j == 1){
                oppositeEdgeVertex1 = 0;
                oppositeEdgeVertex2 = 2;
            }
            if(j == 2){
                oppositeEdgeVertex1 = 0;
                oppositeEdgeVertex2 = 1;
            }
            if(oppositeEdgeVertex1 == -1 || oppositeEdgeVertex2 == -1){
                std::cout << "Problem in weight precomputation" << std::endl;
                return false;
            }
            */

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

    // cot at v1
    cot(0) = (v2 - v1).dot(v3 - v1) / ((v2 - v1).cross(v3 - v1).norm());
    // cot at v2
    cot(1) = (v1 - v2).dot(v3 - v2) / ((v1 - v2).cross(v3 - v2).norm());
    // cot at v3
    cot(2) = (v1 - v3).dot(v2 - v3) / ((v1 - v3).cross(v2 - v3).norm());

    return cot;
}


