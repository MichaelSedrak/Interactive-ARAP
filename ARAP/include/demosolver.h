#ifndef _ARAP_DEMO_ARAPSOLVER_H_
#define _ARAP_DEMO_ARAPSOLVER_H_

#include "Eigen/Dense"
#include "Eigen/Sparse"

namespace arap {
namespace demo {

enum VertexType {
  Fixed,
  Free,
};

// VertexInfo is used to build the vertex mapping from an vertex id to fixed_
// and free_ id. type_ indicates whether this vertex is free or not, and pos_
// is the position of this vertex in either free_ or fixed_.
struct VertexInfo {
  // Default constructor.
  VertexInfo()
    : type(VertexType::Free),
      pos(-1) {}
  VertexInfo(VertexType type_in, int pos_in)
    : type(type_in),
      pos(pos_in) {}
  VertexType type;
  int pos;
};

class DemoArapSolver {
 public:
  DemoArapSolver();

  // Registers data for DemoArapSolver. Once vertices and faces are set then cannot
  // be reset from outside. The only way to reset vertices is to call Solve,
  // which will run ARAP algorithm based on the current state, and update
  // vertices_.
  bool RegisterData(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces,
      const Eigen::VectorXi& fixed, int max_iteration);

  // Pre computes cotangent weight and left matrix in ARAP.
  void Precompute();

  // Solves ARAP problem. |fixed_vertices| is a # of fixed vertices by 3 matrix.
  // Each row in |fixed_vertices| represents the coordinates of a fixed vertex.
  void Solve(const Eigen::MatrixXd& fixed_vertices);

  // Gets vertices_updated_.
  const Eigen::MatrixXd& GetVertexSolution() const;
  // Gets faces_.
  const Eigen::MatrixXi& GetFaces() const { return faces_; }

  // Gets indices of the fixed vertices.
  const Eigen::VectorXi& GetFixedIndices() const { return fixed_; }

 private:
  // Computes the cotangent angle in one face indicated by |face_id|.
  Eigen::Vector3d ComputeCotangent(int face_id) const;

  // Given all the current data members, compute the energy.
  double ComputeEnergy() const;

  // vertices_ is a # of vertices by 3 matrix. Each row in vertices_ represents
  // a vertex's position.
  Eigen::MatrixXd vertices_;
  // vertices_updated_ stores the solution of Solve function. It has the same
  // dimension of vertices_, and for indices in fixed_, vertices_updated_ has
  // the same value as the given parameter in Solve function.
  Eigen::MatrixXd vertices_updated_;
  // faces_ is a # of faces by 3 matrix. Each row contains the indices of this
  // face's three vertices.
  Eigen::MatrixXi faces_;
  // fixed_ is a vector no longer than the # of vertices. It contains the
  // indices of the vertices that we want to fix during the deformation.
  Eigen::VectorXi fixed_;
  // free_ is the indices of free vertices. i.e., the union of fixed_ and
  // free_ is { 1, 2, 3, ..., vertices_.row() - 1 }.
  Eigen::VectorXi free_;
  // vertex_info_ is an array to store information for all vertices. This
  // should be updated in RegisterData function.
  std::vector<VertexInfo> vertex_info_;
  // This is a sparse matrix used to store the cotangent weight. The # of rows
  // and columns equals to the # of vertices, i.e., vertices_.rows(). For any i
  // and j, the definition of cot_weight_(i, j) is as follows:
  // If i != j and (i, j) is not an edge, cot_weight_(i, j) = 0.
  // If i != j and (i, j) is an edge, then cot_weight_(i, j) equals to:
  //   cot_weight_(i, j) = 1/2(cot alpha_ij + cot beta_ij).
  // Note that cot_weight_(i, j) = cot_weight_(j, i), or cot_weight_ is
  // symmetric.
  // If i == j, then cot_weight_(i, i) is defined as the sum of all the -w_ij.
  // The negative sign is intentionally added to help compute lb_operator_.
  Eigen::SparseMatrix<double> cot_weight_;
  // The discrete Laplace-Beltrami operator in the left hand side of the linear
  // system defined in equation (9) in the paper.
  Eigen::SparseMatrix<double> lb_operator_;
  // A vector to store rotations for all the vertices.
  std::vector<Eigen::Matrix3d> rotations_;
  // Sparse linear solver for equation (9) in the paper. Use Cholmod from
  // SuiteSparse.
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver_;
  // Max number of iterations used to solve ARAP.
  int max_iteration_;
};

inline const Eigen::MatrixXd& DemoArapSolver::GetVertexSolution() const {
  return vertices_updated_;
}

}  // namespace demo
}  // namespace arap

#endif  // _ARAP_DEMO_ARAPSOLVER_H_
