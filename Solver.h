#ifndef SOLVER_H
#define SOLVER_H

#include <Eigen/Dense>
#include <vector>

namespace Solver {

    // Computes the (generically unique) fundamental matrix for a V-Umlaut configuration:
    // seven homogeneous correspondences -- five sampled points plus two virtual ones
    // lying on lines through the first correspondence.
    Eigen::Matrix3d v_umlaut(const std::vector<Eigen::Vector3d>& Ps,
                                               const std::vector<Eigen::Vector3d>& Qs);

} // namespace Solver

#endif // SOLVER_H
