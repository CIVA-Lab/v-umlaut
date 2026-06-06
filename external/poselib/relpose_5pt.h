// Copyright (c) 2020, Viktor Larsson
// All rights reserved.
// BSD-3 license. See LICENSE in this directory.

#pragma once

#include <Eigen/Dense>
#include <vector>

namespace poselib {

// Returns essential matrices (up to 10 solutions).
int relpose_5pt(const std::vector<Eigen::Vector3d> &x1, const std::vector<Eigen::Vector3d> &x2,
                std::vector<Eigen::Matrix3d> *essential_matrices);

} // namespace poselib
