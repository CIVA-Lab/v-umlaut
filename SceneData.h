#ifndef SCENEDATA_H
#define SCENEDATA_H

#include <Eigen/Dense>
#include <vector>

struct SceneData {

    std::vector<Eigen::Vector3d> Ps; // 2D points Image 1 (Homogeneous)
    std::vector<Eigen::Vector3d> Qs; // 2D points Image 2 (Homogeneous)

    // 3D world points used to generate the projections
    std::vector<Eigen::Vector3d> world_points;

    Eigen::Matrix3d K1, K2;
    Eigen::Matrix3d R1, R2;
    Eigen::Vector3d t1, t2;
};

#endif // SCENEDATA_H
