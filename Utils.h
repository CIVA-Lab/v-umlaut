#ifndef UTILS_H
#define UTILS_H

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <random>

#include "SceneData.h"

namespace Utils {

    // ---- Geometry ----

    // Convert pixel coordinates to unit bearing vectors via K^{-1} and normalization.
    std::vector<Eigen::Vector3d> to_normalized(const std::vector<Eigen::Vector3d>& pixels,
                                             const Eigen::Matrix3d& K);

    // Project a 3D world point to homogeneous pixel coordinates.
    Eigen::Vector3d project_point(const Eigen::Vector3d& X,
                                  const Eigen::Matrix3d& K,
                                  const Eigen::Matrix3d& R,
                                  const Eigen::Vector3d& t);

    // ---- Error metrics ----

    // Rotation error in degrees (geodesic distance on SO(3)).
    double rotation_error(const Eigen::Matrix3d& R_est, const Eigen::Matrix3d& R_gt);

    // Translation direction error in degrees.
    double translation_error(const Eigen::Vector3d& t_est, const Eigen::Vector3d& t_gt);

    // ---- Pose recovery ----

    // Decompose all E matrices via SVD, pick the (R, t) closest to ground truth.
    // Oracle disambiguation using max(rot_err, trans_err).
    std::pair<Eigen::Matrix3d, Eigen::Vector3d> recover_pose_from_essentials(
        const std::vector<Eigen::Matrix3d>& E_matrices,
        const Eigen::Matrix3d& R_gt,
        const Eigen::Vector3d& t_gt);

    // ---- Statistics ----

    struct Stats {
        double mean, min, q25, median, q75, max;
    };

    double percentile(const std::vector<double>& sorted, double p);
    // Takes v by value: it is sorted internally, so the caller's vector is left unchanged.
    Stats compute_stats(std::vector<double> v);
    void print_stats(const std::string& label, const Stats& s);

    // ---- V-Umlaut specific ----

    // Perturb a point along its line (qa-qb) in cam2 by Gaussian pixel noise.
    // Clamped to [1%, 99%] of the segment.
    void perturb_along_line(Eigen::Vector3d& q, const Eigen::Vector3d& qa,
                          const Eigen::Vector3d& qb, double sigma, std::mt19937& rng);

    // Build V-Umlaut input: first 5 projections + 2 virtual correspondences
    // (3D midpoints of (X1,X2) and (X1,X3), projected into both views).
    // dependent_noise_sigma > 0: perturb Q6/Q7 along their lines in cam2 (pixels).
    void build_v_umlaut_input(const SceneData& data,
                           std::vector<Eigen::Vector3d>& vumlaut_Ps,
                           std::vector<Eigen::Vector3d>& vumlaut_Qs,
                           double dependent_noise_sigma,
                           std::mt19937& rng);

} // namespace Utils

#endif // UTILS_H
