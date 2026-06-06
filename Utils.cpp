#include "Utils.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>

namespace Utils {

    // ---- Geometry ----

    std::vector<Eigen::Vector3d> to_normalized(const std::vector<Eigen::Vector3d>& pixels,
                                              const Eigen::Matrix3d& K) {
        Eigen::Matrix3d K_inv = K.inverse();
        std::vector<Eigen::Vector3d> normalized;
        normalized.reserve(pixels.size());
        for (const Eigen::Vector3d& p : pixels)
            normalized.push_back((K_inv * p).normalized());
        return normalized;
    }

    Eigen::Vector3d project_point(const Eigen::Vector3d& X,
                                 const Eigen::Matrix3d& K,
                                 const Eigen::Matrix3d& R,
                                 const Eigen::Vector3d& t) {
        Eigen::Vector3d X_cam = R * X + t;
        Eigen::Vector3d p = K * X_cam;
        return Eigen::Vector3d(p.x() / p.z(), p.y() / p.z(), 1.0);
    }

    // ---- Error metrics ----

    double rotation_error(const Eigen::Matrix3d& R_est, const Eigen::Matrix3d& R_gt) {
        Eigen::Matrix3d R_diff = R_est * R_gt.transpose();
        double cos_angle = std::clamp((R_diff.trace() - 1.0) / 2.0, -1.0, 1.0);
        return std::acos(cos_angle) * 180.0 / M_PI;
    }

    double translation_error(const Eigen::Vector3d& t_est, const Eigen::Vector3d& t_gt) {
        double n1 = t_est.norm(), n2 = t_gt.norm();
        if (n1 < 1e-12 || n2 < 1e-12) return 180.0;
        double cos_angle = std::clamp(t_est.dot(t_gt) / (n1 * n2), -1.0, 1.0);
        return std::acos(cos_angle) * 180.0 / M_PI;
    }

    // ---- Pose recovery ----
    // R and t decomposition follows PoseLib's implementation.
    // Then, best pose is picked using the ground-truth pose instead of cheirality checks.
    std::pair<Eigen::Matrix3d, Eigen::Vector3d> recover_pose_from_essentials(
                                        const std::vector<Eigen::Matrix3d>& E_matrices,
                                        const Eigen::Matrix3d& R_gt,
                                        const Eigen::Vector3d& t_gt) {
        if (E_matrices.empty())
            return {Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero()};

        Eigen::Matrix3d W;
        W << 0, -1, 0, 1, 0, 0, 0, 0, 1;

        double best_err = std::numeric_limits<double>::max();
        Eigen::Matrix3d best_R = Eigen::Matrix3d::Identity();
        Eigen::Vector3d best_t = Eigen::Vector3d::Zero();

        for (const Eigen::Matrix3d& E : E_matrices) {
            Eigen::JacobiSVD<Eigen::Matrix3d> USV(E, Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::Matrix3d U = USV.matrixU();
            Eigen::Matrix3d Vt = USV.matrixV().transpose();
            if (U.determinant() < 0) U.col(2) *= -1;
            if (Vt.determinant() < 0) Vt.row(2) *= -1;

            Eigen::Matrix3d R_candidates[2] = { U * W * Vt, U * W.transpose() * Vt };
            Eigen::Vector3d t_candidates[2] = { U.col(2), -U.col(2) };

            for (const Eigen::Matrix3d& R : R_candidates) {
                for (const Eigen::Vector3d& t : t_candidates) {
                    double err = std::max(rotation_error(R, R_gt), translation_error(t, t_gt));
                    if (err < best_err) {
                        best_err = err;
                        best_R = R;
                        best_t = t;
                    }
                }
            }
        }
        return {best_R, best_t};
    }

    // ---- Statistics ----

    double percentile(const std::vector<double>& sorted, double p) {
        if (sorted.empty()) return 0.0;
        double idx = p * (sorted.size() - 1);
        size_t lo = static_cast<size_t>(idx);
        size_t hi = std::min(lo + 1, sorted.size() - 1);
        double frac = idx - lo;
        return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
    }

    Stats compute_stats(std::vector<double> v) {
        Stats s{};
        if (v.empty()) return s;
        std::sort(v.begin(), v.end());
        s.mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        s.min = v.front();
        s.q25 = percentile(v, 0.25);
        s.median = percentile(v, 0.50);
        s.q75 = percentile(v, 0.75);
        s.max = v.back();
        return s;
    }

    void print_stats(const std::string& label, const Stats& s) {
        std::cout << label
            << "  mean=" << s.mean
            << "  min=" << s.min
            << "  Q25=" << s.q25
            << "  median=" << s.median
            << "  Q75=" << s.q75
            << "  max=" << s.max
            << std::endl;
    }

    // ---- V-Umlaut specific ----

    void perturb_along_line(Eigen::Vector3d& q, const Eigen::Vector3d& qa,
                          const Eigen::Vector3d& qb, double sigma, std::mt19937& rng) {
        std::normal_distribution<double> noise(0.0, sigma);
        Eigen::Vector2d a = qa.head<2>();
        Eigen::Vector2d b = qb.head<2>();
        Eigen::Vector2d dir = b - a;
        double len = dir.norm();
        if (len < 1e-12) return;
        dir /= len;

        double t_param = (q.head<2>() - a).dot(dir) / len;
        t_param += noise(rng) / len;
        t_param = std::clamp(t_param, 0.01, 0.99);

        Eigen::Vector2d q_new = a + t_param * len * dir;
        q = Eigen::Vector3d(q_new.x(), q_new.y(), 1.0);
    }

    void build_v_umlaut_input(const SceneData& data,
                           std::vector<Eigen::Vector3d>& Ps_vumlaut,
                           std::vector<Eigen::Vector3d>& Qs_vumlaut,
                           double dependent_noise_sigma,
                           std::mt19937& rng) {
        Ps_vumlaut.assign(data.Ps.begin(), data.Ps.begin() + 5);
        Qs_vumlaut.assign(data.Qs.begin(), data.Qs.begin() + 5);

        // Two virtual correspondences: 3D midpoints of (X1,X2) and (X1,X3),
        // projected into both views (epipolar-consistent by construction).
        Eigen::Vector3d X6 = 0.5 * (data.world_points[0] + data.world_points[1]);
        Eigen::Vector3d X7 = 0.5 * (data.world_points[0] + data.world_points[2]);

        Ps_vumlaut.push_back(project_point(X6, data.K1, data.R1, data.t1));
        Qs_vumlaut.push_back(project_point(X6, data.K2, data.R2, data.t2));
        Ps_vumlaut.push_back(project_point(X7, data.K1, data.R1, data.t1));
        Qs_vumlaut.push_back(project_point(X7, data.K2, data.R2, data.t2));

        if (dependent_noise_sigma > 0) {
            perturb_along_line(Qs_vumlaut[5], Qs_vumlaut[0], Qs_vumlaut[1], dependent_noise_sigma, rng);
            perturb_along_line(Qs_vumlaut[6], Qs_vumlaut[0], Qs_vumlaut[2], dependent_noise_sigma, rng);
        }
    }

} // namespace Utils