#include "Solver.h"

#include <cassert>
#include <cmath>

namespace Solver {

    // Per-view projective normalization: builds the homography H = D^{-1} A^{-1}
    // that maps the view's first four points to the canonical coordinates.
    Eigen::Matrix3d normalization_homography(const std::vector<Eigen::Vector3d>& points) {
        assert(points.size() == 7);

        Eigen::Matrix3d A;
        A.col(0) = points[0]; A.col(1) = points[1]; A.col(2) = points[2];

        // D = diag of det(A with column i replaced by points[3]).
        Eigen::Matrix3d D = Eigen::Matrix3d::Zero();
        for (int i = 0; i < 3; ++i) {
            Eigen::Matrix3d M = A;
            M.col(i) = points[3];
            D(i, i) = M.determinant();
        }

        return D.inverse() * A.inverse();
    }

    // Divides by the last coordinate so the point takes the form [u, v, 1].
    Eigen::Vector3d normalize_w(const Eigen::Vector3d& pt) {
        if (std::abs(pt(2)) > 1e-9) {
            return pt / pt(2);
        }
        return pt;
    }

    // Barycentric coordinate s of a virtual point on the canonical edge through e1
    // and the basis vector selected by other_axis (1 is e2, 2 is e3).
    double barycentric_coord(const Eigen::Vector3d& z, int other_axis) {
        double lambda = z(0) / z(other_axis);
        return lambda / (1 + lambda);
    }

    // Closed-form fundamental-matrix entries, generated with Macaulay2.
    Eigen::Matrix3d compute_f_coords(double x1, double x2, double x3,
                                    double y1, double y2, double y3,
                                    double s11, double s12, double s21, double s22) {
        double f21 = 1.0;

        // Pre-calculate powers
        double s11_2 = s11*s11; double s12_2 = s12*s12;
        double s21_2 = s21*s21; double s22_2 = s22*s22;

        // --- f12 ---
        double f12_num = (s11*s12*s21*s22 - s11*s12*s21 - s12*s21*s22 + s12*s21);
        double f12_den = (s11*s12*s21*s22 - s11*s12*s22 - s11*s21*s22 + s11*s22);

        double f12 = -f21 * f12_num / f12_den;

        // --- f20 ---
         double f20_num = (
            x2*y1*s11_2*s12*s21_2 - x1*y2*s11_2*s12*s21_2 - x2*y1*s11_2*s12*s21*s22 + x1*y2*s11_2*s12*s21*s22
            - x3*y2*s11_2*s12*s21*s22 + x2*y3*s11_2*s12*s21*s22 - x2*y1*s11_2*s21_2*s22 + x1*y2*s11_2*s21_2*s22
            + x2*y1*s11*s12*s21_2*s22 - x1*y2*s11*s12*s21_2*s22 + x3*y2*s11*s12*s21_2*s22 - x2*y3*s11*s12*s21_2*s22
            - x2*y1*s11_2*s12*s21 + x3*y2*s11_2*s12*s21 - x2*y1*s11*s12*s21_2 + 2*x1*y2*s11*s12*s21_2
            - x3*y2*s11*s12*s21_2 + x2*y1*s11_2*s12*s22 - x2*y3*s11_2*s12*s22 + 2*x2*y1*s11_2*s21*s22
            - x1*y2*s11_2*s21*s22 - x2*y3*s11_2*s21*s22 - x2*y1*s11*s12*s21*s22 - x1*y2*s11*s12*s21*s22
            + x3*y2*s11*s12*s21*s22 + x2*y3*s11*s12*s21*s22 - x1*y2*s11*s21_2*s22 + x2*y3*s11*s21_2*s22
            + x1*y2*s12*s21_2*s22 - x3*y2*s12*s21_2*s22 + x2*y1*s11*s12*s21 - x3*y2*s11*s12*s21
            - x1*y2*s12*s21_2 + x3*y2*s12*s21_2 - x2*y1*s11_2*s22 + x2*y3*s11_2*s22 + x1*y2*s11*s21*s22
            - x2*y3*s11*s21*s22
        );

        double f20_den = (
            x2*y1*s11_2*s12*s21_2 - x1*y2*s11_2*s12*s21_2 - x3*y1*s11_2*s12*s21*s22 + x1*y3*s11_2*s12*s21*s22
            - x2*y1*s11_2*s21_2*s22 + x1*y2*s11_2*s21_2*s22 + x3*y1*s11*s12*s21_2*s22 - x1*y3*s11*s12*s21_2*s22
            - 2*x2*y1*s11_2*s12*s21 + x3*y1*s11_2*s12*s21 + x1*y2*s11_2*s12*s21 - x3*y1*s11*s12*s21_2
            + x1*y2*s11*s12*s21_2 + x3*y1*s11_2*s12*s22 - x1*y3*s11_2*s12*s22 + 2*x2*y1*s11_2*s21*s22
            - x1*y2*s11_2*s21*s22 - x1*y3*s11_2*s21*s22 - x3*y1*s11*s12*s21*s22 + x1*y3*s11*s12*s21*s22
            - x1*y2*s11*s21_2*s22 + x1*y3*s11*s21_2*s22 + x2*y1*s11_2*s12 - x3*y1*s11_2*s12
            + x3*y1*s11*s12*s21 - x1*y2*s11*s12*s21 - x2*y1*s11_2*s22 + x1*y3*s11_2*s22
            + x1*y2*s11*s21*s22 - x1*y3*s11*s21*s22
        );

        double f20 = -f21 * f20_num / f20_den;

        // --- f10 and f02 ---
        double f10_f02_den = (
            x2*y1*s11_2*s12_2*s21_2*s22 - x1*y2*s11_2*s12_2*s21_2*s22 - x3*y1*s11_2*s12_2*s21*s22_2 + x1*y3*s11_2*s12_2*s21*s22_2
            - x2*y1*s11_2*s12*s21_2*s22_2 + x1*y2*s11_2*s12*s21_2*s22_2 + x3*y1*s11*s12_2*s21_2*s22_2 - x1*y3*s11*s12_2*s21_2*s22_2
            - 2*x2*y1*s11_2*s12_2*s21*s22 + x3*y1*s11_2*s12_2*s21*s22 + x1*y2*s11_2*s12_2*s21*s22 - x2*y1*s11_2*s12*s21_2*s22
            + x1*y2*s11_2*s12*s21_2*s22 - x3*y1*s11*s12_2*s21_2*s22 + x1*y2*s11*s12_2*s21_2*s22 + x3*y1*s11_2*s12_2*s22_2
            - x1*y3*s11_2*s12_2*s22_2 + 2*x2*y1*s11_2*s12*s21*s22_2 + x3*y1*s11_2*s12*s21*s22_2 - x1*y2*s11_2*s12*s21*s22_2
            - 2*x1*y3*s11_2*s12*s21*s22_2 - x3*y1*s11*s12_2*s21*s22_2 + x1*y3*s11*s12_2*s21*s22_2 + x2*y1*s11_2*s21_2*s22_2
            - x1*y2*s11_2*s21_2*s22_2 - x3*y1*s11*s12*s21_2*s22_2 - x1*y2*s11*s12*s21_2*s22_2 + 2*x1*y3*s11*s12*s21_2*s22_2
            + x2*y1*s11_2*s12_2*s22 - x3*y1*s11_2*s12_2*s22 + 2*x2*y1*s11_2*s12*s21*s22 - x3*y1*s11_2*s12*s21*s22
            - x1*y2*s11_2*s12*s21*s22 + x3*y1*s11*s12_2*s21*s22 - x1*y2*s11*s12_2*s21*s22 + x3*y1*s11*s12*s21_2*s22
            - x1*y2*s11*s12*s21_2*s22 - x2*y1*s11_2*s12*s22_2 - x3*y1*s11_2*s12*s22_2 + 2*x1*y3*s11_2*s12*s22_2
            - 2*x2*y1*s11_2*s21*s22_2 + x1*y2*s11_2*s21*s22_2 + x1*y3*s11_2*s21*s22_2 + x3*y1*s11*s12*s21*s22_2
            + x1*y2*s11*s12*s21*s22_2 - 2*x1*y3*s11*s12*s21*s22_2 + x1*y2*s11*s21_2*s22_2 - x1*y3*s11*s21_2*s22_2
            - x2*y1*s11_2*s12*s22 + x3*y1*s11_2*s12*s22 - x3*y1*s11*s12*s21*s22 + x1*y2*s11*s12*s21*s22
            + x2*y1*s11_2*s22_2 - x1*y3*s11_2*s22_2 - x1*y2*s11*s21*s22_2 + x1*y3*s11*s21*s22_2
        );
        
        // Numerator for f10
        double f10_num = (
            -x3*y1*s11_2*s12_2*s21_2*s22 + x3*y2*s11_2*s12_2*s21_2*s22 + x1*y3*s11_2*s12_2*s21_2*s22 - x2*y3*s11_2*s12_2*s21_2*s22
            + x3*y1*s11_2*s12_2*s21*s22_2 - x1*y3*s11_2*s12_2*s21*s22_2 + x3*y1*s11_2*s12*s21_2*s22_2 - x3*y2*s11_2*s12*s21_2*s22_2
            - x1*y3*s11_2*s12*s21_2*s22_2 + x2*y3*s11_2*s12*s21_2*s22_2 - x3*y1*s11*s12_2*s21_2*s22_2 + x1*y3*s11*s12_2*s21_2*s22_2
            + x3*y1*s11_2*s12_2*s21_2 - x3*y2*s11_2*s12_2*s21_2 - x3*y1*s11_2*s12_2*s21*s22 + x2*y3*s11_2*s12_2*s21*s22
            - x3*y1*s11_2*s12*s21_2*s22 + x3*y2*s11_2*s12*s21_2*s22 - x1*y3*s11_2*s12*s21_2*s22 + x2*y3*s11_2*s12*s21_2*s22
            + 3*x3*y1*s11*s12_2*s21_2*s22 - 2*x3*y2*s11*s12_2*s21_2*s22 - 2*x1*y3*s11*s12_2*s21_2*s22 + x2*y3*s11*s12_2*s21_2*s22
            - x3*y1*s11_2*s12*s21*s22_2 + 2*x1*y3*s11_2*s12*s21*s22_2 - x2*y3*s11_2*s12*s21*s22_2 - x3*y1*s11*s12_2*s21*s22_2
            + x1*y3*s11*s12_2*s21*s22_2 + x1*y3*s11_2*s21_2*s22_2 - x2*y3*s11_2*s21_2*s22_2 - x3*y1*s11*s12*s21_2*s22_2
            + 2*x3*y2*s11*s12*s21_2*s22_2 - x2*y3*s11*s12*s21_2*s22_2 + x3*y1*s12_2*s21_2*s22_2 - x1*y3*s12_2*s21_2*s22_2
            - 2*x3*y1*s11*s12_2*s21_2 + 2*x3*y2*s11*s12_2*s21_2 + x3*y1*s11_2*s12*s21*s22 - x2*y3*s11_2*s12*s21*s22
            + x3*y1*s11*s12_2*s21*s22 - x2*y3*s11*s12_2*s21*s22 + x3*y1*s11*s12*s21_2*s22 - 2*x3*y2*s11*s12*s21_2*s22
            + 2*x1*y3*s11*s12*s21_2*s22 - x2*y3*s11*s12*s21_2*s22 - 2*x3*y1*s12_2*s21_2*s22 + x3*y2*s12_2*s21_2*s22
            + x1*y3*s12_2*s21_2*s22 - x1*y3*s11_2*s21*s22_2 + x2*y3*s11_2*s21*s22_2 + x3*y1*s11*s12*s21*s22_2
            - 2*x1*y3*s11*s12*s21*s22_2 + x2*y3*s11*s12*s21*s22_2 - x1*y3*s11*s21_2*s22_2 + x2*y3*s11*s21_2*s22_2
            - x3*y2*s12*s21_2*s22_2 + x1*y3*s12*s21_2*s22_2 + x3*y1*s12_2*s21_2 - x3*y2*s12_2*s21_2
            - x3*y1*s11*s12*s21*s22 + x2*y3*s11*s12*s21*s22 + x3*y2*s12*s21_2*s22 - x1*y3*s12*s21_2*s22
            + x1*y3*s11*s21*s22_2 - x2*y3*s11*s21*s22_2
        );
        double f10 = -f21 * f10_num / f10_f02_den;

        // Numerator for f02
        double f02_num = (
            -x2*y1*s11_2*s12_2*s21_2*s22 + x1*y2*s11_2*s12_2*s21_2*s22 + x2*y1*s11_2*s12_2*s21*s22_2 - x1*y2*s11_2*s12_2*s21*s22_2
            + x3*y2*s11_2*s12_2*s21*s22_2 - x2*y3*s11_2*s12_2*s21*s22_2 + x2*y1*s11_2*s12*s21_2*s22_2 - x1*y2*s11_2*s12*s21_2*s22_2
            - x2*y1*s11*s12_2*s21_2*s22_2 + x1*y2*s11*s12_2*s21_2*s22_2 - x3*y2*s11*s12_2*s21_2*s22_2 + x2*y3*s11*s12_2*s21_2*s22_2
            + x2*y1*s11_2*s12_2*s21_2 - x1*y2*s11_2*s12_2*s21_2 + x1*y2*s11_2*s12_2*s21*s22 - 2*x3*y2*s11_2*s12_2*s21*s22
            + x2*y3*s11_2*s12_2*s21*s22 - x2*y1*s11_2*s12*s21_2*s22 + x1*y2*s11_2*s12*s21_2*s22 + 2*x2*y1*s11*s12_2*s21_2*s22
            - 3*x1*y2*s11*s12_2*s21_2*s22 + 2*x3*y2*s11*s12_2*s21_2*s22 - x2*y3*s11*s12_2*s21_2*s22 - x2*y1*s11_2*s12_2*s22_2
            + x2*y3*s11_2*s12_2*s22_2 - 2*x2*y1*s11_2*s12*s21*s22_2 + x1*y2*s11_2*s12*s21*s22_2 + x2*y3*s11_2*s12*s21*s22_2
            + x2*y1*s11*s12_2*s21*s22_2 + x1*y2*s11*s12_2*s21*s22_2 - x3*y2*s11*s12_2*s21*s22_2 - x2*y3*s11*s12_2*s21*s22_2
            + x1*y2*s11*s12*s21_2*s22_2 - x2*y3*s11*s12*s21_2*s22_2 - x1*y2*s12_2*s21_2*s22_2 + x3*y2*s12_2*s21_2*s22_2
            - x2*y1*s11_2*s12_2*s21 + x3*y2*s11_2*s12_2*s21 - x2*y1*s11*s12_2*s21_2 + 2*x1*y2*s11*s12_2*s21_2
            - x3*y2*s11*s12_2*s21_2 + x2*y1*s11_2*s12_2*s22 - x2*y3*s11_2*s12_2*s22 + 2*x2*y1*s11_2*s12*s21*s22
            - x1*y2*s11_2*s12*s21*s22 - x2*y3*s11_2*s12*s21*s22 - 2*x2*y1*s11*s12_2*s21*s22 - x1*y2*s11*s12_2*s21*s22
            + 2*x3*y2*s11*s12_2*s21*s22 + x2*y3*s11*s12_2*s21*s22 - x1*y2*s11*s12*s21_2*s22 + x2*y3*s11*s12*s21_2*s22
            + 2*x1*y2*s12_2*s21_2*s22 - 2*x3*y2*s12_2*s21_2*s22 + x2*y1*s11_2*s12*s22_2 - x2*y3*s11_2*s12*s22_2
            - x1*y2*s11*s12*s21*s22_2 + x2*y3*s11*s12*s21*s22_2 + x2*y1*s11*s12_2*s21 - x3*y2*s11*s12_2*s21
            - x1*y2*s12_2*s21_2 + x3*y2*s12_2*s21_2 - x2*y1*s11_2*s12*s22 + x2*y3*s11_2*s12*s22
            + x1*y2*s11*s12*s21*s22 - x2*y3*s11*s12*s21*s22
        );
        double f02 = -f21 * f02_num / f10_f02_den;

        // --- f01 ---
        double f01_num = (
            x3*y1*s11*s12_2*s21*s22 - x3*y2*s11*s12_2*s21*s22 - x1*y3*s11*s12_2*s21*s22 + x2*y3*s11*s12_2*s21*s22
            - x3*y1*s11*s12_2*s22_2 + x1*y3*s11*s12_2*s22_2 - x3*y1*s11*s12*s21*s22_2 + x3*y2*s11*s12*s21*s22_2
            + x1*y3*s11*s12*s21*s22_2 - x2*y3*s11*s12*s21*s22_2 + x3*y1*s12_2*s21*s22_2 - x1*y3*s12_2*s21*s22_2
            - x3*y1*s11*s12_2*s21 + x3*y2*s11*s12_2*s21 + x3*y1*s11*s12_2*s22 - x2*y3*s11*s12_2*s22
            + x3*y1*s11*s12*s21*s22 - x3*y2*s11*s12*s21*s22 + x1*y3*s11*s12*s21*s22 - x2*y3*s11*s12*s21*s22
            - 2*x3*y1*s12_2*s21*s22 + x3*y2*s12_2*s21*s22 + x1*y3*s12_2*s21*s22 + x3*y1*s11*s12*s22_2
            - 2*x1*y3*s11*s12*s22_2 + x2*y3*s11*s12*s22_2 - x1*y3*s11*s21*s22_2 + x2*y3*s11*s21*s22_2
            - x3*y2*s12*s21*s22_2 + x1*y3*s12*s21*s22_2 + x3*y1*s12_2*s21 - x3*y2*s12_2*s21
            - x3*y1*s11*s12*s22 + x2*y3*s11*s12*s22 + x3*y2*s12*s21*s22 - x1*y3*s12*s21*s22
            + x1*y3*s11*s22_2 - x2*y3*s11*s22_2
        );

        double f01_den = (
            x2*y1*s11*s12_2*s21*s22 - x1*y2*s11*s12_2*s21*s22 - x3*y1*s11*s12_2*s22_2 + x1*y3*s11*s12_2*s22_2
            - x2*y1*s11*s12*s21*s22_2 + x1*y2*s11*s12*s21*s22_2 + x3*y1*s12_2*s21*s22_2 - x1*y3*s12_2*s21*s22_2
            - x2*y1*s11*s12_2*s22 + x3*y1*s11*s12_2*s22 - x2*y1*s11*s12*s21*s22 + x1*y2*s11*s12*s21*s22
            - x3*y1*s12_2*s21*s22 + x1*y2*s12_2*s21*s22 + x2*y1*s11*s12*s22_2 + x3*y1*s11*s12*s22_2
            - 2*x1*y3*s11*s12*s22_2 + x2*y1*s11*s21*s22_2 - x1*y2*s11*s21*s22_2 - x3*y1*s12*s21*s22_2
            - x1*y2*s12*s21*s22_2 + 2*x1*y3*s12*s21*s22_2 + x2*y1*s11*s12*s22 - x3*y1*s11*s12*s22
            + x3*y1*s12*s21*s22 - x1*y2*s12*s21*s22 - x2*y1*s11*s22_2 + x1*y3*s11*s22_2
            + x1*y2*s21*s22_2 - x1*y3*s21*s22_2
        );

        double f01 = -f21 * f01_num / f01_den;

        Eigen::Matrix3d F;
        F(0, 0) = 0.0; F(1, 1) = 0.0; F(2, 2) = 0.0;
        F(0, 1) = f01; F(0, 2) = f02;
        F(1, 0) = f10; F(1, 2) = f12;
        F(2, 0) = f20; F(2, 1) = f21;

        return F;
    }

    Eigen::Matrix3d v_umlaut(const std::vector<Eigen::Vector3d>& Ps, const std::vector<Eigen::Vector3d>& Qs) {
        assert(Ps.size() == 7 && Qs.size() == 7);

        const Eigen::Matrix3d H1 = normalization_homography(Ps);
        const Eigen::Matrix3d H2 = normalization_homography(Qs);

        // 5th correspondence in the canonical frame.
        const Eigen::Vector3d z5 = normalize_w(H1 * Ps[4]);
        const Eigen::Vector3d w5 = normalize_w(H2 * Qs[4]);

        // Barycentric coordinates of the two virtual correspondences.
        const double s11 = barycentric_coord(H1 * Ps[5], 1);
        const double s21 = barycentric_coord(H2 * Qs[5], 1);
        const double s12 = barycentric_coord(H1 * Ps[6], 2);
        const double s22 = barycentric_coord(H2 * Qs[6], 2);

        const Eigen::Matrix3d F_coords =
            compute_f_coords(z5(0), z5(1), z5(2), w5(0), w5(1), w5(2), s11, s12, s21, s22);

        // Map F back to the original coordinate system.
        Eigen::Matrix3d F = H2.transpose() * F_coords * H1;

        return F;

    }

} // namespace Solver
