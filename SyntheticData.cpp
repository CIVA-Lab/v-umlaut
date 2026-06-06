#include "SyntheticData.h"
#include "Utils.h"

#include <iostream>
#include <cmath>

// World coordinate system (ENU):
//   +X = East,  +Y = North,  +Z = Up
//
// Camera coordinate system (OpenCV):
//   +X = image right,  +Y = image down,  +Z = forward
//
// Extrinsic convention:
//   X_cam = R * X_world + t,  with  t = -R * C
// where C is the camera center in world coordinates.
//
// Since R maps world coordinates into camera coordinates, its rows are the
// camera-frame axes used by the world-to-camera transform:
//   row 0 = camera +X (right)
//   row 1 = camera +Y (cam_y)
//   row 2 = camera +Z (forward)
static Eigen::Matrix3d look_at(const Eigen::Vector3d& from,
                              const Eigen::Vector3d& to) {
    Eigen::Vector3d forward = (to - from).normalized();

    // world_up is a roll hint, fall back to +Y if forward is near-parallel to +Z.
    Eigen::Vector3d world_up(0.0, 0.0, 1.0);
    if (std::abs(forward.dot(world_up)) > 0.99)
        world_up = Eigen::Vector3d(0.0, 1.0, 0.0);

    // Camera +X. For a level camera looking north in ENU, this gives East.
    Eigen::Vector3d right = forward.cross(world_up).normalized();

    // Camera +Y. Following the OpenCV convention, this corresponds to image-down.
    Eigen::Vector3d cam_y = forward.cross(right).normalized();

    Eigen::Matrix3d R;
    R.row(0) = right.transpose();
    R.row(1) = cam_y.transpose();
    R.row(2) = forward.transpose();
    return R;
}

// Check if a 3D point projects inside the image bounds for a given camera.
static bool is_visible(const Eigen::Vector3d& X,
                      const Eigen::Matrix3d& K,
                      const Eigen::Matrix3d& R,
                      const Eigen::Vector3d& t,
                      int width, int height) {
    if ((R * X + t).z() <= 0.0) return false; // must be in front of the camera

    Eigen::Vector3d p = Utils::project_point(X, K, R, t);
    return p.x() >= 0 && p.x() < width && p.y() >= 0 && p.y() < height;
}

// Add independent Gaussian noise to the x,y pixel coordinates of projected points.
static void add_pixel_noise(SceneData& data, double sigma, std::mt19937& rng) {
    std::normal_distribution<double> noise(0.0, sigma);
    for (Eigen::Vector3d& p : data.Ps)
        p = Eigen::Vector3d(p.x() + noise(rng), p.y() + noise(rng), 1.0);
    for (Eigen::Vector3d& q : data.Qs)
        q = Eigen::Vector3d(q.x() + noise(rng), q.y() + noise(rng), 1.0);
}

namespace SyntheticData {

    // Scene/camera parameters per mode (the two configurations used in the paper).
    struct SceneConfig {
        Eigen::Matrix3d K;
        Eigen::Vector3d box_min, box_max, scene_center;
        double cylinder_radius, cam_z_low, cam_z_high;
        int image_width, image_height;
    };

    static SceneConfig scene_config(Mode mode) {
        SceneConfig c;
        if (mode == Mode::WAMI) {
            // WAMI: long focal length, flat ground scene, cameras high above inside a cylinder.
            c.K << 17000, 0, 3300,
                   0, 17000, 2200,
                   0, 0, 1;
            c.box_min = Eigen::Vector3d(-1000, -1000, 0);
            c.box_max = Eigen::Vector3d(1000, 1000, 100);
            c.scene_center = Eigen::Vector3d(0, 0, 50);
            c.cylinder_radius = 1500.0;
            c.cam_z_low = 1500.0;
            c.cam_z_high = 2000.0;
            c.image_width = 6600;
            c.image_height = 4400;
        } else {
            // General: VGA intrinsics, small scene, cameras inside an encompassing cylinder.
            c.K << 800, 0, 320,
                   0, 800, 240,
                   0, 0, 1;
            c.box_min = Eigen::Vector3d(-2, -2, 0);
            c.box_max = Eigen::Vector3d(2, 2, 4);
            c.scene_center = Eigen::Vector3d(0, 0, 2);
            c.cylinder_radius = 8.0;
            c.cam_z_low = 0.0;
            c.cam_z_high = 8.0;
            c.image_width = 640;
            c.image_height = 480;
        }
        return c;
    }

    void generate(SceneData& data, std::mt19937& rng, Mode mode, double independent_noise_sigma) {
        data.Ps.clear();
        data.Qs.clear();
        data.world_points.clear();

        const SceneConfig cfg = scene_config(mode);
        const Eigen::Matrix3d& K = cfg.K; // aliased: used for both views and every projection

        data.K1 = K;
        data.K2 = K;

        // Distributions for 3D point sampling
        std::uniform_real_distribution<double> x_distrib(cfg.box_min.x(), cfg.box_max.x());
        std::uniform_real_distribution<double> y_distrib(cfg.box_min.y(), cfg.box_max.y());
        std::uniform_real_distribution<double> z_distrib(cfg.box_min.z(), cfg.box_max.z());

        // Sample two cameras uniformly inside a solid cylinder.
        std::uniform_real_distribution<double> angle_distrib(0.0, 2.0 * M_PI);
        std::uniform_real_distribution<double> unit_distrib(0.0, 1.0);
        std::uniform_real_distribution<double> cam_z_distrib(cfg.cam_z_low, cfg.cam_z_high);

        // sqrt corrects for the disk area element so points are uniform.
        double r1 = cfg.cylinder_radius * std::sqrt(unit_distrib(rng));
        double r2 = cfg.cylinder_radius * std::sqrt(unit_distrib(rng));
        double theta1 = angle_distrib(rng);
        double theta2 = angle_distrib(rng);
        Eigen::Vector3d C1(r1 * std::cos(theta1),
                           r1 * std::sin(theta1),
                           cam_z_distrib(rng));
        Eigen::Vector3d C2(r2 * std::cos(theta2),
                           r2 * std::sin(theta2),
                           cam_z_distrib(rng));

        // Both cameras look at the scene center
        Eigen::Matrix3d R1 = look_at(C1, cfg.scene_center);
        Eigen::Matrix3d R2 = look_at(C2, cfg.scene_center);
        Eigen::Vector3d t1 = -R1 * C1;
        Eigen::Vector3d t2 = -R2 * C2;

        data.R1 = R1;  data.R2 = R2;
        data.t1 = t1;  data.t2 = t2;

        // Sample 7 points visible in both cameras
        for (int i = 0; i < 7; ++i) {
            for (int attempt = 0; ; ++attempt) {
                if (attempt > 100000) {
                    std::cerr << "Error: Could not find visible point after 100000 attempts." << std::endl;
                    exit(EXIT_FAILURE);
                }
                Eigen::Vector3d X(x_distrib(rng), y_distrib(rng), z_distrib(rng));
                if (is_visible(X, K, R1, t1, cfg.image_width, cfg.image_height) &&
                    is_visible(X, K, R2, t2, cfg.image_width, cfg.image_height)) {
                    data.world_points.push_back(X);
                    break;
                }
            }
        }

        // Project all 7 points to both images
        for (const Eigen::Vector3d& X : data.world_points) {
            data.Ps.push_back(Utils::project_point(X, K, R1, t1));
            data.Qs.push_back(Utils::project_point(X, K, R2, t2));
        }

        // Add pixel noise to all 7 projected points
        add_pixel_noise(data, independent_noise_sigma, rng);
    }

} // namespace SyntheticData
