#include <iostream>
#include <iomanip>
#include <string>
#include <random>
#include <map>

#include "SceneData.h"
#include "SyntheticData.h"
#include "Solver.h"
#include "Utils.h"

#include "external/poselib/relpose_5pt.h"
#include "external/poselib/relpose_7pt.h"

void print_usage() {
    std::cout << "Usage:\n";
    std::cout << "  ./v_umlaut [wami|general] [num_iterations]"
        << " [--independent-point-noise sigma] [--dependent-point-noise sigma]\n";
    std::cout << "    --independent-point-noise: 2D pixel noise on the 5 independent points, both views\n";
    std::cout << "    --dependent-point-noise:   1D pixel noise on the 2 dependent points, cam2 only\n";
}

bool find_option(int argc, char* argv[], const std::string& option, std::string& value) {
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == option) {
            if (i + 1 < argc) {
                value = argv[i + 1];
                return true;
            } else {
                std::cerr << "Error: Option " << option << " requires a value." << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    return false;
}

int main(int argc, char* argv[]) {

    if (argc > 1 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
        print_usage();
        return EXIT_SUCCESS;
    }

    SyntheticData::Mode synth_mode = SyntheticData::Mode::WAMI;
    std::string mode_name = "wami";
    int num_iterations = 1000;
    int next_arg = 1;

    if (next_arg < argc) {
        std::string arg = argv[next_arg];
        if (arg == "wami" || arg == "general") {
            mode_name = arg;
            synth_mode = (arg == "general") ? SyntheticData::Mode::General
                : SyntheticData::Mode::WAMI;
            next_arg++;
        }
    }
    if (next_arg < argc && std::string(argv[next_arg]) != "--independent-point-noise"
        && std::string(argv[next_arg]) != "--dependent-point-noise")
        num_iterations = std::stoi(argv[next_arg++]);

    double independent_noise = 0.0;
    std::string independent_noise_val;
    if (find_option(argc, argv, "--independent-point-noise", independent_noise_val))
        independent_noise = std::stod(independent_noise_val);

    double dependent_noise = 0.0;
    std::string dependent_noise_val;
    if (find_option(argc, argv, "--dependent-point-noise", dependent_noise_val))
        dependent_noise = std::stod(dependent_noise_val);

    std::string dependent_info;
    if (dependent_noise > 0) {
        dependent_info = ", dependent-point-noise=" + std::to_string(dependent_noise);
    }

    std::cout << "Running " << num_iterations << " synthetic " << mode_name
        << " iterations (independent-point-noise=" << independent_noise
        << dependent_info
        << ")...\n" << std::endl;

    // ---- Run experiment ----

    std::mt19937 rng(std::random_device{}());
    SceneData data;

    struct MethodErrors {
        std::vector<double> rot;
        std::vector<double> trans;
    };
    std::map<std::string, MethodErrors> methods;

    for (int i = 0; i < num_iterations; ++i) {
        SyntheticData::generate(data, rng, synth_mode, independent_noise);

        Eigen::Matrix3d R_gt = data.R2 * data.R1.transpose();
        Eigen::Vector3d t_gt = data.t2 - R_gt * data.t1;

        std::vector<Eigen::Vector3d> Ps_normalized = Utils::to_normalized(data.Ps, data.K1);
        std::vector<Eigen::Vector3d> Qs_normalized = Utils::to_normalized(data.Qs, data.K2);
        std::vector<Eigen::Vector3d> Ps_normalized_5(Ps_normalized.begin(), Ps_normalized.begin() + 5);
        std::vector<Eigen::Vector3d> Qs_normalized_5(Qs_normalized.begin(), Qs_normalized.begin() + 5);

        // V-Umlaut (5 general + 2 midpoints)
        {
            std::vector<Eigen::Vector3d> Ps_vumlaut, Qs_vumlaut;
            Utils::build_v_umlaut_input(data, Ps_vumlaut, Qs_vumlaut, dependent_noise, rng);
            std::vector<Eigen::Vector3d> Ps_vumlaut_normalized = Utils::to_normalized(Ps_vumlaut, data.K1);
            std::vector<Eigen::Vector3d> Qs_vumlaut_normalized = Utils::to_normalized(Qs_vumlaut, data.K2);

            Eigen::Matrix3d E_est = Solver::v_umlaut(Ps_vumlaut_normalized, Qs_vumlaut_normalized);
            auto [R_est, t_est] = Utils::recover_pose_from_essentials({E_est}, R_gt, t_gt);

            methods["V-Umlaut"].rot.push_back(Utils::rotation_error(R_est, R_gt));
            methods["V-Umlaut"].trans.push_back(Utils::translation_error(t_est, t_gt));
        }

        // PoseLib 7pt (7 general points)
        {
            std::vector<Eigen::Matrix3d> F_matrices;
            poselib::relpose_7pt(Ps_normalized, Qs_normalized, &F_matrices);

            auto [R_est, t_est] = Utils::recover_pose_from_essentials(F_matrices, R_gt, t_gt);

            methods["PoseLib-7pt"].rot.push_back(Utils::rotation_error(R_est, R_gt));
            methods["PoseLib-7pt"].trans.push_back(Utils::translation_error(t_est, t_gt));
        }

        // PoseLib 5pt (first 5 general points)
        {
            std::vector<Eigen::Matrix3d> E_matrices;
            poselib::relpose_5pt(Ps_normalized_5, Qs_normalized_5, &E_matrices);

            auto [R_est, t_est] = Utils::recover_pose_from_essentials(E_matrices, R_gt, t_gt);

            methods["PoseLib-5pt"].rot.push_back(Utils::rotation_error(R_est, R_gt));
            methods["PoseLib-5pt"].trans.push_back(Utils::translation_error(t_est, t_gt));
        }

        if ((i + 1) % 100 == 0 || i == num_iterations - 1)
            std::cout << "  Completed " << (i + 1) << "/" << num_iterations << std::endl;
    }

    // ---- Report results ----

    std::cout << "\n--- Results over " << num_iterations << " iterations ---" << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    for (auto& [name, errs] : methods) {
        std::cout << "\n[" << name << "]" << std::endl;
        Utils::print_stats("  Rotation error (deg):   ", Utils::compute_stats(errs.rot));
        Utils::print_stats("  Translation error (deg):", Utils::compute_stats(errs.trans));
    }

    return EXIT_SUCCESS;
}