#ifndef SYNTHETICDATA_H
#define SYNTHETICDATA_H

#include "SceneData.h"
#include <random>

namespace SyntheticData {

    enum class Mode { WAMI, General };

    // Generate a synthetic scene with 7 world points projected to two cameras.
    // Ps/Qs contain 7 projected pixel coordinates (homogeneous).
    // world_points contains the 7 3D points.
    // independent_noise_sigma: Gaussian pixel noise applied to all projections.
    void generate(SceneData& data, std::mt19937& rng, Mode mode = Mode::WAMI,
                  double independent_noise_sigma = 0.0);
}

#endif // SYNTHETICDATA_H
