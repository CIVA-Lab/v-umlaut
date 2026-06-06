Minimal extraction from PoseLib (https://github.com/PoseLib/PoseLib).
BSD-3 license, full text in the `LICENSE` file in this directory.
Copyright (c) 2020-2021, Viktor Larsson.

Only the following components are extracted:
- 5-point relative pose solver (Nister)
- 7-point relative pose solver
- Supporting utilities (Sturm bracketing, univariate polynomial solvers)

## Changes from the original PoseLib source

### Include path adjustments
All files originally used PoseLib's internal include paths (e.g.,
`#include "PoseLib/misc/essential.h"`). These have been changed to flat
relative includes within this folder:
- `relpose_5pt.cc`: `#include "PoseLib/misc/essential.h"` and
  `#include "PoseLib/misc/sturm.h"` replaced with `#include "sturm.h"`.
  The essential.h include was removed entirely (only needed by the
  CameraPose overload which was removed).
- `relpose_7pt.cc`: `#include "PoseLib/misc/univariate.h"` replaced with
  `#include "univariate.h"`.
- `camera_pose.h`: `#include "PoseLib/misc/quaternion.h"` replaced with
  `#include "quaternion.h"`. Includes for `"alignment.h"` and
  `"misc/camera_models.h"` removed (not needed).

### relpose_5pt.cc — second overload removed
The original file has two overloads of relpose_5pt:
1. Returns `std::vector<Eigen::Matrix3d>` (essential matrices) — kept.
2. Returns `std::vector<CameraPose>` (calls motion_from_essential internally) — removed.
