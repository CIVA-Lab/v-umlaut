# Linear Fundamental Matrix Estimation from 7 or 5 Points

Repository for our CVPR 2026 paper: [Linear Fundamental Matrix Estimation from 7 or 5 Points](https://openaccess.thecvf.com/content/CVPR2026/html/Kucukpinar_Linear_Fundamental_Matrix_Estimation_from_7_or_5_Points_CVPR_2026_paper.html)

This repository contains a C++ implementation of the **V-Umlaut** fundamental-matrix
solver, without a RANSAC pipeline, using synthetic experiments to understand how
different noise models affect the solver's results.

The solver's linear formula is derived with [Macaulay2](https://macaulay2.com/) in
[`derivation/v_umlaut.m2`](derivation/v_umlaut.m2), and implemented in
[`Solver.cpp`](Solver.cpp).

To use the V-Umlaut solver inside [PoseLib](https://github.com/PoseLib/PoseLib)'s RANSAC
pipeline, as a drop-in minimal fundamental matrix solver or for calibrated relative pose
(optionally with early non-minimal refitting, ENM),
[see the v-umlaut branch in our PoseLib fork](https://github.com/CIVA-Lab/PoseLib/tree/v-umlaut).

Our real-data experiments in the paper follow the evaluation pipeline of this
[repository](https://github.com/kocurvik/threeview) for the paper
[Practical Solutions to the Relative Pose of Three Calibrated Cameras](https://openaccess.thecvf.com/content/CVPR2025/papers/Tzamos_Practical_Solutions_to_the_Relative_Pose_of_Three_Calibrated_Cameras_CVPR_2025_paper.pdf)
by Tzamos et al., CVPR 2025. [Our fork of this pipeline](https://github.com/CIVA-Lab/threeview)
reproduces the real-data figures of our paper (Figure 5 and supplementary Figure 6) for
`5pF-V-Umlaut` and `5pF-V-Umlaut + ENM (5p)`.

## Synthetic Experiments

These experiments reproduce **Section 4 / Figure 4** of our paper: the relative-pose
accuracy of `5pF-V-Umlaut` against the `5pE` (5-point essential) and `7pF` (7-point
fundamental) solvers, under two pixel-noise models, on two camera/scene configurations,
without RANSAC. The synthetic data generation and the two noise models are described in
detail in **Section 4** of our paper.

### Dependencies
- A C++17 compiler and CMake ≥ 3.10
- Eigen 3

The two baseline solvers used for comparison, Nister's 5-point relative-pose solver (`5pE`)
and the 7-point fundamental-matrix solver (`7pF`), are minimally extracted from their
[PoseLib](https://github.com/PoseLib/PoseLib) implementations and included in this repository
under `external/poselib/`. PoseLib itself is therefore not a dependency.

### Build
```bash
# Linux / macOS
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
#   -> build/v_umlaut

# Windows
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -DEigen3_DIR=<prefix>/share/eigen3/cmake
cmake --build build --config Release
#   -> build/Release/v_umlaut.exe
```

### Reproducing Figure 4
After building, run (Python 3 with matplotlib):

```bash
python plot_figure4.py    # --binary <path> if the binary is not in build/
```
It runs the binary for every data point of Figure 4 (both scenes, both noise models,
σ = 0 to 2 px, 10,000 trials each) in parallel, saves the median and interquartile range of
the rotation and translation errors to `build/figure4.csv`, and draws the figure to
`build/figure4.pdf` and `build/figure4.png`. Runs are not seeded, so the medians vary by a
few percent between runs.

### Reproducing a Figure 4 data point
Each run performs many trials and prints the median rotation/translation error for all
three solvers. One *data point* in Figure 4 corresponds to a single (scene, noise-model, σ)
setting. For the `Standard` scene at σ = 1 px over 10,000 trials:

```bash
# Figure 4a — independent-point noise
./v_umlaut general 10000 --independent-point-noise 1.0

# Figure 4b — dependent-point noise
./v_umlaut general 10000 --dependent-point-noise 1.0
```
Swap `general` → `wami` for the WAMI configuration, and vary σ (e.g. `0.5`, `1.5`, `2.0`)
to move along the x-axis of Figure 4. The three methods are printed as `PoseLib-5pt` (5pE),
`PoseLib-7pt` (7pF), and `V-Umlaut` (5pF-V-Umlaut).

## Citation

```bibtex
@InProceedings{V-Umlaut_2026_CVPR,
    author    = {Kucukpinar, Taci Ata and Mogollon, Juan and Fraser, Joshua and Duff, Timothy and Palaniappan, Kannappan},
    title     = {Linear Fundamental Matrix Estimation from 7 or 5 Points},
    booktitle = {Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR)},
    month     = {June},
    year      = {2026},
    pages     = {21464-21473}
}
```

If you cite our paper, please consider also citing:
- [PLMP – Point-Line Minimal Problems for Projective SfM](https://openaccess.thecvf.com/content/ICCV2025/html/Kiehn_PLMP_-_Point-Line_Minimal_Problems_for_Projective_SfM_ICCV_2025_paper.html) (ICCV 2025), whose Table 1 lists the V-Umlaut configuration and classifies it as a degree-1 minimal problem.
- [Practical Solutions to the Relative Pose of Three Calibrated Cameras](https://openaccess.thecvf.com/content/CVPR2025/html/Tzamos_Practical_Solutions_to_the_Relative_Pose_of_Three_Calibrated_Cameras_CVPR_2025_paper.html) (CVPR 2025), whose evaluation pipeline our real-data experiments use.

```bibtex
@InProceedings{Kiehn_2025_ICCV,
    author    = {Kiehn, Kim and Ahlb\"ack, Albin and Kohn, Kathl\'en},
    title     = {PLMP - Point-Line Minimal Problems for Projective SfM},
    booktitle = {Proceedings of the IEEE/CVF International Conference on Computer Vision (ICCV)},
    month     = {October},
    year      = {2025},
    pages     = {8558-8567}
}

@InProceedings{Tzamos_2025_CVPR,
    author    = {Tzamos, Charalambos and Kocur, Viktor and Ding, Yaqing and Barath, Daniel and Haladova, Zuzana Berger and Sattler, Torsten and Kukelova, Zuzana},
    title     = {Practical Solutions to the Relative Pose of Three Calibrated Cameras},
    booktitle = {Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR)},
    month     = {June},
    year      = {2025},
    pages     = {21913-21923}
}
```
