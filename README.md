# FlexRotor — compliant tine disc coupling of a power harrow rotor

FlexRotor designs the compliant coupling between the shaft and the tine disc of a power
harrow rotor from the rotor geometry, calibrated on soil bin records. In the rotor studied, the tine disc is
mounted on the shaft on bearings and connected to it by a preloaded torsion spring with
end stops. The program works on the raw records of the shaft torque and the shaft and
disc encoders, and it gives the quantities a designer needs to choose the spring:
the threshold torque at which the disc yields, the coupling characteristic under
operation in soil, the dry friction in the coupling, and the stiffness and threshold torque
that keep the disc in compliant operation in a target soil. From these it dimensions the
torsion spring (EN 13906-3), tests the disc in a virtual soil bin and computes the tine paths
and the tine clearance of adjacent rotors of a row.

It is a small C++17 program and library with no external dependencies.

What it does:

* **Soil resistance torque from raw records.** It reads the logger records (2 ms), unwraps
  the 32-bit timestamps, converts the encoder counts to the disc lag angle and the shaft
  angle, detects the soil phase (disc lifted from the initial stop, or the middle part of
  the carriage travel) and subtracts the drive-train loss torque Ms(ω) measured in air.
  A zero-phase fourth-order Butterworth filter (identical to `scipy.signal.filtfilt`)
  separates the slow soil load from the tine-entry components.
* **Dynamic model of the disc.** J α̈ = Mg − k(α + β) − μα̇ − Mc tanh(α̇/v0) − Mo(α, α̇), with a
  preloaded spring, viscous damping, dry friction and penalty end stops, integrated by RK4
  with sub-steps. A quasi-static play model is included for comparison.
* **Identification and identifiability.** Bounded Levenberg–Marquardt fit of the spring
  stiffness, preload, damping and dry friction to the measured disc angle of several runs,
  with the parameter correlation matrix. Strongly correlated stiffness and preload show
  that the records determine the threshold torque M* = kβ rather than k and β separately.
  The viscous-only model is fitted for comparison.
* **Predictive capability.** Leave-one-run-out and leave-one-speed-group-out validation.
* **Order analysis.** Angle-synchronous Fourier coefficients of the shaft torque and the
  disc angle at the tine-entry order, the soil torque component corrected for disc inertia,
  and the response of the linearised coupling.
* **Tine kinematics.** Change of the disc angle over one tine period, the resulting bite
  length factor of the trochoidal tine path and the power dissipated in the coupling.
* **Coupling characteristic in soil.** Torque–angle relation from 5° angle bins of the
  smoothed records, M = M0 + kα α.
* **Mean-load model and selection rule.** The mean soil torque is described as a component
  independent of the bite length plus the work of cutting a unit volume of soil,
  M = Ma + d_setting + K s_k with s_k = 2πu/(zω); the specific work of cutting is e = 2πK/(zhb).
  Combined with the characteristic, the model predicts the disc angle from the kinematic
  parameters alone (leave-one-out). The selection rule gives the minimum stiffness and the
  threshold torque that keep the disc between the stops over the operating range of a
  target soil.
* **Spring selection map.** Fractions of time on the initial stop, in compliant operation
  and on the end stop over a grid of stiffness and threshold torque, simulated with the
  measured soil torque.
* **Rotor pair (prediction).** Two adjacent counter-rotating rotors on a common gear train:
  relative disc angle, tine-order torque of the pair relative to the sum of the two rotors and,
  for a given rotor spacing and tip radius, the minimum tip-to-tip distance. Rotor 2 is the same
  record shifted in time; the interaction of the rotors through the soil is not represented, so
  the results are predictions to be verified on a rotor pair.
* **Spring dimensioning (`spring`).** Torsion spring after EN 13906-3: rate k = E d⁴/(3667 D n),
  static bending stress 32M/(πd³) at the end stop against 0.7 Rm of EN 10270-1 wire, stress
  fatigue check by the Gerber criterion with the corrected stresses (strength of a repeated stress 0.50 Rm at 10⁶ and
  0.53 Rm at 10⁵ cycles for torsion springs of patented wire, not shot peened; slow fluctuation and whole operating
  range), coil diameter wound up by the preload
  and the rotation range (clearance on the arbor) and body length. Checks a given spring and
  returns, for every standard wire diameter, the most compact spring that meets the selection
  rule within the limits Do_max, Di_min and L_max of the configuration file, the largest rotation
  range of the disc that a spring within these limits allows, and the length of the shortest spring
  over the rotation range (spring_range.csv).
* **Virtual soil bin test (`virtual`).** The soil load is the mean torque of the load model plus
  an Ornstein–Uhlenbeck fluctuation with the standard deviation and correlation time estimated
  from the records; the disc model gives the disc angle and the time at the stops. Each run is
  predicted without its own record (validation), or any spring is tested (`--params`,
  `--soil-scale` for a harder soil). The random numbers are generated in the program, so the
  results do not depend on the compiler.
* **Rotor row (`row`).** Two adjacent counter-rotating rotors with synthetic loads (correlation
  `load_corr`): tine tip paths, minimum tine clearance per revolution, its dependence on the
  axis distance and on the tine phase offset, and the largest soil element not crossed by any
  tip path (largest empty circle on a 1 mm raster, row repeated laterally).
* **Figures.** SVG plots with thin lines, axes that end on a labelled tick and the legend
  above the plot.

## Build

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build
```

Tested with GCC 13 and Clang 18 on Linux; the code uses only the C++17 standard library and
the repository includes a continuous-integration workflow for Linux, Windows and macOS. The tests check the
signal processing against SciPy reference values, the static equilibrium of the disc model,
the recovery of known parameters from synthetic data and the published results.

## Use

```
flexrotor summary    examples/runs.csv --out out
flexrotor identify   examples/runs.csv --setting N2 --out out
flexrotor cv         examples/runs.csv --setting N2
flexrotor orders     examples/runs.csv --setting N2 --params out/params.txt
flexrotor kinematics examples/runs.csv --setting N2 --params out/params.txt
flexrotor map        examples/runs.csv --setting N2 --params out/params.txt --out out
flexrotor loadmodel  examples/runs.csv --out out
flexrotor pair       examples/runs.csv --setting N2 --out out
flexrotor spring     examples/runs.csv --config examples/rig_sggw.cfg --setting N2 --params out/params.txt --out out
flexrotor virtual    examples/runs.csv --config examples/rig_sggw.cfg --setting N2 --params out/params.txt --out out
flexrotor row        examples/runs.csv --config examples/rig_sggw.cfg --setting N2 --params out/params.txt --out out
```

`run_examples.sh` (or `.bat`) runs the whole example, including the virtual soil bin test of
the designed springs.

Every command accepts `--config rig.cfg`, a plain `key = value` description of the logger
channels and scaling, the drive-train loss polynomial, the filter and soil-phase thresholds,
the disc inertia and rotation range, the rotor geometry (tip radius, tines, axis distance,
phase offset), the coupling spring and the design limits. `examples/rig_sggw.cfg` lists all
keys with the values of the SGGW soil bin rotor (the built-in defaults); copy it to process
records of another rotor, coupling or logger without recompiling. Unknown keys are rejected.

The run table lists the records: `name,file,setting,variant,rpm,u`. Runs with the same
setting and variant are averaged by `loadmodel`; `cv` groups the runs by the nominal shaft
speed. Records exported by the logger to Excel are converted to the CSV read by the program
with `tools/myrio_xlsx2csv.py` (columns: timestamp, shaft encoder, disc encoder, carriage
encoder, torque ×10, speed ×10). Channel scaling, the drive-train loss polynomial and the
disc inertia are set in the configuration file (or in `fr::RigConfig` when the library is
called from C++).

## Example data

`examples/data` contains the 19 soil bin records of the compliant rotor at a working depth
of 14 cm in sandy loam: setting N1 (25 September 2024) and setting N2 (27 September 2024),
each at three shaft speeds (150, 200 and 250 rpm) and three forward speeds
(0.23–0.32 m/s). Running the commands above reproduces the published results:

| quantity | FlexRotor | published |
|---|---|---|
| threshold torque M*, setting N2 | 11.80 N·m | 11.8 N·m |
| dry friction Mc | 1.42 N·m | 1.4 N·m |
| disc angle RMSE (R²), dry friction model | 3.93° (0.957) | 3.9° (0.96) |
| disc angle RMSE, viscous-only model | 5.55° | 5.6° |
| correlation of k and β | −0.99 | −0.99 |
| order-2 disc response | 0.03–0.23° | 0.03–0.23° |
| bite length factor, 2.5–97.5% band over all runs / maximum | 0.95–1.13 / 1.87 | typically 3–6% / 87% |
| characteristic in soil, N1 | M0 3.97 N·m, kα 256 N·mm/° | 3.6 N·m, 255 N·mm/° |
| characteristic in soil, N2 | M0 12.65 N·m, kα 116 N·mm/° | 12.7 N·m, 114 N·mm/° |
| specific work of cutting | 17.2 kJ/m³ | 17.3–17.9 kJ/m³ |
| minimum stiffness from the selection rule | 134 N·mm/° | 134 N·mm/° |

Design results for the same rotor (not published before):

| quantity | FlexRotor |
|---|---|
| virtual soil bin test, mean disc angle of unseen runs | RMSE 5.8° (R² 0.88) |
| prototype spring (d 5 mm, D 46.5 mm, n 5.5): bending stress at the measured torque | 1404–1967 N/mm² (0.7 Rm = 1162 N/mm²) |
| space of the prototype (drawings): shaft 35 mm, body length 40 mm, outer diameter 70 mm | no spring meets the rule for 90° |
| designed spring, tested soil, prototype space | d 6.5 mm, D 63.2 mm, n 4.5, k 353 N·mm/°, β 43°, rotation range 36° |
| fatigue safety factor, fluctuation (10⁶) / operating range (10⁵): designed / prototype spring | 1.29 / 1.12; 0.64 / 0.55 |
| designed spring, tested soil, full 90° range | d 6.5 mm, D 63.5 mm, n 11.75, body length 86 mm |
| prototype space, 1.5× higher mean torque | rotation range 21° |
| tine clearance at a = 217 mm, φ = 90°: rigid / compliant | 153 mm / 121 mm |
| phase offset for a clearance ≥ 20 mm: rigid / compliant | 11–169° / 30–147° |

Small differences in the characteristic and the mean-load model come from the bin-count
rule and the averaging of the shaft speed; they are within the run repeatability.

## Support and contributing

Questions and bug reports: GitHub issues, or karol_tucki@sggw.edu.pl. See `CONTRIBUTING.md`.

## Licence and citation

MIT licence (see `LICENSE.txt`). If you use FlexRotor, please cite the software
(`CITATION.cff`; version 1.0.0: https://doi.org/10.5281/zenodo.24419296) and the studies in which the methods were developed:

* K. Tucki, Ł. Gruz, L. Mieszkalski, Response of a power harrow rotor with a compliantly
  mounted tine disc to soil resistance: soil bin identification, threshold torque and
  friction-limited compliance, Journal of Terramechanics (under review).
* K. Tucki, Ł. Gruz, L. Mieszkalski, A load model and a selection rule for a torsionally
  compliant tine disc coupling of a rotary harrow rotor: soil bin tests of two settings,
  Soil & Tillage Research (under review).
