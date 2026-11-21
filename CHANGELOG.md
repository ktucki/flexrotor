# Changelog

## 1.0.0 (2026-11-21)
First release: reading soil bin records, soil resistance torque, disc model with dry
friction and end stops, identification with parameter correlation, cross-validation, order
analysis, tine kinematics, coupling characteristic in soil, mean-load model, selection rule
and spring selection maps; SVG figures; tests against SciPy and the published results.
Rig and rotor description in a configuration file (`--config`, `examples/rig_sggw.cfg`), so
that records of another logger, rotor or coupling are processed without recompiling.
Rotor pair prediction (`pair`): relative disc angle, tine-order torque and tip distance of two
adjacent rotors from records of a single rotor.
Per-revolution minimum tip distance of a rotor pair (pair_distance.csv); specific work of cutting for any number of tines.
Design path: torsion spring dimensioning after EN 13906-3 (`spring`), virtual soil bin test with
a synthetic soil load (`virtual`) and tine paths, clearance and uncut element of a rotor row (`row`).
Fatigue check of the spring by the Gerber criterion (keys `fatigue_r6`, `fatigue_r5`).
