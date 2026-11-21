# Contributing to FlexRotor

Questions, bug reports and suggestions are welcome as GitHub issues
(https://github.com/ktucki/flexrotor/issues). For a bug, please attach the run table, the
configuration file and, if possible, one record that shows the problem.

Pull requests should:

* build without warnings with `-Wall -Wextra -Wpedantic` (GCC, Clang) or `/W4` (MSVC);
* keep the library free of external dependencies (C++17 standard library only);
* add a check to `tests/test_main.cpp` for new functions and keep all existing checks passing
  (`ctest --test-dir build`);
* keep the example results in `README.md` unchanged, or explain in the pull request why
  they change.

To use FlexRotor with another rotor or logger, start from `examples/rig_sggw.cfg`.
