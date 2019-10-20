# Testing Facility

This part of repo defines optional tools for testing and debugging.

* Common testing routines are the subject of dedicated shared
library (see `common/`).
* Utilities for plotting and demo processing are located at `utils/` directory.
* Unit testing facility is located at `units/` directory.

Unit tests may also provide code coverage tests.

## Note about Code Coverage

Important technical information about unit testing coverage is obtained with a
specific tools: the unit tests are accomodated with some initial scuffolding
for code-coverage done with [gcov](https://gcc.gnu.org/onlinedocs/gcc/Gcov.html).
To cope the CMake script and code coverage tools, the
[third-party module](https://github.com/bilke/cmake-modules/blob/master/CodeCoverage.cmake)
is applied. User may perform debug build with `CODECOV` CMake option enabled
to gain the access to special `dph-na64-ut-tests-dbg_coverage` target.

