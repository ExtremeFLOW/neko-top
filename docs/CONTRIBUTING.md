# Contributing to Neko-TOP {#contributing}

This document describes the guidelines for contributing to Neko-TOP. Neko-TOP is
a topology optimization code based on the Neko library. As such, the guidelines
for contributing to Neko-TOP are the same as those for contributing to Neko.
Please refer to the 
for more details.

## Checks for the Pull Request

When a pull request is submitted, a series of continuous integration tests will
be run. A pull request will not be accepted nor merged into `develop` until it
passes the test suite.

- The pull request should be based on the `develop` branch.
- The pull request should pass the test suite.  
  The test suite can be run locally by executing `./setup.sh -t` in the root
  directory of the repository.
- The pull request should not introduce any new warnings or errors from the
  compilers.  
  This can be checked by compiling the code in `Debug` mode. This can be
  done by executing `CMAKE_BUILD_TYPE=Debug ./setup.sh` in the root
  directory of the repository.
- The pull request should adhere to the code style guidelines.  
  The flint tool can be used to check the code style. This can be done
  by executing `flint -R flinter_rc.yml [sources/ | examples/ | tests/]` in the root directory of the
  repository. Please note that the flint tool is not included in the
  repository. It can be installed by executing `pip install flint` in the
  terminal.
