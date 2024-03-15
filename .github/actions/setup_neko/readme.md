# Installation and caching of Neko

This action uses caching to accelerate the building of Neko. The action will
download the Neko source code and build it using the provided compiler. The
action will also cache the build to speed up future builds.

The version of Neko to build can be specified using the `version` input. The
default version is `v0.7.2`. If the version is set to `develop` the action will
build the latest commit on the `develop` branch and append the date to the 
hash. This ensures we update the cache regularly.

## Inputs

| Name               | Optional | Description                                                                                    | Default                         |
| ------------------ | -------- | ---------------------------------------------------------------------------------------------- | ------------------------------- |
| `install-dir`      | Yes      | The directory to install Neko to.                                                              | `/home/runner/pkg/neko`         |
| `working-dir`      | Yes      | The directory to work in.                                                                      | `/home/runner/tmp/neko`         |
| `os`               | Yes      | The operating system to use for building Neko. Which should allow the use of matrix workflows. | `runner.os`                     |
| `compiler`         | Yes      | The compiler to use for building Neko. The compiler should be available in the PATH.           | `gfortran`.                     |
| `compiler-options` | Yes      | The compiler options to use for building Neko.                                                 | `-O3`                           |
| `build-options`    | Yes      | The build options to use for building Neko.                                                    | `-j $(nproc)`.                  |
| `version`          | Yes      | The version of Neko to build.                                                                  | `v0.7.2`.                       |
| `json-fortran-dir` | Yes      | The directory where JSON-Fortran was installed.                                                | `/home/runner/pkg/json-fortran` |

## Outputs

| Name          | Description                             |
| ------------- | --------------------------------------- |
| `install-dir` | The directory where Neko was installed. |

## Example usage

The following example uses the Neko action to build Neko using the `gfortran`
compiler with the `-O3` optimization level and the `--parallel 4` cmake build
option.

Additionally the next step will capture the install location and print where the
Neko was installed.

```yaml
name: Build Neko

on: [push, pull_request]

jobs:
  build:
    runs-on: ubuntu-latest
    steps:

    - name: Setup Neko
      uses: ./.github/actions/setup_neko
      with:
        compiler: 'gfortran'
        compiler-options: '-O3'
        build-options: '--parallel=4'

    - name: Echo install directory
      env: 
        NEKO_DIR: ${{ steps.setup-neko.outputs.install-dir }}
      run: echo "neko was installed to $NEKO_DIR"
```