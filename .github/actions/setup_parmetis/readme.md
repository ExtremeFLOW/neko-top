# Installation and caching of ParMetis

This action uses caching to accelerate the building of ParMetis. The action will
download the ParMetis source code and build it using the provided compiler. The
action will also cache the build to speed up future builds.

It should be noted that the download is done using wget from a github reference
mentioned by the Neko documentation.

## Inputs

| Name             | Optional | Description                                    | Default                     |
| ---------------- | -------- | ---------------------------------------------- | --------------------------- |
| `install-dir`    | Yes      | The directory to install the ParMetis library. | `/home/runner/pkg/parmetis` |
| `working-dir`    | Yes      | The directory to work in.                      | `/home/runner/tmp/parmetis` |
| `os`             | Yes      | The operating system to use.                   | `${{ runner.os }}`          |
| `compiler`       | Yes      | The compiler to use.                           | `mpif90`                    |
| `compiler-flags` | Yes      | The compiler flag to use.                      | `-O3`                       |
| `build-options`  | Yes      | The build option to use.                       | `--parallel $(nproc)`       |
| `version`        | Yes      | The version of the ParMetis library to use.    | `4.0.3`                     |

## Outputs

| Name          | Description                                 |
| ------------- | ------------------------------------------- |
| `install-dir` | The directory where ParMetis was installed. |

## Example usage

The following example uses the ParMetis action to build ParMetis using the `gfortran`
compiler with the `-O3` optimization level and the `--parallel 4` cmake build
option.

Additionally the next step will capture the install location and print where the
ParMetis was installed.

```yaml
name: Build ParMetis

on: [push, pull_request]

jobs:
  build:
    runs-on: ubuntu-latest
    steps:

    - name: Setup ParMetis
      id: setup-parmetis
      uses: ./.github/actions/setup_parmetis
      with:
        compiler: 'gfortran'
        compiler-options: '-O3'
        build-options: '--parallel=4'

    - name: Echo install directory
      env: 
        PARMETIS_DIR: ${{ steps.setup-parmetis.outputs.install-dir }}
      run: echo "parmetis was installed to $PARMETIS_DIR"
```