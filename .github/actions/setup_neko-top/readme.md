# Installation and caching of Neko-TOP

This action uses caching to accelerate the building of Neko-TOP. The action
will download the Neko-TOP source code and build it using the provided
compiler. The action will also cache the build to speed up future builds.

## Inputs

| Name               | Optional | Description                                                                                        | Default                     |
| ------------------ | -------- | -------------------------------------------------------------------------------------------------- | --------------------------- |
| `install-dir`      | Yes      | The directory to install Neko-TOP to.                                                              | `/home/runner/pkg/neko-top` |
| `working-dir`      | Yes      | The directory to work in.                                                                          | `/home/runner/tmp/neko-top` |
| `os`               | Yes      | The operating system to use for building Neko-TOP. Which should allow the use of matrix workflows. | `runner.os`                 |
| `compiler`         | Yes      | The compiler to use for building Neko-TOP. The compiler should be available in the PATH.           | `gfortran`.                 |
| `compiler-options` | Yes      | The compiler options to use for building Neko-TOP.                                                 | `-O3`                       |
| `build-options`    | Yes      | The build options to use for building Neko-TOP.                                                    | `--parallel=$(nproc)`.      |
| `version`          | Yes      | The version of Neko-TOP to build.                                                                  | `develop`.                  |

## Outputs

| Name          | Description                                 |
| ------------- | ------------------------------------------- |
| `install-dir` | The directory where Neko-TOP was installed. |

## Example usage

The following example uses the Neko-TOP action to build Neko-TOP using
the `gfortran` compiler with the `-O3` optimization level and the `--parallel 4`
cmake build option.

Additionally the next step will capture the install location and print where the
Neko-TOP was installed.

```yaml
name: Build Neko-TOP

on: [push, pull_request]

jobs:
  build:
    runs-on: ubuntu-latest
    steps:

    - name: Setup Neko-TOP
      uses: ./.github/actions/setup_neko-top
      with:
        compiler: 'gfortran'
        compiler-options: '-O3'
        build-options: '--parallel=4'

    - name: Echo install directory
      env: 
        NEKO_TOP_DIR: ${{ steps.setup-neko-top.outputs.install-dir }}
      run: echo "Neko-TOP was installed to $NEKO_TOP_DIR"
```