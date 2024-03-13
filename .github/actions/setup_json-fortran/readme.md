# Setup of the JSON-Fortran action

This action uses caching to accelerate the building of JSON-Fortran. The action
will download the JSON-Fortran source code and build it using the provided
compiler. The action will also cache the build to speed up future builds.

## Inputs

- `install-dir`  
**Optional** The directory to install JSON-Fortran to.  
*Default* `/home/runner/pkg/json-fortran`
- `working-dir`  
**Optional** The directory to work in.  
*Default* `/home/runner/tmp/json-fortran`
- `os`  
**Optional** The operating system to use for building JSON-Fortran. Which should
allow the use of matrix workflows.  
*Default* is to use the runner's operating system.
- `compiler`  
**Optional** The compiler to use for building JSON-Fortran. The compiler should
be available in the PATH.  
*Default* `gfortran`.
- `compiler-options`  
**Optional** The compiler options to use for building JSON-Fortran.  
*Default* `-O3`
- `build-options`  
**Optional** The build options to use for building JSON-Fortran.  
*Default* `--parallel=$(nproc)`.

## Example usage

```yaml
uses: ./.github/actions/setup_json-fortran
with:
  compiler: 'gfortran'
  compiler-options: '-O3'
  build-options: '-j 4'
```
