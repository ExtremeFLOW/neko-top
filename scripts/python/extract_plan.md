# Plan: VTKHDF extract / split tool for Neko & Neko-TOP output

## Context

Neko and Neko-TOP write temporal VTKHDF output as two parts:
- a small main file, `<stem>.vtkhdf`, holding the static mesh and the step metadata;
- one HDF5 file per output step, `<stem>.data/<step>.h5`.

The main file joins the per-step files together through a printf-style virtual dataset (VDS) mapping `<stem>.data/%b.h5`. Production runs reach many GB to TB, so these files cannot be loaded into memory, and post-processing machines cannot hold a whole run. `scripts/python/vtkhdf_extract.py` is currently only a comment listing the wanted operations. We will add a CLI (`main()` + argparse subcommands) and a backend module that use h5py with bounded memory and never modify the original data.

Decisions confirmed with the user:
- **`extract`** writes one step as a **static snapshot**. The output is a standalone, non-temporal VTKHDF with no `Steps` group and no VDS: all data is materialised, and the step's time is kept as `FieldData/Time`.
- **`split`** uses **reference originals**. Each chunk `<stem>_steps_<a>-<b>.vtkhdf` is written next to the input. It holds:
  - a copy of the mesh;
  - `Steps` metadata for the chunk;
  - per-field VDS with one explicit mapping per step, pointing at the original `<stem>.data/<i>.h5` by relative path. No field data is copied.

  A manifest `<stem>_steps_<a>-<b>.files` lists the chunk file and its `.h5` files, relative to the input directory. `rsync -a --files-from=<manifest> <input_dir>/ host:dst/` therefore ships exactly one self-contained chunk. The footprint of a chunk is the chunk file plus its referenced `.h5` files.
- Tests use the local fixtures only: `vtkhdf_test_data/` and `external/neko/tests/unit/`. `logs/` is a live run directory and must not be used.

## Format facts the design relies on

These were verified against Neko's writer and reader (`external/neko/src/io/vtkhdf_file.F90`), VTK 9.3–9.7.1 `vtkHDFReader` source, the HDF5 2.0 source in `external/hdf5/src`, and read-only checks on the fixtures.

**Mesh and metadata as Neko writes them:**
- `/VTKHDF` carries `Type` (fixed 16-char NULLTERM ASCII) and `Version` (`int32[2] = [2,6]`).
- The mesh is static and written once:

  | Dataset        | Type and shape                               |
  | -------------- | -------------------------------------------- |
  | `Points`       | f64 `(N,3)`                                  |
  | `Connectivity` | i32                                          |
  | `Offsets`      | i32, `C+P` entries (each part restarts at 0) |
  | `Types`        | u8                                           |

- `NumberOfPoints`, `NumberOfCells` and `NumberOfConnectivityIds` are i64 with `NSteps·P` entries.
- `Steps` has `NSteps` (int32 attribute), `Values` (f64) and `NumberOfParts` (= P). `PartOffsets`, `PointOffsets`, `CellOffsets` and `ConnectivityIdOffsets` are all 0. `PointDataOffsets/<f>` holds `s·N` for step s.

**PointData:**
- Each field is a VDS with one printf mapping: start 0, stride = block = N rows, count UNLIMITED, source `<stem>.data/%b.h5`, dataset `<f>`, selection ALL.
- Scalars have shape `(·,)`; `Velocity` has shape `(·,3)`.
- The dataset attribute `Attribute` is `Scalars` or `Vectors`.
- The source-space extent is stored as `()`, so N must be derived from `NumberOfPoints`.

**VTK reader (`vtkHDFReader`):**
- Pieces = `len(NumberOfPoints)/NSteps`. A file with no `Steps` group has exactly P entries.
- A file is temporal only if `NSteps > 1`, so a 1-step file is read as static and its `Steps` group is ignored.
- `Steps/Values` must be strictly increasing and exactly NSteps long.
- Every PointData array of a temporal file needs `PointDataOffsets/<f>`. Without it, VTK silently shows step 0.
- `FieldData` datasets must be at least 1-D, and must not appear in temporal files (VTK adds its own `Time` there).
- VTK 9.7 (ParaView 6.2) keeps only the last array tagged `Attribute=Scalars`.

**Neko reader:**
- For an explicit (non-`%b`) VDS, Neko's reader uses mapping index == step counter. It opens the source file read-only and reads the root dataset named like the field.
- Non-VDS PointData is read directly.

**HDF5 and h5py pitfalls:**

| Pitfall                                                              | Consequence                                                                                   |
| -------------------------------------------------------------------- | --------------------------------------------------------------------------------------------- |
| `.shape`, `len()`, `.ndim`, `repr()` or slicing on Neko's printf VDS | Probes and opens every per-step file and keeps them open                                      |
| VDS sources take the access mode of the file holding the VDS         | Reading through a file open for writing opens the originals read-write                        |
| Missing VDS sources                                                  | Read back silently as the fill value                                                          |
| h5py high-level attribute writes                                     | Change NULLTERM to NULLPAD or variable-length, so attributes need a byte-exact low-level copy |
| `libver='latest'` / `'v200'`                                         | Output unreadable by ParaView 6.2's HDF5 1.14.2                                               |
| `%` in VDS names                                                     | Printf-parsed (only `%b` and `%%` are valid)                                                  |

## Files

| File                                     | Change                                                                                                                                                                                                                                                                                                                                                   |
| ---------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `scripts/python/vtkhdf_extract.py`       | Replace the comment stub with a thin driver: module docstring with usage, `parse_args()` with subcommands `info`, `extract`, `split`, and `main()` that dispatches, prints a summary and turns `VtkhdfError` into `exit 1` with the message. Style as `scripts/python/gif_maker.py`: argparse, pathlib, PEP 604 hints, small functions. About 120 lines. |
| `scripts/python/vtkhdf_backend.py` (new) | All HDF5 logic. About 350 lines.                                                                                                                                                                                                                                                                                                                         |
| `scripts/python/requirements.txt`        | Add `numpy`, which is imported directly (next to `h5py`).                                                                                                                                                                                                                                                                                                |

## CLI

```
vtkhdf_extract.py info    FILE
vtkhdf_extract.py extract FILE (--step S | --time T) [-o OUT] [--overwrite]
vtkhdf_extract.py split   FILE (--steps-per-chunk N | --time-window W | --max-size SIZE)
                          [--dry-run] [--overwrite]
```

**`info`:**
- Prints the type and version, and whether the file is temporal.
- Prints NSteps, the time range and whether time is strictly increasing.
- Prints the mesh sizes: N, cells, P and cell types.
- Prints each field with its dtype, components and storage (printf VDS, explicit VDS or inline).
- Prints the main file size, the per-step and total size of the referenced files, and any missing per-step files (`stat` only).

**`extract`:**
- `--step` accepts negative values (`-1` is the last step).
- `--time` picks the nearest step. A tie is an error that names the candidate steps.
- The default output is `<input_dir>/<stem>_step_<s>.vtkhdf`, with `s` zero-padded to the width of `NSteps-1`.

**`split`:**
- Writes next to the input. For each chunk it prints the file name, the steps, the time range and the footprint, then the `rsync` command.
- `--time-window W` uses absolute windows `[mW, (m+1)W)`, with index `floor(t/W + 1e-9)`; the epsilon absorbs float noise. Empty windows are skipped, and strictly increasing times are required.
- `--max-size` takes bytes with an optional `K/M/G/T` suffix (powers of 1024).
- `--dry-run` validates and prints the plan but writes nothing.

## Backend design (`vtkhdf_backend.py`)

**Types:**
- `BLOCK_BYTES = 64 MiB` is the maximum buffer size for any bulk copy.
- `VtkhdfError(Exception)`.
- `FieldInfo(name, dtype, row_shape, kind, names, dset)`:
  - `kind` is one of `printf`, `explicit` or `inline`;
  - `names` holds the stored VDS file name or names.
- `FileInfo(path, temporal, n_steps, times, n_parts, n_points, n_cells, n_conn, part_counts, fields)`.

**Functions:**
- `open_input(path)`: `h5py.File(path, "r", locking="best-effort")`. This is the only way inputs and source files are opened.
- `inspect(f, path) -> FileInfo`: reads and validates the metadata (see below).
- `step_source(info, field, step) -> (Path, dset, first_row)`: finds where a step's rows live.
  - `printf`: substitute `%b → step` and `%% → %`.
  - `explicit`: use mapping number `step`.
  - `inline`: rows `step·N` of the main-file dataset.
  - Relative names resolve against the input's directory, which is HDF5's own rule.
- `check_sources(info, steps) -> list[(step, files, nbytes)]`:
  - opens each referenced file `"r"` one at a time;
  - checks each field's dataset holds at least `first_row + N` rows of the right row shape;
  - collects the unique files and their on-disk bytes.
- `find_step(info, step=None, time=None) -> int`.
- `plan_by_count(n_steps, k)`, `plan_by_time(times, width)` and `plan_by_size(step_bytes, fixed_bytes, limit)`, each returning `list[range]`.
- `write_extract(f, info, step, out, overwrite)`.
- `write_chunk(f, info, steps, out, overwrite) -> list[Path]` (the referenced files).
- `write_manifest(path, entries, overwrite)`.
- Helpers:
  - `copy_attrs(src, dst, skip=...)`: byte-exact copy using `h5a.create` with the source `TypeID` and `mtype=ftype`.
  - `copy_rows(src, row0, nrows, dst)`: one reusable buffer of at most `BLOCK_BYTES`, moved with `read_direct` and `write_direct`.
  - `copy_mesh(src_vtkhdf, dst_vtkhdf)`.
  - `atomic_output(path, overwrite)`, a context manager.
  - `parse_size`, `format_size`.

**`inspect` rules:**
- Only `[:NSteps]` and `[:P]` slices of the metadata arrays are read. `NumberOf*` can hold `NSteps·P` entries, so they are never read whole.
- `.shape`, `.ndim`, `len()`, `repr()` and slices of an input PointData VDS are never used. Only `ds.dtype`, `ds.is_virtual` and `ds.id.get_create_plist()` are allowed.

Each failure below raises a `VtkhdfError` that states why:
1. `/VTKHDF` must exist with `Type == "UnstructuredGrid"`.
2. `temporal` means a `Steps` group exists.
   - `NSteps` is authoritative, and `Values[:NSteps]` gives the times.
   - `NumberOfParts[:NSteps]` must be constant; that value is P.
   - The part, point, cell and connectivity offsets must all be 0. Otherwise the error is "time-varying mesh/partitioning not supported".
3. N, C and K are the sums of `NumberOfPoints`, `NumberOfCells` and `NumberOfConnectivityIds` over `[:P]`. Then:
   - `Points` must be `(N,3)`;
   - `Types` must be `(C,)`;
   - `Offsets` must be `(C+P,)`;
   - `Connectivity` must be `(K,)`.
4. `CellData` containing datasets is an error (Neko never writes it). An input `FieldData` group is ignored, with a warning unless it only holds `Time`.
5. For each PointData field, `PointDataOffsets/<f>[:NSteps]` must equal `arange·N`. The storage is classified from the dcpl:
   - **printf**: one mapping whose name contains `%b`, a regular hyperslab with start 0, stride = block = N rows, count `h5s.UNLIMITED`, and source selection ALL.
   - **explicit**: NSteps mappings, where mapping i has select bounds of exactly rows `[iN,(i+1)N)` and source selection ALL.
   - **inline**: a non-VDS dataset, as written by Neko before 2026-03.
   - Anything else is unsupported.

   `row_shape` comes from `get_virtual_vspace(0).shape[1:]`, which needs no probing. For inline datasets it comes from the dataset shape.
6. Names containing `%` other than `%b` are refused. If `HDF5_VDS_PREFIX` is set, the tool warns, because HDF5 would then resolve sources differently.

### `extract` output

The file is non-temporal and contains no VDS.

| Object                                       | Content                                                                                                               |
| -------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| `/VTKHDF` attrs                              | Byte-exact copies of the input `/VTKHDF` attributes (`Type`, `Version`)                                               |
| `Points`, `Connectivity`, `Offsets`, `Types` | `copy_rows` streams of the input (same dtype, contiguous)                                                             |
| `NumberOfPoints/Cells/ConnectivityIds`       | The P per-part entries, exactly P (VTK counts pieces by length)                                                       |
| `PointData/<f>`                              | Shape `(N,)+row_shape` with the input dtype, streamed from `step_source(f, s)`. Attributes copied per the tag policy. |
| `FieldData/Time`                             | f64 of shape `(1,)` holding `times[s]`. It must not be scalar, because VTK rejects scalar dataspaces.                 |
| `/` attrs                                    | Provenance: the source file name and step index. Readers ignore these.                                                |

### `split` chunk output (`<stem>_steps_<a>-<b>.vtkhdf`, k steps, zero-padded indices)

| Object                | Content                                                                                                                                                                                                                                                                                                                     |
| --------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `/VTKHDF` attrs, mesh | Same as `extract`                                                                                                                                                                                                                                                                                                           |
| `NumberOf*`           | `k·P` entries: the per-part entries repeated for each step, as Neko does. Written per step, so never as one huge array.                                                                                                                                                                                                     |
| `Steps`               | `NSteps` int32 = k. `Values` f64 holds the original times, not rebased. `NumberOfParts` = P. The other offsets are 0, and `PointDataOffsets/<f>` = `arange(k)·N`. All int64, length k.                                                                                                                                      |
| `PointData/<f>`       | A VDS with fixed shape `(k·N,)+row_shape`, the input dtype, and fill value NaN for floats, so a missing source shows up. Mapping j covers rows `[jN,(j+1)N)` and points to `relpath(source(a+j), input_dir)` : `dset`, select ALL. Mappings are assigned in ascending j, because Neko's reader needs mapping index == step. |
| `FieldData/Time`      | Only when k == 1: VTK reads 1-step files as static and would otherwise lose the time.                                                                                                                                                                                                                                       |
| `/` attrs             | Provenance: the source file and steps a..b                                                                                                                                                                                                                                                                                  |

- The manifest `<stem>_steps_<a>-<b>.files` contains the chunk file, then each referenced file in step order, one per line, relative to the input directory.
- Every source must lie inside the input directory, i.e. its `relpath` has no `..`. Otherwise `--files-from` cannot ship it, and split refuses.
- `split` also refuses inline layouts: the data lives in the main file, so there is nothing to reference.

### Chunk planning

- **count**: consecutive runs of N steps.
- **time**: group the steps by window index.
- **size**: greedy packing, with footprint estimated as `mesh payload + 1 MiB + Σ step-file bytes + 256 B per step·field` (the VDS mapping is about 90 B).
  - If one step alone exceeds the limit, split refuses.
  - After writing, the tool prints the actual footprint (`getsize(chunk)` plus the files) and warns if it is over the limit.
- **All modes:**
  - Validate every source with `check_sources` before writing anything. HDF5 would otherwise silently fill missing data.
  - Warn, without failing, if the times are not strictly increasing. Neko-TOP optimizer restarts can write `t=0` mid-file, and VTK reads the original the same way.

### Safety, memory, compatibility

**Originals are never modified:**
- No code path opens an input or a referenced file other than through `open_input("r")`.
- Data is never read through a file that is open for writing.

**Outputs:**
- If a planned output exists, the tool refuses unless `--overwrite` is given.
- Even with `--overwrite`, it refuses any output that resolves (`realpath`) to the input or to a referenced source.
- Each file is written to `.<name>.tmp-<pid>` in its destination directory (h5py mode `"x"`), closed, then moved with `os.replace`. The temp file is deleted on error. Each manifest is written after its chunk, in the same way.

**Memory:**
- Bulk data moves only through `copy_rows`, so peak memory is about `BLOCK_BYTES`.
- Metadata reads are O(NSteps + P).
- `split` copies no field data. It adds about 90 B of metadata per step per field.
- Per-step files are opened one at a time and closed, so file descriptors do not accumulate.

**Compatibility:**
- `libver` keeps the h5py default, so outputs stay readable by HDF5 1.10–1.14, including ParaView 6.2.
- Chunk files are readable by VTK 9.3–9.7 (the same temporal semantics) and by Neko's own reader.
- Extracts can be read by Neko's Brinkman `file` source term, e.g. to start a run from a design.

### `Attribute` tag policy (pending the user's answer)

Recommendation: in outputs, drop `Attribute="Scalars"` and keep the other tags (`Vectors`). Every array is then visible in ParaView 6.2. Older VTK ignores the tags, and Neko never reads them. The originals are untouched either way.

## Verification (during implementation)

Scratch work dir: `/tmp/claude-1000/-home-tife-Projects-neko-top/25140904-a3a4-4939-8327-c225ddb71c3d/scratchpad/vtkhdf_check`.

Copy into it with `cp -a`:
- `vtkhdf_test_data/{design_0,forward_fields_0}` (`.vtkhdf` + `.data`);
- `external/neko/tests/unit/test_write_time_0{.vtkhdf,.data}`, which has P=2 and f8 data.

The tool must never be run inside `vtkhdf_test_data/`.

0. Record `sha256sum` and mtimes of everything in `vtkhdf_test_data/` and the copies. Re-check at the end: they must be identical.
1. Run `info` on all three. It must report 120 / 121 / 2 steps, P = 1/1/2, the fields, and no missing files.
2. **`extract`:**
   - Runs: `design_0 --step 5`, `--step -1`, `--time 57.4` (→ 57); `forward_fields_0 --step 120`; `test_write_time_0 --step 1`.
   - h5py checks, with the file opened `"r"`:
     - no dataset is virtual;
     - each PointData array is `np.array_equal` to the per-step file;
     - the mesh equals the original;
     - `NumberOf*` has length P;
     - `FieldData/Time` is correct;
     - there is no `Steps` group.
3. **`split`** (in scratch):
   - Runs:
     - `design_0 --steps-per-chunk 25` gives 0-24 … 100-119;
     - `--time-window 30 --overwrite`;
     - `--max-size 40M --dry-run`, then for real;
     - `forward_fields_0 --steps-per-chunk 60` gives a final 1-step chunk 120-120;
     - `test_write_time_0 --steps-per-chunk 1` gives P=2 chunks of 1 step each.
   - Checks on the files, reopened `"r"`:
     - each chunk's VDS rows equal the per-step files;
     - mapping i points to `<stem>.data/<a+i>.h5` with no `%b`;
     - the `Steps` arrays match the table above;
     - the manifest lists exactly the chunk plus its files;
     - the chunk size is about the mesh size.
4. **Explicit-mapping path:** re-split a chunk with `--steps-per-chunk 5`, and extract a step from a chunk. Both must equal the original data.
5. **Transfer:**
   - Run `rsync -a --files-from=<manifest> $S/ $S/remote/`, then read the chunk from `remote/` with h5py, using a different cwd. The data must match.
   - Delete one `.h5` file in `remote/`. That step must read as NaN, and `info` must list it as missing.
6. **VTK end-to-end:**
   - Use `pvpython`: `/usr/local/paraview/bin/pvpython` if it has been reinstalled (it is currently empty), otherwise `/mnt/c/Program Files/ParaView 6.2.0/bin/pvpython.exe` with `wslpath -w` paths and `WSLENV=HDF5_USE_FILE_LOCKING/u` set to `FALSE`.
   - Create a fresh `vtkHDFReader` per file, because its cache goes stale across files.
   - Compare the original, the chunks and the extracts: the time values, and the per-step array sums for each array. For the originals, enable one array at a time because of the tag issue.
   - P=2 must give 2 partitions, and a 1-step chunk must read as static with FieldData `Time`.
7. **Bounded memory:** import the backend, set `BLOCK_BYTES = 4096`, and rerun 2 and 3. The outputs must be byte-identical in data, which exercises the block boundaries. Also record `/usr/bin/time -v` max RSS for a normal run.
8. **Error paths:** each must fail with a clear message, without a traceback, and leave no temp file behind:
   - a missing per-step file (deleted in a scratch copy);
   - an existing output without `--overwrite`;
   - `-o` pointing at the input, even with `--overwrite`;
   - a non-temporal input (`external/neko/tests/unit/test_write_0.vtkhdf`);
   - `--max-size` smaller than one step.

## Related findings (not part of this change)

- **Neko writer tags.** Neko tags every scalar as `Scalars` (`vtkhdf_file.F90:1046-1058`), so in ParaView 6.2 `design_0` shows only `brinkman_amplitude` and `forward_fields_0` loses `Pressure`. This affects `visualisation/mixer/design_evolution.py`. The proper fix is upstream in Neko.
- **Neko-TOP sensitivity stage names.** `mapping_handler_set_stage_names` (`sources/mapping_functions/mapping_handler.f90:178`) runs before the sensitivity stages exist (`:201-204`), so the stages are named `Field` and VTKHDF's duplicate-skip drops two of them.
