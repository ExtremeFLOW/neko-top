#!/usr/bin/env python3
"""Extrude a *3D one-layer* Nek5000/Neko XUP fld to multiple z layers.

The source is already three-dimensional (lz > 1), with an element layer
spanning z0..z1. We preserve the GLL dimensions and write an ordinary 3D
Neko/Nek5000 field file; this does NOT generate the solver's connectivity mesh.

Usage:
    python extrude_fld.py baseflow0.f00000 baseflow3d0.f00000 \
        --length 4 --layers 4

Default geometry is 'flat': x,y come from the lower z face and are constant
along z. This is appropriate for an extrusion of a nominally 2D geometry.
The 'periodic' mode instead preserves the original interior mesh warping,
while removing the linear mismatch between the two spanwise faces.

Field ordering: global IDs, all elements' X, all elements' U, all elements'
P, then per-element min/max for X, U, P. This supports single-file XUP only.
"""
from __future__ import annotations

import argparse
import math
import re
import struct
from pathlib import Path

import numpy as np

HEADER_SIZE = 132
TEST_PATTERN = 6.54321


def get_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--length', type=float, default=4.0, help='output spanwise length (default 4)')
    parser.add_argument('--layers', type=int, default=4, help='output element layers (default 4)')
    parser.add_argument('--z-min', type=float, default=0.0, help='output starting z (default 0)')
    parser.add_argument('--geometry', choices=('flat', 'periodic', 'raw'), default='flat',
                        help='flat: straight extrusion (default); periodic: keep interior warping while matching spanwise faces; raw: copy original x,y without modification')
    parser.add_argument('--keep-w', action='store_true', help='keep original spanwise velocity; default w=0 as in the original script')
    parser.add_argument('--seam-tolerance', type=float, default=1.0e-5, help='tolerance for warnings about mismatched spanwise faces')
    parser.add_argument('--no-meta', action='store_true',
                        help='do not write the accompanying .nek5000 series metadata file')
    return parser.parse_args()


def check_header(raw: bytes):
    if len(raw) != HEADER_SIZE or raw[:4] != b'#std':
        raise ValueError('Expected a standard 132-byte Nek5000/Neko header')
    parts = raw.decode('ascii').split()
    if len(parts) < 12:
        raise ValueError(f'Incomplete header: {parts!r}')
    word_size, lx, ly, lz, nel, nel_global = map(int, parts[1:7])
    if word_size not in (4, 8) or min(lx, ly) < 2 or lz < 2 or nel < 1:
        raise ValueError('Expected a 3D input with lz>=2 and word size 4 or 8')
    if parts[11] != 'XUP':
        raise ValueError(f'Only XUP is supported; input has {parts[11]!r}')
    if nel != nel_global or int(parts[9]) != 1 or int(parts[10]) != 1:
        raise ValueError('Only single-file fld output is supported')
    return word_size, lx, ly, lz, nel


def max_seam_difference(array: np.ndarray) -> float:
    """Max top-bottom discrepancy across each element for array (...,lz,ly,lx)."""
    return float(np.max(np.abs(array[..., -1, :, :] - array[..., 0, :, :])))


def minmax(array: np.ndarray, endian: str) -> np.ndarray:
    """Return per-element, per-component min,max pairs in Nek metadata format."""
    # Vector: (nel,3,lz,ly,lx); pressure is promoted to (nel,1,lz,ly,lx).
    f = array.reshape(array.shape[0], array.shape[1], -1)
    result = np.stack((f.min(axis=2), f.max(axis=2)), axis=2)
    return result.astype(np.dtype(endian + 'f4'), copy=False)


def series_metadata_path(output: Path):
    """Return Neko sidecar path and file-series info for `name0.f00000`."""
    match = re.fullmatch(r'(.+)0\.f(\d{5})', output.name)
    if match is None:
        raise ValueError('Output filename must have the Neko single-file pattern '
                         '`name0.f00000` to create a series metadata file. '
                         'Use --no-meta to disable metadata generation.')
    basename, counter = match.group(1), int(match.group(2))
    return output.with_name(f'{basename}0.nek5000'), basename, counter


def write_series_metadata(output: Path):
    """Write the text .nek5000 sidecar Neko uses to resolve fld names."""
    meta_path, basename, counter = series_metadata_path(output)
    # For example: baseflow3d0.f00000 -> baseflow3d0.nek5000,
    # with filetemplate: baseflow3d%01d.f%05d
    contents = (f'filetemplate: {basename}%01d.f%05d\n'
                f'firsttimestep: {counter}\n'
                'numtimesteps: 1\n')
    meta_path.write_text(contents, encoding='ascii')
    return meta_path


def main():
    args = get_args()
    if args.input.resolve() == args.output.resolve():
        raise ValueError('Input and output paths must differ')
    if not args.no_meta:
        series_metadata_path(args.output)  # Validate before creating output
    if args.layers < 1 or not math.isfinite(args.length) or args.length <= 0:
        raise ValueError('Need --layers >= 1 and --length > 0')
    if not math.isfinite(args.z_min):
        raise ValueError('--z-min must be finite')

    with args.input.open('rb') as f:
        header = f.read(HEADER_SIZE)
        testbytes = f.read(4)
    word_size, lx, ly, lz, nel = check_header(header)
    if len(testbytes) != 4:
        raise ValueError('Missing endian test pattern')
    if math.isclose(struct.unpack('<f', testbytes)[0], TEST_PATTERN, abs_tol=1e-5):
        endian = '<'
    elif math.isclose(struct.unpack('>f', testbytes)[0], TEST_PATTERN, abs_tol=1e-5):
        endian = '>'
    else:
        raise ValueError('Invalid endian test pattern')

    dtype = np.dtype(endian + ('f4' if word_size == 4 else 'f8'))
    lxyz = lx * ly * lz
    ids_start = HEADER_SIZE + 4
    fields_start = ids_start + 4 * nel
    xyz_bytes = nel * 3 * lxyz * word_size
    p_bytes = nel * lxyz * word_size
    metadata_start = fields_start + 2 * xyz_bytes + p_bytes
    expected_size = metadata_start + nel * 14 * 4
    if args.input.stat().st_size != expected_size:
        raise ValueError(f'Unexpected input size: got {args.input.stat().st_size}, expected {expected_size}')

    source_ids = np.memmap(args.input, dtype=np.dtype(endian + 'i4'), mode='r',
                           offset=ids_start, shape=(nel,))
    if not np.array_equal(np.sort(source_ids), np.arange(1, nel + 1)):
        raise ValueError('Input element IDs are not a permutation of 1..nel')
    mesh = np.memmap(args.input, dtype=dtype, mode='r', offset=fields_start,
                     shape=(nel, 3, lz, ly, lx))
    vel = np.memmap(args.input, dtype=dtype, mode='r', offset=fields_start + xyz_bytes,
                    shape=(nel, 3, lz, ly, lx))
    pres = np.memmap(args.input, dtype=dtype, mode='r',
                     offset=fields_start + 2 * xyz_bytes, shape=(nel, 1, lz, ly, lx))
    in_meta = np.memmap(args.input, dtype=np.dtype(endian + 'f4'), mode='r',
                        offset=metadata_start, shape=(nel * 14,))
    if not np.all(np.isfinite(mesh)) or not np.all(np.isfinite(vel)) or not np.all(np.isfinite(pres)):
        raise ValueError('Input field has NaNs/Infs')
    expected_in_meta = np.concatenate([
        minmax(mesh, endian).ravel(), minmax(vel, endian).ravel(), minmax(pres, endian).ravel()
    ])
    if not np.allclose(in_meta, expected_in_meta, atol=1.e-5, rtol=1.e-5):
        raise ValueError('Input metadata does not agree with stored field blocks')

    z_source = mesh[:, 2]
    z0 = float(np.min(z_source))
    z1 = float(np.max(z_source))
    if z1 <= z0:
        raise ValueError('Source has zero spanwise thickness')
    if np.max(np.abs(z_source[:, 0, :, :] - z0)) > args.seam_tolerance or \
       np.max(np.abs(z_source[:, -1, :, :] - z1)) > args.seam_tolerance:
        raise ValueError('Input must consist of a single layer with common bottom and top z faces')

    # Input XYZ and field blocks may have arbitrary input element-ID order.
    # Canonicalise to element IDs 1..nel before creating the layer-major output,
    # so the generated file maps to 1..(nel*layers) in n2to3 ordering.
    src_order = np.argsort(source_ids)
    mesh = np.asarray(mesh[src_order])
    vel = np.asarray(vel[src_order])
    pres = np.asarray(pres[src_order])

    geometry_seam = max(max_seam_difference(mesh[:, 0]), max_seam_difference(mesh[:, 1]))
    flow_seam = max(max_seam_difference(vel[:, c]) for c in range(2))
    pressure_seam = max_seam_difference(pres[:, 0])
    print(f'Input: {nel} elements, {lx}x{ly}x{lz} GLL points, source z=[{z0:g},{z1:g}]')
    print(f'Input face mismatch: XY={geometry_seam:.7g}, UV={flow_seam:.7g}, P={pressure_seam:.7g}')
    if args.geometry == 'raw' and geometry_seam > args.seam_tolerance:
        print('WARNING: raw geometry has discontinuous XY faces between replicated layers!')

    # Construct one normalised geometry block, later shift/rescale its z.
    coords = mesh.copy()
    if args.geometry == 'flat':
        coords[:, :2] = mesh[:, :2, 0:1, :, :]
    elif args.geometry == 'periodic':
        t = (mesh[:, 2] - z0) / (z1 - z0)
        delta = mesh[:, :2, -1, :, :] - mesh[:, :2, 0, :, :]
        coords[:, :2] -= delta[:, :, None, :, :] * t[:, None, :, :, :]
    # raw: preserve XY as-is.
    coords[:, 2] = (mesh[:, 2] - z0) / (z1 - z0)

    out_uv = vel.copy()
    if not args.keep_w:
        out_uv[:, 2, :, :, :] = 0
    out_p = pres.copy()

    # Reuse original header bytes and patch only the global/local element counts.
    # This avoids modifying other fixed-width fields such as time and counter.
    new_nel = nel * args.layers
    if new_nel > 2_147_483_647:
        raise ValueError('Too many elements for 32-bit global IDs')
    hdr = bytearray(header)
    hdr[16:26] = f'{new_nel:10d}'.encode('ascii')
    hdr[27:37] = f'{new_nel:10d}'.encode('ascii')
    check_header(bytes(hdr))

    f32 = np.dtype(endian + 'f4')
    i32 = np.dtype(endian + 'i4')
    dz = args.length / args.layers
    with args.output.open('wb') as f:
        f.write(hdr)
        f.write(testbytes)
        np.arange(1, new_nel + 1, dtype=i32).tofile(f)
        # X: one layer at a time, each with ALL its elements.
        for layer in range(args.layers):
            xyz = coords.copy()
            xyz[:, 2] = args.z_min + layer * dz + coords[:, 2] * dz
            xyz.tofile(f)
        # U block across all layers.
        for layer in range(args.layers):
            out_uv.tofile(f)
        # P block across all layers.
        for layer in range(args.layers):
            out_p.tofile(f)
        # Min/max: X block, then U block, then P block, all in element-major order.
        for layer in range(args.layers):
            xyz = coords.copy()
            xyz[:, 2] = args.z_min + layer * dz + coords[:, 2] * dz
            minmax(xyz, endian).tofile(f)
        meta_u = minmax(out_uv, endian)
        for layer in range(args.layers):
            meta_u.tofile(f)
        meta_p = minmax(out_p, endian)
        for layer in range(args.layers):
            meta_p.tofile(f)

    size = args.output.stat().st_size
    expected_output_size = 136 + 4 * new_nel + new_nel * (7 * lxyz * word_size + 14 * 4)
    if size != expected_output_size:
        raise RuntimeError(f'Incorrect output length: {size} vs {expected_output_size}')

    output_geometry_seam = max(max_seam_difference(coords[:, 0]), max_seam_difference(coords[:, 1]))
    print(f'Wrote: {args.output} ({size:,} bytes), {new_nel} elements, '
          f'z=[{args.z_min:g},{args.z_min + args.length:g}]')
    print(f'Replicated geometry face mismatch: {output_geometry_seam:.7g}')
    if output_geometry_seam > args.seam_tolerance:
        print('WARNING: output field geometry has mismatched adjacent spanwise faces.')
    if max(flow_seam, pressure_seam) > args.seam_tolerance:
        print('WARNING: repeating these field values also creates discontinuities in velocity/pressure.')
    if not args.no_meta:
        meta_path = write_series_metadata(args.output)
        print(f'Wrote field-series metadata: {meta_path}')
    print('IMPORTANT: the solver mesh must have matching geometry and layer-major element ordering.')


if __name__ == '__main__':
    main()
