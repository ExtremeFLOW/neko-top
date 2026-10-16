#!/usr/bin/env python3
"""Repeat a one-element-thick Nek5000/Neko fld slab in the spanwise direction."""

from __future__ import annotations

import argparse
import math
import struct
from pathlib import Path


HEADER_SIZE = 132
TEST_PATTERN = 6.54321


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a spanwise-extruded 3D baseflow field by repeating a slab "
            "fld file. The input is expected to contain X, U, and P."
        )
    )
    parser.add_argument("input", type=Path, help="input slab field, e.g. baseflow0.f00000")
    parser.add_argument("output", type=Path, help="output 3D field, e.g. baseflow3d0.f00000")
    parser.add_argument("--length", type=float, default=4.0, help="spanwise length")
    parser.add_argument("--layers", type=int, default=8, help="number of spanwise element layers")
    return parser.parse_args()


def parse_header(raw: bytes) -> dict[str, object]:
    text = raw.decode("ascii")
    parts = text.split()
    if len(parts) < 12 or parts[0] != "#std":
        raise ValueError("Unsupported fld header")
    return {
        "word_size": int(parts[1]),
        "lx": int(parts[2]),
        "ly": int(parts[3]),
        "lz": int(parts[4]),
        "nel": int(parts[5]),
        "time": float(parts[7].replace("D", "E")),
        "counter": int(parts[8]),
        "rdcode": parts[11],
    }


def write_header(info: dict[str, object], nel: int) -> bytes:
    rdcode = str(info["rdcode"])
    header = (
        f"#std {int(info['word_size']):1d} "
        f"{int(info['lx']):2d} {int(info['ly']):2d} {int(info['lz']):2d} "
        f"{nel:10d} {nel:10d} "
        f"{float(info['time']):20.13E} "
        f"{int(info['counter']):9d} {1:6d} {1:6d} "
        f"{rdcode:<10s}"
    )
    return header.ljust(HEADER_SIZE).encode("ascii")


def unpack_values(data: bytes, endian: str, word_size: int) -> tuple[float, ...]:
    code = "f" if word_size == 4 else "d"
    return struct.unpack(endian + str(len(data) // word_size) + code, data)


def pack_values(values: list[float], endian: str, word_size: int) -> bytes:
    code = "f" if word_size == 4 else "d"
    return struct.pack(endian + str(len(values)) + code, *values)


def minmax_metadata(
    chunks: list[tuple[tuple[float, ...], ...]],
    endian: str,
) -> bytes:
    values: list[float] = []
    for fields in chunks:
        for field in fields:
            values.extend([float(min(field)), float(max(field))])
    return struct.pack(endian + str(len(values)) + "f", *values)


def main() -> None:
    args = parse_args()
    if args.layers <= 0:
        raise ValueError("--layers must be positive")
    if args.length <= 0.0:
        raise ValueError("--length must be positive")

    raw = args.input.read_bytes()
    info = parse_header(raw[:HEADER_SIZE])
    word_size = int(info["word_size"])
    if word_size not in (4, 8):
        raise ValueError(f"Unsupported fld word size: {word_size}")

    little = struct.unpack("<f", raw[HEADER_SIZE:HEADER_SIZE + 4])[0]
    big = struct.unpack(">f", raw[HEADER_SIZE:HEADER_SIZE + 4])[0]
    if math.isclose(little, TEST_PATTERN, rel_tol=0.0, abs_tol=1.0e-5):
        endian = "<"
    elif math.isclose(big, TEST_PATTERN, rel_tol=0.0, abs_tol=1.0e-5):
        endian = ">"
    else:
        raise ValueError("Could not identify fld endian test pattern")

    rdcode = str(info["rdcode"])
    if rdcode != "XUP":
        raise ValueError(f"Expected rdcode XUP, got {rdcode!r}")

    lx = int(info["lx"])
    ly = int(info["ly"])
    lz = int(info["lz"])
    old_nel = int(info["nel"])
    lxyz = lx * ly * lz
    new_nel = old_nel * args.layers

    offset = HEADER_SIZE + 4
    idx_size = old_nel * 4
    data_offset = offset + idx_size

    vector_bytes_per_el = 3 * lxyz * word_size
    scalar_bytes_per_el = lxyz * word_size
    x_offset = data_offset
    u_offset = x_offset + old_nel * vector_bytes_per_el
    p_offset = u_offset + old_nel * vector_bytes_per_el

    dz = args.length / args.layers
    out = bytearray()
    out += write_header(info, new_nel)
    out += raw[HEADER_SIZE:HEADER_SIZE + 4]
    out += struct.pack(endian + str(new_nel) + "i", *range(1, new_nel + 1))

    x_chunks: list[tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...]]] = []
    u_chunks: list[tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...]]] = []
    p_chunks: list[tuple[tuple[float, ...]]] = []

    for layer in range(args.layers):
        for el in range(old_nel):
            start = x_offset + el * vector_bytes_per_el
            vals = unpack_values(raw[start:start + vector_bytes_per_el], endian, word_size)
            x = vals[:lxyz]
            y = vals[lxyz:2 * lxyz]
            z_old = vals[2 * lxyz:]
            z_min = min(z_old)
            z_max = max(z_old)
            scale = dz / (z_max - z_min) if z_max > z_min else 0.0
            z = tuple(layer * dz + (zz - z_min) * scale for zz in z_old)
            x_chunks.append((x, y, z))
            out += pack_values(list(x) + list(y) + list(z), endian, word_size)

    zeros = [0.0] * lxyz
    for layer in range(args.layers):
        for el in range(old_nel):
            start = u_offset + el * vector_bytes_per_el
            vals = unpack_values(raw[start:start + vector_bytes_per_el], endian, word_size)
            u = vals[:lxyz]
            v = vals[lxyz:2 * lxyz]
            w = tuple(zeros)
            u_chunks.append((u, v, w))
            out += pack_values(list(u) + list(v) + list(w), endian, word_size)

    for layer in range(args.layers):
        for el in range(old_nel):
            start = p_offset + el * scalar_bytes_per_el
            p = unpack_values(raw[start:start + scalar_bytes_per_el], endian, word_size)
            p_chunks.append((p,))
            out += pack_values(list(p), endian, word_size)

    out += minmax_metadata(x_chunks, endian)
    out += minmax_metadata(u_chunks, endian)
    out += minmax_metadata(p_chunks, endian)

    args.output.write_bytes(out)
    print(
        f"Wrote {args.output} with {new_nel} elements "
        f"({args.layers} layers over length {args.length:g})."
    )


if __name__ == "__main__":
    main()
