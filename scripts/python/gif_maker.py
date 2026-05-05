import re
from argparse import ArgumentParser
from collections import defaultdict
from pathlib import Path

from PIL import Image

DEFAULT_FPS = 10
TRANSPARENCY_INDEX = 255

QUALITY_PRESETS = {
    "low": {
        "colors": 256,
        "dither": Image.Dither.NONE,
        "optimize": False,
    },
    "medium": {
        "colors": 256,
        "dither": Image.Dither.FLOYDSTEINBERG,
        "optimize": True,
    },
    "high": {
        "colors": 256,
        "dither": Image.Dither.FLOYDSTEINBERG,
        "optimize": True,
    },
}


def png_files_in_order(folder: Path) -> list[Path]:
    return sorted(
        [
            p for p in folder.iterdir()
            if p.is_file() and p.suffix.lower() == ".png"
        ],
        key=lambda p: p.name.lower(),
    )


def numbered_image_files_by_pattern(
    folder: Path,
    pattern: str,
) -> dict[str, list[Path]]:
    try:
        filename_regex = re.compile(pattern, re.IGNORECASE)
    except re.error as error:
        raise ValueError(f"invalid regex pattern: {error}") from error

    grouped_files: dict[str, list[tuple[int, str, Path]]] = defaultdict(list)
    for file_path in folder.iterdir():
        if not file_path.is_file():
            continue
        if filename_regex.fullmatch(file_path.name) is None:
            continue

        number_match = re.search(r"(\d+)$", file_path.stem)
        if number_match is None:
            continue
        frame_number = int(number_match.group(1))
        prefix = file_path.stem[:number_match.start()]
        grouped_files[prefix].append(
            (frame_number, file_path.name.lower(), file_path))

    result: dict[str, list[Path]] = {}
    for prefix, files_with_numbers in grouped_files.items():
        files_with_numbers.sort(key=lambda item: (item[0], item[1]))
        result[prefix] = [item[2] for item in files_with_numbers]

    return dict(sorted(result.items(), key=lambda item: item[0].lower()))


def collapse_repeated_underscores(value: str) -> str:
    return re.sub(r"_+", "_", value)


def build_global_palette(frames: list[Image.Image],
                         colors: int) -> Image.Image:
    sample_size = (64, 64)
    strip = Image.new("RGB", (sample_size[0], sample_size[1] * len(frames)))

    for index, frame in enumerate(frames):
        sample = frame.convert("RGB")
        sample.thumbnail(sample_size, Image.Resampling.LANCZOS)
        y_offset = index * sample_size[1]
        strip.paste(sample, (0, y_offset))

    return strip.convert("P", palette=Image.Palette.ADAPTIVE, colors=colors)


def has_any_transparency(frames: list[Image.Image]) -> bool:
    for frame in frames:
        if "A" in frame.getbands():
            alpha_min, _ = frame.getchannel("A").getextrema()
            if alpha_min < 255:
                return True
        if frame.info.get("transparency") is not None:
            return True
    return False


def to_gif_palette(
    frame: Image.Image,
    palette_image: Image.Image,
    dither: int,
    transparency_index: int | None,
) -> Image.Image:
    paletted = frame.convert("RGB").quantize(palette=palette_image,
                                             dither=dither)
    if transparency_index is None:
        return paletted

    alpha = frame.convert("RGBA").getchannel("A")
    transparent_mask = alpha.point(lambda value: 255 if value < 128 else 0)
    paletted.paste(transparency_index, mask=transparent_mask)
    paletted.info["transparency"] = transparency_index
    return paletted


def make_gif_from_folder(
    folder: Path,
    output_dir: Path,
    fps: int = DEFAULT_FPS,
    colors: int = 256,
    dither: int = Image.Dither.FLOYDSTEINBERG,
    optimize: bool = True,
    output_name: str | None = None,
) -> None:
    png_files = png_files_in_order(folder)
    if not png_files:
        print(f"Skipping {folder.name}: no PNG files found")
        return

    frame_duration_ms = int(1000 / fps)
    frames = [Image.open(png_file) for png_file in png_files]
    gif_name = output_name if output_name is not None else f"{folder.name}.gif"
    output_path = output_dir / gif_name

    palette_frames: list[Image.Image] = []
    shared_palette: Image.Image | None = None
    try:
        use_transparency = has_any_transparency(frames)
        palette_colors = min(colors, 255) if use_transparency else colors
        transparency_index = TRANSPARENCY_INDEX if use_transparency else None

        shared_palette = build_global_palette(frames, palette_colors)
        palette_frames = [
            to_gif_palette(frame, shared_palette, dither, transparency_index)
            for frame in frames
        ]
        first, rest = palette_frames[0], palette_frames[1:]
        save_options = {
            "save_all": True,
            "append_images": rest,
            "duration": frame_duration_ms,
            "loop": 0,
            "optimize": optimize,
            # disposal=2 restores to background before the next frame.
            "disposal": 2,
        }
        if transparency_index is not None:
            save_options["transparency"] = transparency_index

        first.save(output_path, **save_options)
    finally:
        for frame in frames:
            frame.close()
        for frame in palette_frames:
            frame.close()
        if shared_palette is not None:
            shared_palette.close()

    print(f"Saved {output_path.name} from {len(png_files)} frame(s)")


def make_gif_from_files(
    frame_files: list[Path],
    output_path: Path,
    fps: int,
    colors: int,
    dither: int,
    optimize: bool,
) -> None:
    if not frame_files:
        print(f"Skipping {output_path.name}: no matching frame files found")
        return

    frame_duration_ms = int(1000 / fps)
    frames = [Image.open(frame_file) for frame_file in frame_files]

    palette_frames: list[Image.Image] = []
    shared_palette: Image.Image | None = None
    try:
        use_transparency = has_any_transparency(frames)
        palette_colors = min(colors, 255) if use_transparency else colors
        transparency_index = TRANSPARENCY_INDEX if use_transparency else None

        shared_palette = build_global_palette(frames, palette_colors)
        palette_frames = [
            to_gif_palette(frame, shared_palette, dither, transparency_index)
            for frame in frames
        ]
        first, rest = palette_frames[0], palette_frames[1:]
        save_options = {
            "save_all": True,
            "append_images": rest,
            "duration": frame_duration_ms,
            "loop": 0,
            "optimize": optimize,
            # disposal=2 restores to background before the next frame.
            "disposal": 2,
        }
        if transparency_index is not None:
            save_options["transparency"] = transparency_index

        first.save(output_path, **save_options)
    finally:
        for frame in frames:
            frame.close()
        for frame in palette_frames:
            frame.close()
        if shared_palette is not None:
            shared_palette.close()

    print(f"Saved {output_path} from {len(frame_files)} frame(s)")


def colors_value(value: str) -> int:
    colors = int(value)
    if not 2 <= colors <= 256:
        raise ValueError("colors must be between 2 and 256")
    return colors


def parse_args() -> tuple[Path, Path, str, int, int, int, bool]:
    parser = ArgumentParser(
        description="Create a GIF from PNG frames in a folder")
    parser.add_argument(
        "input_folder",
        help="Folder containing frame images to match",
    )
    parser.add_argument(
        "output_folder",
        nargs="?",
        default=None,
        help="Folder where generated GIF files are written",
    )
    parser.add_argument("--fps", type=int, default=DEFAULT_FPS)
    parser.add_argument(
        "--quality",
        choices=sorted(QUALITY_PRESETS),
        default="medium",
        help="Preset controlling colors, dithering, and optimization",
    )
    parser.add_argument(
        "--colors",
        type=colors_value,
        default=None,
        help=
        "Override palette size (2-256). Higher usually means better quality and bigger files.",
    )
    parser.add_argument(
        "--no-optimize",
        action="store_true",
        help="Disable GIF frame optimization",
    )
    parser.add_argument(
        "--pattern",
        required=True,
        help=("Regex pattern for filenames inside input_folder, for example: "
              "Temperature_[0-9]*.png"),
    )

    args = parser.parse_args()
    if args.fps <= 0:
        raise ValueError("fps must be greater than 0")

    input_folder = Path(args.input_folder)
    input_folder = input_folder.resolve()

    if not input_folder.exists() or not input_folder.is_dir():
        parser.error(
            f"input_folder does not exist or is not a directory: {input_folder}"
        )

    if args.output_folder is None:
        output_folder = input_folder.parent
    else:
        output_folder = Path(args.output_folder)
        output_folder = output_folder.resolve()
    output_folder.mkdir(parents=True, exist_ok=True)

    fps = args.fps
    preset = QUALITY_PRESETS[args.quality]
    colors = args.colors if args.colors is not None else preset["colors"]
    optimize = preset["optimize"] and not args.no_optimize

    return input_folder, output_folder, args.pattern, fps, colors, preset[
        "dither"], optimize


def main() -> None:
    input_folder, output_folder, pattern, fps, colors, dither, optimize = parse_args(
    )

    frames_by_prefix = numbered_image_files_by_pattern(input_folder, pattern)
    if not frames_by_prefix:
        print(
            f"Skipping {input_folder.name}: no matching JPG/PNG files for pattern '{pattern}'"
        )
        return

    for prefix, frame_files in frames_by_prefix.items():
        output_stem = collapse_repeated_underscores(
            f"{prefix}_{input_folder.name}")
        output_name = f"{output_stem}.gif"
        output_path = output_folder / output_name
        make_gif_from_files(
            frame_files,
            output_path,
            fps=fps,
            colors=colors,
            dither=dither,
            optimize=optimize,
        )


if __name__ == "__main__":
    main()
