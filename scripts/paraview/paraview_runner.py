#!/usr/bin/env pvpython
"""Load a ParaView state file (.pvsm) and export an animation of all frames.

The output frames are stored in a subdirectory of the data directory, named 
after the state file.

Usage:
    pvpython paraview_runner.py state.pvsm data_dir

Additional arguments:
    --resolution H: Set the height of the output images (default: 360)
    --frames START END: Export only frames from START to END (inclusive)
    --overwrite: Overwrite existing files in the output directory
"""

import argparse
import os

from paraview.simple import (
    GetActiveViewOrCreate,
    GetAnimationScene,
    LoadState,
    LoadPlugin,
    SaveAnimation,
)


def main():
    parser = argparse.ArgumentParser(description="Run a ParaView animation export.")
    parser.add_argument("state", help="Path to the .pvsm state file")
    parser.add_argument("data_dir", help="Path to the data directory")
    parser.add_argument(
        "--resolution", type=int, default=360, help="Height of the output images (default: 360)"
    )
    parser.add_argument(
        "--frames", nargs=2, type=int, default=None, metavar=("START", "END")
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing files"
    )
    args = parser.parse_args()

    if not os.path.isfile(args.state):
        raise FileNotFoundError(args.state)
    if not os.path.isdir(args.data_dir):
        raise NotADirectoryError(args.data_dir)
    
    LoadPlugin("VTKmFilters", remote=True)

    LoadState(
        args.state,
        data_directory=args.data_dir,
        restrict_to_data_directory=True,
    )

    view = GetActiveViewOrCreate("RenderView")
    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    scene.PlayMode = "Snap To TimeSteps"

    original_size = list(view.ViewSize)
    original_ppi = view.PPI

    scale = args.resolution / original_size[1]
    output_resolution = [int(original_size[0] * scale), args.resolution]

    view.ViewSize = output_resolution
    view.PPI = round(original_ppi * scale)
    view.UpdateVTKObjects()

    # Set the frame index range to export.
    frame_count = len(scene.TimeKeeper.TimestepValues)
    if args.frames:
        start_frame, end_frame = args.frames
    else:
        start_frame, end_frame = 0, frame_count - 1

    # Determine the output location and ensure the directory exists
    state_name = os.path.splitext(os.path.basename(args.state))[0]
    output_dir = os.path.join(args.data_dir, "visualisations", state_name)
    os.makedirs(output_dir, exist_ok=True)

    # Define the output file path for the animation frames
    animation_file = os.path.join(output_dir, "frame.png")

    if not args.overwrite:
        frame_list = list(range(start_frame, end_frame + 1))
        frame_windows = []
        window_start = frame_list[0]
        for frame_index in frame_list:
            if os.path.isfile(os.path.join(output_dir, f"frame_{frame_index:04d}.png")):
                if window_start < frame_index:
                    frame_windows.append((window_start, frame_index - 1))
                window_start = frame_index + 1
        if window_start <= frame_list[-1]:
            frame_windows.append((window_start, frame_list[-1]))
    else:
        frame_windows = [(start_frame, end_frame)]

    # Export the animation frames for each frame window
    for frame_window in frame_windows:
        SaveAnimation(
            animation_file,
            view,
            FrameWindow=(frame_window[0], frame_window[1]),
            ImageResolution=output_resolution,
            SuffixFormat="_{:04d}",
        )


if __name__ == "__main__":
    main()
