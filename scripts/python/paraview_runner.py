#!/usr/bin/env pvpython
"""Load a ParaView state file (.pvsm) and export an animation of all frames.

Usage:
    pvpython paraview_runner.py state.pvsm data_dir output.png [--resolution H] [--frames START END]
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
        "--resolution", nargs=1, type=int, default=1080, metavar=("H")
    )
    parser.add_argument(
        "--frames", nargs=2, type=int, default=None, metavar=("START", "END")
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

    # Set the view size based on the resolution argument retaining the aspect ratio of the original view
    original_size = view.ViewSize
    aspect_ratio = original_size[0] / original_size[1]
    new_height = args.resolution[0]
    new_width = int(new_height * aspect_ratio)
    view.ViewSize = [new_width, new_height]
    output_resolution = view.ViewSize

    # Set the animation bounds
    if args.frames:
        start_frame, end_frame = args.frames
        scene.StartTime = start_frame
        scene.EndTime = end_frame
    else:
        scene.StartTime = scene.TimeKeeper.TimestepValues[0]
        scene.EndTime = scene.TimeKeeper.TimestepValues[-1]


    # Determine the output location and ensure the directory exists
    state_name = os.path.splitext(os.path.basename(args.state))[0]
    output_dir = os.path.join(args.data_dir, "visualisations", state_name)
    os.makedirs(output_dir, exist_ok=True)

    # Define the output file path for the animation frames
    animation_file = os.path.join(output_dir, "frame.png")

    # Export all frames, start to end.
    SaveAnimation(
        animation_file,
        view,
        ImageResolution=output_resolution,
        SuffixFormat="_{:04d}",
    )


if __name__ == "__main__":
    main()
