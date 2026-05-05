import os
import sys
import paraview
from paraview.simple import *


def _check_mpi_compatibility():
    """Guard against MPI launcher/runtime mismatches.

    pvbatch bundles its own MPICH. If launched with an incompatible MPI
    (e.g. OpenMPI mpirun), each process runs in an isolated single-process
    session and all internal rank APIs silently return 0.

    Two independent strategies are used:

    1. mpi4py ABI check — if mpi4py is installed but raises an ImportError
       (not ModuleNotFoundError), its ABI is incompatible with pvbatch's
       bundled MPI, which is a reliable indicator of a launcher mismatch.
    2. Controller size check — compare the world size reported by the MPI
       launcher environment against what ParaView's process controller sees.
       A discrepancy means processes are running in isolated sessions.
    """
    _hint = ('Use a launcher that matches this ParaView build (MPICH):\n'
             '  mpiexec.mpich -n N pvbatch visualise.py')

    # Strategy 1: mpi4py ABI mismatch.
    try:
        import mpi4py  # noqa: F401
    except ModuleNotFoundError:
        pass  # mpi4py not installed — skip this check
    except ImportError as exc:
        print(
            f'ERROR: MPI stack mismatch detected (mpi4py ABI conflict: {exc}).\n{_hint}',
            file=sys.stderr,
        )
        sys.exit(1)

    # Strategy 2: compare launcher-reported world size with controller size.
    _launcher_size = None
    for _key in ('OMPI_COMM_WORLD_SIZE', 'PMI_SIZE', 'PMIX_SIZE',
                 'MV2_COMM_WORLD_SIZE', 'SLURM_NTASKS'):
        _v = os.environ.get(_key)
        if _v is not None:
            try:
                _launcher_size = int(_v)
            except ValueError:
                pass
            break

    if _launcher_size is not None and _launcher_size > 1:
        from paraview import servermanager as _sm
        _ctrl = _sm.vtkProcessModule.GetProcessModule().GetGlobalController()
        _ctrl_size = _ctrl.GetNumberOfProcesses() if _ctrl else 1
        if _ctrl_size < _launcher_size:
            print(
                f'ERROR: MPI launcher/runtime mismatch: launcher reports '
                f'{_launcher_size} processes but ParaView controller sees '
                f'{_ctrl_size}. Processes are running in isolated sessions.\n{_hint}',
                file=sys.stderr,
            )
            sys.exit(1)


_check_mpi_compatibility()
del _check_mpi_compatibility

# ----------------------------------------------------------------
# Parse command-line arguments
# ----------------------------------------------------------------
import argparse as _argparse


def _parse_args():
    p = _argparse.ArgumentParser(
        description='Visualise Neko static-mixer results with ParaView.', )
    p.add_argument(
        'input_dir',
        help='Directory containing the simulation output files.',
    )
    p.add_argument(
        '--stride',
        type=int,
        default=1,
        metavar='N',
        help='Process every Nth timestep (default: 1).',
    )
    p.add_argument(
        '--brinkman-clip',
        type=float,
        default=0.5,
        metavar='VALUE',
        help='Clip value for the Brinkman indicator field (default: 0.5).',
    )
    p.add_argument(
        '--text-color',
        choices=['white', 'black'],
        default='black',
        help='Color of all text and outline elements in the output '
        '(default: black). Each color variant is written to its own '
        'subdirectory inside the output directory.',
    )
    p.add_argument(
        '--overwrite',
        action='store_true',
        default=False,
        help='Overwrite existing output frames (default: skip existing frames, '
        'allowing interrupted runs to be resumed).',
    )
    p.add_argument(
        '--output-dir',
        default='visualisation',
        metavar='DIR',
        help='Output directory relative to input_dir (default: visualisation). '
        'The chosen text color is appended as a subdirectory so that '
        'black and white variants never overwrite each other.',
    )
    return p.parse_args()


_args = _parse_args()
del _parse_args

input_dir = os.path.realpath(_args.input_dir)
output_dir = os.path.join(input_dir, _args.output_dir, _args.text_color)
stride = _args.stride
brinkman_clip_value = _args.brinkman_clip
text_color = [1.0, 1.0, 1.0
              ] if _args.text_color == 'white' else [0.0, 0.0, 0.0]
overwrite = _args.overwrite
del _args

# state file generated using paraview version 6.1.0

paraview.compatibility.major = 6
paraview.compatibility.minor = 1

#### import the simple module from the paraview

# Load VTKm accelerated filters plugin
LoadDistributedPlugin('VTKmFilters', ns=globals())

# Ensure extractor output directory exists
os.makedirs(output_dir, exist_ok=True)

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# Create a new 'Render View'
velocityView = CreateView('RenderView')
velocityView.Set(
    ViewSize=[1920, 1080],
    AxesGrid='Grid Axes 3D Actor',
    OrientationAxesVisibility=0,
    CenterOfRotation=[2.0, 0.5, 0.5],
    CameraPosition=[4.77493228495059, -2.6157057373401145, 1.6192436657419564],
    CameraFocalPoint=[
        -0.8086255538276163, 4.860015290263143, -1.7409940566481248
    ],
    CameraViewUp=[-0.2112642288118684, 0.2650839449526201, 0.9407964326850324],
)

# Create a secondary 'Render View' — identical camera, temperature slices only
temperatureView = CreateView('RenderView')
temperatureView.Set(
    ViewSize=[1920, 1080],
    AxesGrid='Grid Axes 3D Actor',
    OrientationAxesVisibility=0,
    CenterOfRotation=[2.0, 0.5, 0.5],
    CameraPosition=[4.77493228495059, -2.6157057373401145, 1.6192436657419564],
    CameraFocalPoint=[
        -0.8086255538276163, 4.860015290263143, -1.7409940566481248
    ],
    CameraViewUp=[-0.2112642288118684, 0.2650839449526201, 0.9407964326850324],
)

SetActiveView(None)

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Annotate Time'
annotateTime1 = AnnotateTime(registrationName='AnnotateTime1')
annotateTime1.Format = 'Time: {time:.2f}'

# create a new 'VTKHDF Reader'
brinkman_field = VTKHDFReader(registrationName='Brinkman',
                              FileName=[
                                  input_dir + '/brinkman_0.vtkhdf',
                              ])
brinkman_field.PointArrayStatus = ['brinkman_indicator']

# create a new 'Merge Blocks'
brinkman_block = MergeBlocks(registrationName='BrinkmanBlock',
                             Input=brinkman_field)

# create a new 'VTKm Clip'
BrinkmanClip = VTKmClip(registrationName='BrinkmanClip', Input=brinkman_block)
BrinkmanClip.Set(
    Scalars=['POINTS', 'brinkman_indicator'],
    ClipValue=brinkman_clip_value,
)

# create a new 'VTKHDF Reader'
physics_field = VTKHDFReader(registrationName='Physics Field',
                             FileName=[input_dir + '/fields/field_0.vtkhdf'])
physics_field.PointArrayStatus = ['Velocity', 'temperature']

# create a new 'Merge Blocks'
physics_block = MergeBlocks(registrationName='PhysicsBlock',
                            Input=physics_field)

# create a new 'Slice'
Physics0 = Slice(registrationName='Physics x: 0', Input=physics_block)
Physics0.Set(
    SliceType='Plane',
    SliceOffsetValues=[0.0],
    PointMergeMethod='Uniform Binning',
)
Physics0.SliceType.Origin = [1e-10, 0.5, 0.5]
Physics0.HyperTreeGridSlicer.Origin = [2.0, 0.5, 0.5]

# create a new 'Slice'
Physics4 = Slice(registrationName='Physics x: 4', Input=physics_block)
Physics4.Set(
    SliceType='Plane',
    SliceOffsetValues=[0.0],
    PointMergeMethod='Uniform Binning',
)
Physics4.SliceType.Origin = [4.0, 0.5, 0.5]
Physics4.HyperTreeGridSlicer.Origin = [2.0, 0.5, 0.5]

# create a new 'Clip'
PhysicsClip = Clip(registrationName='PhysicsClip', Input=physics_block)
PhysicsClip.Set(
    ClipType='Plane',
    Invert=0,
    Crinkleclip=1,
)

# init the 'Plane' selected for 'ClipType'
PhysicsClip.ClipType.Set(
    Origin=[2.0, 0.5, 0.5],
    Normal=[0.0, 1.0, 0.0],
)

# init the 'Plane' selected for 'HyperTreeGridClipper'
PhysicsClip.HyperTreeGridClipper.Origin = [2.0, 0.5, 0.5]

# ----------------------------------------------------------------
# setup color maps and opacity maps before display to avoid
# redundant pipeline updates when Show() is first called
# ----------------------------------------------------------------

# get color transfer function/color map for 'temperature'
temperatureLUT = GetColorTransferFunction('temperature')
temperatureLUT.Set(
    AutomaticRescaleRangeMode='Never',
    RGBPoints=[
        # scalar, red, green, blue
        0.0,
        0.0,
        0.0,
        0.34902,
        0.03125,
        0.039216,
        0.062745,
        0.380392,
        0.0625,
        0.062745,
        0.117647,
        0.411765,
        0.09375,
        0.090196,
        0.184314,
        0.45098,
        0.125,
        0.12549,
        0.262745,
        0.501961,
        0.15625,
        0.160784,
        0.337255,
        0.541176,
        0.1875,
        0.2,
        0.396078,
        0.568627,
        0.21875,
        0.239216,
        0.454902,
        0.6,
        0.25,
        0.286275,
        0.521569,
        0.65098,
        0.28125,
        0.337255,
        0.592157,
        0.701961,
        0.3125,
        0.388235,
        0.654902,
        0.74902,
        0.34375,
        0.466667,
        0.737255,
        0.819608,
        0.375,
        0.572549,
        0.819608,
        0.878431,
        0.40625,
        0.654902,
        0.866667,
        0.909804,
        0.4375,
        0.752941,
        0.917647,
        0.941176,
        0.46875,
        0.823529,
        0.956863,
        0.968627,
        0.5,
        0.941176,
        0.984314,
        0.988235,
        0.5,
        0.988235,
        0.960784,
        0.901961,
        0.52,
        0.988235,
        0.945098,
        0.85098,
        0.54,
        0.980392,
        0.898039,
        0.784314,
        0.5625,
        0.968627,
        0.835294,
        0.698039,
        0.59375,
        0.94902,
        0.733333,
        0.588235,
        0.625,
        0.929412,
        0.65098,
        0.509804,
        0.65625,
        0.909804,
        0.564706,
        0.435294,
        0.6875,
        0.878431,
        0.458824,
        0.352941,
        0.71875,
        0.839216,
        0.388235,
        0.286275,
        0.75,
        0.760784,
        0.294118,
        0.211765,
        0.78125,
        0.701961,
        0.211765,
        0.168627,
        0.8125,
        0.65098,
        0.156863,
        0.129412,
        0.84375,
        0.6,
        0.094118,
        0.094118,
        0.875,
        0.54902,
        0.066667,
        0.098039,
        0.90625,
        0.501961,
        0.05098,
        0.12549,
        0.9375,
        0.45098,
        0.054902,
        0.172549,
        0.96875,
        0.4,
        0.054902,
        0.192157,
        1.0,
        0.34902,
        0.070588,
        0.211765,
    ],
    NanColor=[0.25, 0.0, 0.0],
    ScalarRangeInitialized=1.0,
)

# get opacity transfer function/opacity map for 'temperature'
temperaturePWF = GetOpacityTransferFunction('temperature')
temperaturePWF.ScalarRangeInitialized = 1

# get color transfer function/color map for 'brinkman_indicator'
brinkman_indicatorLUT = GetColorTransferFunction('brinkman_indicator')
brinkman_indicatorLUT.Set(
    RGBPoints=GenerateRGBPoints(
        preset_name='X Ray',
        range_min=0.5,
        range_max=0.993492066860199,
    ),
    ColorSpace='RGB',
    NanColor=[1.0, 0.0, 0.0],
    ScalarRangeInitialized=1.0,
)

# get opacity transfer function/opacity map for 'brinkman_indicator'
brinkman_indicatorPWF = GetOpacityTransferFunction('brinkman_indicator')
brinkman_indicatorPWF.Set(
    Points=[0.5, 0.0, 0.5, 0.0, 0.993492066860199, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# get color transfer function/color map for 'Velocity'
velocityLUT = GetColorTransferFunction('Velocity')
velocityLUT.Set(
    RGBPoints=GenerateRGBPoints(
        preset_name='Linear Green (Gr4L)',
        range_min=0.0,
        range_max=3.5,
    ),
    NanColor=[0.25, 0.0, 0.0],
    ScalarRangeInitialized=1.0,
)

# get opacity transfer function/opacity map for 'Velocity'
velocityPWF = GetOpacityTransferFunction('Velocity')
velocityPWF.Set(
    Points=[0.0, 0.0, 0.5, 0.0, 3.5, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# ----------------------------------------------------------------
# setup the visualization in view 'velocityView'
# (only objects that remain visible are shown)
# ----------------------------------------------------------------

# show data from physics_block
mergeBlocks1Display_vel = Show(physics_block, velocityView,
                               'UnstructuredGridRepresentation')

# trace defaults for the display properties.
mergeBlocks1Display_vel.Set(
    Representation='Outline',
    AmbientColor=[0.0, 0.0, 0.0],
    ColorArrayName=[None, ''],
    DiffuseColor=[0.0, 0.0, 0.0],
    LineWidth=3.0,
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
mergeBlocks1Display_vel.ScaleTransferFunction.Points = [
    0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0
]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
mergeBlocks1Display_vel.OpacityTransferFunction.Points = [
    0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0
]

# show data from BrinkmanClip
BrinkmanClipDisplay_vel = Show(BrinkmanClip, velocityView,
                               'UnstructuredGridRepresentation')

# trace defaults for the display properties.
BrinkmanClipDisplay_vel.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', ''],
    ComputePointNormals=1,
    FeatureAngle=45.0,
    Scale=[0.99999, 0.99999, 0.99999],
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
BrinkmanClipDisplay_vel.ScaleTransferFunction.Points = [
    0.13141125440597534, 0.0, 0.5, 0.0, 0.9962556958198547, 1.0, 0.5, 0.0
]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
BrinkmanClipDisplay_vel.OpacityTransferFunction.Points = [
    0.13141125440597534, 0.0, 0.5, 0.0, 0.9962556958198547, 1.0, 0.5, 0.0
]

# init the 'Polar Axes Representation' selected for 'PolarAxes'
BrinkmanClipDisplay_vel.PolarAxes.Scale = [0.99999, 0.99999, 0.99999]

# show data from BrinkmanClip
BrinkmanClip_1Display_vel = Show(OutputPort(BrinkmanClip, 1), velocityView,
                                 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
BrinkmanClip_1Display_vel.Set(
    Representation='Surface',
    ColorArrayName=[None, ''],
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
)

# show data from PhysicsClip
clip1Display_vel = Show(PhysicsClip, velocityView,
                        'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip1Display_vel.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Velocity'],
    LookupTable=velocityLUT,
    NonlinearSubdivisionLevel=2,
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
clip1Display_vel.ScaleTransferFunction.Points = [
    0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0
]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
clip1Display_vel.OpacityTransferFunction.Points = [
    0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0
]

# show data from Physics4
slice1Display_vel = Show(Physics4, velocityView, 'GeometryRepresentation')

# trace defaults for the display properties.
slice1Display_vel.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Velocity'],
    LookupTable=velocityLUT,
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
slice1Display_vel.ScaleTransferFunction.Points = [
    0.0, 0.0, 0.5, 0.0, 3.5, 1.0, 0.5, 0.0
]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
slice1Display_vel.OpacityTransferFunction.Points = [
    0.0, 0.0, 0.5, 0.0, 3.5, 1.0, 0.5, 0.0
]

# show data from Physics0
slice2Display_vel = Show(Physics0, velocityView, 'GeometryRepresentation')

# trace defaults for the display properties.
slice2Display_vel.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Velocity'],
    LookupTable=velocityLUT,
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
)

# show data from annotateTime1
annotateTime1Display_vel = Show(annotateTime1, velocityView,
                                'TextSourceRepresentation')

# trace defaults for the display properties.
annotateTime1Display_vel.Set(
    WindowLocation='Any Location',
    Position=[0.8899484270242393, 0.0215],
    Bold=1,
    Shadow=1,
    FontSize=36,
    Color=text_color,
)

# ----------------------------------------------------------------
# setup color legend bars for each visible scalar field
# ----------------------------------------------------------------

# get color legend/bar for brinkman_indicatorLUT in view velocityView
brinkman_indicatorLUTColorBar_vel = GetScalarBar(brinkman_indicatorLUT,
                                                 velocityView)
brinkman_indicatorLUTColorBar_vel.Visibility = 0

# get color legend/bar for temperatureLUT in view velocityView
temperatureLUTColorBar_vel = GetScalarBar(temperatureLUT, velocityView)
temperatureLUTColorBar_vel.Visibility = 0

# get color legend/bar for velocityLUT in view velocityView
velocityLUTColorBar_vel = GetScalarBar(velocityLUT, velocityView)
velocityLUTColorBar_vel.Set(
    AutoOrient=0,
    Orientation='Horizontal',
    Title='Velocity',
    ComponentTitle='',
    TitleBold=1,
    TitleShadow=1,
    TitleFontSize=24,
    TitleColor=text_color,
    LabelBold=1,
    LabelShadow=1,
    LabelFontSize=18,
    LabelColor=text_color,
    ScalarBarThickness=20,
    DrawScalarBarOutline=1,
    ScalarBarOutlineColor=[0.0, 0.0, 0.0],
    ScalarBarOutlineThickness=2,
    AutomaticLabelFormat=0,
    LabelFormat='{:<#6.1f}',
    DrawTickMarks=0,
    DrawTickLabels=0,
    AddRangeLabels=0,
)

# set color bar visibility
velocityLUTColorBar_vel.Visibility = 1

# show color legend
clip1Display_vel.SetScalarBarVisibility(velocityView, True)
slice1Display_vel.SetScalarBarVisibility(velocityView, True)
slice2Display_vel.SetScalarBarVisibility(velocityView, True)

# ----------------------------------------------------------------
# setup the visualization in view 'temperatureView'
# (temperature slices only — mirrors velocityView camera)
# ----------------------------------------------------------------

# show domain outline
mergeBlocks1Display_temp = Show(physics_block, temperatureView,
                                'UnstructuredGridRepresentation')
mergeBlocks1Display_temp.Set(
    Representation='Outline',
    AmbientColor=[0.0, 0.0, 0.0],
    ColorArrayName=[None, ''],
    DiffuseColor=[0.0, 0.0, 0.0],
    LineWidth=3.0,
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
)

# show Brinkman structure
BrinkmanClipDisplay_temp = Show(BrinkmanClip, temperatureView,
                                'UnstructuredGridRepresentation')
BrinkmanClipDisplay_temp.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', ''],
    ComputePointNormals=1,
    FeatureAngle=45.0,
    Scale=[0.99999, 0.99999, 0.99999],
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
)
BrinkmanClipDisplay_temp.ScaleTransferFunction.Points = [
    0.13141125440597534, 0.0, 0.5, 0.0, 0.9962556958198547, 1.0, 0.5, 0.0
]
BrinkmanClipDisplay_temp.OpacityTransferFunction.Points = [
    0.13141125440597534, 0.0, 0.5, 0.0, 0.9962556958198547, 1.0, 0.5, 0.0
]
BrinkmanClipDisplay_temp.PolarAxes.Scale = [0.99999, 0.99999, 0.99999]

# show secondary output of Brinkman clip
BrinkmanClip_1Display_temp = Show(OutputPort(BrinkmanClip, 1), temperatureView,
                                  'UnstructuredGridRepresentation')
BrinkmanClip_1Display_temp.Set(
    Representation='Surface',
    ColorArrayName=[None, ''],
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
)

# show clip colored by temperature
clip1Display_temp = Show(PhysicsClip, temperatureView,
                         'UnstructuredGridRepresentation')
clip1Display_temp.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'temperature'],
    LookupTable=temperatureLUT,
    NonlinearSubdivisionLevel=2,
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
)
clip1Display_temp.ScaleTransferFunction.Points = [
    0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0
]
clip1Display_temp.OpacityTransferFunction.Points = [
    0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0
]

# show Physics4 colored by temperature
slice1Display_temp = Show(Physics4, temperatureView, 'GeometryRepresentation')
slice1Display_temp.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'temperature'],
    LookupTable=temperatureLUT,
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
)
slice1Display_temp.ScaleTransferFunction.Points = [
    -0.04006486013531685, 0.0, 0.5, 0.0, 0.0, 1.0, 0.5, 0.0
]
slice1Display_temp.OpacityTransferFunction.Points = [
    -0.04006486013531685, 0.0, 0.5, 0.0, 0.0, 1.0, 0.5, 0.0
]

# show Physics0 colored by temperature
slice2Display_temp = Show(Physics0, temperatureView, 'GeometryRepresentation')
slice2Display_temp.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'temperature'],
    LookupTable=temperatureLUT,
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
)

# show time annotation
annotateTime1Display_temp = Show(annotateTime1, temperatureView,
                                 'TextSourceRepresentation')
annotateTime1Display_temp.Set(
    WindowLocation='Any Location',
    Position=[0.8899484270242393, 0.0215],
    Bold=1,
    Shadow=1,
    FontSize=36,
    Color=text_color,
)

# temperature color bar for temperatureView
temperatureLUTColorBar_temp = GetScalarBar(temperatureLUT, temperatureView)
temperatureLUTColorBar_temp.Set(
    AutoOrient=0,
    Orientation='Horizontal',
    Title='Temperature',
    ComponentTitle='',
    TitleBold=1,
    TitleShadow=1,
    TitleFontSize=24,
    TitleColor=text_color,
    LabelBold=1,
    LabelShadow=1,
    LabelFontSize=18,
    LabelColor=text_color,
    ScalarBarThickness=20,
    DrawScalarBarOutline=1,
    ScalarBarOutlineColor=[0.0, 0.0, 0.0],
    ScalarBarOutlineThickness=2,
    AutomaticLabelFormat=0,
    LabelFormat='{:<#6.1f}',
    DrawTickMarks=0,
    DrawTickLabels=0,
    AddRangeLabels=0,
    DrawAnnotations=0,
)
temperatureLUTColorBar_temp.Visibility = 1

# show color legends
clip1Display_temp.SetScalarBarVisibility(temperatureView, True)
slice1Display_temp.SetScalarBarVisibility(temperatureView, True)
slice2Display_temp.SetScalarBarVisibility(temperatureView, True)

# ----------------------------------------------------------------
# setup animation scene, tracks and keyframes
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get time animation track
timeAnimationCue = GetTimeTrack()
timeKeeper = GetTimeKeeper()
timeKeeper.SuppressedTimeSources = brinkman_field
timesteps = list(timeKeeper.TimestepValues)
max_animation_time = timesteps[-1] if timesteps else 0.0
animationScene = GetAnimationScene()

# initialize the animation scene
animationScene.Set(
    ViewModules=[velocityView, temperatureView],
    Cues=timeAnimationCue,
    AnimationTime=max_animation_time,
    PlayMode='Snap To TimeSteps',
)

# ----------------------------------------------------------------
# restore active source
SetActiveSource(None)

# ----------------------------------------------------------------
# Save PNG frames — only visit every Nth timestep so no work is wasted on
# intermediate steps that would never be written.

# Use ParaView's own process module to obtain the MPI rank. This requires
# pvbatch to be launched with a matching MPI implementation (i.e. mpiexec.mpich
# for this ParaView build which is linked against MPICH).
from paraview import servermanager as _sm

mpi_rank = _sm.vtkProcessModule.GetProcessModule().GetPartitionId()

for step_idx, t in enumerate(timesteps[::stride]):
    # Compute the original (un-strided) timestep index used for file names so
    # that Velocity_000002.png corresponds to the second raw timestep when
    # stride=2, making it easy to correlate frames with simulation output.
    raw_idx = step_idx * stride
    velocity_output = f'{output_dir}/Velocity_{raw_idx:06d}.png'
    temperature_output = f'{output_dir}/Temperature_{raw_idx:06d}.png'

    if not overwrite and os.path.exists(velocity_output) and os.path.exists(
            temperature_output):
        if mpi_rank == 0:
            print(f'Skipping step {raw_idx} (output files already exist; '
                  f'use --overwrite to force regeneration)')
        continue
    elif mpi_rank == 0:
        print(f'Processing step {raw_idx} (t={t:.2f})...')

    animationScene.AnimationTime = t
    SaveScreenshot(
        velocity_output,
        velocityView,
        ImageResolution=velocityView.ViewSize,
        TransparentBackground=1,
        CompressionLevel=0,
    )
    SaveScreenshot(
        temperature_output,
        temperatureView,
        ImageResolution=temperatureView.ViewSize,
        TransparentBackground=1,
        CompressionLevel=0,
    )
