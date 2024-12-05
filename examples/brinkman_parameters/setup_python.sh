# Clone git repository pyNekTools

ROOT_FOLDER=$(realpath $(dirname $0)/../..)
CURRENT_DIR=$(pwd)

# Activate the virtual environment
if [ ! -d $ROOT_FOLDER/.venv ]; then
    python3 -m venv $ROOT_FOLDER/.venv
fi
source $ROOT_FOLDER/.venv/bin/activate
python3 -m pip install h5py tqdm

# Clone pyNekTools
rm -fr $ROOT_FOLDER/external/pyNekTools
git clone --depth=1 git@github.com:ExtremeFLOW/pyNekTools.git $ROOT_FOLDER/external/pyNekTools

# Install pyNekTools
cd $ROOT_FOLDER/external/pyNekTools
python3 -m pip install -r requirements.txt
python3 -m pip install .
cd $CURRENT_DIR
