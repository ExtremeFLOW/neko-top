# Paraview for working with Neko-TOP {#paraview}
\tableofcontents

Paraview is one of the main systems we use to inspect and visualise our work.
However it can be tricky to get the integration between python and Paraview to
work.

This document will give a short set of notes on how to set it up if needed.

## TLDR: Quick copy-paste set of commands for Ubuntu

```shell
sudo apt-get install git cmake build-essential libgl1-mesa-dev libxt-dev libqt5x11extras5-dev libqt5help5 qttools5-dev qtxmlpatterns5-dev-tools libqt5svg5-dev python3-dev python3-numpy libopenmpi-dev libtbb-dev ninja-build qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools libboost-dev
cd ~/
git clone --recursive https://gitlab.kitware.com/paraview/paraview.git
cmake -S paraview/ -B paraview/build -DPARAVIEW_USE_PYTHON=ON -DPARAVIEW_USE_MPI=ON -DVTK_SMP_IMPLEMENTATION_TYPE=TBB -DCMAKE_BUILD_TYPE=Release -DPARAVIEW_ENABLE_VISITBRIDGE=ON
cmake --build paraview/build --parallel
sudo cmake --install paraview/build
```

## Paraview setup

We managed to make it work by compiling paraview from scratch. Here is a couple
of useful links:

- https://kitware.github.io/paraview-docs/latest/cxx/md__builds_gitlab-kitware-sciviz-ci_Documentation_dev_build.html
- https://www.paraview.org/Wiki/ParaView/Python_Scripting

The following dependencies were required:
```
sudo apt-get install git cmake build-essential libgl1-mesa-dev libxt-dev libqt5x11extras5-dev libqt5help5 qttools5-dev qtxmlpatterns5-dev-tools libqt5svg5-dev python3-dev python3-numpy libopenmpi-dev libtbb-dev ninja-build qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools libboost-dev
```

We have managed to get Paraview up and running as the rendering engine. by
conducting the following:
- Clone paraview and build manually.
- Remember to build with the following options:
  - `PARAVIEW_USE_PYTHON`=ON
  - `PARAVIEW_USE_MPI`=ON
  - `PARAVIEW_ENABLE_VISITBRIDGE`=ON
  - `VTK_SMP_IMPLEMENTATION_TYPE`=TBB

## Paraview in python

In order to use paraview with python, one must make sure the libraries and
packages can be located by python.

Setting the following environment variables should be the way to go:
- `PYTHONPATH`
- `LD_LIBRARY_PATH`

However we do recommend to complete the install of Paraview by using the 
`cmake --install` command. That should prevent the need for setting the
`LD_LIBRARY_PATH`.

The Jupyter notebook `visualisation.ipynb` utilizes paraview, so practical
examples can be found there.