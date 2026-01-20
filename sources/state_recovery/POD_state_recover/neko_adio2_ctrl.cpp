#include <adios2.h>
#include <mpi.h>
#include <iostream>
#include <thread>
#include <chrono>

// Global ADIOS instance & rank/size come from your field streamer TU
extern adios2::ADIOS adios;
extern int rank, size;

static adios2::IO io_ctrl;
static adios2::Engine ctrl_writer;
static adios2::Engine ctrl_reader;
static adios2::ADIOS *adios_ctrl = nullptr;

static adios2::Variable<int>    v_mode, v_phase, v_step;
static adios2::Variable<double> v_time;

static adios2::Variable<int> v_mode_cmd, v_phase_cmd;

static bool ctrl_inited = false;

extern "C" void adios2_ctrl_initialize_(const int *comm_int)
{
    if (ctrl_inited) return;

    MPI_Comm comm = MPI_Comm_f2c(*comm_int);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (rank == 0)
    {
        if (!adios_ctrl)
        {
            adios_ctrl = new adios2::ADIOS(MPI_COMM_SELF);
        }

        io_ctrl = adios_ctrl->DeclareIO("ctrlIO");
        io_ctrl.SetEngine("SST");

        v_mode  = io_ctrl.DefineVariable<int>("mode");
        v_phase = io_ctrl.DefineVariable<int>("phase");
        v_step  = io_ctrl.DefineVariable<int>("step");
        v_time  = io_ctrl.DefineVariable<double>("time");

        ctrl_writer = io_ctrl.Open("neko_ctrl_state", adios2::Mode::Write);
        ctrl_reader = io_ctrl.Open("neko_ctrl_cmd",   adios2::Mode::Read);

        std::cout << "neko_ctrl: initialized SST control streams (rank 0)\n";
    }

    ctrl_inited = true;
}

extern "C" void adios2_ctrl_finalize_()
{
    if (!ctrl_inited) return;

    if (rank == 0)
    {
        std::cout << "neko_ctrl: closing control streams (rank 0)\n";
        ctrl_writer.Close();
        ctrl_reader.Close();
        delete adios_ctrl;
        adios_ctrl = nullptr;
    }

    ctrl_inited = false;
}

extern "C" void adios2_ctrl_put_state_(
    const int    *mode,
    const int    *phase,
    const int    *step,
    const double *time)
{
    if (!ctrl_inited) return;
    if (rank != 0) return;

    ctrl_writer.BeginStep();
    ctrl_writer.Put(v_mode,  *mode);
    ctrl_writer.Put(v_phase, *phase);
    ctrl_writer.Put(v_step,  *step);
    ctrl_writer.Put(v_time,  *time);
    ctrl_writer.EndStep();
}

extern "C" void adios2_ctrl_wait_cmd_(int *mode_cmd, int *phase_cmd)
{
    if (!ctrl_inited) return;
    if (rank != 0) return;

    while (true)
    {
        auto status = ctrl_reader.BeginStep();

        if (status == adios2::StepStatus::OK)
        {
            v_mode_cmd  = io_ctrl.InquireVariable<int>("mode_cmd");
            v_phase_cmd = io_ctrl.InquireVariable<int>("phase_cmd");

            int m = *mode_cmd;
            int p = *phase_cmd;

            if (v_mode_cmd)  ctrl_reader.Get(v_mode_cmd,  m, adios2::Mode::Sync);
            if (v_phase_cmd) ctrl_reader.Get(v_phase_cmd, p, adios2::Mode::Sync);

            ctrl_reader.EndStep();

            *mode_cmd  = m;
            *phase_cmd = p;
            break;
        }
        else if (status == adios2::StepStatus::NotReady)
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
            continue;
        }
        else
        {
            *mode_cmd  = 9; // MODE_STOP
            *phase_cmd = 0;
            break;
        }
    }
}
