// Same 1D cantilever beam problem as tests/regression/mma in neko-top, driven by
// the UNMODIFIED MMA.cc from topopt_in_petsc (Aage).
#include <petsc.h>
#include <cstring>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include "MMA.h"

static const PetscInt n = 221184, m = 11, ncon = 10;
static const double L_total = 2.0, b = 0.02, rho = 7800.0, E = 210.0e9, P = 1000.0;
static const double h_min = 0.005, h_max = 0.05, u_tip_max = 0.25, sigma_max = 250e6;

// Port of fill_constraint_indices from neko-top tests/regression/mma/driver.f90 (1-based)
static void fill_idx(std::vector<int>& idx, int nc, int np, int ds) {
    int cpp = nc / np, rc = nc % np, bs = ds / np, rs = ds % np;
    int k = 0, cum = 1;
    for (int i = 1; i <= np; i++) {
        int pc = (i <= rc) ? cpp + 1 : cpp;
        int ps = (i <= rs) ? bs + 1 : bs;
        int st = cum, en = st + ps - 1;
        for (int j = 1; j <= pc; j++) {
            if (k >= nc) break;
            idx[k++] = std::min(st + (j - 1), en);
        }
        cum = en + 1;
    }
}

static void evaluate(const PetscScalar* x, double& f0, PetscScalar* df0, double* g, PetscScalar** dg,
                     const std::vector<double>& Delta, const std::vector<int>& idx, double Le) {
    std::vector<double> h(n);
    for (PetscInt k = 0; k < n; k++) h[k] = h_min + (h_max - h_min) * x[k];
    double s = 0.0;
    for (PetscInt k = 0; k < n; k++) s += h[k];
    f0 = rho * b * Le * s;
    for (PetscInt k = 0; k < n; k++) df0[k] = rho * b * Le * (h_max - h_min);
    s = 0.0;
    for (PetscInt k = 0; k < n; k++) s += Delta[k] / (b * h[k] * h[k] * h[k] / 12.0) * (P / E);
    g[0] = s / u_tip_max - 1.0;
    for (PetscInt i = 0; i < m; i++)
        for (PetscInt k = 0; k < n; k++) dg[i][k] = 0.0;
    for (PetscInt k = 0; k < n; k++)
        dg[0][k] = Delta[k] * (P * (-36.0) * (h_max - h_min) / (E * b)) / (h[k] * h[k] * h[k] * h[k]) / u_tip_max;
    for (int i = 0; i < ncon; i++) {
        int j = idx[i] - 1;
        double xe = Le * (double)(idx[i] - 1);
        double Ie = b * h[j] * h[j] * h[j] / 12.0, ce = h[j] / 2.0, Me = P * (L_total - xe);
        g[1 + i] = (Me * ce / Ie) / sigma_max - 1.0;
        dg[1 + i][j] = Me * ((1.0 / (2.0 * Ie)) - (ce * 3.0 * b * h[j] * h[j] / 12.0) / (Ie * Ie)) * (h_max - h_min) / sigma_max;
    }
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);
    std::string out = (argc > 1) ? argv[1] : "ref";
    int nit = (argc > 2) ? atoi(argv[2]) : 20;
    double movlim = (argc > 3) ? atof(argv[3]) : -1.0;

    std::vector<int> idx(ncon);
    fill_idx(idx, ncon, ncon, n);
    double Le = L_total / (double)n;
    std::vector<double> Delta(n);
    for (PetscInt k = 0; k < n; k++)
        Delta[k] = (std::pow(L_total - Le * (double)k, 3) - std::pow(L_total - Le * (double)(k + 1), 3)) / 3.0;

    Vec x, xold, xmin, xmax, xminL, xmaxL, dfdx, *dgdx;
    VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, n, &x);
    VecDuplicate(x, &xold); VecDuplicate(x, &xmin); VecDuplicate(x, &xmax); VecDuplicate(x, &dfdx);
    VecDuplicateVecs(x, m, &dgdx);
    VecSet(x, 0.5); VecSet(xmin, 0.0); VecSet(xmax, 1.0);
    VecDuplicate(x, &xminL); VecDuplicate(x, &xmaxL); VecCopy(xmin, xminL); VecCopy(xmax, xmaxL);
    PetscScalar a[m], c[m], d[m], gx[m];
    for (int i = 0; i < m; i++) { a[i] = 0.0; c[i] = 1000.0; d[i] = 1.0; }
    MMA* mma = new MMA(n, m, x, a, c, d);   // default asymptotes 0.5 / 0.7 / 1.2

    double f0;
    auto eval = [&]() {
        PetscScalar *xv, *dfv, **dgv;
        VecGetArray(x, &xv); VecGetArray(dfdx, &dfv); VecGetArrays(dgdx, m, &dgv);
        evaluate(xv, f0, dfv, gx, dgv, Delta, idx, Le);
        VecRestoreArray(x, &xv); VecRestoreArray(dfdx, &dfv); VecRestoreArrays(dgdx, m, &dgv);
    };
    eval();
    FILE* fc = fopen((out + ".csv").c_str(), "w");
    FILE* fx = fopen((out + ".x").c_str(), "wb");
    fprintf(fc, "iter,f0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,kktmax,kktnorm,maxdx\n");
    fprintf(fc, "0,%.16e", f0);
    for (int i = 0; i < m; i++) fprintf(fc, ",%.16e", gx[i]);
    fprintf(fc, ",0,0,0\n");
    for (int itr = 1; itr <= nit; itr++) {
        VecCopy(x, xold);
        if (movlim > 0.0) mma->SetOuterMovelimit(0.0, 1.0, movlim, x, xminL, xmaxL);  // as main.cc
        mma->Update(x, dfdx, gx, dgdx, xminL, xmaxL);
        eval();
        PetscScalar n2, ninf;
        mma->KKTresidual(x, dfdx, gx, dgdx, xmin, xmax, &n2, &ninf);
        Vec dx; VecDuplicate(x, &dx); VecWAXPY(dx, -1.0, xold, x);
        PetscReal mdx; VecNorm(dx, NORM_INFINITY, &mdx); VecDestroy(&dx);
        fprintf(fc, "%d,%.16e", itr, f0);
        for (int i = 0; i < m; i++) fprintf(fc, ",%.16e", gx[i]);
        fprintf(fc, ",%.16e,%.16e,%.16e\n", ninf, n2, mdx);
        fflush(fc);
        const PetscScalar* xv; VecGetArrayRead(x, &xv); fwrite(xv, sizeof(double), n, fx); VecRestoreArrayRead(x, &xv);
    }
    fclose(fc); fclose(fx);
    delete mma;
    PetscFinalize();
    return 0;
}
