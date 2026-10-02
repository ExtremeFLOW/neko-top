// Minimal PETSc stand-in: just the Vec API used by topopt_in_petsc MMA.cc and
// the cross-check driver beam_ref.cc. MPI is real (mpicxx); Vec is a plain
// local array with PETSC_DECIDE layout as in PETSc. Arithmetic follows PETSc's
// VecAXPBYPCZ_Seq / VecAXPY_Seq / VecWAXPY_Seq special-casing.
#ifndef PETSC_SHIM_H
#define PETSC_SHIM_H
#include <mpi.h>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>

typedef int    PetscErrorCode;
typedef int    PetscInt;
typedef double PetscScalar;
typedef double PetscReal;
typedef enum { PETSC_FALSE = 0, PETSC_TRUE = 1 } PetscBool;
typedef enum { NORM_1 = 0, NORM_2 = 1, NORM_FROBENIUS = 2, NORM_INFINITY = 3 } NormType;

#define PETSC_COMM_WORLD MPI_COMM_WORLD
#define MPIU_SCALAR MPI_DOUBLE
#define PETSC_DECIDE (-1)
#define PetscMax(a, b) (((a) < (b)) ? (b) : (a))
#define PetscAbsReal(a) (((a) < 0) ? -(a) : (a))

struct _p_Vec {
    PetscInt     n, N;
    PetscScalar* a;
};
typedef _p_Vec* Vec;

inline PetscErrorCode PetscInitialize(int* argc, char*** argv, const char*, const char*) {
    MPI_Init(argc, argv);
    return 0;
}
inline PetscErrorCode PetscFinalize() {
    MPI_Finalize();
    return 0;
}
inline PetscErrorCode PetscPrintf(MPI_Comm comm, const char* fmt, ...) {
    int r;
    MPI_Comm_rank(comm, &r);
    if (r == 0) {
        va_list ap;
        va_start(ap, fmt);
        vprintf(fmt, ap);
        va_end(ap);
    }
    return 0;
}
inline PetscErrorCode PetscErrorPrintf(const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    return 0;
}

inline Vec shim_vec_new(PetscInt n, PetscInt N) {
    Vec v = new _p_Vec;
    v->n  = n;
    v->N  = N;
    v->a  = (PetscScalar*)calloc(n > 0 ? n : 1, sizeof(PetscScalar));
    return v;
}
inline PetscErrorCode VecCreateMPI(MPI_Comm comm, PetscInt n, PetscInt N, Vec* v) {
    int r, s;
    MPI_Comm_rank(comm, &r);
    MPI_Comm_size(comm, &s);
    if (n == PETSC_DECIDE) n = N / s + ((N % s) > r ? 1 : 0);
    *v = shim_vec_new(n, N);
    return 0;
}
inline PetscErrorCode VecDuplicate(Vec x, Vec* y) {
    *y = shim_vec_new(x->n, x->N);
    return 0;
}
inline PetscErrorCode VecDuplicateVecs(Vec x, PetscInt m, Vec** y) {
    *y = new Vec[m > 0 ? m : 1];
    for (PetscInt i = 0; i < m; i++) (*y)[i] = shim_vec_new(x->n, x->N);
    return 0;
}
inline PetscErrorCode VecDestroy(Vec* v) {
    if (*v) {
        free((*v)->a);
        delete *v;
        *v = nullptr;
    }
    return 0;
}
inline PetscErrorCode VecDestroyVecs(PetscInt m, Vec** v) {
    for (PetscInt i = 0; i < m; i++) VecDestroy(&(*v)[i]);
    delete[] *v;
    *v = nullptr;
    return 0;
}
inline PetscErrorCode VecGetLocalSize(Vec x, PetscInt* n) {
    *n = x->n;
    return 0;
}
inline PetscErrorCode VecGetArray(Vec x, PetscScalar** a) {
    *a = x->a;
    return 0;
}
inline PetscErrorCode VecRestoreArray(Vec, PetscScalar** a) {
    *a = nullptr;
    return 0;
}
inline PetscErrorCode VecGetArrayRead(Vec x, const PetscScalar** a) {
    *a = x->a;
    return 0;
}
inline PetscErrorCode VecRestoreArrayRead(Vec, const PetscScalar** a) {
    *a = nullptr;
    return 0;
}
inline PetscErrorCode VecGetArrays(const Vec* x, PetscInt m, PetscScalar*** a) {
    *a = new PetscScalar*[m > 0 ? m : 1];
    for (PetscInt i = 0; i < m; i++) (*a)[i] = x[i]->a;
    return 0;
}
inline PetscErrorCode VecRestoreArrays(const Vec*, PetscInt, PetscScalar*** a) {
    delete[] *a;
    *a = nullptr;
    return 0;
}
inline PetscErrorCode VecSet(Vec x, PetscScalar s) {
    for (PetscInt i = 0; i < x->n; i++) x->a[i] = s;
    return 0;
}
// y <- x
inline PetscErrorCode VecCopy(Vec x, Vec y) {
    if (x != y) memcpy(y->a, x->a, x->n * sizeof(PetscScalar));
    return 0;
}
// y <- y + alpha x
inline PetscErrorCode VecAXPY(Vec y, PetscScalar alpha, Vec x) {
    for (PetscInt i = 0; i < y->n; i++) y->a[i] += alpha * x->a[i];
    return 0;
}
// w <- alpha x + y
inline PetscErrorCode VecWAXPY(Vec w, PetscScalar alpha, Vec x, Vec y) {
    for (PetscInt i = 0; i < w->n; i++) w->a[i] = alpha * x->a[i] + y->a[i];
    return 0;
}
// z <- alpha x + beta y + gamma z   (VecAXPBYPCZ_Seq special cases)
inline PetscErrorCode VecAXPBYPCZ(Vec z, PetscScalar alpha, PetscScalar beta, PetscScalar gamma, Vec x, Vec y) {
    PetscScalar *zz = z->a, *xx = x->a, *yy = y->a;
    PetscInt     n = z->n;
    if (alpha == 1.0) {
        for (PetscInt i = 0; i < n; i++) zz[i] = xx[i] + beta * yy[i] + gamma * zz[i];
    } else if (gamma == 1.0) {
        for (PetscInt i = 0; i < n; i++) zz[i] = alpha * xx[i] + beta * yy[i] + zz[i];
    } else if (gamma == 0.0) {
        for (PetscInt i = 0; i < n; i++) zz[i] = alpha * xx[i] + beta * yy[i];
    } else {
        for (PetscInt i = 0; i < n; i++) zz[i] = alpha * xx[i] + beta * yy[i] + gamma * zz[i];
    }
    return 0;
}
inline PetscErrorCode VecNorm(Vec x, NormType t, PetscReal* r) {
    if (t != NORM_INFINITY) {
        PetscErrorPrintf("petsc shim: only NORM_INFINITY implemented\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    PetscReal loc = 0.0;
    for (PetscInt i = 0; i < x->n; i++) loc = PetscMax(loc, PetscAbsReal(x->a[i]));
    MPI_Allreduce(&loc, r, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return 0;
}
#endif
