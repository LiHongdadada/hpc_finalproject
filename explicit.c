static char help[] = "Solves a 10000x10000 linear system.\n\n";

#include <petscksp.h>
#include <petscmath.h>
#include <petscsys.h>
#include <petscviewerhdf5.h>
#include <math.h>

#define PI 3.14159265
#define FILE "explicit.h5"
int main(int argc, char **args)
{
    PetscInitialize(&argc, &args, (char *)0, help);
    Mat G_B, G_A;
    Vec x, T, temp_vec, times, G_q, G_Q;
    PetscErrorCode ierr;
    KSP ksp;
    PC pc;
    MPI_Comm comm;
    PetscInt maxit, r, row, rank;
    PetscReal h, dt, t = 0, val, temp_B[4][4] = {{0.1111111111, 0.0555555556, 0.0277777778, 0.0555555556}, {0.0555555556, 0.1111111111, 0.0555555556, 0.0277777778}, {0.0277777778, 0.0555555556, 0.1111111111, 0.0555555556}, {0.0555555556, 0.0277777778, 0.0555555556, 0.1111111111}};
    PetscReal temp_A1[4][4] = {{0.6666666667, -0.1666666667, -0.3333333333, -0.1666666667}, {-0.1666666667, 0.6666666667, -0.1666666667, -0.3333333333}, {-0.3333333333, -0.1666666667, 0.6666666667, -0.1666666667}, {-0.1666666667, -0.3333333333, -0.1666666667, 0.6666666667}};
    comm = MPI_COMM_WORLD;
    /*get parameters from command line */
    PetscOptionsGetReal(NULL, NULL, "-h", &h, NULL);
    PetscOptionsGetReal(NULL, NULL, "-dt", &dt, NULL);
    PetscOptionsGetInt(NULL, NULL, "-maxit", &maxit, NULL);
    PetscOptionsGetInt(NULL, NULL, "-restart", &r, NULL);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    CHKERRQ(ierr);
    PetscReal xi[4] = {-0.5773, 0.5773, 0.5773, -0.5773}, eta[4] = {-0.5773, -0.5773, 0.5773, 0.5773};
    PetscInt n = 1 / h;
    PetscInt num_of_nodes = (n + 1) * (n + 1), num_of_elements = n * n;
    PetscInt iter = 0;
    PetscScalar temp_val, norm = 1.0, tempNorm = 0.0, tempV = 1.0;
    PetscInt index = 0;
    PetscReal t_v, Q[4];
    PetscReal nodes[num_of_nodes][3],B[4][4],A[4][4];
    PetscScalar data[3];
    PetscViewer h5;
    PetscInt elements[num_of_elements][5];

    /*create nodes and elements */
    for (int i = 0; i < num_of_nodes; i++)
    {
        int a = i / (n + 1);
        nodes[i][0] = i;
        nodes[i][1] = (i - a * (n + 1)) * h;
        nodes[i][2] = a * h;
    }
    for (int i = 0; i < num_of_elements; i++)
    {
        int a = i / n;
        elements[i][0] = i;
        elements[i][1] = i + a;
        elements[i][2] = i + a + 1;
        elements[i][3] = i + a + n + 2;
        elements[i][4] = i + a + n + 1;
    }

    /* create each matrix. B is the matrix before time. */
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            B[i][j]= temp_B[i][j] * (h * h / dt);
        }
    }


    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            A[i][j] = B[i][j] - temp_A1[i][j];
        }
    }
    /* G_* is the global matrix */
    ierr = MatCreateAIJ(comm, PETSC_DECIDE, PETSC_DECIDE, num_of_nodes, num_of_nodes, num_of_nodes, PETSC_NULL, num_of_nodes, PETSC_NULL, &G_A);
    CHKERRQ(ierr);
    ierr = MatSetUp(G_A);
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(G_A);
    CHKERRQ(ierr);

    for (int n = 0; n < num_of_elements; n++)
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                PetscInt ii = elements[n][i + 1], jj = elements[n][j + 1];
                temp_val=A[i][j];
                ierr = MatSetValues(G_A, 1, &ii, 1, &jj, &temp_val, ADD_VALUES);
                CHKERRQ(ierr);
            }
        }
    }
    ierr = MatAssemblyBegin(G_A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(G_A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    ierr = MatCreateAIJ(comm, PETSC_DECIDE, PETSC_DECIDE, num_of_nodes, num_of_nodes, num_of_nodes, PETSC_NULL, num_of_nodes, PETSC_NULL, &G_B);
    CHKERRQ(ierr);
    ierr = MatSetUp(G_B);
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(G_B);
    CHKERRQ(ierr);
    for (int n = 0; n < num_of_elements; n++)
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                PetscInt ii = elements[n][i + 1], jj = elements[n][j + 1];
                temp_val=B[i][j];
                ierr = MatSetValues(G_B, 1, &ii, 1, &jj, &temp_val, ADD_VALUES);
                CHKERRQ(ierr);
            }
        }
    }
    ierr = MatAssemblyBegin(G_B, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(G_B, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    /*add boundary condition*/
    for (int i = 0; i < n + 1; i++)
    {
        ierr = MatZeroRows(G_B, 1, &i, 1, 0, 0);
        CHKERRQ(ierr);
        ierr = MatZeroRows(G_A, 1, &i, 1, 0, 0);
        CHKERRQ(ierr);
    }
    for (int i = n + 1; i < n * (n + 1) + 1; i += n + 1)
    {
        ierr = MatZeroRows(G_B, 1, &i, 1, 0, 0);
        CHKERRQ(ierr);
        ierr = MatZeroRows(G_A, 1, &i, 1, 0, 0);
        CHKERRQ(ierr);
    }
    for (int i = n * (n + 1) + 1; i < num_of_nodes; i++)
    {
        ierr = MatZeroRows(G_B, 1, &i, 1, 0, 0);
        CHKERRQ(ierr);
        ierr = MatZeroRows(G_A, 1, &i, 1, 0, 0);
        CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(G_B, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(G_B, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyBegin(G_A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(G_A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    /* G_q is the heat flux */
    ierr = VecCreate(PETSC_COMM_WORLD, &G_q);
    CHKERRQ(ierr);
    ierr = VecSetSizes(G_q, PETSC_DECIDE, num_of_nodes);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(G_q);
    CHKERRQ(ierr);
    ierr = VecSet(G_q, 0);
    CHKERRQ(ierr);

    for (int i = 2 * n + 1; i < num_of_nodes - 1; i += n + 1)
    {
        val = h;
        ierr = VecSetValues(G_q, 1, &i, &val, INSERT_VALUES);
        CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(G_q);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(G_q);
    CHKERRQ(ierr);
    /*G_Q is the heat source */
    ierr = VecCreate(PETSC_COMM_WORLD, &G_Q);
    CHKERRQ(ierr);
    ierr = VecSetSizes(G_Q, PETSC_DECIDE, num_of_nodes);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(G_Q);
    CHKERRQ(ierr);

    for (int k = 0; k < num_of_elements; k++)
    {
        for (int tt = 0; tt < 4; tt++)
        {
            Q[tt] = 0.0;
        }
        for (int jj = 0; jj < 4; jj++)
        {
            Q[0] += (h * h / 32.0) * (sin(PI * (h / 2.0 * xi[jj] + nodes[elements[k][1]][1] + h / 2.0)) + sin(PI * (h / 2.0 * eta[jj] + nodes[elements[k][1]][2] + h / 2.0))) * (1 - xi[jj]) * (1 - eta[jj]);
            Q[1] += (h * h / 32.0) * (sin(PI * (h / 2.0 * xi[jj] + nodes[elements[k][1]][1] + h / 2.0)) + sin(PI * (h / 2.0 * eta[jj] + nodes[elements[k][1]][2] + h / 2.0))) * (1 + xi[jj]) * (1 - eta[jj]);
            Q[2] += (h * h / 32.0) * (sin(PI * (h / 2.0 * xi[jj] + nodes[elements[k][1]][1] + h / 2.0)) + sin(PI * (h / 2.0 * eta[jj] + nodes[elements[k][1]][2] + h / 2.0))) * (1 + xi[jj]) * (1 + eta[jj]);
            Q[3] += (h * h / 32.0) * (sin(PI * (h / 2.0 * xi[jj] + nodes[elements[k][1]][1] + h / 2.0)) + sin(PI * (h / 2.0 * eta[jj] + nodes[elements[k][1]][2] + h / 2.0))) * (1 - xi[jj]) * (1 + eta[jj]);
        }
        for (int i = 0; i < 4; i++)
        {
            row = elements[k][i + 1];
            val = Q[i];
            VecSetValues(G_Q, 1, &row, &val, ADD_VALUES);
            CHKERRQ(ierr);
        }
    }
    ierr = VecAssemblyBegin(G_Q);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(G_Q);
    CHKERRQ(ierr);

    /*merge matrixes in the right of the equation */
    ierr = VecAXPY(G_Q, -1, G_q);
    CHKERRQ(ierr);
    /* create result temperature vector and record iteration times */
    ierr = VecCreate(PETSC_COMM_WORLD, &T);
    CHKERRQ(ierr);
    ierr = VecSetSizes(T, PETSC_DECIDE, num_of_nodes);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(T);
    CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &times);
    CHKERRQ(ierr);
    ierr = VecSetSizes(times, PETSC_DECIDE, 3);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(times);
    CHKERRQ(ierr);

    /* judge wether the program need restart */
    if (r == 1)
    {
        ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, FILE, FILE_MODE_READ, &h5);
        CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)T, "T");
        CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)times, "iteration_times");
        CHKERRQ(ierr);
        ierr = VecLoad(T, h5);
        CHKERRQ(ierr);
        ierr = VecLoad(times, h5);
        CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&h5);
        CHKERRQ(ierr);
        ierr = VecGetValues(times, 1, &index, &t);
        CHKERRQ(ierr);
    }
    else
    {
        /************温度初值************/
        if (rank == 0)
        {
            ierr = VecSet(T, 0);
            CHKERRQ(ierr);
            for (int i = 0; i < n + 1; i++)
            {
                val = 1;
                ierr = VecSetValues(T, 1, &i, &val, INSERT_VALUES);
                CHKERRQ(ierr);
            }
            for (int i = n + 1; i < n * (n + 1) + 1; i += n + 1)
            {
                val = 1;
                ierr = VecSetValues(T, 1, &i, &val, INSERT_VALUES);
                CHKERRQ(ierr);
            }
            for (int i = n * (n + 1) + 1; i < num_of_nodes; i++)
            {
                val = 1;
                ierr = VecSetValues(T, 1, &i, &val, INSERT_VALUES);
                CHKERRQ(ierr);
            }
        }

        ierr = VecAssemblyBegin(T);
        CHKERRQ(ierr);
        ierr = VecAssemblyEnd(T);
        CHKERRQ(ierr);
    }

    /* create some temporary vector */
    ierr = VecCreate(PETSC_COMM_WORLD, &temp_vec);
    CHKERRQ(ierr);
    ierr = VecSetSizes(temp_vec, PETSC_DECIDE, num_of_nodes);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(temp_vec);
    CHKERRQ(ierr);
    ierr = VecSet(temp_vec, 0);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(temp_vec);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(temp_vec);
    CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &x);
    CHKERRQ(ierr);
    ierr = VecSetSizes(x, PETSC_DECIDE, num_of_nodes);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);
    CHKERRQ(ierr);
    ierr = VecSet(x, 0);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(x);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x);
    CHKERRQ(ierr);

    /******** ksp *********/
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
    CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, G_B, G_B);
    CHKERRQ(ierr);
    ierr = KSPGetPC(ksp, &pc);
    CHKERRQ(ierr);
    ierr = PCSetType(pc, PCJACOBI);
    CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, 1.e-10, 1.e-50, PETSC_DEFAULT, PETSC_DEFAULT);
    CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);
    CHKERRQ(ierr);
    while (PetscAbsReal(t) < maxit * dt)
    {
        /*tempV = PetscAbsScalar(tempNorm - norm);
 *         tempNorm = norm;*/
        t += dt;
        ierr = MatMult(G_A, T, temp_vec);
        CHKERRQ(ierr);
        ierr = VecAXPY(temp_vec, 1, G_Q);
        CHKERRQ(ierr);
        ierr = KSPSolve(ksp, temp_vec, x);
        CHKERRQ(ierr); /* solve equation */
        ierr = VecCopy(x, T);
        CHKERRQ(ierr); /* transfer it to result vector */
        if (rank == 0)
        {
            for (int i = 0; i < n + 1; i++)
            {
                val = 1;
                ierr = VecSetValues(T, 1, &i, &val, INSERT_VALUES);
                CHKERRQ(ierr);
            }
            for (int i = n + 1; i < n * (n + 1) + 1; i += n + 1)
            {
                val = 1;
                ierr = VecSetValues(T, 1, &i, &val, INSERT_VALUES);
                CHKERRQ(ierr);
            }
            for (int i = n * (n + 1) + 1; i < num_of_nodes; i++)
            {
                val = 1;
                ierr = VecSetValues(T, 1, &i, &val, INSERT_VALUES);
                CHKERRQ(ierr);
            }
        }

        ierr = VecAssemblyBegin(T);
        CHKERRQ(ierr);
        ierr = VecAssemblyEnd(T);
        CHKERRQ(ierr);
        /*VecNorm(T, NORM_2, &norm);*/
        /***********每十步迭代记录一次数据***************/
        iter += 1;
        if ((iter % 10) == 0)
        {
            data[0] = t;
            data[1] = h;
            data[2] = dt;
            ierr = VecSet(times, 0);
            CHKERRQ(ierr);
            for (index = 0; index < 3; index++)
            {
                t_v = data[index];
                ierr = VecSetValues(times, 1, &index, &t_v, INSERT_VALUES);
                CHKERRQ(ierr);
            }

            ierr = VecAssemblyBegin(times);
            CHKERRQ(ierr);
            ierr = VecAssemblyEnd(times);
            CHKERRQ(ierr);
            ierr = PetscViewerCreate(PETSC_COMM_WORLD, &h5);
            CHKERRQ(ierr);
            ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, FILE, FILE_MODE_WRITE, &h5);
            CHKERRQ(ierr);
            ierr = PetscObjectSetName((PetscObject)T, "T");
            CHKERRQ(ierr);
            ierr = PetscObjectSetName((PetscObject)times, "iteration_times");
            CHKERRQ(ierr);
            ierr = VecView(times, h5);
            CHKERRQ(ierr);
            ierr = VecView(T, h5);
            CHKERRQ(ierr);
            ierr = PetscViewerDestroy(&h5);
            CHKERRQ(ierr);
        }
    }
    /*show results */
    ierr = PetscPrintf(PETSC_COMM_WORLD, "solution\n");
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "iter:%d\n", iter);
    CHKERRQ(ierr);
    ierr = VecView(T, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ(ierr);

    ierr = MatDestroy(&G_A);
    CHKERRQ(ierr);
    ierr = MatDestroy(&G_B);
    CHKERRQ(ierr);
    ierr = VecDestroy(&temp_vec);
    CHKERRQ(ierr);
    ierr = VecDestroy(&T);
    CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    ierr = VecDestroy(&times);
    CHKERRQ(ierr);
    ierr = VecDestroy(&G_q);
    CHKERRQ(ierr);
    ierr = PetscFinalize();
    CHKERRQ(ierr);
    return ierr;
}

