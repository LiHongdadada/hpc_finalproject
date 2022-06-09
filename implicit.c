static char help[] = "Solves a 10000x10000 linear system.\n\n";

#include <petscksp.h>
#include <petscmath.h>
#include <petscsys.h>
#include <petscviewerhdf5.h>

#include <math.h>

#define PI 3.14159265
#define FILE "implicit.h5"
int main(int argc, char **args)
{
    PetscInitialize(&argc, &args, (char *)0, help);
    Mat G_B, G_A, B, A;
    Vec x, T, temp_vec, times, G_q, G_Q;
    PetscErrorCode ierr;
    KSP ksp;
    PC pc;
    MPI_Comm comm;
    PetscInt maxit, r, row;
    PetscReal h, dt, t = 0, val, temp_B[4][4] = {{0.1111111111, 0.0555555556, 0.0277777778, 0.0555555556}, {0.0555555556, 0.1111111111, 0.0555555556, 0.0277777778}, {0.0277777778, 0.0555555556, 0.1111111111, 0.0555555556}, {0.0555555556, 0.0277777778, 0.0555555556, 0.1111111111}};
    PetscReal temp_A1[4][4] = {{0.6666666667, -0.1666666667, -0.3333333333, -0.1666666667}, {-0.1666666667, 0.6666666667, -0.1666666667, -0.3333333333}, {-0.3333333333, -0.1666666667, 0.6666666667, -0.1666666667}, {-0.1666666667, -0.3333333333, -0.1666666667, 0.6666666667}};
    comm = MPI_COMM_WORLD;
    PetscOptionsGetReal(NULL, NULL, "-h", &h, NULL);
    PetscOptionsGetReal(NULL, NULL, "-dt", &dt, NULL);
    PetscOptionsGetInt(NULL, NULL, "-maxit", &maxit, NULL);
    PetscOptionsGetInt(NULL, NULL, "-restart", &r, NULL);
    PetscReal xi[4] = {-0.5773, 0.5773, 0.5773, -0.5773}, eta[4] = {-0.5773, -0.5773, 0.5773, 0.5773};
    PetscInt n = 1 / h;
    PetscInt num_of_nodes = (n + 1) * (n + 1), num_of_elements = n * n;
    PetscInt iter = 0;
    PetscScalar temp_val;
    PetscInt index = 0;
    PetscReal t_v, Q[4];
    PetscReal nodes[num_of_nodes][3];
    PetscScalar data[3];
    PetscViewer h5; /*创建输出*/

    PetscInt elements[num_of_elements][5];
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

    ierr = MatCreateAIJ(comm, PETSC_DECIDE, PETSC_DECIDE, 4, 4, 4, PETSC_NULL, 4, PETSC_NULL, &B);
    CHKERRQ(ierr);
    ierr = MatSetUp(B);
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(B);
    CHKERRQ(ierr);
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            val = temp_B[i][j] * (h * h / dt);
            ierr = MatSetValues(B, 1, &i, 1, &j, &val, INSERT_VALUES);
            CHKERRQ(ierr);
        }
    }
    ierr = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    /*	PetscPrintf(PETSC_COMM_WORLD,"B\n");
            MatView(B, PETSC_VIEWER_STDOUT_WORLD);
    */
    ierr = MatCreateAIJ(comm, PETSC_DECIDE, PETSC_DECIDE, 4, 4, 4, PETSC_NULL, 4, PETSC_NULL, &A);
    CHKERRQ(ierr);
    ierr = MatSetUp(A);
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);
    CHKERRQ(ierr);
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {

            ierr = MatGetValue(B, i, j, &temp_val);
            CHKERRQ(ierr);
            val = temp_val + temp_A1[i][j];
            ierr = MatSetValues(A, 1, &i, 1, &j, &val, INSERT_VALUES);
            CHKERRQ(ierr);
        }
    }
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    /*        PetscPrintf(PETSC_COMM_WORLD,"A\n");
            MatView(A, PETSC_VIEWER_STDOUT_WORLD);
    */
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
                ierr = MatGetValue(A, i, j, &temp_val);
                CHKERRQ(ierr);
                ierr = MatSetValues(G_A, 1, &ii, 1, &jj, &temp_val, ADD_VALUES);
                CHKERRQ(ierr);
            }
        }
    }
    ierr = MatAssemblyBegin(G_A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(G_A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    /*       PetscPrintf(PETSC_COMM_WORLD,"before G_A\n");

           MatView(G_A, PETSC_VIEWER_STDOUT_WORLD);
   */
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
                ierr = MatGetValue(B, i, j, &temp_val);
                CHKERRQ(ierr);
                ierr = MatSetValues(G_B, 1, &ii, 1, &jj, &temp_val, ADD_VALUES);
                CHKERRQ(ierr);
            }
        }
    }
    ierr = MatAssemblyBegin(G_B, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(G_B, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    /*      PetscPrintf(PETSC_COMM_WORLD,"before G_B\n");

          MatView(G_B, PETSC_VIEWER_STDOUT_WORLD);
  */
    /******** 置一划零法 **********/
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
    /*     PetscPrintf(PETSC_COMM_WORLD,"AFTER G_B\n");

         MatView(G_B, PETSC_VIEWER_STDOUT_WORLD);
 */
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
            ierr = VecSetValues(G_Q, 1, &row, &val, ADD_VALUES);
            CHKERRQ(ierr);
        }
    }
    ierr = VecAssemblyBegin(G_Q);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(G_Q);
    CHKERRQ(ierr);
    ierr = VecAXPY(G_Q, -1, G_q);
    CHKERRQ(ierr);

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
        ierr = PetscPrintf(PETSC_COMM_WORLD, "t:%f\n", t);
        CHKERRQ(ierr);
    }
    else
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
        for (int i = 2 * n + 1; i < num_of_nodes - 1; i += n + 1)
        {
            val = 1;
            ierr = VecSetValues(T, 1, &i, &val, INSERT_VALUES);
            CHKERRQ(ierr);
        }
        ierr = VecAssemblyBegin(T);
        CHKERRQ(ierr);
        ierr = VecAssemblyEnd(T);
        CHKERRQ(ierr);
        t = 0;
    }
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
    ierr = KSPSetOperators(ksp, G_A, G_A);
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

        t += dt;
        ierr = MatMult(G_B, T, temp_vec);
        CHKERRQ(ierr);
        ierr = VecAXPY(temp_vec, 1, G_Q);
        CHKERRQ(ierr);

        ierr = KSPSolve(ksp, temp_vec, x);
        CHKERRQ(ierr); /* 求解出下一时间 */
        ierr = VecCopy(x, T);
        CHKERRQ(ierr); /* 求解后的结果传到T作为初始值 */
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

        ierr = VecAssemblyBegin(T);
        CHKERRQ(ierr);
        ierr = VecAssemblyEnd(T);
        CHKERRQ(ierr);
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
    ierr = PetscPrintf(PETSC_COMM_WORLD, "solution\n");
    CHKERRQ(ierr);

    ierr = VecView(T, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ(ierr);

    ierr = MatDestroy(&B);
    CHKERRQ(ierr);
    ierr = MatDestroy(&G_A);
    CHKERRQ(ierr);
    ierr = MatDestroy(&G_B);
    CHKERRQ(ierr);
    ierr = MatDestroy(&A);
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
