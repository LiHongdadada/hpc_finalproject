static char help[] = "Solves a 10000x10000 linear system";

#include <petscksp.h>
#include <petscmath.h>
#include <stdio.h>
#include <math.h>
#define PI 3.14159265

int main(int argc, char **args)
{
    PetscInitialize(&argc, &args, (char *)0, help);
    Mat G_B, G_A, B, A;
    Vec G_Q;
    MPI_Comm comm;
    PetscReal h, dt, val, temp_B[][4] = {{0.1111, 0.0556, 0.0278, 0.0556}, {0.0556, 0.1111, 0.0556, 0.0278}, {0.0278, 0.0556, 0.1111, 0.0556}, {0.0556, 0.0278, 0.0556, 0.1111}};
    PetscReal temp_A1[4][4] = {{0.6667, -0.1667, -0.3333, -0.1667}, {-0.1667, 0.6667, -0.1667, -0.3333}, {-0.3333, -0.1667, 0.6667, -0.1667}, {-0.1667, -0.3333, -0.1667, 0.6667}};
    comm = MPI_COMM_WORLD;
    PetscReal xi[4] = {-0.5773, 0.5773, 0.5773, -0.5773}, eta[4] = {-0.5773, -0.5773, 0.5773, 0.5773};
    PetscOptionsGetReal(NULL, NULL, "-h", &h, NULL);
    PetscOptionsGetReal(NULL, NULL, "-dt", &dt, NULL);
    PetscInt n = 1 / h;
    PetscInt num_of_nodes = (n + 1) * (n + 1), num_of_elements = n * n;
    PetscScalar temp_val;
    PetscReal nodes[num_of_nodes][3], Q[4];
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
        int a = i / (n + 1);
        elements[i][0] = i;
        elements[i][1] = i + a;
        elements[i][2] = i + a + 1;
        elements[i][3] = i + a + n + 1;
        elements[i][4] = i + a + n + 2;
    }
    PetscPrintf(PETSC_COMM_WORLD, "nodes:%f\n", nodes[3][1]);
    PetscPrintf(PETSC_COMM_WORLD, "e:%d\n", elements[0][3]);

    MatCreateAIJ(comm, PETSC_DECIDE, PETSC_DECIDE, 4, 4, 4, PETSC_NULL, 4, PETSC_NULL, &B);
    MatSetUp(B);
    MatSetFromOptions(B);
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            val = temp_B[i][j] * (h * h / dt);
            MatSetValues(B, 1, &i, 1, &j, &val, INSERT_VALUES);
        }
    }
    MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
    MatView(B, PETSC_VIEWER_STDOUT_WORLD);

    MatCreateAIJ(comm, PETSC_DECIDE, PETSC_DECIDE, 4, 4, 4, PETSC_NULL, 4, PETSC_NULL, &A);
    MatSetUp(A);
    MatSetFromOptions(A);
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {

            MatGetValue(B, i, j, &temp_val);
            val = temp_val + temp_A1[i][j];
            MatSetValues(A, 1, &i, 1, &j, &val, INSERT_VALUES);
        }
    }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    MatView(A, PETSC_VIEWER_STDOUT_WORLD);

    MatCreateAIJ(comm, PETSC_DECIDE, PETSC_DECIDE, num_of_nodes, num_of_nodes, num_of_nodes, PETSC_NULL, num_of_nodes, PETSC_NULL, &G_A);
    MatSetUp(G_A);
    MatSetFromOptions(G_A);
    for (int n = 0; n < num_of_elements; n++)
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {

                PetscInt ii = elements[n][i + 1], jj = elements[n][j + 1];
                MatGetValue(A, i, j, &temp_val);
                MatSetValues(G_A, 1, &ii, 1, &jj, &temp_val, ADD_VALUES);
            }
        }
    }
    MatAssemblyBegin(G_A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(G_A, MAT_FINAL_ASSEMBLY);
    MatView(G_A, PETSC_VIEWER_STDOUT_WORLD);

    MatCreateAIJ(comm, PETSC_DECIDE, PETSC_DECIDE, num_of_nodes, num_of_nodes, num_of_nodes, PETSC_NULL, num_of_nodes, PETSC_NULL, &G_B);
    MatSetUp(G_B);
    MatSetFromOptions(G_B);
    for (int n = 0; n < num_of_elements; n++)
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {

                PetscInt ii = elements[n][i + 1], jj = elements[n][j + 1];
                MatGetValue(B, i, j, &temp_val);
                MatSetValues(G_B, 1, &ii, 1, &jj, &temp_val, ADD_VALUES);
            }
        }
    }
    MatAssemblyBegin(G_B, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(G_B, MAT_FINAL_ASSEMBLY);
    MatView(G_B, PETSC_VIEWER_STDOUT_WORLD);

    VecCreate(PETSC_COMM_WORLD, &G_Q);
    VecSetSizes(G_Q, PETSC_DECIDE, num_of_nodes);
    VecSetFromOptions(G_Q);
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
        for (int i = 0; i < num_of_nodes; i++)
        {

            val = Q[i];
            VecSetValues(G_Q, 1, &i, &val, ADD_VALUES);
        }
    }
    VecAssemblyBegin(G_Q);
    VecAssemblyEnd(G_Q);
    VecView(G_Q, PETSC_VIEWER_STDOUT_WORLD);

    MatDestroy(&B);
    MatDestroy(&G_A);
    MatDestroy(&G_B);
    MatDestroy(&A);
    VecDestroy(&G_Q);
    PetscFinalize();
    return 0;
}

void assemble_G_q(double *G_q, int num_of_nodes, int n, double h);
void huayifa_G_Qaq(double *G_Q, double *G_q, int n, int num_of_nodes);
void huayifa_G_A(double **G_A, int num_of_elements, int num_of_nodes, int n, double h);
void huayifa_G_B(double **G_B, int num_of_elements, int num_of_nodes, int n, double h);

void assemble_G_q(double *G_q, int num_of_nodes, int n, double h)
{
    // double G_q[num_of_nodes]={0};
    int h1 = 1; // heat flux
    int i = 0;
    for (int i = 0; i < num_of_nodes; i++)
    {
        G_q[i] = 0;
    }

    for (i = n * (n + 1); i < num_of_nodes; ++i)
    {
        if (i == n * (n + 1))
        {
            G_q[i] += h * h1 / 2;
        }
        else if (i != n * (n + 1) && i != num_of_nodes - 1)
        {
            G_q[i] += h * h1;
        }
        else if (i == num_of_nodes - 1)
        {
            G_q[i] += h * h1 / 2;
        }
    }
    int j = 0;
    for (j = n; j < num_of_nodes; j += n + 1)
    {
        if (j == n)
        {
            G_q[j] += h * h1 / 2;
        }
        else if (j != n && j != num_of_nodes - 1)
        {
            G_q[j] += h * h1;
        }
        else if (j == num_of_nodes - 1)
        {
            G_q[j] += h * h1 / 2;
        }
    }
}

void huayifa_G_Qaq(double *G_Q, double *G_q, int n, int num_of_nodes)
{
    int i = 0;
    int j = 0;
    for (i = 0; i < num_of_nodes; i++)
    {
        G_Q[i] = G_Q[i] - G_q[i];
    }
    for (i = 0; i < n; i++)
    {
        G_Q[i] = 0;
    }
    for (j = 0; j < n * (n + 1) + 1; j += n + 1)
    {
        G_Q[j] = 0;
    }
}

void huayifa_G_A(double **G_A, int num_of_elements, int num_of_nodes, int n, double h)
{
    int i = 0;
    int j = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < num_of_nodes; j++)
        {
            if (j == i)
            {
                G_A[i][j] = 1;
            }
            else
            {
                G_A[i][j] = 0;
            }
        }
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < num_of_nodes; j++)
        {
            if (j == i)
            {
                G_A[j][i] = 1;
            }
            else
            {
                G_A[j][i] = 0;
            }
        }
    }
    for (i = 0; i < n * (n + 1) + 1; i += n + 1)
    {
        for (j = 0; j < num_of_nodes; j++)
        {
            if (j == i)
            {
                G_A[i][j] = 1;
            }
            else
            {
                G_A[i][j] = 0;
            }
        }
    }
    for (i = 0; i < n * (n + 1) + 1; i += n + 1)
    {
        for (j = 0; j < num_of_nodes; j++)
        {
            if (j == i)
            {
                G_A[j][i] = 1;
            }
            else
            {
                G_A[j][i] = 0;
            }
        }
    }
}

void huayifa_G_B(double **G_B, int num_of_elements, int num_of_nodes, int n, double h)
{
    int i = 0;
    int j = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < num_of_nodes; j++)
        {
            if (j == i)
            {
                G_B[i][j] = 1;
            }
            else
            {
                G_B[i][j] = 0;
            }
        }
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < num_of_nodes; j++)
        {
            if (j == i)
            {
                G_B[j][i] = 1;
            }
            else
            {
                G_B[j][i] = 0;
            }
        }
    }
    for (i = 0; i < n * (n + 1) + 1; i += n + 1)
    {
        for (j = 0; j < num_of_nodes; j++)
        {
            if (j == i)
            {
                G_B[i][j] = 1;
            }
            else
            {
                G_B[i][j] = 0;
            }
        }
    }
    for (i = 0; i < n * (n + 1) + 1; i += n + 1)
    {
        for (j = 0; j < num_of_nodes; j++)
        {
            if (j == i)
            {
                G_B[j][i] = 1;
            }
            else
            {
                G_B[j][i] = 0;
            }
        }
    }
}
