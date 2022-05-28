static char help[] = "Solves a 10000x10000 linear system";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <petscksp.h>
#include <petscmath.h>

#define PI 3.14159265

void q_matrix(double q[4], double h);
double *Q_matrix(double h, double nodes[][3], int elements[][5], int i);
void B_matrix(double B[][4], double h);
void elements_matrix(int elements[][5], int num_of_elements, int n);
void nodes_matrix(double nodes[][3], int num_of_nodes, int n, double h);
void A_matrix(double A[][4], double B[][4], double A1[][4]);
void assemble_G_A(double **G_A, double A[][4], int num_of_elements, int num_of_nodes, int elments[][5]);
void assemble_G_B(double **G_B, double B[][4], int num_of_elements, int num_of_nodes, int elments[][5]);
void assembel_G_Q(double *G_Q, double h, double nodes[][3], int num_of_elements, int num_of_nodes, int elments[][5]);
void assemble_G_q(double *G_q, int num_of_nodes, int n, double h);
void huayifa_G_Qaq(double *G_Q, double *G_q, int n, int num_of_nodes);
void huayifa_G_A(double **G_A, int num_of_elements, int num_of_nodes, int n, double h);
void huayifa_G_B(double **G_B, int num_of_elements, int num_of_nodes, int n, double h);
double **mallocMatrix(int row, int col);
void freeMatrix(double **a);

int main(int argc, char **args)
{
    /******************** petsc *********************/
    Vec T, QQ, b, x, w, t;
    Mat AA, BB;
    KSP ksp; /*设置求解方法*/
    PC pc;   /*设置求解参数*/
    PetscErrorCode ierr;
    PetscInt n = 100, num_of_nodes = (n + 1) * (n + 1); /*这是将区域分成n*n块*/
    PetscInt i = 0, ii = 0, j = 0, col4[4], col6[6], col8[8], col9[9], rank;
    /*其中i,ii是矩阵和向量的角标，rank为程序并行化所需参数*/
    PetscReal h = 0.1, q0 = 0;
    PetscReal k = 50; /* time steps */
    PetscScalar value[9];
    /******************** petsc *********************/
    // double h=1;
    // int i = 0, j = 0;
    int length = 1;
    // int n = length / h;
    int num_of_elements = n * n;
    // int num_of_nodes = (n + 1) * (n + 1);
    double nodes[num_of_nodes][3];
    double q[4];
    // double T[num_of_nodes];
    // double Tdt[num_of_nodes];
    double B[4][4] = {{0.1111, 0.0556, 0.0278, 0.0556}, {0.0556, 0.1111, 0.0556, 0.0278}, {0.0278, 0.0556, 0.1111, 0.0556}, {0.0556, 0.0278, 0.0556, 0.1111}};
    int elements[num_of_elements][5];
    double A1[4][4] = {{0.6667, -0.1667, -0.3333, -0.1667}, {-0.1667, 0.6667, -0.1667, -0.3333}, {-0.3333, -0.1667, 0.6667, -0.1667}, {-0.1667, -0.3333, -0.1667, 0.6667}};
    double A[4][4] = {0};
    double *Q;
    // double **G_A = mallocMatrix(num_of_nodes, num_of_nodes);
    // double **G_B = mallocMatrix(num_of_nodes, num_of_nodes);
    // // double **G_Q = mallocMatrix(num_of_nodes, num_of_nodes);
    // // double **G_q = mallocMatrix(num_of_nodes, num_of_nodes);
    // double *G_Q = (double *)malloc(sizeof(double) * num_of_nodes);
    // double *G_q = (double *)malloc(sizeof(double) * num_of_nodes);

    // inititalize matrices.
    // q_matrix(q, h);
    // nodes_matrix(nodes, num_of_nodes, n, h);
    // elements_matrix(elements, num_of_elements, n);
    // Q = Q_matrix(h, nodes, elements, 0);
    // B_matrix(B, h);
    // A_matrix(A, B, A1);

    // /***********global***************/
    // assemble_G_A(G_A, A, num_of_elements, num_of_nodes, elements);
    // assembel_G_Q(G_Q, h, nodes, num_of_elements, num_of_nodes, elements);
    // assemble_G_q(G_q, num_of_nodes, n, h);
    // assemble_G_B(G_B, B, num_of_elements, num_of_nodes, elements);
    // huayifa_G_Qaq(G_Q, G_q, n, num_of_nodes);
    // huayifa_G_A(G_A, num_of_elements, num_of_nodes, n, h);
    // huayifa_G_B(G_B, num_of_elements, num_of_nodes, n, h);
    /******************** petsc *********************/

    ierr = PetscInitialize(&argc, &args, (char *)0, help);
    if (ierr)
        return ierr; /*初始化Petsc*/
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    CHKERRQ(ierr); /*设置并行MPI参数*/
    ierr = VecCreate(PETSC_COMM_WORLD, &T);
    CHKERRQ(ierr); /*创建一个并行空间*/
    ierr = VecSetSizes(T, PETSC_DECIDE, num_of_nodes);
    CHKERRQ(ierr); /*创建一个长度(n+1)*(n+1)的向量*/
    ierr = VecSetFromOptions(T);
    CHKERRQ(ierr); /*从选项数据库中配置向量*/
    ierr = VecSet(T, 0);
    CHKERRQ(ierr);
    ierr = VecDuplicate(T, &QQ);
    CHKERRQ(ierr); /*将T的格式赋给Q*/
    ierr = VecDuplicate(T, &b);
    CHKERRQ(ierr); /*将T的格式赋给b*/
    ierr = VecDuplicate(T, &x);
    CHKERRQ(ierr); /*将T的格式赋给x*/
    ierr = VecDuplicate(T, &t);
    CHKERRQ(ierr); /*将T的格式赋给t*/
    ierr = VecDuplicate(T, &w);
    CHKERRQ(ierr); /*将T的格式赋给t*/
    ierr = VecAssemblyBegin(T);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(T);
    CHKERRQ(ierr);

    ierr = VecCopy(T, b);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(b);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);
    CHKERRQ(ierr);

    ierr = VecCopy(T, x);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(x);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x);
    CHKERRQ(ierr);

    ierr = VecCopy(T, t);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(t);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(t);
    CHKERRQ(ierr);

    ierr = VecCopy(T, w);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(w);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(w);
    CHKERRQ(ierr);

    ierr = VecSet(QQ, 0);
    CHKERRQ(ierr);
    if (rank == 0)
    {
        for (i = 0; i < num_of_nodes; i++)
        {
            q0 = G_Q[i];
            ierr = VecSetValues(QQ, 1, &i, &q0, INSERT_VALUES);
            CHKERRQ(ierr); /*将向量的对应位置的值进行修改*/
        }
    }
    ierr = VecAssemblyBegin(QQ);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(QQ);
    CHKERRQ(ierr);
    ierr = VecView(QQ, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ(ierr);

    ierr = MatCreate(PETSC_COMM_WORLD, &AA);
    CHKERRQ(ierr); /*在并行空间创建一个矩阵A*/
    ierr = MatSetSizes(AA, PETSC_DECIDE, PETSC_DECIDE, num_of_nodes, num_of_nodes);
    CHKERRQ(ierr); /*设置矩阵的行数和列数*/
    ierr = MatSetFromOptions(AA);
    CHKERRQ(ierr); /*从选项数据库中配置矩阵*/
    ierr = MatSetUp(AA);
    CHKERRQ(ierr); /*开始建立矩阵*/
    for (ii = 0; ii < num_of_nodes; ii++)
    {
        if (ii == 0)
        {
            value[0] = G_A[ii][ii];
            value[1] = G_A[ii][ii + 1];
            value[2] = G_A[ii][ii + n + 1];
            value[3] = G_A[ii][ii + n + 2];
            col4[0] = ii;
            col4[1] = ii + 1;
            col4[2] = ii + n + 1;
            col4[3] = ii + n + 2;
            ierr = MatSetValues(AA, 1, &ii, 4, col4, value, INSERT_VALUES);
            CHKERRQ(ierr);
        }
        else if (ii > 0 && ii < n + 1)
        {
            value[0] = G_A[ii][ii - 1];
            value[1] = G_A[ii][ii];
            value[2] = G_A[ii][ii + 1];
            value[3] = G_A[ii][ii + n];
            value[4] = G_A[ii][ii + n + 1];
            value[5] = G_A[ii][ii + n + 2];
            col6[0] = ii - 1;
            col6[1] = ii;
            col6[2] = ii + 1;
            col6[3] = ii + n;
            col6[4] = ii + n + 1;
            col6[5] = ii + n + 2;
            ierr = MatSetValues(AA, 1, &ii, 6, col6, value, INSERT_VALUES);
            CHKERRQ(ierr);
        }
        else if (ii == n + 1)
        {
            value[0] = G_A[ii][ii - n - 1];
            value[1] = G_A[ii][ii - n];
            value[2] = G_A[ii][ii - 1];
            value[3] = G_A[ii][ii];
            value[4] = G_A[ii][ii + 1];
            value[5] = G_A[ii][ii + n];
            value[6] = G_A[ii][ii + n + 1];
            value[7] = G_A[ii][ii + n + 2];
            col8[0] = ii - n - 1;
            col8[1] = ii - n;
            col8[2] = ii - 1;
            col8[3] = ii;
            col8[4] = ii + 1;
            col8[5] = ii + n;
            col8[6] = ii + n + 1;
            col8[7] = ii + n + 2;
            ierr = MatSetValues(AA, 1, &ii, 8, col8, value, INSERT_VALUES);
            CHKERRQ(ierr);
        }
        else if (ii > n + 1 && ii < (n + 1) * n - 1)
        {
            value[0] = G_A[ii][ii - n - 2];
            value[1] = G_A[ii][ii - n - 1];
            value[2] = G_A[ii][ii - n];
            value[3] = G_A[ii][ii - 1];
            value[4] = G_A[ii][ii];
            value[5] = G_A[ii][ii + 1];
            value[6] = G_A[ii][ii + n];
            value[7] = G_A[ii][ii + n + 1];
            value[8] = G_A[ii][ii + n + 2];
            col9[0] = ii - n - 2;
            col9[1] = ii - n - 1;
            col9[2] = ii - n;
            col9[3] = ii - 1;
            col9[4] = ii;
            col9[5] = ii + 1;
            col9[6] = ii + n;
            col9[7] = ii + n + 1;
            col9[8] = ii + n + 2;
            ierr = MatSetValues(AA, 1, &ii, 9, col9, value, INSERT_VALUES);
            CHKERRQ(ierr);
        }
        else if (ii == (n + 1) * n - 1)
        {
            value[0] = G_A[ii][ii - n - 2];
            value[1] = G_A[ii][ii - n - 1];
            value[2] = G_A[ii][ii - n];
            value[3] = G_A[ii][ii - 1];
            value[4] = G_A[ii][ii];
            value[5] = G_A[ii][ii + 1];
            value[6] = G_A[ii][ii + n];
            value[7] = G_A[ii][ii + n + 1];
            col8[0] = ii - n - 2;
            col8[1] = ii - n - 1;
            col8[2] = ii - n;
            col8[3] = ii - 1;
            col8[4] = ii;
            col8[5] = ii + 1;
            col8[6] = ii + n;
            col8[7] = ii + n + 1;
            ierr = MatSetValues(AA, 1, &ii, 8, col8, value, INSERT_VALUES);
            CHKERRQ(ierr);
        }
        else if (ii > (n + 1) * n - 1 && ii < (n + 1) * (n + 1) - 1)
        {
            value[0] = G_A[ii][ii - n - 2];
            value[1] = G_A[ii][ii - n - 1];
            value[2] = G_A[ii][ii - n];
            value[3] = G_A[ii][ii - 1];
            value[4] = G_A[ii][ii];
            value[5] = G_A[ii][ii + 1];
            col6[0] = ii - n - 2;
            col6[1] = ii - n - 1;
            col6[2] = ii - n;
            col6[3] = ii - 1;
            col6[4] = ii;
            col6[5] = ii + 1;
            ierr = MatSetValues(AA, 1, &ii, 6, col6, value, INSERT_VALUES);
            CHKERRQ(ierr);
        }
        else if (ii == (n + 1) * (n + 1) - 1)
        {
            value[0] = G_A[ii][ii - n - 2];
            value[1] = G_A[ii][ii - n - 1];
            value[2] = G_A[ii][ii - 1];
            value[3] = G_A[ii][ii];
            col4[0] = ii - n - 2;
            col4[1] = ii - n - 1;
            col4[2] = ii - 1;
            col4[3] = ii;
            ierr = MatSetValues(AA, 1, &ii, 4, col4, value, INSERT_VALUES);
            CHKERRQ(ierr);
        }
    }
    /* Assemble the A matrix */
    ierr = MatAssemblyBegin(AA, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr); /*通知其余并行块将矩阵统一*/
    ierr = MatAssemblyEnd(AA, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr); /*结束通知*/
    ierr = MatView(AA, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ(ierr); /*打印矩阵，检查是否出错*/

    ierr = MatCreate(PETSC_COMM_WORLD, &BB);
    CHKERRQ(ierr); /*在并行空间创建一个矩阵B*/
    ierr = MatSetSizes(BB, PETSC_DECIDE, PETSC_DECIDE, num_of_nodes, num_of_nodes);
    CHKERRQ(ierr); /*设置矩阵的行数和列数*/
    ierr = MatSetFromOptions(BB);
    CHKERRQ(ierr); /*从选项数据库中配置矩阵*/
    ierr = MatSetUp(BB);
    CHKERRQ(ierr); /*开始建立矩阵*/
    for (ii = 0; ii < num_of_nodes; ii++)
    {
        if (ii == 0)
        {
            value[0] = G_B[ii][ii];
            value[1] = G_B[ii][ii + 1];
            value[2] = G_B[ii][ii + n + 1];
            value[3] = G_B[ii][ii + n + 2];
            col4[0] = ii;
            col4[1] = ii + 1;
            col4[2] = ii + n + 1;
            col4[3] = ii + n + 2;
            ierr = MatSetValues(BB, 1, &ii, 4, col4, value, INSERT_VALUES);
            CHKERRQ(ierr);
        }
        else if (ii > 0 && ii < n + 1)
        {
            value[0] = G_B[ii][ii - 1];
            value[1] = G_B[ii][ii];
            value[2] = G_B[ii][ii + 1];
            value[3] = G_B[ii][ii + n];
            value[4] = G_B[ii][ii + n + 1];
            value[5] = G_B[ii][ii + n + 2];
            col6[0] = ii - 1;
            col6[1] = ii;
            col6[2] = ii + 1;
            col6[3] = ii + n;
            col6[4] = ii + n + 1;
            col6[5] = ii + n + 2;
            ierr = MatSetValues(BB, 1, &ii, 6, col6, value, INSERT_VALUES);
            CHKERRQ(ierr);
        }
        else if (ii == n + 1)
        {
            value[0] = G_B[ii][ii - n - 1];
            value[1] = G_B[ii][ii - n];
            value[2] = G_B[ii][ii - 1];
            value[3] = G_B[ii][ii];
            value[4] = G_B[ii][ii + 1];
            value[5] = G_B[ii][ii + n];
            value[6] = G_B[ii][ii + n + 1];
            value[7] = G_B[ii][ii + n + 2];
            col8[0] = ii - n - 1;
            col8[1] = ii - n;
            col8[2] = ii - 1;
            col8[3] = ii;
            col8[4] = ii + 1;
            col8[5] = ii + n;
            col8[6] = ii + n + 1;
            col8[7] = ii + n + 2;
            ierr = MatSetValues(BB, 1, &ii, 8, col8, value, INSERT_VALUES);
            CHKERRQ(ierr);
        }
        else if (ii > n + 1 && ii < (n + 1) * n - 1)
        {
            value[0] = G_B[ii][ii - n - 2];
            value[1] = G_B[ii][ii - n - 1];
            value[2] = G_B[ii][ii - n];
            value[3] = G_B[ii][ii - 1];
            value[4] = G_B[ii][ii];
            value[5] = G_B[ii][ii + 1];
            value[6] = G_B[ii][ii + n];
            value[7] = G_B[ii][ii + n + 1];
            value[8] = G_B[ii][ii + n + 2];
            col9[0] = ii - n - 2;
            col9[1] = ii - n - 1;
            col9[2] = ii - n;
            col9[3] = ii - 1;
            col9[4] = ii;
            col9[5] = ii + 1;
            col9[6] = ii + n;
            col9[7] = ii + n + 1;
            col9[8] = ii + n + 2;
            ierr = MatSetValues(BB, 1, &ii, 9, col9, value, INSERT_VALUES);
            CHKERRQ(ierr);
        }
        else if (ii == (n + 1) * n - 1)
        {
            value[0] = G_B[ii][ii - n - 2];
            value[1] = G_B[ii][ii - n - 1];
            value[2] = G_B[ii][ii - n];
            value[3] = G_B[ii][ii - 1];
            value[4] = G_B[ii][ii];
            value[5] = G_B[ii][ii + 1];
            value[6] = G_B[ii][ii + n];
            value[7] = G_B[ii][ii + n + 1];
            col8[0] = ii - n - 2;
            col8[1] = ii - n - 1;
            col8[2] = ii - n;
            col8[3] = ii - 1;
            col8[4] = ii;
            col8[5] = ii + 1;
            col8[6] = ii + n;
            col8[7] = ii + n + 1;
            ierr = MatSetValues(BB, 1, &ii, 8, col8, value, INSERT_VALUES);
            CHKERRQ(ierr);
        }
        else if (ii > (n + 1) * n - 1 && ii < (n + 1) * (n + 1) - 1)
        {
            value[0] = G_B[ii][ii - n - 2];
            value[1] = G_B[ii][ii - n - 1];
            value[2] = G_B[ii][ii - n];
            value[3] = G_B[ii][ii - 1];
            value[4] = G_B[ii][ii];
            value[5] = G_B[ii][ii + 1];
            col6[0] = ii - n - 2;
            col6[1] = ii - n - 1;
            col6[2] = ii - n;
            col6[3] = ii - 1;
            col6[4] = ii;
            col6[5] = ii + 1;
            ierr = MatSetValues(BB, 1, &ii, 6, col6, value, INSERT_VALUES);
            CHKERRQ(ierr);
        }
        else if (ii == (n + 1) * (n + 1) - 1)
        {
            value[0] = G_B[ii][ii - n - 2];
            value[1] = G_B[ii][ii - n - 1];
            value[2] = G_B[ii][ii - 1];
            value[3] = G_B[ii][ii];
            col4[0] = ii - n - 2;
            col4[1] = ii - n - 1;
            col4[2] = ii - 1;
            col4[3] = ii;
            ierr = MatSetValues(BB, 1, &ii, 4, col4, value, INSERT_VALUES);
            CHKERRQ(ierr);
        }
    }
    /* Assemble the B matrix */
    ierr = MatAssemblyBegin(BB, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr); /*通知其余并行块将矩阵统一*/
    ierr = MatAssemblyEnd(BB, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr); /*结束通知*/
    ierr = MatView(BB, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ(ierr); /*打印矩阵，检查是否出错*/
    /* Q-A*T */
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
    CHKERRQ(ierr); /*创建ksp解空间*/
    ierr = KSPSetOperators(ksp, BB, BB);
    CHKERRQ(ierr); /*设置方程左侧的系数*/
    ierr = KSPGetPC(ksp, &pc);
    CHKERRQ(ierr); /*设置矩阵求解的相关系数*/
    ierr = PCSetType(pc, PCJACOBI);
    CHKERRQ(ierr); /*设置pc的默认参数*/
    ierr = KSPSetTolerances(ksp, 1.e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    CHKERRQ(ierr); /*设置各种误差值*/
    ierr = KSPSetFromOptions(ksp);
    CHKERRQ(ierr); /*从选项数据库中配置ksp解空间*/

    while (PetscAbsReal(k) < 51)
    {
        k += 1;
        ierr = MatMult(AA, t, b);
        ierr = VecWAXPY(w, -1, t, QQ);
        ierr = KSPSolve(ksp, w, T);
        CHKERRQ(ierr);
        ierr = VecCopy(T, t);
        CHKERRQ(ierr); /*将T的值赋给t*/
    }
    ierr = VecView(T, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ(ierr); /*打印向量*/

    ierr = VecDestroy(&T);
    CHKERRQ(ierr); /*关闭向量T*/
    ierr = VecDestroy(&QQ);
    CHKERRQ(ierr); /*关闭向量Q*/
    ierr = VecDestroy(&x);
    CHKERRQ(ierr); /*关闭向量x*/
    ierr = VecDestroy(&t);
    CHKERRQ(ierr); /*关闭向量t*/
    ierr = VecDestroy(&w);
    CHKERRQ(ierr); /*关闭向量w*/
    ierr = VecDestroy(&b);
    CHKERRQ(ierr); /*关闭向量b*/
    ierr = MatDestroy(&AA);
    CHKERRQ(ierr); /*关闭矩阵A*/
    ierr = MatDestroy(&BB);
    CHKERRQ(ierr); /*关闭矩阵B*/

    ierr = PetscFinalize(); /*结束并行*/
    /*************************************** petsc end ***************************************/
    // free matrices
    freeMatrix(G_A);
    freeMatrix(G_B);
    free(G_Q);
    return ierr;
}

void nodes_matrix(double nodes[][3], int num_of_nodes, int n, double h)
{
    for (int i = 0; i < num_of_nodes; i++)
    {
        int a = i / (n + 1);
        nodes[i][0] = i;
        nodes[i][1] = (i - a * (n + 1)) * h;
        nodes[i][2] = a * h;
    }
}

void elements_matrix(int elements[][5], int num_of_elements, int n)
{
    for (int i = 0; i < num_of_elements; i++)
    {
        int a = i / (n + 1);
        elements[i][0] = i;
        elements[i][1] = i + a;
        elements[i][2] = i + a + 1;
        elements[i][3] = i + a + n + 1;
        elements[i][4] = i + a + n + 2;
    }
}

void q_matrix(double q[4], double h)
{
    q[0] = 0;
    q[1] = 1 * h / 2;
    q[2] = 1 * h / 2;
    q[3] = 0;
}

double *Q_matrix(double h, double nodes[][3], int elements[][5], int k)
{
    double xi[4], eta[4];
    double *Q = (double *)malloc(sizeof(double) * 4);
    xi[0] = -0.5773;
    xi[1] = 0.5773;
    xi[2] = 0.5773;
    xi[3] = -0.5773;

    eta[0] = -0.5773;
    eta[1] = -0.5773;
    eta[2] = 0.5773;
    eta[3] = 0.5773;

    for (int i = 0; i < 4; i++)
    {
        Q[i] = 0;
    }

    for (int jj = 0; jj < 4; jj++)
    {
        Q[0] += (h * h / 32.0) * (sin(PI * (h / 2.0 * xi[jj] + nodes[elements[k][1]][1] + h / 2.0)) + sin(PI * (h / 2.0 * eta[jj] + nodes[elements[k][1]][2] + h / 2.0))) * (1 - xi[jj]) * (1 - eta[jj]);
        Q[1] += (h * h / 32.0) * (sin(PI * (h / 2.0 * xi[jj] + nodes[elements[k][1]][1] + h / 2.0)) + sin(PI * (h / 2.0 * eta[jj] + nodes[elements[k][1]][2] + h / 2.0))) * (1 + xi[jj]) * (1 - eta[jj]);
        Q[2] += (h * h / 32.0) * (sin(PI * (h / 2.0 * xi[jj] + nodes[elements[k][1]][1] + h / 2.0)) + sin(PI * (h / 2.0 * eta[jj] + nodes[elements[k][1]][2] + h / 2.0))) * (1 + xi[jj]) * (1 + eta[jj]);
        Q[3] += (h * h / 32.0) * (sin(PI * (h / 2.0 * xi[jj] + nodes[elements[k][1]][1] + h / 2.0)) + sin(PI * (h / 2.0 * eta[jj] + nodes[elements[k][1]][2] + h / 2.0))) * (1 - xi[jj]) * (1 + eta[jj]);
    }
    return Q;
}

void B_matrix(double B[][4], double h, double dt)
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            B[i][j] *= (h * h / dt);
        }
    }
}

void A_matrix(double A[][4], double B[][4], double A1[][4])
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            A[i][j] = B[i][j] + A1[i][j];
        }
    }
}

void assemble_G_A(double **G_A, double A[][4], int num_of_elements, int num_of_nodes, int elments[][5])
{
    for (int i = 0; i < num_of_nodes; i++)
    {
        for (int j = 0; j < num_of_nodes; j++)
        {
            G_A[i][j] = 0.0;
        }
    }
    for (int n = 0; n < num_of_elements; n++)
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                G_A[elments[n][i + 1]][elments[n][j + 1]] += A[i][j];
            }
        }
    }
}

void assemble_G_B(double **G_B, double B[][4], int num_of_elements, int num_of_nodes, int elments[][5])
{
    for (int i = 0; i < num_of_nodes; i++)
    {
        for (int j = 0; j < num_of_nodes; j++)
        {
            G_B[i][j] = 0.0;
        }
    }
    for (int n = 0; n < num_of_elements; n++)
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                G_B[elments[n][i + 1]][elments[n][j + 1]] += B[i][j];
            }
        }
    }
}

void assembel_G_Q(double *G_Q, double h, double nodes[][3], int num_of_elements, int num_of_nodes, int elments[][5])
{
    double *Q;
    for (int i = 0; i < num_of_nodes; i++)
    {
        G_Q[i] = 0;
    }
    for (int n = 0; n < num_of_elements; n++)
    {
        Q = Q_matrix(h, nodes, elments, n);
        for (int i = 0; i < num_of_nodes; i++)
        {
            G_Q[elments[n][i + 1]] += Q[i];
        }
    }
    free(Q);
}

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

double **mallocMatrix(int row, int col)
{
    double **ret = (double **)malloc(sizeof(double *) * row); // a row pointer
    double *p = (double *)malloc(sizeof(double) * row * col); // pointer of whole array

    int i = 0;
    if (ret && p) // check whether they are both null.
    {
        for (i = 0; i < row; i++)
        {
            ret[i] = (p + i * col); // ret[i] is a row pointer that points to each row.
        }
    }
    else
    {
        free(ret);
        free(p);
        ret = NULL;
        p = NULL;
    }

    return ret;
}

void freeMatrix(double **a)
{
    free(a[0]); // free array first.
    free(a);    // free row pointer.
}
