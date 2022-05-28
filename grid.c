#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <hdf5.h>
#include <assert.h>
#define PI 3.14159265
#define FILE "results.h5"

void q_matrix(double q[4], double h);
double *Q_matrix(double h, double nodes[][3], int elements[][5], int i);
void B_matrix(double B[][4], double h, double dt);
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

int main(int argc, char *argv[])
{
    double h = 0.5;
    double dt = 0.5;
    // assert(argc == 4);
    // h = atof(argv[1]);
    // dt = atof(argv[2]);

    int length = 1;
    int n = length / h;
    int num_of_elements = n * n;
    int num_of_nodes = (n + 1) * (n + 1);
    double nodes[num_of_nodes][3];
    double q[4];
    double T[num_of_nodes];
    double Tdt[num_of_nodes];
    double B[4][4] = {{0.1111, 0.0556, 0.0278, 0.0556}, {0.0556, 0.1111, 0.0556, 0.0278}, {0.0278, 0.0556, 0.1111, 0.0556}, {0.0556, 0.0278, 0.0556, 0.1111}};
    int elements[num_of_elements][5];
    double A1[4][4] = {{0.6667, -0.1667, -0.3333, -0.1667}, {-0.1667, 0.6667, -0.1667, -0.3333}, {-0.3333, -0.1667, 0.6667, -0.1667}, {-0.1667, -0.3333, -0.1667, 0.6667}};
    double A[4][4] = {0};
    double *Q;
//    double **G_A = mallocMatrix(num_of_nodes, num_of_nodes);
    double **G_B = mallocMatrix(num_of_nodes, num_of_nodes);
    // double **G_Q = mallocMatrix(num_of_nodes, num_of_nodes);
    // double **G_q = mallocMatrix(num_of_nodes, num_of_nodes);
    double *G_Q = (double *)malloc(sizeof(double) * num_of_nodes);
    double *G_q = (double *)malloc(sizeof(double) * num_of_nodes);
	double *G_A =(double *)malloc(sizeof(double) * num_of_nodes* num_of_nodes);
    hid_t file_id, group_id, dataspace_G_B_id, dataspace_G_A_id, dataspace_G_Q_id, dataspace_G_q_id, dataspace_T_id, dataspace_Tdt_id, dataset_G_B_id, dataset_G_A_id, dataset_G_Q_id, dataset_G_q_id, dataset_T_id, dataset_Tdt_id;
    herr_t status;
    hsize_t G_B_dims[2], G_A_dims[2], G_Q_dims = num_of_nodes, G_q_dims = num_of_nodes, T_dims = num_of_nodes, Tdt_dims = num_of_nodes; // see whether dims can be only one number,i.e. for q and Q.
    G_B_dims[0] = num_of_nodes;
    G_B_dims[1] = num_of_nodes;
    G_A_dims[0] = num_of_nodes;
    G_A_dims[1] = num_of_nodes;

    // if (*argv[3] == 'r')
    //{
        file_id = H5Fopen(FILE, H5F_ACC_RDWR, H5P_DEFAULT);
        dataset_G_B_id=H5Dopen(file_id,"/results/G_B",H5P_DEFAULT);
        dataset_G_A_id=H5Dopen(file_id,"/results/G_A",H5P_DEFAULT);
        dataset_G_Q_id=H5Dopen(file_id,"/results/G_Q",H5P_DEFAULT);
        dataset_G_q_id=H5Dopen(file_id,"/results/G_q",H5P_DEFAULT);

    //    status = H5Dread(dataset_G_B_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, G_B);
        status = H5Dread(dataset_G_A_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, G_A);
      //  status = H5Dread(dataset_G_Q_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, G_Q);
       // status = H5Dread(dataset_G_q_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, G_q);
  /*
        status = H5Dclose(dataset_G_B_id);
        status = H5Dclose(dataset_G_A_id);
        status = H5Dclose(dataset_G_Q_id);
        status = H5Dclose(dataset_G_q_id);
        status = H5Fclose(file_id);
  */ // }
  //  else
   // {

/*        file_id = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        group_id = H5Gcreate(file_id, "/results", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        dataspace_G_B_id = H5Screate_simple(2, G_B_dims, NULL);
        dataset_G_B_id = H5Dcreate2(file_id, "/results/G_B", H5T_NATIVE_DOUBLE, dataspace_G_B_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        dataspace_G_A_id = H5Screate_simple(2, G_A_dims, NULL);
        dataset_G_A_id = H5Dcreate2(file_id, "/results/G_A", H5T_NATIVE_DOUBLE, dataspace_G_A_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        dataspace_G_Q_id = H5Screate_simple(1, &G_Q_dims, NULL);
        dataset_G_Q_id = H5Dcreate2(file_id, "/results/G_Q", H5T_NATIVE_DOUBLE, dataspace_G_Q_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        dataspace_G_q_id = H5Screate_simple(1, &G_q_dims, NULL);
        dataset_G_q_id = H5Dcreate2(file_id, "/results/G_q", H5T_NATIVE_DOUBLE, dataspace_G_q_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // inititalize matrices.
        q_matrix(q, h);
        nodes_matrix(nodes, num_of_nodes, n, h);
        elements_matrix(elements, num_of_elements, n);
        Q = Q_matrix(h, nodes, elements, 0);
        B_matrix(B, h, dt);
        A_matrix(A, B, A1);
        assemble_G_A(G_A, A, num_of_elements, num_of_nodes, elements);
        assembel_G_Q(G_Q, h, nodes, num_of_elements, num_of_nodes, elements);
        assemble_G_q(G_q, num_of_nodes, n, h);
        assemble_G_B(G_B, B, num_of_elements, num_of_nodes, elements);
        huayifa_G_Qaq(G_Q, G_q, n, num_of_nodes);
        huayifa_G_A(G_A, num_of_elements, num_of_nodes, n, h);
        huayifa_G_B(G_B, num_of_elements, num_of_nodes, n, h);

        /************single*****************/
/*
        status = H5Dwrite(dataset_G_B_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, G_B);
        status = H5Dwrite(dataset_G_A_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, G_A);
        status = H5Dwrite(dataset_G_Q_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, G_Q);
        status = H5Dwrite(dataset_G_q_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, G_q);

        status = H5Dclose(dataset_G_B_id);
        status = H5Sclose(dataspace_G_B_id);
        status = H5Dclose(dataset_G_A_id);
        status = H5Sclose(dataspace_G_A_id);
        status = H5Dclose(dataset_G_Q_id);
        status = H5Sclose(dataspace_G_Q_id);
        status = H5Dclose(dataset_G_q_id);
        status = H5Sclose(dataspace_G_q_id);

        status = H5Gclose(group_id);
        status = H5Fclose(file_id);
        
 */ //  }


        }
    */

    /*  printf("single A matrix:\n");
      for (int i = 0; i < 4; i++)
      {
          for (int j = 0; j < 4; j++)
          {
              printf("%lf ", A[i][j]);
          }
          printf("\n");
      }

      printf("single Q matrix:\n");

      for (int j = 0; j < 4; j++)
      {
          printf("%f \n", Q[j]);
      }
      printf("\n");

      printf("single B matrix:\n");
      for (int i = 0; i < 4; i++)
      {
          for (int j = 0; j < 4; j++)
          {
              printf("%lf ", B[i][j]);
          }
          printf("\n");
      }
  */
    /***********global***************/

    printf("global A matrix:\n");
printf("%lf ", G_A[3]);
/*
    for (int i = 0; i < num_of_nodes; i++)
    {
        for (int j = 0; j < num_of_nodes; j++)
        {
            printf("%lf ", G_A[i][j]);
        }
        printf("\n");
<<<<<<< HEAD
    }
    /*
        printf("global Q matrix:\n");
=======
    }*/
/*
    printf("global Q matrix:\n");

    for (int j = 0; j < num_of_nodes; j++)
    {
        printf("%lf ", G_Q[j]);
    }
    printf("\n");
>>>>>>> 3a89d2e03beb3d0c7bcb57df7db016a68c5af7d0

        for (int j = 0; j < num_of_nodes; j++)
        {
            printf("%lf ", G_Q[j]);
        }
        printf("\n");

        printf("global q matrix:\n");
        for (int j = 0; j < num_of_nodes; j++)
        {
            printf("%lf ", G_q[j]);
        }
        printf("\n");

        printf("global B matrix:\n");
        for (int i = 0; i < num_of_nodes; i++)
        {
            for (int j = 0; j < num_of_nodes; j++)
            {
                printf("%lf ", G_B[i][j]);
            }
            printf("\n");
        }

        printf("huayifa Q matrix:\n");

        for (int j = 0; j < num_of_nodes; j++)
        {
            printf("%lf ", G_Q[j]);
        }
        printf("\n");

        printf("huayifa A matrix:\n");

        for (int i = 0; i < num_of_nodes; i++)
        {
            for (int j = 0; j < num_of_nodes; j++)
            {
                printf("%lf ", G_A[i][j]);
            }
            printf("\n");
        }

        printf("huayifa B matrix:\n");

        for (int i = 0; i < num_of_nodes; i++)
        {
            for (int j = 0; j < num_of_nodes; j++)
            {
                printf("%lf ", G_B[i][j]);
            }
            printf("\n");
        }*/
    status = H5Dclose(dataset_G_B_id);
<<<<<<< HEAD
    status = H5Dclose(dataset_G_A_id);
    status = H5Dclose(dataset_G_Q_id);
    status = H5Dclose(dataset_G_q_id);
    status = H5Fclose(file_id);
    // free matrices
    freeMatrix(G_A);
=======
        status = H5Dclose(dataset_G_A_id);
        status = H5Dclose(dataset_G_Q_id);
        status = H5Dclose(dataset_G_q_id);
        status = H5Fclose(file_id);
  // free matrices
   // freeMatrix(G_A);
>>>>>>> 3a89d2e03beb3d0c7bcb57df7db016a68c5af7d0
    freeMatrix(G_B);
    free(G_Q);
    free(G_q);
    return 0;
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
