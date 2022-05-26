#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265

void q_matrix(double q[4], double h);
double* Q_matrix(double h, double nodes[][3], int elements[][5], int i);
void B_matrix(double B[][4], double h);
void elements_matrix(int elements[][5], int num_of_elements, int n);
void nodes_matrix(double nodes[][3], int num_of_nodes, int n, double h);
void A_matrix(double A[][4], double B[][4], double A1[][4]);
void assemble_G_A(double **G_A, double A[][4], int num_of_elements, int num_of_nodes, int elments[][5]);
void assemble_G_B(double **G_B, double B[][4], int num_of_elements, int num_of_nodes, int elments[][5]);
void assembel_G_Q(double *G_Q, double h,double nodes[][3], int num_of_elements, int num_of_nodes, int elments[][5]);
double **mallocMatrix(int row, int col);
void freeMatrix(double **a);

int main(int argc, char *argv[])
{
	double h;
	if (argc == 2)
	{
		h = atof(argv[1]);
	}
	int i = 0, j = 0;
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
	double **G_A = mallocMatrix(num_of_nodes, num_of_nodes);
	double **G_B = mallocMatrix(num_of_nodes, num_of_nodes);
	//double **G_Q = mallocMatrix(num_of_nodes, num_of_nodes);
	//double **G_q = mallocMatrix(num_of_nodes, num_of_nodes);
	double *G_Q=(double *)malloc(sizeof(double)*num_of_nodes);

	// inititalize matrices.
	q_matrix(q, h);
	nodes_matrix(nodes, num_of_nodes, n, h);
	elements_matrix(elements, num_of_elements, n);
	Q_matrix( h, nodes, elements, i);
	B_matrix(B, h);
	A_matrix(A, B, A1);



	//free matrices
	freeMatrix(G_A);
	freeMatrix(G_B);
	freeMatrix(G_Q);
	freeMatrix(G_q);
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

double* Q_matrix(double h, double nodes[][3], int elements[][5], int n)
{
	double xi[4], eta[4];
	double* Q=(double*)malloc(sizeof(double)*4);
	xi[0] = -0.5773;
	xi[1] = 0.5773;
	xi[2] = 0.5773;
	xi[3] = -0.5773;

	eta[0] = -0.5773;
	eta[1] = -0.5773;
	eta[2] = 0.5773;
	eta[3] = 0.5773;

	for (int jj = 0; jj < 4; jj++)
	{
		Q[0] += (h * h / 32.0) * (sin(PI * (h / 2.0 * xi[jj] + nodes[elements[n][1]][1] + h / 2.0)) + sin(PI * (h / 2.0 * eta[jj] + nodes[elements[n][1]][2] + h / 2.0))) * (1 - xi[jj]) * (1 - eta[jj]);
		Q[1] += (h * h / 32.0) * (sin(PI * (h / 2.0 * xi[jj] + nodes[elements[n][1]][1] + h / 2.0)) + sin(PI * (h / 2.0 * eta[jj] + nodes[elements[n][1]][2] + h / 2.0))) * (1 + xi[jj]) * (1 - eta[jj]);
		Q[2] += (h * h / 32.0) * (sin(PI * (h / 2.0 * xi[jj] + nodes[elements[n][1]][1] + h / 2.0)) + sin(PI * (h / 2.0 * eta[jj] + nodes[elements[n][1]][2] + h / 2.0))) * (1 + xi[jj]) * (1 + eta[jj]);
		Q[3] += (h * h / 32.0) * (sin(PI * (h / 2.0 * xi[jj] + nodes[elements[n][1]][1] + h / 2.0)) + sin(PI * (h / 2.0 * eta[jj] + nodes[elements[n][1]][2] + h / 2.0))) * (1 - xi[jj]) * (1 + eta[jj]);
	}
	return Q;
}

void B_matrix(double B[][4], double h)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			B[i][j] *= h;
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
				G_A[elments[n][i+1]][elments[n][j+1]] += A[i][j];
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
				G_B[elments[n][i+1]][elments[n][j+1]] += B[i][j];
			}
		}
	}
}


void assembel_G_Q(double *G_Q, double h,double nodes[][3], int num_of_elements, int num_of_nodes, int elments[][5])
{
	double* Q;
	for (int i = 0; i < num_of_nodes; i++)
	{
		G_Q[i]=0;
	}
	for (int n = 0; n < num_of_elements; n++)
	{
		Q=Q_matrix(h,nodes,elments,n);
		for (int i = 0; i < num_of_nodes; i++)
		{
			G_Q[elments[n][i+1]]+=Q[i];
		}	
	}
	free(Q);
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
	free(a);	// free row pointer.
}
