#include<stdio.h>
#include<stdlib.h>
#include<math.h>

# define PI 3.14159265

void q_matrix(float q[4] ,float h);
void Q_matrix(float Q[4],float h,float nodes[][3],int elements[][5],int i);
void B_matrix(float B[][4], float h);
void elements_matrix(int elements[][5],int num_of_elements,int n);
void nodes_matrix(float nodes[][3],int num_of_nodes,int n,float h);
void A_matrix(float A[][4], float B[][4], float A1[][4]);

int main(int argc,char* argv[]){
	float h;	
	if(argc==2)
	{
		h=atof(argv[1]);
	}	
	int i=0,j=0;
	int length=1;
	int n=length/h;
	int num_of_elements = n*n;
	int num_of_nodes=(n+1)*(n+1);
	float nodes[num_of_nodes][3];
	float q[4];
	float Q[4]={0};
	float T[num_of_nodes];
	float Tdt[num_of_nodes];
	float B[4][4]={{0.1111,0.0556,0.0278,0.0556},{0.0556,0.1111,0.0556,0.0278},{0.0278,0.0556,0.1111,0.0556},{0.0556,0.0278,0.0556,0.1111}};
	int elements[num_of_elements][5];
	float A1[4][4]={{0.6667,-0.1667,-0.3333,-0.1667},{-0.1667,0.6667 , 0 , -0.3333},{-0.3333,-0.1667,0.6667,-0.1667},{-0.1667,-0.3333,-0.1667,0.6667}};
	float A[4][4]={0};
	float G_A[num_of_nodes][num_of_nodes];
	float G_B[num_of_nodes][num_of_nodes];
	float G_Q[num_of_nodes][num_of_nodes];
	float G_q[num_of_nodes][num_of_nodes];

	//inititalize matrices.
	q_matrix(q,h);
	nodes_matrix(nodes,num_of_nodes, n,h);
	elements_matrix(elements,num_of_elements,n);
	Q_matrix(Q,h,nodes,elements,i);
	B_matrix(B,h);
	A_matrix(A,B,A1);



	return 0;
}

void nodes_matrix(float nodes[][3],int num_of_nodes,int n,float h)
{
	for (int i=0;i<num_of_nodes;i++)
	{
		int a=i/(n+1);
		nodes[i][0]=i;
		nodes[i][1]=(i-a*(n+1))*h;
		nodes[i][2]=a*h;
	}
}

void elements_matrix(int elements[][5],int num_of_elements,int n)
{
	for (int i=0;i<num_of_elements;i++)
	{
		int a=i/(n+1);
		elements[i][0]=i;
		elements[i][1]=i+a;
		elements[i][2]=i+a+1;
		elements[i][3]=i+a+n+1;
		elements[i][4]=i+a+n+2;
		
	}
}


void  q_matrix(float q[4] ,float h)
{
	q[0]=0;
	q[1]=1*h/2;
	q[2]=1*h/2;
	q[3]=0;
}


void Q_matrix(float Q[4],float h,float nodes[][3],int elements[][5],int i)
{
	float xi[4],eta[4];
	xi[0]=-0.5773;
	xi[1]=0.5773;
	xi[2]=0.5773;
	xi[3]=-0.5773;

	eta[0]=-0.5773;
	eta[1]=-0.5773;
	eta[2]=0.5773;
	eta[3]=0.5773;

	
	for (int jj = 0; jj < 4; jj++)
	{
		Q[0]+=(h*h/32.0)*(sin(PI*(h/2.0*xi[jj]+nodes[elements[i][1]][1]+h/2.0))+sin(PI*(h/2.0*eta[jj]+nodes[elements[i][1]][2]+h/2.0)))*(1-xi[jj])*(1-eta[jj]);
		Q[1]+=(h*h/32.0)*(sin(PI*(h/2.0*xi[jj]+nodes[elements[i][1]][1]+h/2.0))+sin(PI*(h/2.0*eta[jj]+nodes[elements[i][1]][2]+h/2.0)))*(1+xi[jj])*(1-eta[jj]);
		Q[2]+=(h*h/32.0)*(sin(PI*(h/2.0*xi[jj]+nodes[elements[i][1]][1]+h/2.0))+sin(PI*(h/2.0*eta[jj]+nodes[elements[i][1]][2]+h/2.0)))*(1+xi[jj])*(1+eta[jj]);
		Q[3]+=(h*h/32.0)*(sin(PI*(h/2.0*xi[jj]+nodes[elements[i][1]][1]+h/2.0))+sin(PI*(h/2.0*eta[jj]+nodes[elements[i][1]][2]+h/2.0)))*(1-xi[jj])*(1+eta[jj]);
	}
}


void B_matrix(float B[][4], float h)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			B[i][j]*=h;
		}
		
	}
	
}


void A_matrix(float A[][4], float B[][4], float A1[][4])
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			A[i][j]=B[i][j]+A1[i][j];
		}
		
	}
	
}
/*void assemble_G_A()
{
	
}
*/
