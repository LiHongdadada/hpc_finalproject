#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(int argc,char* argv[]){
	float h;	
	if(argc==2)
	{
		h=atof(argv[1]);
	}	
	int i,j;//,k;
	int length=1;
	int n=length/h;
	int num_of_elements = n*n;
	int num_of_nodes=(n+1)*(n+1);
	//nodes information
	float nodes[num_of_nodes][3];
	for (i=0;i<num_of_nodes;i++)
	{
		int a=i/(n+1);
		nodes[i][0]=i;
		nodes[i][1]=(i-a*(n+1))*h;
		nodes[i][2]=a*h;
	}
	//elements information
	int elements[num_of_elements][5];
	for (i=0;i<num_of_elements;i++)
	{
		int a=i/(n+1);
		elements[i][0]=i;
		elements[i][1]=i+a;
		elements[i][2]=i+a+1;
		elements[i][3]=i+a+n+1;
		elements[i][4]=i+a+n+2;
		
	}
	printf("nodes[%d][%d] is : \n",n+1,n+1);
	for(i = 0;i<num_of_nodes;i++)
	{
		for (j=0;j<3;j++)
		{
			printf("%f\t",nodes[i][j]);
		}
		printf("\n");
	}
	printf("elements[%d][%d] is : \n",n,n);
	for(i = 0;i<num_of_elements;i++)
	{
		for (j=0;j<5;j++)
		{
			printf("%d\t",elements[i][j]);
		}
		printf("\n");
	}
	return 0;
}

