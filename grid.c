#include<stdio.h>
#include<stdlib.h>

int main(int argc,char* argv[]){
	float h;	
	if(argc==2)
	{
		h=atof(argv[1]);
	}	
	int i,j;//,k;
	int length=1;
	int number=length/h;
	float x[number+1][number+1];
	float y[number+1][number+1];
	for(i=0;i<number+1;i++)
	{
		for(j=0;j<number+1;j++)
		{
			x[i][j]=0;			
			x[i][j]=j*h;
			y[i][j]=0;
			y[i][j]=i*h;
		}
	}
	printf("x[%d][%d] is : \n",number+1,number+1);
	for(i = 0;i<number+1;i++)
	{
		for (j=0;j<number+1;j++)
		{
			printf("%f\t",x[i][j]);
		}
		printf("\n");
	}
	printf("y[%d][%d] is : \n",number+1,number+1);
	for(i = 0;i<number+1;i++)
	{
		for (j=0;j<number+1;j++)
		{
			printf("%f\t",y[i][j]);
		}
		printf("\n");
	}
	return 0;
}

