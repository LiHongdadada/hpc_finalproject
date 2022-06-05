static char help[] = "Solves a 10000x10000 linear system.\n\n";

#include <petscksp.h>
#include <petscmath.h>
#include <petscsys.h>
#include <petscviewerhdf5.h>

#include<math.h>

#define PI 3.14159265

int main(int argc, char **args)
{
    	PetscInitialize(&argc, &args, (char *)0, help);
    	Mat G_B, G_A, B, A;
    	Vec x, T, temp_vec, tem, G_q;
    	KSP ksp;
    	PC  pc;
    	MPI_Comm comm;
	PetscInt maxit;
    	PetscReal h, dt, t=0, val, temp_B[4][4] = {{0.1111111111, 0.0555555556, 0.0277777778, 0.0555555556}, {0.0555555556, 0.1111111111, 0.0555555556, 0.0277777778}, {0.0277777778, 0.0555555556, 0.1111111111, 0.0555555556}, {0.0555555556, 0.0277777778, 0.0555555556, 0.1111111111}};
    	PetscReal temp_A1[4][4] = {{0.6666666667, -0.1666666667, -0.3333333333, -0.1666666667}, {-0.1666666667, 0.6666666667, -0.1666666667, -0.3333333333}, {-0.3333333333, -0.1666666667, 0.6666666667, -0.1666666667}, {-0.1666666667, -0.3333333333, -0.1666666667, 0.6666666667}};
    	comm = MPI_COMM_WORLD;
    	PetscOptionsGetReal(NULL, NULL, "-h", &h, NULL);
    	PetscOptionsGetReal(NULL, NULL, "-dt", &dt, NULL);
	    PetscOptionsGetInt(NULL, NULL, "-maxit", &maxit, NULL);
    	PetscInt n = 1 / h;
    	PetscInt num_of_nodes = (n + 1) * (n + 1), num_of_elements = n * n;
    	PetscInt iter=0;
	    PetscScalar temp_val;
        PetscInt    index;
        PetscReal   t_v;
        PetscScalar data[2];
        PetscViewer    h5;    /*创建输出*/

    	PetscInt elements[num_of_elements][5];
  	for (int i = 0; i < num_of_elements; i++)
  	{
  	    int a = i / n ;
  	    elements[i][0] = i;
  	    elements[i][1] = i + a;
  	    elements[i][2] = i + a + 1;
  	    elements[i][3] = i + a + n + 2;
  	    elements[i][4] = i + a + n + 1;
  	}
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
/*	PetscPrintf(PETSC_COMM_WORLD,"B\n");
    	MatView(B, PETSC_VIEWER_STDOUT_WORLD);
*/
    	MatCreateAIJ(comm, PETSC_DECIDE, PETSC_DECIDE, 4, 4, 4, PETSC_NULL, 4, PETSC_NULL, &A);
    	MatSetUp(A);
    	MatSetFromOptions(A);
    	for (int i = 0; i < 4; i++)
    	{
    	    for (int j = 0; j < 4; j++)
    	    {
	
    	        MatGetValue(B, i, j, &temp_val);
    	        val = temp_val- temp_A1[i][j];
    	        MatSetValues(A, 1, &i, 1, &j, &val, INSERT_VALUES);
    	    }
    	}
    	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
/*        PetscPrintf(PETSC_COMM_WORLD,"A\n");
    	MatView(A, PETSC_VIEWER_STDOUT_WORLD);
*/
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
 /*       PetscPrintf(PETSC_COMM_WORLD,"before G_A\n");

    	MatView(G_A, PETSC_VIEWER_STDOUT_WORLD);
*/
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
  /*      PetscPrintf(PETSC_COMM_WORLD,"before G_B\n");

        MatView(G_B, PETSC_VIEWER_STDOUT_WORLD);
*/
        /******** 置一划零法 **********/
        for (int i = 0; i < n+1; i++)
        {
                MatZeroRows(G_B, 1, &i, 1, 0, 0);
                MatZeroRows(G_A, 1, &i, 1, 0, 0);
        }
        for (int i = n+1; i < n * (n + 1) + 1; i += n + 1)
        {
                MatZeroRows(G_B, 1, &i, 1, 0, 0);
                MatZeroRows(G_A, 1, &i, 1, 0, 0);
        }
        for (int i = n * (n + 1) + 1; i < num_of_nodes; i ++)
        {
                MatZeroRows(G_B, 1, &i, 1, 0, 0);
                MatZeroRows(G_A, 1, &i, 1, 0, 0);
        }
    	MatAssemblyBegin(G_B, MAT_FINAL_ASSEMBLY);
    	MatAssemblyEnd(G_B, MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(G_A, MAT_FINAL_ASSEMBLY);
    	MatAssemblyEnd(G_A, MAT_FINAL_ASSEMBLY);
   /*     PetscPrintf(PETSC_COMM_WORLD,"AFTER G_B\n");

        MatView(G_B, PETSC_VIEWER_STDOUT_WORLD);
*/
	VecCreate(PETSC_COMM_WORLD, &G_q);
	VecSetSizes(G_q, PETSC_DECIDE, num_of_nodes);
    VecSetFromOptions(G_q);
	VecSet(G_q,0);
	for (int i = 2*n+1; i < num_of_nodes-1; i += n + 1)
	{
	        val = h;
	        VecSetValues(G_q, 1, &i, &val, INSERT_VALUES);
	}
	VecAssemblyBegin(G_q);
   	 VecAssemblyEnd(G_q);
//	 PetscPrintf(PETSC_COMM_WORLD,"after G_q\n");
//	VecView(G_q ,PETSC_VIEWER_STDOUT_WORLD);
	

	VecCreate(PETSC_COMM_WORLD, &temp_vec);
	VecSetSizes(temp_vec, PETSC_DECIDE, num_of_nodes);
    VecSetFromOptions(temp_vec);
	VecSet(temp_vec,0);
	VecAssemblyBegin(temp_vec);
    VecAssemblyEnd(temp_vec);
        
        VecCreate(PETSC_COMM_WORLD, &x);
	VecSetSizes(x, PETSC_DECIDE, num_of_nodes);
    	VecSetFromOptions(x);
	VecSet(x,0);
	VecAssemblyBegin(x);
    VecAssemblyEnd(x);
      
    VecCreate(PETSC_COMM_WORLD, &T);
	VecSetSizes(T, PETSC_DECIDE, num_of_nodes);
    VecSetFromOptions(T);
	VecSet(T, 0);
            for (int i = 0; i < n+1; i++)
        {
                val = 1;
                VecSetValues(T, 1, &i, &val, INSERT_VALUES);
        }
        for (int i = n+1; i < n * (n + 1) + 1; i += n + 1)
        {
                val = 1;
                VecSetValues(T, 1, &i, &val, INSERT_VALUES);
        }
        for (int i = n * (n + 1) + 1; i < num_of_nodes; i ++)
        {
                val = 1;
                VecSetValues(T, 1, &i, &val, INSERT_VALUES);
        }

	    VecAssemblyBegin(T);
        VecAssemblyEnd(T);
       
	VecCreate(PETSC_COMM_WORLD, &tem);
	VecSetSizes(tem, PETSC_DECIDE, 2);
    VecSetFromOptions(tem); 

    /******** ksp *********/
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, G_B, G_B);
	KSPGetPC(ksp, &pc);
	PCSetType(pc, PCJACOBI);
	KSPSetTolerances(ksp, 1.e-10, 1.e-50, PETSC_DEFAULT, PETSC_DEFAULT);
	KSPSetFromOptions(ksp);

	while(PetscAbsReal(t)<maxit*dt){

		t += dt;
		MatMult(G_A, T, temp_vec);
		VecAXPY(temp_vec, 1, G_q);
	//	PetscPrintf(PETSC_COMM_WORLD,"add G_q\n");
	//        VecView(G_q ,PETSC_VIEWER_STDOUT_WORLD);

	    	KSPSolve(ksp, temp_vec, x); /* 求解出下一时间 */
		VecCopy(x, T); /* 求解后的结果传到T作为初始值 */
	        for (int i = 0; i < n+1; i++)
			{
        	        val = 1;
       	        	VecSetValues(T, 1, &i, &val, INSERT_VALUES);
      	 	}
     		for (int i = n+1; i < n * (n + 1) + 1; i += n + 1)
       		{
               		val = 1;
             		VecSetValues(T, 1, &i, &val, INSERT_VALUES);
	        }
	        for (int i = n * (n + 1) + 1; i < num_of_nodes; i ++)
	        {
	                val = 1;
	                VecSetValues(T, 1, &i, &val, INSERT_VALUES);
	        }

	        VecAssemblyBegin(T);
      		VecAssemblyEnd(T);
        iter +=1;
        if ((iter%10)==0)
        {
            data[0] = dt; data[1] = t;
            VecSet(tem,0);
            for (index=0;index<2;index++)
            {
                t_v=data[index];
                VecSetValues(tem,1,&index,&t_v,INSERT_VALUES);
            }
            VecAssemblyBegin(tem);
            VecAssemblyEnd(tem);

            PetscViewerCreate(PETSC_COMM_WORLD,&h5);
            PetscViewerHDF5Open(PETSC_COMM_WORLD,"explicit.h5", FILE_MODE_WRITE, &h5);
            PetscObjectSetName((PetscObject) T, "explicit-vector");
            PetscObjectSetName((PetscObject) tem, "explicit-necess-data");
            VecView(tem, h5);
            VecView(T, h5);
            PetscViewerDestroy(&h5);
        }

	}
	PetscPrintf(PETSC_COMM_WORLD,"solution\n");

	VecView(T ,PETSC_VIEWER_STDOUT_WORLD);
	
    	MatDestroy(&B);
    	MatDestroy(&G_A);
    	MatDestroy(&G_B);
    	MatDestroy(&A);
    	VecDestroy(&temp_vec);
    	VecDestroy(&T);
    	VecDestroy(&x);
        VecDestroy(&tem);
	VecDestroy(&G_q);
    	PetscFinalize();
    	return 0;
}

