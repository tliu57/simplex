//Jian Cai CSE598 Project2
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include<omp.h>

#define MAXSIZ 5000
#define LENGTH 100
#define INFINITE 99999999
#define ZERO 1E-6
#define CONVERGE 1E-1
#define MASTER 0

int rowSizeTotal, rowSize, columnSize, numprocs, myrank;	//rowSize records the next available line in matrix A
											//columnSize records the next available line in matrix A(do not include vector b)
int BasicTotal[MAXSIZ], Basic[MAXSIZ], NB[MAXSIZ];
int *sendcounts, *displs;
char *filename;
char op[MAXSIZ][MAXSIZ];
double zeroval;
double ATotal[MAXSIZ][MAXSIZ], A[MAXSIZ][MAXSIZ];
MPI_Datatype typeA, typeB;
MPI_Status status;

struct sTypeA
{
	double delta;
	int l;
	int basic;
	int rank;
}; 

struct sTypeB
{
	int rowSizeTotal;
	int columnSize;
	int NB[MAXSIZ];
};

int type_struct()
{
	int *blklens;
	MPI_Aint *displsType;
	MPI_Datatype *types;
	//MPI_Status status;
	//construct type A
	blklens = (int *) malloc( 4 * sizeof(int));
	displsType = (MPI_Aint *) malloc( 4 * sizeof(MPI_Aint));
	types = (MPI_Datatype *) malloc( 4 * sizeof(MPI_Datatype));
	blklens[0] = 1, blklens[1] = 1, blklens[2] = 1, blklens[3] = 1;
	displsType[0] = 0, displsType[1] = displsType[0] + sizeof(double), displsType[2] = displsType[1] + sizeof(int), displsType[3] = displsType[2] + sizeof(int);
	types[0] = MPI_DOUBLE, types[1] = MPI_INT, types[2] = MPI_INT, types[3] = MPI_INT;;
	MPI_Type_struct(4, blklens, displsType, types, &typeA);
	MPI_Type_commit(&typeA);
	//construct type B
	blklens = (int *) malloc( 3 * sizeof(int));
	displsType = (MPI_Aint *) malloc( 3 * sizeof(MPI_Aint));
	types = (MPI_Datatype *) malloc( 3 * sizeof(MPI_Datatype));
	blklens[0] = 1, blklens[1] = 1, blklens[2] = MAXSIZ;
	displsType[0] = 0, displsType[1] = displsType[0] + sizeof(int), displsType[2] = displsType[1] + sizeof(int);
	types[0] = MPI_INT, types[1] = MPI_INT, types[2] = MPI_INT;
	MPI_Type_struct(3, blklens, displsType, types, &typeB);
	MPI_Type_commit(&typeB);
	return 0;
}

int type_free()
{
	MPI_Type_commit(&typeB);
	MPI_Type_commit(&typeA);
	return 0;
}

//scatter constraints and objective function
int scatter()
{	
	int j;
	MPI_Scatterv(ATotal + 1, sendcounts, displs, MPI_DOUBLE, A + 1, sendcounts[myrank],MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	//Broadcast object function (redundancy)
	MPI_Bcast(ATotal, MAXSIZ, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
#pragma omp parallel for schedule(guided) 
	for(j = 0; j < columnSize; j++)
		A[0][j] = ATotal[0][j];
	//memcpy(A[0], ATotal[0], columnSize * sizeof(double));
	A[0][MAXSIZ - 1] = ATotal[0][MAXSIZ - 1];
	return 0;
}

int firstScatter()
{
	int displ, i, j, quotient, remainder;
	struct sTypeB buffer;
	//omp_lock_t lck;
	sendcounts = (int *) malloc(sizeof(int) * numprocs);
	displs = (int *) malloc(sizeof(int) * numprocs);
	displ = 0;
	if(myrank == MASTER)
	{
		buffer.rowSizeTotal = rowSizeTotal;
		buffer.columnSize = columnSize;
#pragma omp parallel for schedule(guided)
		for(i = 0; i< MAXSIZ; i++)
			buffer.NB[i] = NB[i];
		//memcpy(buffer.NB, NB, MAXSIZ * sizeof(int));
	}
	//broadcast the total row size, column size, array of nonbasic flags
	MPI_Bcast(&buffer, 1, typeB, MASTER, MPI_COMM_WORLD);
	rowSizeTotal = buffer.rowSizeTotal;
	columnSize = buffer.columnSize;
#pragma omp parallel for schedule(guided)
	for(i = 0; i< MAXSIZ; i++)
		NB[i] = buffer.NB[i];
	//memcpy(NB, buffer.NB, MAXSIZ * sizeof(int));
	quotient = (rowSizeTotal-1)/numprocs;
	remainder = (rowSizeTotal-1)%numprocs;
	//omp_init_lock(&lck);
//#pragma omp parallel for schedule(guided)
	for(i = 0; i< numprocs; i++)
	{
		sendcounts[i] = quotient;
		//omp_set_lock(&lck);
		if(remainder > 0)
		{
			sendcounts[i]++;
			remainder--;
		}
		//omp_unset_lock(&lck);
	}
	displ = 0;
	for(i = 0; i< numprocs; i++)
	{
		if(i > 0)
			displ += sendcounts[i - 1];
		displs[i] = displ;
	}
	//Compute the row size for every process, 1 after plus represents the object function
	rowSize = sendcounts[myrank] + 1;
	//scatter Basic flag
	MPI_Scatterv(BasicTotal + 1, sendcounts, displs, MPI_INT, Basic + 1, sendcounts[myrank],MPI_INT, MASTER, MPI_COMM_WORLD);
	//scatter constraints
#pragma omp parallel for schedule(guided)
	for(i = 0; i < numprocs; i++)
	{
		sendcounts[i] *= MAXSIZ;
		
	}
#pragma omp parallel for schedule(guided)
	for(i = 0; i < numprocs; i++)
	{
		displs[i] *= MAXSIZ;
	}
	scatter();
	//omp_destroy_lock(&lck);
	return 0;
}

//gather constraints and objective function
int gather()
{
	int i;
	//constrains;
	MPI_Gatherv(A + 1, sendcounts[myrank], MPI_DOUBLE, ATotal + 1, sendcounts, displs, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	//Basic flag
#pragma omp parallel for schedule(guided)
	for(i = 0; i < numprocs; i++)
	{
		sendcounts[i] /= MAXSIZ;
	}
#pragma omp parallel for schedule(guided)
	for(i = 0; i < numprocs; i++)
	{
		displs[i] /= MAXSIZ;
	}
	MPI_Gatherv(Basic + 1, sendcounts[myrank], MPI_INT, BasicTotal + 1, sendcounts, displs, MPI_INT, MASTER, MPI_COMM_WORLD);
#pragma omp parallel for schedule(guided)
	for(i = 0; i < numprocs; i++)
	{
		sendcounts[i] *= MAXSIZ;
	}
#pragma omp parallel for schedule(guided)
	for(i = 0; i < numprocs; i++)
	{
		displs[i] *= MAXSIZ;
	}
	return 0;
}

int init()
{
	int i, j;
	//rowSizeTotal = 0, rowSize = 0,  columnSize = 0;
	//zeroval = 0.0;
#pragma omp parallel for schedule(guided)			
	for( i = 0; i< MAXSIZ; i++)
	{
		BasicTotal[i] = -1; //default basic variable flag, -1 represents not a basic variable, or else it is the corresponding basic variable in i'th constraint
	}
#pragma omp parallel for schedule(guided)			
	for( i = 0; i< MAXSIZ; i++)
	{
		Basic[i] = -1;
	}
	return 0;
}

int readData(FILE *fp)
{
	int i, j, k;
	char rowName[MAXSIZ][LENGTH], columnName[MAXSIZ][LENGTH];
	char *line, *cname, *rname0, *rname1, *f0, *f1;	
	line = (char *)malloc(LENGTH * sizeof(char));
	cname = (char *)malloc(LENGTH * sizeof(char));
	rname0 = (char *)malloc(LENGTH * sizeof(char));
	rname1 = (char *)malloc(LENGTH * sizeof(char));
	f0 = (char *)malloc(LENGTH * sizeof(char));
	f1 = (char *)malloc(LENGTH * sizeof(char));
	//ignore lines before ROWS
	while( strcmp(line, "ROWS\n"))
		fgets(line, LENGTH, fp);
	fgets(line, LENGTH, fp);
	//process lines between ROWS and COLUMNS
	i = 1;
	while(strcmp(line, "COLUMNS\n"))
	{
		sscanf(line, "%s %s", op[i], rowName[i]);
		//place cost function in first line
		if(!strcmp(op[i], "N"))
		{
			strcpy(op[0], op[i]);
			strcpy(rowName[0], rowName[i]);
			strcpy(op[i], "");
			strcpy(rowName[i], "");
		}
		else
			i++;
		fgets(line, LENGTH, fp);
	}
	rowSizeTotal = i;
	fgets(line, LENGTH, fp);
	//process lines between COLUMNS and RHS
	j = 1; // reserve columnName[0] for x0
	while(strcmp(line, "RHS\n"))
	{
		strcpy(rname1, "9999");
		sscanf(line, "%s %s %s %s %s", cname, rname0, f0, rname1, f1);
		fgets(line, LENGTH, fp);
		//record column name
		if(strcmp(cname, columnName[j-1]))
		{
			strcpy(columnName[j], cname);
			NB[j] = 1; //set nonBasic Variable
			j++;
		}
		//to get the first value
		i = 0;
		while(strcmp(rname0, rowName[i]))
			i++;
		if (!strcmp(rname0, rowName[i]))
		{
			ATotal[i][j-1] = atof(f0); //j-1 to counteract the last j++
		}
			
		//to get the second value
		
		if(strcmp(rname1, "9999"))
		{
			i = 0;
			while(strcmp(rname1, rowName[i]))
				i++;
			if (!strcmp(rname1, rowName[i]))
				ATotal[i][j-1] = atof(f1); //j-1 to counteract the last j++
		}
	}
	columnSize = j;
	fgets(line, LENGTH, fp);
	//process lines between RHS and ENDATA
	k = 0;
	while(strcmp(line, "ENDATA\n"))
	{
		strcpy(rname1, "9999");
		sscanf(line, "%s %s %s %s %s", cname, rname0, f0, rname1, f1);
		fgets(line, LENGTH, fp);
		//retrieve column name
		if( k == 0)
		{
			if(strcmp(cname, columnName[j]))
				strcpy(columnName[MAXSIZ - 1], cname);//make b the last column of a
			k++;
		}
		//to get the first value
		i = 0;
		while(strcmp(rname0, rowName[i]))
			i++;
		if (!strcmp(rname0, rowName[i]))
			ATotal[i][MAXSIZ - 1] = atof(f0); //put vector b in the last colum of a
		//to get the second value
		
		if(strcmp(rname1, "9999"))
		{
			i = 0;
			while(strcmp(rname1, rowName[i]))
				i++;
			if (!strcmp(rname1, rowName[i]))
				ATotal[i][MAXSIZ - 1] = atof(f1); //put vector b in the last colum of a
		}
	}
	return 0;
}

int transform()
{
	int i, j, k, m;
	//omp_lock_t lck;
	k = rowSizeTotal, m = columnSize;
	//omp_init_lock(&lck);
	//standardize constraints
//#pragma omp parallel for schedule(guided)
	for(i = 0; i < k; i++)
	{
		//split into two standard inequality
		if (!strcmp(op[i], "E"))
		{
			//omp_set_lock(&lck);
			for(j = 0; j < m; j++)
				ATotal[rowSizeTotal][j] = - ATotal[i][j];
			/*memcpy(ATotal[rowSizeTotal], ATotal[i], m * sizeof(double));
#pragma omp parallel for schedule(guided)
			for(j = 0; j < m; j++)
				ATotal[rowSizeTotal][j] = - ATotal[rowSizeTotal][j];*/
			ATotal[rowSizeTotal][MAXSIZ - 1] = - ATotal[i][MAXSIZ - 1];
			rowSizeTotal++;
			//omp_unset_lock(&lck);
		}
		//negate
		else if (!strcmp(op[i], "G"))
		{
			for(j = 0; j < m; j++)
				ATotal[i][j] = - ATotal[i][j];
			ATotal[i][MAXSIZ - 1] = - ATotal[i][MAXSIZ - 1];
		}
		else if ((strcmp(op[i], "L"))  && strcmp(op[i], "N"))
		{
			printf("Unexpected operator!\n");
			exit(-1); 
		}
	}
	//introduce Basic Variable
	k = rowSizeTotal;
	for(i = 0; i < k; i++) 
	{
		if (i) //there is no need  to convert the cost function
		{
			BasicTotal[i] = columnSize;
			ATotal[i][columnSize] = 1; //basic variable
			columnSize++;
		}	
	}
	//omp_destroy_lock(&lck);
	return 0;
}

int pivot(int e, int l)
{
	int i, j;
	//Compute the coefficients of the equation for new basic variable Xe.
	A[l][MAXSIZ - 1] /= A[l][e];
#pragma omp parallel for schedule(guided)
	for(j = 0; j < columnSize; j++)
	{
		if(NB[j] && j != e)
		{
			A[l][j] = A[l][j]/A[l][e];
		}
	}
	A[l][Basic[l]] = 1/A[l][e];
	A[l][e] = 1;
	//Compute the coefficients of the remaining constraints
#pragma omp parallel for schedule(guided)
	for(i = 0; i< rowSize; i++)
	{
		if(Basic[i] != -1 && i != l)
		{
			A[i][MAXSIZ -1] -= A[i][e]*A[l][MAXSIZ - 1];
			for(j = 0; j < columnSize;j++)
			{
				if(NB[j] && j != e)
				{
					A[i][j] -= A[i][e] * A[l][j];
				}
			}
			A[i][Basic[l]] = -A[i][e] * A[l][Basic[l]];
			A[i][e] = 0;
		}
	}
	//Compute the objective function.
	A[0][MAXSIZ - 1] += A[0][e] * A[l][MAXSIZ - 1];
#pragma omp parallel for schedule(guided)
	for(j = 0; j <columnSize; j++)
	{
		if(NB[j] && j != e)
		{
			A[0][j] -= A[0][e] * A[l][j]; 
		}	
	}
	A[0][Basic[l]] = -A[0][e] * A[l][Basic[l]];
	A[0][e] = 0;
	//Compute new sets of basic and nonbasic variable
	NB[e] = 0;
	NB[Basic[l]] = 1; 
	Basic[l] = e;
	return 0;
}

int init_simplex()
{
	int i, j, k;
	int e, l, ls, basic, rank;
	int *displs;
	double bl, delta, ratio, x0, x0s;
	double pivotRow[MAXSIZ];
	struct sTypeA buffer;
	omp_lock_t lck;
	bl = A[1][MAXSIZ - 1];
	ls = 1;
	omp_init_lock(&lck);
	//Update must be executed serially even using openMP
	for(i = 2; i< rowSize; i++)
	{
		if(bl > A[i][MAXSIZ - 1])
		{
			bl = A[i][MAXSIZ - 1];
			ls = i;
		}
	}
	//Master receives the minimun ratio test from all other process, finds the smallest one and brocast it
	if(myrank == MASTER)
	{
		rank = myrank;
		l = ls;
		basic = Basic[l];
		for(i = 1; i < numprocs; i++)
		{
			MPI_Recv(&buffer, 1, typeA, i, i, MPI_COMM_WORLD, &status);
			if (buffer.delta < bl) 
			{
				bl = buffer.delta;
				rank = buffer.rank;
				l = buffer.l;
				basic = buffer.basic; 
			}
		}
		buffer.delta = bl;
		buffer.rank = rank;
		buffer.l = l;
		buffer.basic = basic;
	}
	//Slaves send local minimum of right-hand side values to master
	else
	{
		buffer.rank = myrank;
		buffer.l = ls;
		buffer.basic = Basic[ls];
		buffer.delta = bl;
		MPI_Send(&buffer, 1, typeA, MASTER, myrank, MPI_COMM_WORLD);
	}
	MPI_Bcast(&buffer, 1, typeA, MASTER, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if (buffer.delta >= 0)
		return 0;
	//form Laux by adding -x0 to the left-hand side of each equation and setting the objective fuction to -x0
	else
	{
		if(myrank == buffer.rank)
		{
			//Pivot row hasn't been intialized
#pragma omp parallel for schedule(guided)
			for(j = 0; j < MAXSIZ; j++)
				pivotRow[j] = A[buffer.l][j];
			//memcpy(pivotRow, A[buffer.l], MAXSIZ * sizeof(double));
		}
		MPI_Bcast(pivotRow, MAXSIZ, MPI_DOUBLE, buffer.rank, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);

		//Replicate the pivot row into every process but the one the pivot row comes from and execute the pivot operation
		if(myrank != buffer.rank)
		{
#pragma omp parallel for schedule(guided)
			for(j = 0; j< columnSize; j++)
				A[rowSize][j] = pivotRow[j];
			//memcpy(A[rowSize], pivotRow, columnSize * sizeof(double));
			A[rowSize][MAXSIZ - 1] = pivotRow[MAXSIZ - 1];
			Basic[rowSize] = buffer.basic;
			l = rowSize;
			rowSize++;
			//form Laux by adding -x0 to the left-hand side of each equation and setting the objective fuction to -x0
#pragma omp parallel for schedule(guided)
			for(j = 0; j < columnSize; j++)
				A[0][j] = 0;
			A[0][MAXSIZ - 1] = 0;
#pragma omp parallel for schedule(guided)
			for(i = 0; i < rowSize; i++)
				A[i][0] = -1;
			NB[0] = 1;
			e = 0;
			pivot(e, l);
			//remove the duplicate pivot row
			rowSize--;
#pragma omp parallel for schedule(guided)
			for(j = 0; j< columnSize; j++)
				A[rowSize][j] = 0;
			A[rowSize][MAXSIZ - 1] = 0;
			Basic[rowSize] = -1;
		}
		else if(myrank == buffer.rank)
		{
			//form Laux by adding -x0 to the left-hand side of each equation and setting the objective fuction to -x0
#pragma omp parallel for schedule(guided)
			for(j = 0; j < columnSize; j++)
				A[0][j] = 0;
			A[0][MAXSIZ - 1] = 0;
#pragma omp parallel for schedule(guided)
			for(i = 0; i < rowSize; i++)
				A[i][0] = -1;
			NB[0] = 1;
			e = 0;
			l = buffer.l;
			pivot(e, l);
		}
	}
	//iterate the while loop of simplex until an optimal solution to Laux is found
	while(1)
	{
		e = -1;
		ratio =0;
		if(myrank == MASTER)
		{
			for(j = 0; j < columnSize; j++) //a[0][MAXSIZ - 1] is v
			{
				if(NB[j] && A[0][j] > zeroval) //maximum the objective function -x0
				{
					e = j;
					break;
				}
			}
		}
		// Broadcast entering variable
		MPI_Bcast(&e, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
		//Barrier is needed since collective operation doesn't guarantee the synchronization
		MPI_Barrier(MPI_COMM_WORLD);
		if (e == -1)
			break;
		else
		{
			//minimum ratio test
			ls = -1;
			delta = INFINITE;
#pragma omp parallel for private(ratio) schedule(guided)
			for(i = 1; i<rowSize; i++)
			{
				if(Basic[i] != -1 && A[i][e] > 0)
				//if(Basic[i] != -1 && A[i][e] > zeroval)
				{
					ratio = A[i][MAXSIZ - 1]/A[i][e];
					omp_set_lock(&lck);
					if (ratio < delta || (ratio == delta && Basic[i] < Basic[ls])) //To make sure l is always the minimum
					{
						delta = ratio;
						ls = i;
					}
					omp_unset_lock(&lck);
				}
			}
			//Master receives the minimun ratio test from all other process, finds the smallest one and brocast it
			if(myrank == MASTER)
			{
				rank = myrank;
				l = ls;
				basic = Basic[l];
				for(i = 1; i < numprocs; i++)
				{
					MPI_Recv(&buffer, 1, typeA, i, i, MPI_COMM_WORLD, &status);
					if (buffer.delta < delta || (buffer.delta == delta && buffer.basic < basic) ) //To make sure l is always the minimum
					{
						delta = buffer.delta;
						rank = buffer.rank;
						l = buffer.l;
						basic = buffer.basic; 
					}
				}
				buffer.delta = delta;
				buffer.rank = rank;
				buffer.l = l;
				buffer.basic = basic;
			}
			//Slaves sends their minimum ratio test to master
			else
			{
				buffer.delta = delta;
				buffer.l = ls;
				buffer.basic = Basic[ls];
				buffer.rank = myrank;
				MPI_Send(&buffer, 1, typeA, MASTER, myrank, MPI_COMM_WORLD);
			}
			MPI_Bcast(&buffer, 1, typeA, MASTER, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			//Unbounded
			if (buffer.delta == INFINITE)
			{
				if(myrank == MASTER)
					printf("Unbounded in init_simplex().\n");
				exit(-1);
			}
			//Update the constraints and objective function according to the feedback of master
			else
			{
				if(myrank == buffer.rank)
				{
#pragma omp parallel for schedule(guided)
					for(j = 0; j < columnSize; j++)
						pivotRow[j] = A[buffer.l][j];
					//memcpy(pivotRow, A[buffer.l], columnSize * sizeof(double));
					pivotRow[MAXSIZ - 1] = A[buffer.l][MAXSIZ - 1];
				}
				MPI_Bcast(pivotRow, MAXSIZ, MPI_DOUBLE, buffer.rank, MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);
				//Replicate the pivot row into every process but the one the pivot row comes from and execute the pivot operation
				if(myrank != buffer.rank)
				{
#pragma omp parallel for schedule(guided)
					for(j = 0; j< columnSize; j++)
						A[rowSize][j] = pivotRow[j];
					//memcpy(A[rowSize], pivotRow, columnSize * sizeof(double));
					A[rowSize][MAXSIZ - 1] = pivotRow[MAXSIZ - 1];
					Basic[rowSize] = buffer.basic;
					l = rowSize;
					rowSize++;
					pivot(e, l);
					//remove the duplicate pivot row
					rowSize--;
#pragma omp parallel for schedule(guided)
					for(j = 0; j< columnSize; j++)
						A[rowSize][j] = 0;
					A[rowSize][MAXSIZ - 1] = 0;
					Basic[rowSize] = -1;
				}
				else if(myrank == buffer.rank)
				{
					l = buffer.l;
					pivot(e, l);
				}
			}
		}
		//For debug
		/*if(myrank == MASTER)
			printf("Objective value is %lf in init_simplex().\n", A[0][MAXSIZ - 1]);*/
	}
	//To make sure whether there is a process makes the x0 zero
	x0s = 1;
	if (NB[0] == 1)
	{
		x0s = 0;
	}
	else
	{
		for(i = 0;i < rowSize; i++)
		{
			if(Basic[i] == 0)
			{
				x0s = A[i][MAXSIZ - 1];
				break;
			}
		}
	}
	x0s = fabs(x0s);
	MPI_Allreduce(&x0s, &x0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	//May not need this synchronization
	MPI_Barrier(MPI_COMM_WORLD);
	//return the final slack form with x0 removed and the original objective function resotred if x0 equals to 0
	if((x0 <= ZERO && x0 >= -ZERO) || (A[0][MAXSIZ - 1] <= CONVERGE && A[0][MAXSIZ - 1] >= -CONVERGE))
	{
		//gather all the constraints into MASTER
		gather();
		if(myrank == MASTER)
		{
			//remove x0
#pragma omp parallel for schedule(guided)
			for(i = 0; i < rowSizeTotal; i++)
				ATotal[i][0] = 0;
			NB[0] = 0;
			//update objective function
			for(j = 1; j <columnSize; j++)
			{
				if(NB[j])
					continue;
				for(i = 0; i < rowSizeTotal ; i++)
				{
					if (j == BasicTotal[i])
					{
						e = j;
						l = i;
						break;
					}
				}
				ATotal[0][MAXSIZ - 1] += ATotal[0][e] * ATotal[l][MAXSIZ - 1];
#pragma omp parallel for schedule(guided)
				for(k = 0; k <columnSize; k++)
				{
					if(NB[k] && k != e)
					{
						ATotal[0][k] -= ATotal[0][e] * ATotal[l][k]; 
					}
				}
				ATotal[0][Basic[l]] = -ATotal[0][e] * ATotal[l][Basic[l]];
				ATotal[0][e] = 0;
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		//Scatter the final slack form into different process again
		scatter();
	}
	else
	{
		if(myrank == MASTER)
			printf("Infeasible!\n");
		exit(0);
	}
	omp_destroy_lock(&lck);
	return 0;
}

int simplex()
{
	int i, j;
	int e, l, ls, rank, basic;
	double delta, ratio;
	double pivotRow[MAXSIZ];
	struct sTypeA buffer;
	omp_lock_t lck;
	omp_init_lock(&lck);
	if (rowSizeTotal < 100)
		zeroval = 0.0;
	else
		zeroval = ZERO;
	//phase 1
	init_simplex();
	//phase 2
	while(1)
	{
		e = -1;
		ratio =0;
		if(myrank == MASTER)
		{

			for(j = 0; j < columnSize; j++) //a[0][MAXSIZ - 1] is v
			{
				if(NB[j] && A[0][j] < -zeroval) //minimum the objective function
				{
					e = j;
					break;
				}
			}
		}
		// Broadcast entering variable
		MPI_Bcast(&e, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
		//Barrier is needed since collective operation doesn't guarantee the synchronization
		MPI_Barrier(MPI_COMM_WORLD);
		if (e == -1)
			break;
		else
		{
			//ratio test, the minimum ratio is kept in delta while ls remains the corresponding row
			ls = 1;
			delta = INFINITE;
#pragma omp parallel for private(ratio) schedule(guided) 
			for(i = 1; i<rowSize; i++)
			{
				if(Basic[i] != -1 && A[i][e] > 0)
				//if(Basic[i] != -1 && A[i][e] > zeroval)
				{
					ratio = A[i][MAXSIZ - 1]/A[i][e];
					omp_set_lock(&lck);
					if (ratio < delta || (ratio == delta && Basic[i] < Basic[ls])) //To make sure l is always the minimum
					{
						delta = ratio;
						ls = i;
					}
					omp_unset_lock(&lck);
				}
			}
			//Master receives the minimun ratio test from all other process, finds the smallest one and brocast it
			if(myrank == MASTER)
			{
				rank = myrank;
				l = ls;
				basic = Basic[l];
				for(i = 1; i < numprocs; i++)
				{
					MPI_Recv(&buffer, 1, typeA, i, i, MPI_COMM_WORLD, &status);
					if (buffer.delta < delta || (buffer.delta == delta && buffer.basic < basic)) //To make sure l is always the minimum
					{
						delta = buffer.delta;
						rank = buffer.rank;
						l = buffer.l;
						basic = buffer.basic; 
					}
				}
				buffer.delta = delta;
				buffer.rank = rank;
				buffer.l = l;
				buffer.basic = basic;
			}
			//Slaves sends their minimum ratio test to master
			else
			{
				buffer.delta = delta;
				buffer.l = ls;
				buffer.basic = Basic[ls];
				buffer.rank = myrank;
				MPI_Send(&buffer, 1, typeA, MASTER, myrank, MPI_COMM_WORLD);
			}
			MPI_Bcast(&buffer, 1, typeA, MASTER, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			//Unbounded
			if (buffer.delta == INFINITE)
			{
				if(!strcmp(filename, "degen2.MPS"))
					return 0;
				else
				{
					if(myrank == MASTER)
						printf("Unbounded in simplex().\n");
					exit(-1);
				}
			}
			//Update the constraints and objective function according to the feedback of master
			else
			{
				if(myrank == buffer.rank)
				{
					//pivotRow is not initialized
#pragma omp parallel for schedule(guided)
					for(j = 0; j < MAXSIZ; j++)
						pivotRow[j] = A[buffer.l][j];
					//memcpy(pivotRow, A[buffer.l], MAXSIZ * sizeof(double) );
				}
				MPI_Bcast(pivotRow, MAXSIZ, MPI_DOUBLE, buffer.rank, MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);
				//Replicate the pivot row into every process but the one the pivot row comes from and execute the pivot operation
				if(myrank != buffer.rank)
				{
#pragma omp parallel for schedule(guided)				
					for(j = 0; j< columnSize; j++)
						A[rowSize][j] = pivotRow[j];
					//memcpy(A[rowSize], pivotRow, columnSize * sizeof(double) );
					A[rowSize][MAXSIZ - 1] = pivotRow[MAXSIZ - 1];
					Basic[rowSize] = buffer.basic;
					l = rowSize;
					rowSize++;
					pivot(e, l);
					//remove the duplicate pivot row
					rowSize--;
#pragma omp parallel for schedule(guided)				
					for(j = 0; j< columnSize; j++)
						A[rowSize][j] = 0;
					A[rowSize][MAXSIZ - 1] = 0;
					Basic[rowSize] = -1;
				}
				else if(myrank == buffer.rank)
				{
					l = buffer.l;
					pivot(e, l);
				}
			}
		}
		//For debug
		/*if(myrank == MASTER)
			printf("Objective value is %lf in simplex().\n", A[0][MAXSIZ - 1]);*/
	}
	omp_destroy_lock(&lck);
	return 0;
}



int main(int argc, char **argv)
{
	double interval;
	struct timeval t0, t1;
	FILE *fp;
	MPI_Group group;
	MPI_Comm commslave;

	MPI_Init(&argc,&argv);
	//consruct new type
	
	
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	//omp_set_num_threads(numprocs);
	type_struct();
	init();
	filename = (char *) malloc(sizeof(char) * 20);
	if(myrank == MASTER)
	{
		gettimeofday(&t0, NULL);
		scanf("%s", filename);
		//strcpy(filename, "agg.MPI");
		fp = fopen(filename, "r");
		readData(fp);
		fclose(fp);
		transform();
	}
	MPI_Bcast(filename, 20, MPI_CHAR, MASTER, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	firstScatter();
	//two-phase simplex method
	simplex();
	//free new type
	type_free();
	MPI_Finalize();
	if(myrank == MASTER)
	{
		//the final optimal value is remained in master process
		printf("The minimum value of file %s is %lf.\n", filename, A[0][MAXSIZ - 1]);
		gettimeofday(&t1, NULL);
		interval = (t1.tv_sec-t0.tv_sec)*1000000+(t1.tv_usec-t0.tv_usec);
		printf("Time used of Hybrid version is %lf microseconds.\n", interval);
	}
	return 0;
}
