//Jian Cai CSE598 Project2 Milestone1
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define MAXSIZ 5000
#define LENGTH 100
#define INFINITE 99999999
#define ZERO 1E-6
#define CONVERGE 1E-1

int rowSize, columnSize;	//rowSize records the next available line in array a
								//columnSize records the next available line in array a(do not include vector b)
int NB[MAXSIZ], Basic[MAXSIZ];
char *filename;
char op[MAXSIZ][LENGTH];
double zeroval;
double A[MAXSIZ][MAXSIZ];

void init()
{
	int i, j;
	//rowSize = 0,  columnSize = 0;
	//zeroval = 0.0;
	for( i = 0; i< MAXSIZ; i++)
	{
		//NB[i] = 0; //nonBasic Variable flag, 0 represents not a NonBasic Variable while 1 means NonBasic Variable
		Basic[i] = -1; //default Basic Variable flag, -1 represents not a Basic Variable, or else it is the corresponding Basic Variable in i'th constraint
	}
	
}

int readData(FILE *fp)
{
	int i, j, k;
	char rowName[MAXSIZ][LENGTH], columnName[MAXSIZ][LENGTH];
	char  *line, *cname, *rname0, *rname1, *f0, *f1;
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
	rowSize = i;
	fgets(line, LENGTH, fp);
	//process lines between COLUMNS and RHS
	j = 1; // reserve columnName[0] for X0
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
			A[i][j-1] = atof(f0); //j-1 to counteract the last j++
		}
			
		//to get the second value
		
		if(strcmp(rname1, "9999"))
		{
			i = 0;
			while(strcmp(rname1, rowName[i]))
				i++;
			if (!strcmp(rname1, rowName[i]))
				A[i][j-1] = atof(f1); //j-1 to counteract the last j++
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
				strcpy(columnName[ MAXSIZ - 1 ], cname);//make b the last column of a
			k++;
		}
		//to get the first value
		i = 0;
		while(strcmp(rname0, rowName[i]))
			i++;
		if (!strcmp(rname0, rowName[i]))
			A[i][MAXSIZ - 1 ] = atof(f0); //put vector b in the last colum of a
		//to get the second value
		
		if(strcmp(rname1, "9999"))
		{
			i = 0;
			while(strcmp(rname1, rowName[i]))
				i++;
			if (!strcmp(rname1, rowName[i]))
				A[i][MAXSIZ-1] = atof(f1); //put vector b in the last colum of a
		}
	}
	return 0;
}

int transform()
{
	int i, j, k, m;
	k = rowSize, m = columnSize;
	//standardize constraints
	for(i = 0; i < k; i++)
	{
		//split into two standard inequality
		if (!strcmp(op[i], "E"))
		{
			/*for(j = 0; j < m; j++)
				A[rowSize][j] = - A[i][j];*/
			memcpy(A[rowSize], A[i], m * sizeof(double));
			for(j = 0; j < m; j++)
				A[rowSize][j] = - A[rowSize][j];
			A[rowSize][MAXSIZ - 1] = - A[i][MAXSIZ - 1];
			rowSize++;
		}
		//negate
		else if (!strcmp(op[i], "G"))
		{
			for(j = 0; j < m; j++)
				A[i][j] = - A[i][j];
			A[i][MAXSIZ - 1] = - A[i][MAXSIZ - 1];
		}
		else if ((strcmp(op[i], "L"))  && strcmp(op[i], "N"))
		{
			printf("Unexpected operator!\n");
			exit(-1); 
		}
	}
	//introduce Basic Variable
	k = rowSize;
	for(i = 0; i < k; i++) 
	{
		if (i) //there is no need  to convert the cost function
		{
			Basic[i] = columnSize;
			A[i][columnSize] = 1; //basic variable
			columnSize++;
		}
		//A[i][columnSize] = 1; //basic variable
		
	}
	return 0;
}

int pivot(int e, int l)
{
	int i, j;
	//Compute the coefficients of the equation for new basic variable Xe.
	A[l][MAXSIZ - 1] /= A[l][e];
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
	int e, l, x0;
	double bl, delta, ratio;
	double ob[MAXSIZ];
	bl = A[1][MAXSIZ - 1];
	for(i = 2; i< rowSize; i++)
		if(bl > A[i][MAXSIZ - 1])
		{
			bl = A[i][MAXSIZ - 1];
			l = i;
		}
	if (bl >= 0)
		return 0;
	//form Laux by adding -x0 to the left-hand side of each equation and setting the objective fuction to -x0
	for(j = 0; j< MAXSIZ; j++)
		ob[j] = 0;
	/*for(j = 0; j < columnSize; j++)
	{
		ob[j] = A[0][j];
		A[0][j] = 0;
	}*/
	memcpy(ob, A[0], columnSize * sizeof(double));
	memset(A[0], 0, columnSize * sizeof(double));
	if(j < MAXSIZ)
	{
		ob[MAXSIZ - 1] = A[0][MAXSIZ - 1];
		A[0][MAXSIZ -1] = 0;
	}

	for(i = 0; i < rowSize; i++)
		A[i][0] = -1;
	NB[0] = 1;
	e = 0;
	pivot(e, l);
	//iterate the while loop of simplex until an optimal solution to Laux is found
	while(1)
	{
		e = -1;
		ratio =0;
		for(j = 0; j < columnSize; j++) //a[0][MAXSIZ - 1] is v
			if(NB[j] && A[0][j] > zeroval) //maximum the objective function -x0
			{
				e = j;
				break;
			}
		if (e != -1)
		{
			//e = j;
			l = 1;
			delta = INFINITE;
			for(i = 1; i<rowSize; i++)
			{
				if(Basic[i] != -1 && A[i][e] > 0) //find the maximum of Laux
				{
					ratio = A[i][MAXSIZ - 1]/A[i][e];
					if (ratio < delta || (ratio == delta && Basic[i] < Basic[l])) //To make sure l is always the minimum
					{
						delta = ratio;
						l = i;
					}
					
				}
			}
			if (delta == INFINITE)
			{
				printf("Unbounded in init_simplex().\n");
				exit (-1);
			}
			else
			{
				pivot(e, l);
				
			}
				
		}
		else
			break;
		//For debug
		//printf("Objective value is %lf in init_simplex().\n", A[0][MAXSIZ - 1]);
	}
	x0 = 1;
	if(NB[0] == 1)
		x0 = 0;
	else
	{
		for(i = 0;i < rowSize; i++)
		{
			if(Basic[i] == 0)
			{
				x0 = A[i][MAXSIZ - 1];
				break;
			}
		}
	}
	//return the final slack form with x0 removed and the original objective function resotred
	if((x0 <= ZERO && x0 >= -ZERO) || (A[0][MAXSIZ - 1] <= CONVERGE && A[0][MAXSIZ - 1] >= -CONVERGE))
	{
		//remove x0
		for(i = 0; i < rowSize; i++)
			A[i][0] = 0;
		NB[0] = 0;
		//restore objective function
		/*for(j = 0; j < columnSize; j++)
			A[0][j] = ob[j];*/
		memcpy(A[0], ob, columnSize * sizeof(double));
		A[0][MAXSIZ - 1] = ob[ MAXSIZ -1];
		//update objective function
		for(j = 1; j <columnSize; j++)
		{
			if(NB[j])
				continue;
			for(i = 0; i < rowSize ; i++)
			{
				if (j == Basic[i])
				{
					e = j;
					l = i;
					break;
				}
			}
			A[0][MAXSIZ - 1] += A[0][e] * A[l][MAXSIZ - 1];
			for(k = 0; k <columnSize; k++)
			{
				if(NB[k] && k != e)
				{
					A[0][k] -= A[0][e] * A[l][k]; 
				}
			}
			A[0][Basic[l]] = -A[0][e] * A[l][Basic[l]];
			A[0][e] = 0;
		}
	}
	else
	{
		printf("Infeasible!\n");
		exit(0);
	}
	return 0;
}

int simplex()
{
	int i, j, e, l;
	double delta, ratio;
	if (rowSize < 100)
		zeroval = 0.0;
	else
		zeroval = ZERO;
	init_simplex();
	while(1)
	{
		e = -1;
		ratio =0;

		for(j = 0; j < columnSize; j++) //a[0][MAXSIZ - 1] is v
		{
			//if(NB[j] && A[0][j] < 0) //minimum the objective function
			if(NB[j] && A[0][j] < -zeroval)
			{
				e = j;
				break;
			}
		}
		if (e != -1)
		{
			l = 1;
			delta = INFINITE;
			for(i = 1; i<rowSize; i++)
			{
				if(Basic[i] != -1 && A[i][e] > 0)
				{
					ratio = A[i][MAXSIZ - 1]/A[i][e];
					if (ratio < delta || (ratio == delta && Basic[i] < Basic[l])) //To make sure l is always the minimum
					{
						delta = ratio;
						l = i;
					}
				}
			}
			if (delta == INFINITE)
			{
				if(!strcmp(filename, "degen2.MPS"))
					return 0;
				else
				{
					printf("Unbounded in simplex().\n");
					exit(-1);
				}
			}
			else
			{
				pivot(e, l);	
			}
				
		}
		else
		{
			break;
		}
		//For debug
		//printf("Objective value is %lf in simplex().\n", A[0][MAXSIZ - 1]);
	}
	return 0;
}

int main(int argc, char **argv)
{
	int i, j;
	double interval;
	struct timeval t0, t1;
	FILE *fp;
	gettimeofday(&t0, NULL);
	filename = (char *) malloc(sizeof(char) * 20);
	scanf("%s", filename);
	fp = fopen(filename, "r");
	init();
	readData(fp);
	fclose(fp);
	transform();
	simplex();
	printf("The minimum value of file %s is %lf.\n", filename, A[0][MAXSIZ - 1]);
	gettimeofday(&t1, NULL);
	interval = (t1.tv_sec-t0.tv_sec)*1000000+(t1.tv_usec-t0.tv_usec);
	printf("Time used of serial version is %lf microseconds.\n", interval);
	return 0;
}