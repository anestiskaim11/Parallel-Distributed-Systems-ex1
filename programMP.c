/*
****************************************************************
authors: Anestis Kaimakamidis 9627 Anastasios Gerontopoulos 9682
webmail: anestisk@ece.auth.gr, ganastas@ece.auth.gr
****************************************************************
OPENMP IMPLEMENTATION
*/
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "mmio.h"
#include <sys/time.h>
#include <sys/times.h>
struct timeval startwtime, endwtime;
static clock_t st_time;
static clock_t en_time;
static struct tms st_cpu;
static struct tms en_cpu;

double seq_time;
double p_time;

//thread variables
int  maxThreads = 16;
omp_lock_t writelock;

typedef struct s_sparse{
	int *v;
	int *row_ind;
	int *col_ptr;
    int n;
	int m;
} sparse;

/*
    Open file(sparse *a): takes the pointer to an empty sparse matrix and fills its values reading an mtx file
*/
void openFile(sparse *a){
    FILE *f;
    //read values
    int M, N, nz;
    int i, *I, *J;
    long int *val;
    int *b;
    mm_read_unsymmetric_sparse("com-Youtube.mtx", &M, &N, &nz, &val, &I, &J);
    //initialise given matrix A
    a->row_ind = (int *) malloc(2*nz * sizeof(int));
    a->col_ptr = (int *) malloc((N+1) * sizeof(int));
    a->n = N;
    a->m = 2 *nz;
    //b: helper array explained in report
    b = (int*)malloc((N+1)*sizeof(int));

    //initialize col_ptr
    for(int i = 0; i < N+1; i++){
        a->col_ptr[i] = 0;
    }
    //fill col_ptr with the number of non-zero values of each column
    for(int i = 0; i < nz; i++){
       a->col_ptr[J[i]+1]++;
       a->col_ptr[I[i]+1]++;
       a->row_ind[2*i] = -1;
       a->row_ind[2*i+1] = -1;
    }
    //create col_ptr final form
    //initialise b
    for(int i = 1; i < N+1; i++){
       int tmp = a->col_ptr[i];
       a->col_ptr[i] = a->col_ptr[i-1] + tmp;
       b[i] = a->col_ptr[i];
    }
    a->col_ptr[N] = 2*nz;
    b[0] = 0;
    //fill row_ind using b
    for(int i = 0; i < nz; i++){
        a->row_ind[b[J[i]]] = I[i];
        a->row_ind[b[I[i]]] = J[i];
        b[I[i]]++;
        b[J[i]]++;
    }
    printf("READ OK\n");
    free(b);
    free(J);
    free(I);
}

/*
    int mult(sparse *a, int row, int col): function that takes as arguments a sparse matrix a row and a column
    returns: (int) the result of the multiplication A*A for this cell
*/
int mult(sparse *a, int row, int col){
	int result = 0, changed = 0, lasti = 0, lastj = 0;;
    int i = a->col_ptr[col], j = a->col_ptr[row];
    if(a->col_ptr[col] == a->col_ptr[col+1]) return 0;
    if(a->col_ptr[row] == a->col_ptr[row+1]) return 0;
    while(1){
       lasti = i;
       lastj = j;
       if(a->row_ind[i] == a->row_ind[j]){
           result++;
           if(j < a->col_ptr[row+1] - 1) j++;
           if(i < a->col_ptr[col+1] - 1) i++;
       }
       else if(a->row_ind[i] > a->row_ind[j]){
        	if(j < a->col_ptr[row+1] - 1) j++;
        	else if(i < a->col_ptr[col+1] - 1)i++;
       }
       else if(a->row_ind[i] < a->row_ind[j]){
            if(i < a->col_ptr[col+1] - 1) i++;
            else if(j < a->col_ptr[row+1] - 1) j++;
       }
       if(i == a->col_ptr[col+1] - 1 && j == a->col_ptr[row+1] - 1 && (lasti < i || lastj < j)) changed = 1;
       if((i >= a->col_ptr[col+1] - 1) && (j >= a->col_ptr[row+1] - 1)){
            if(a->row_ind[i] == a->row_ind[j] && changed) result++;
            break;
       }
    }
    return result;
}

/*
    Thread_calc(sparse *a, sparse *c, int id, int *table): function that makes calculations for each thread given id and working table
    writes back data to the sparse pointer c
*/
void Thread_calc(sparse *a, sparse *c, int id, int *table){
    for(int i = table[id]; i < table[id+1]; i++){
    	if(a->col_ptr[i] > a->col_ptr[i-1]){
    		for(int j = a->col_ptr[i-1]; j < a->col_ptr[i]; j++){
    			int result = mult(a, a->row_ind[j] , i - 1);
    			//c.v = realloc(c.v, (c.m + 1)*sizeof(long int));
    			//c.row_ind = realloc(c.row_ind, (c.m + 1)*sizeof(long int));
                if(result != 0){
                    omp_set_lock(&writelock);
                    c->m++;
                    omp_unset_lock(&writelock);
                }
    			c->v[j] = result;
    			c->row_ind[j] = a->row_ind[j];
            }
        }
    }
}
/*
   DoCalc(sparse *a, sparse *c): initialises threads and controls the work given to each thread with the table array
*/
void DoCalc(sparse *a, sparse *c, int *table){
    omp_init_lock(&writelock);
    int counter = 0, index = 1;
    table[0] = 1;
    for(int i = 1; i < a->n + 1; i++){
        counter += (a->col_ptr[i] - a->col_ptr[i-1]);
        if(counter >= a->m/maxThreads){
            table[index] = i;
            index++;
            counter = 0;
        }
        if(i == a->n) table[maxThreads] = i + 1;
    }
    #pragma omp parallel num_threads(maxThreads)
    {
	    Thread_calc(a, c, omp_get_thread_num(), table);
	}
}
int main(int argc, char *argv[]){
	sparse a;
	openFile(&a);
	/*printf("v: ");
    	for(int i = 0; i < a.m; i++){
    		printf("%ld ", a.v[i]);
    	}
    	printf("\n");
    	printf("rowind: ");
    	for(int i = 0; i < a.m; i++){
    		printf("%d ", a.row_ind[i]);
    	}
    	printf("\n");
    	printf("colptr: ");
    	for(int i = 0; i < a.n + 1; i++){
    		printf("%d ", a.col_ptr[i]);
        }*/
    sparse c;
    c.n = a.n;
    c.m = 0;
    c.v = (int*)malloc(a.m*sizeof(int));
    c.row_ind = (int*)malloc(a.m*sizeof(int));
    c.col_ptr = (int*)malloc((a.n+1)*sizeof(int));
    c.col_ptr[0] = 0;
    for(int i = 1; i < a.n+1; i++) c.col_ptr[i] = a.col_ptr[i];
    int *table = (int*)malloc((maxThreads+1)*sizeof(int));
    double sum = 0;
    for(int i = 0; i < 8; i++){
        gettimeofday (&startwtime, NULL);
        st_time = times(&st_cpu);
        DoCalc(&a, &c, table);
        if(i == 0) printf("nonzeros: %d\n", c.m);
        en_time = times(&en_cpu);

        gettimeofday (&endwtime, NULL);
         p_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
                      + endwtime.tv_sec - startwtime.tv_sec);
        printf("function returned in: %f\n", p_time);
        if(i > 2) sum += p_time;
    }
    printf("average: %f\n", sum/5);
    //printf("\n");
    //printf("nonzeros: %d\n", c.m);
	/*printf("v: ");
	for(int i = 0; i < a.m; i++){
		printf("%d ", c.v[i]);
	}
	printf("\nrowind: ");
	for(int i = 0; i < a.m; i++){
		printf("%d ", c.row_ind[i]);
	}
	printf("\ncolptr: ");
	for(int i = 0; i < c.n + 1; i++){
		printf("%d ", c.col_ptr[i]);
	}*/

	//free(a.v);
	free(table);
	free(a.row_ind);
	free(a.col_ptr);
	free(c.v);
    free(c.row_ind);
    free(c.col_ptr);
	return 0;
}
