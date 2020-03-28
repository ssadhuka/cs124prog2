#include <iostream>
#include <climits>
#include <cmath>
#include <time.h>
using namespace std;

#define CROSS_OVER 128

void standard_mat_mul(int** a, int** b, int** result, int dim)
{
	for(int i = 0; i < dim; i++)
	{
		//res->rows[i] = (int*)malloc(sizeof(int) * res->col_len); 
		/*for(int m = 0; m < res->col_len; m++)
		{
			res->rows[i][m] = 0; 
		}*/
		for(int j = 0; j < dim; j++)
		{
			for(int k = 0; k < dim; k++)
			{
				result[i][k] += a[i][j] * b[j][k]; 
			}
		}
	}
}
void mat_add(int** a, int** b, int** result, int dim)
{
	for(int i = 0; i < dim; i++)
	{
		for(int j = 0; j < dim; j++)
		{
			result[i][j] = a[i][j] + b[i][j];
		}
	}
}
void mat_sub(int** a, int** b, int** result, int dim)
{
	for(int i = 0; i < dim; i++)
	{
		for(int j = 0; j < dim; j++)
		{
			result[i][j] = a[i][j] - b[i][j];
		}
	}
}
void add_res(int** a, int** result, int r_offset, int c_offset, int r, int c)
{
	for(int i = 0; i < r; i++)
	{
		for(int j = 0; j < c; j++)
		{
			result[r_offset + i][c_offset + j] += a[i][j];
		}
	}
}
void sub_res(int** a, int** result, int r_offset, int c_offset, int r, int c)
{
	for(int i = 0; i < r; i++)
	{
		for(int j = 0; j < c; j++)
		{
			result[r_offset + i][c_offset + j] -= a[i][j];
		}
	}
}
void zero_buf(int** buf, int dim)
{
	for(int i = 0; i < dim; i++)
	{
		for(int j = 0; j < dim; j++)
		{
			buf[i][j] = 0; 
		}
	}
}
void strassen_mat_mul(int** m1, int** m2, int** result, int dim)
{
	if(dim == 1)
	{
		standard_mat_mul(m1,m2,result,1); 
	}
	int n_over_2; 
	int** a; 
	int** b; 
	int** c; 
	int** d; 
	int** e; 
	int** f; 
	int** g; 
	int** h; 
	int** add_buf_1; 
	int** add_buf_2; 
	int** mult_buf; 
	int** ae_bg; 
	int** af_bh; 
	int** cd_dg; 
	int** cf_dh; 
	if(dim % 2 == 0)
	{
		n_over_2 = dim/2;
		a = m1;
		b = (int**)malloc(sizeof(int*)*n_over_2);
		c = &m1[n_over_2];
		d = (int**)malloc(sizeof(int*)*n_over_2);
		e = m2;  
		f = (int**)malloc(sizeof(int*)*n_over_2);
		g = &m2[n_over_2];
		h = (int**)malloc(sizeof(int*)*n_over_2);
		add_buf_1 = (int**)malloc(sizeof(int*)*n_over_2);
		add_buf_2 = (int**)malloc(sizeof(int*)*n_over_2);
		mult_buf = (int**)malloc(sizeof(int*)*n_over_2);
		for(int i = 0; i < n_over_2; i++)
		{
			b[i] = &m1[i][n_over_2]; 
			d[i] = &m1[n_over_2 + i][n_over_2];  
			f[i] = &m2[i][n_over_2]; 
			h[i] = &m2[n_over_2 + i][n_over_2]; 
			add_buf_1[i] = (int*)malloc(sizeof(int)*n_over_2);
			add_buf_2[i] = (int*)malloc(sizeof(int)*n_over_2);
			mult_buf[i] = (int*)malloc(sizeof(int)*n_over_2);
			for(int j = 0; j < n_over_2; j++)
			{
				add_buf_1[i][j] = 0; 
				add_buf_2[i][j] = 0; 
				mult_buf[i][j] = 0; 
			}
		}
		//Compute P1
		mat_sub(f,h,add_buf_1, n_over_2); 
		if(dim > CROSS_OVER)
		{
			strassen_mat_mul(a, add_buf_1, mult_buf, n_over_2); 
		}
		else
		{
			standard_mat_mul(a, add_buf_1, mult_buf, n_over_2); 
		}
		//add to af + bh
		add_res(mult_buf, result, 0, n_over_2, n_over_2, n_over_2); 
		//add to cf + dh
		add_res(mult_buf, result, n_over_2, n_over_2, n_over_2, n_over_2); 
		zero_buf(mult_buf, n_over_2); 
		//Compute P2 
		mat_add(a,b,add_buf_1,n_over_2);
		if(dim > CROSS_OVER)
		{ 
			strassen_mat_mul(add_buf_1, h, mult_buf, n_over_2); 
		}
		else
		{
			standard_mat_mul(add_buf_1, h, mult_buf, n_over_2); 
		}
		//add to af + bh
		add_res(mult_buf, result, 0, n_over_2, n_over_2, n_over_2); 
		//subtract from ae + bg
		sub_res(mult_buf, result, 0, 0, n_over_2, n_over_2); 
		zero_buf(mult_buf, n_over_2); 
		
		//Compute P3
		mat_add(c,d,add_buf_1,n_over_2); 
		if(dim > CROSS_OVER)
		{ 
			strassen_mat_mul(add_buf_1, e, mult_buf, n_over_2);
		} 
		else
		{
			standard_mat_mul(add_buf_1, e, mult_buf, n_over_2);
		}
		//add to ce + dg
		add_res(mult_buf, result, n_over_2, 0, n_over_2, n_over_2); 
		//subtract from cf + dh
		sub_res(mult_buf, result, n_over_2, n_over_2, n_over_2, n_over_2); 
		zero_buf(mult_buf, n_over_2); 
		//Compute P4
		mat_sub(g,e,add_buf_1,n_over_2); 
		if(dim > CROSS_OVER)
		{ 
			strassen_mat_mul(d, add_buf_1, mult_buf, n_over_2);
		} 
		else
		{
			standard_mat_mul(d, add_buf_1, mult_buf, n_over_2);
		}
		//add to ce + dg
		add_res(mult_buf, result, n_over_2, 0, n_over_2, n_over_2); 
		//add to ae + bg
		add_res(mult_buf, result, 0, 0, n_over_2, n_over_2); 
		zero_buf(mult_buf, n_over_2); 
		//Compute P5
		mat_add(a,d,add_buf_1,n_over_2); 
		mat_add(e,h,add_buf_2,n_over_2); 
		if(dim > CROSS_OVER)
		{ 
			strassen_mat_mul(add_buf_1, add_buf_2, mult_buf, n_over_2); 
		}
		else
		{
			standard_mat_mul(add_buf_1, add_buf_2, mult_buf, n_over_2); 
		}
		//add to cf + dh
		add_res(mult_buf, result, n_over_2, n_over_2, n_over_2, n_over_2); 
		//add to ae + bg
		add_res(mult_buf, result, 0, 0, n_over_2, n_over_2); 
		zero_buf(mult_buf, n_over_2); 
		//Compute P6
		mat_sub(b,d,add_buf_1,n_over_2); 
		mat_add(g,h,add_buf_2,n_over_2); 
		if(dim > CROSS_OVER)
		{
			strassen_mat_mul(add_buf_1, add_buf_2, mult_buf, n_over_2); 
		}
		else
		{
			standard_mat_mul(add_buf_1, add_buf_2, mult_buf, n_over_2); 
		}
		//add to ae + bg
		add_res(mult_buf, result, 0, 0, n_over_2, n_over_2); 
		zero_buf(mult_buf, n_over_2); 
		//Compute P7
		mat_sub(a,c,add_buf_1,n_over_2); 
		mat_add(e,f,add_buf_2,n_over_2);
		if(dim > CROSS_OVER)
		{ 
			strassen_mat_mul(add_buf_1, add_buf_2, mult_buf, n_over_2); 
		}
		else
		{
			standard_mat_mul(add_buf_1, add_buf_2, mult_buf, n_over_2); 
		}
		//subtract from cf + dh
		sub_res(mult_buf, result, n_over_2, n_over_2, n_over_2, n_over_2); 
		//zero_buf(mult_buf, n_over_2);
		free(b);  
		free(d);  
		free(f);
		free(h);   
		
		for(int i = 0; i < n_over_2; i++)
		{
			free(add_buf_1[i]);
			free(add_buf_2[i]);
			free(mult_buf[i]);
		}
		free(add_buf_1);
		free(add_buf_2);
		free(mult_buf);
		


	}
	else
	{
		n_over_2 = dim/2 + 1; 
	}

}
bool matrix_equal(int** a, int** b, int dim)
{
	 for(int i = 0; i < dim; i++)
     {
     	for(int j = 0; j < dim; j++)
     	{
     		if(a[i][j] != b[i][j])
     		{
     			return false; 
     		} 
     	}
     }
     return true; 
}
int main(int argc, char **argv) 
{

   if(0/*argc != 4*/)
   {
   	cout << "Invalid Argument Number" << endl; 
   }
   else
   {


   	 //read from command line - TODO file processing
     /*int version = atoi(argv[1]); 
     int dimension = atoi(argv[2]); 
     int filename = atoi(argv[3]); */

   	 int test_dim = 1024; 

     int** c = (int**)malloc(sizeof(int*) * test_dim); 
     int** d = (int**)malloc(sizeof(int*) * test_dim); 
     int** result = (int**)malloc(sizeof(int*) * test_dim); 
     int** result_std = (int**)malloc(sizeof(int*) * test_dim); 
     for(int i = 0; i < test_dim; i++)
     {	
     	result[i] = (int*)malloc(sizeof(int) * test_dim);
     	result_std[i] = (int*)malloc(sizeof(int) * test_dim);
     	d[i] = (int*)malloc(sizeof(int) * test_dim);
     	c[i] = (int*)malloc(sizeof(int) * test_dim);
     	for(int j = 0; j < test_dim; j++)
     	{
     		result[i][j] = 0; 
     		result_std[i][j] = 0; 
     		//Fill up matrices with some arbitrary values - you can change this
     		d[i][j] = (i + 1) * (j + 1); 
     		c[i][j] = 2 * (j + 1) * (i + 1);
     	}
     }
     /*for(int i = 0; i < test_dim; i++)
     {
     	for(int j = 0; j < test_dim; j++)
     	{
     		cout << c[i][j] << " "; 
     	}
     	cout << endl; 
     }
     cout << endl; 
     for(int i = 0; i < test_dim; i++)
     {
     	for(int j = 0; j < test_dim; j++)
     	{
     		cout << d[i][j] << " "; 
     	}
     	cout << endl; 
     }
     cout << endl;*/ 
     srand(time(NULL));
	 clock_t begin_time = clock();
     strassen_mat_mul(c,d,result,test_dim);
     cout << "Time: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
     /*for(int i = 0; i < test_dim; i++)
     {
     	for(int j = 0; j < test_dim; j++)
     	{
     		cout << result[i][j] << " "; 
     	}
     	cout << endl; 
     }
     cout << endl;*/
     //zero_buf(result, test_dim); 
     begin_time = clock();
     standard_mat_mul(c,d, result_std, test_dim);
     cout << "Time: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
     /*for(int i = 0; i < test_dim; i++)
     {
     	for(int j = 0; j < test_dim; j++)
     	{
     		cout << result[i][j] << " "; 
     	}
     	cout << endl; 
     }*/
     //Check if std and strassen results are equal
     cout << matrix_equal(result, result_std, test_dim) << endl;

    }
    return 0; 
}