#include <iostream>
#include <climits>
#include <cmath>
#include <time.h>
#include <fstream>
using namespace std;

int CROSS_OVER = 80; 

void standard_mat_mul(int** a, int** b, int** result, int dim)
{
	for(int i = 0; i < dim; i++)
	{
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
	}
	else
	{
		n_over_2 = dim/2 + 1; 
		a = m1;
		b = (int**)malloc(sizeof(int*)*n_over_2);
		c = (int**)malloc(sizeof(int*)*n_over_2);
		d = (int**)malloc(sizeof(int*)*n_over_2);
		e = m2;  
		f = (int**)malloc(sizeof(int*)*n_over_2);
		g = (int**)malloc(sizeof(int*)*n_over_2);
		h = (int**)malloc(sizeof(int*)*n_over_2);
		add_buf_1 = (int**)malloc(sizeof(int*)*n_over_2);
		add_buf_2 = (int**)malloc(sizeof(int*)*n_over_2);
		mult_buf = (int**)malloc(sizeof(int*)*n_over_2);
		for(int i = 0; i < n_over_2; i++)
		{
			b[i] = (int*)malloc(sizeof(int)*n_over_2);
			c[i] = (int*)malloc(sizeof(int)*n_over_2);
			d[i] = (int*)malloc(sizeof(int)*n_over_2);
			f[i] = (int*)malloc(sizeof(int)*n_over_2);
			g[i] = (int*)malloc(sizeof(int)*n_over_2);
			h[i] = (int*)malloc(sizeof(int)*n_over_2);
			add_buf_1[i] = (int*)malloc(sizeof(int)*n_over_2);
			add_buf_2[i] = (int*)malloc(sizeof(int)*n_over_2);
			mult_buf[i] = (int*)malloc(sizeof(int)*n_over_2);
			if(i == n_over_2 - 1)
			{
				for(int j = 0; j < n_over_2; j++)
				{
					c[i][j] = 0;
					d[i][j] = 0; 
					g[i][j] = 0;
					h[i][j] = 0;
					add_buf_1[i][j] = 0; 
					add_buf_2[i][j] = 0; 
					mult_buf[i][j] = 0; 
					if(j == n_over_2 - 1)
					{
						b[i][j] = 0; 
						f[i][j] = 0;
					}
					else
					{
						b[i][j] = m1[i][n_over_2+j];
						f[i][j] = m2[i][n_over_2+j];
					}
				}
			}
			else
			{
				for(int j = 0; j < n_over_2; j++)
				{
					if(j == n_over_2 - 1)
					{
						b[i][j] = 0; 
						d[i][j] = 0; 
						f[i][j] = 0; 
						h[i][j] = 0; 
					}
					else
					{
						b[i][j] = m1[i][n_over_2+j];
						d[i][j] = m1[i+n_over_2][n_over_2+j]; 
						f[i][j] = m2[i][n_over_2+j];
						h[i][j] = m2[i+n_over_2][n_over_2+j];  
					}
					c[i][j] = m1[i+n_over_2][j];
					g[i][j] = m2[i+n_over_2][j];
					add_buf_1[i][j] = 0; 
					add_buf_2[i][j] = 0; 
					mult_buf[i][j] = 0; 
				}
			}
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
	//add to af + bh and then cf + dh
	if(dim % 2 == 0)
	{
		add_res(mult_buf, result, 0, n_over_2, n_over_2, n_over_2); 
		add_res(mult_buf, result, n_over_2, n_over_2, n_over_2, n_over_2); 
	}
	else
	{
		add_res(mult_buf, result, 0, n_over_2, n_over_2, n_over_2-1); 
		add_res(mult_buf, result, n_over_2, n_over_2, n_over_2-1, n_over_2-1); 
	}
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
	//add to af + bh then subtract from ae + bg
	if(dim % 2 == 0)
	{
		add_res(mult_buf, result, 0, n_over_2, n_over_2, n_over_2); 
		sub_res(mult_buf, result, 0, 0, n_over_2, n_over_2); 
	}
	else
	{
		add_res(mult_buf, result, 0, n_over_2, n_over_2, n_over_2-1); 
		sub_res(mult_buf, result, 0, 0, n_over_2, n_over_2); 
	}
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
	if(dim % 2 == 0)
	{
		//add to ce + dg
		add_res(mult_buf, result, n_over_2, 0, n_over_2, n_over_2); 
		//subtract from cf + dh
		sub_res(mult_buf, result, n_over_2, n_over_2, n_over_2, n_over_2); 
	}
	else
	{
		//add to ce + dg
		add_res(mult_buf, result, n_over_2, 0, n_over_2-1, n_over_2); 
		//subtract from cf + dh
		sub_res(mult_buf, result, n_over_2, n_over_2, n_over_2-1, n_over_2-1);
	}
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
	if(dim % 2 == 0)
	{
		//add to ce + dg
		add_res(mult_buf, result, n_over_2, 0, n_over_2, n_over_2); 
		//add to ae + bg
		add_res(mult_buf, result, 0, 0, n_over_2, n_over_2); 
	}
	else
	{
		//add to ce + dg
		add_res(mult_buf, result, n_over_2, 0, n_over_2-1, n_over_2); 
		//add to ae + bg
		add_res(mult_buf, result, 0, 0, n_over_2, n_over_2); 
	}	
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
	if(dim % 2 == 0)
	{
		//add to cf + dh
		add_res(mult_buf, result, n_over_2, n_over_2, n_over_2, n_over_2); 
		//add to ae + bg
		add_res(mult_buf, result, 0, 0, n_over_2, n_over_2);
	}
	else
	{
		//add to cf + dh
		add_res(mult_buf, result, n_over_2, n_over_2, n_over_2-1, n_over_2-1); 
		//add to ae + bg
		add_res(mult_buf, result, 0, 0, n_over_2, n_over_2);
	} 
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
	//add to ae + bg - same for both even and odd cases
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
	if(dim % 2 == 0)
	{
		//subtract from cf + dh
		sub_res(mult_buf, result, n_over_2, n_over_2, n_over_2, n_over_2); 
	}
	else
	{
		//subtract from cf + dh
		sub_res(mult_buf, result, n_over_2, n_over_2, n_over_2-1, n_over_2-1);
	}
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

   if(argc > 4)
   {
   	cout << "Too many arguments" << endl; 
   }
   else
   {
   	 //read from command line - TODO file processing
     int version = atoi(argv[1]); 
     srand(time(NULL));
     int test_dim;
     int count; 
     clock_t begin_time;
     int** c; 
     int** d; 
     int** result; 
     int** result_std; 
     int** result_3; 
     int** result_std_3; 
     double rand_n; 
     char* filename;
     if(argc == 4)
     {
     	test_dim = atoi(argv[2]); 
     	filename = argv[3]; 
     }
     int cur; 
     if(version == 0) //grading output
     {
     	ifstream test_input; 
     	test_input.open(filename);
     	c = (int**)malloc(sizeof(int*) * test_dim); 
    	d = (int**)malloc(sizeof(int*) * test_dim); 
     	result = (int**)malloc(sizeof(int*) * test_dim);
     	for(int i = 0; i < test_dim; i++)
     	{
     		result[i] = (int*)malloc(sizeof(int) * test_dim);
     		c[i] = (int*)malloc(sizeof(int) * test_dim);
     		for(int j = 0; j < test_dim; j++)
     		{
     			result[i][j] = 0; 
     			test_input >> c[i][j];
     		}
     	} 
     	for(int i = 0; i < test_dim; i++)
     	{
     		d[i] = (int*)malloc(sizeof(int) * test_dim);
     		for(int j = 0; j < test_dim; j++)
     		{
     			test_input >> d[i][j];
     		}
     	} 
     	strassen_mat_mul(c,d,result,test_dim);
     	for(int i = 0; i < test_dim; i++)
     	{
     		cout << result[i][i] << "\n";
     	}
     	test_input.close(); 
     }
     else if(version == 1) //For experimental tests
     {
     	//Set Cross over to large value so that it is ignored
     	//Therefore, only 1 step of strassen is automatically executed
     	CROSS_OVER = 10000; 
     	float strass_time; 
     	float std_time; 
     	float strass; 
    	float std; 
     	for(test_dim = 5; test_dim < 150; test_dim++)
     	{
     		c = (int**)malloc(sizeof(int*) * test_dim); 
    		d = (int**)malloc(sizeof(int*) * test_dim); 
     		result = (int**)malloc(sizeof(int*) * test_dim); 
     		result_std = (int**)malloc(sizeof(int*) * test_dim); 
     		for(int i = 0; i < test_dim; i++)
     		{
     			d[i] = (int*)malloc(sizeof(int) * test_dim);
     			c[i] = (int*)malloc(sizeof(int) * test_dim);
     		}
     		for(int i = 0; i < test_dim; i++)
     		{	
     			result[i] = (int*)malloc(sizeof(int) * test_dim);
     			result_std[i] = (int*)malloc(sizeof(int) * test_dim);
     			for(int j = 0; j < test_dim; j++)
     			{
     				result[i][j] = 0; 
     				result_std[i][j] = 0; 
     				rand_n = ((float)(rand()))/((float)(RAND_MAX)); 
     				if(rand_n < 0.5)
     				{
     					d[i][j] = 0; 
     				}
     				else
     				{
     					d[i][j] = 1; 
     				}
     				rand_n = ((float)(rand()))/((float)(RAND_MAX)); 
     				if(rand_n < 0.5)
     				{
     					c[i][j] = 0; 
     				}
     				else
     				{
     					c[i][j] = 1; 
     				}
     		}
     		}
     		begin_time = clock();
     		strassen_mat_mul(d,d,result,test_dim);
     		strass_time = float(clock () - begin_time)/CLOCKS_PER_SEC; 
     		begin_time = clock();
     		standard_mat_mul(d,d,result_std,test_dim);
     		std_time = float(clock () - begin_time)/CLOCKS_PER_SEC;
     		//Print out dimensions where 1 step strassen beats standard
     		if(strass_time - std_time < 0)
     		{
     			cout << test_dim << endl; 
     		}
     	}
     }
     else //Random Triangles
     {
     	int trials = 1; 
     	test_dim = 1024; 
   	 	double p = 0.01;  
     	c = (int**)malloc(sizeof(int*) * test_dim); 
     	d = (int**)malloc(sizeof(int*) * test_dim); 
     	result = (int**)malloc(sizeof(int*) * test_dim); 
     	result_3 = (int**)malloc(sizeof(int*) * test_dim); 
     	for(int i = 0; i < test_dim; i++)
     	{
     		d[i] = (int*)malloc(sizeof(int) * test_dim);
     		c[i] = (int*)malloc(sizeof(int) * test_dim);
     		result[i] = (int*)malloc(sizeof(int) * test_dim);
     		result_3[i] = (int*)malloc(sizeof(int) * test_dim);
     	}
     	for(int m = 0; m < 5; m++)
     	{
     		for(int i = 0; i < test_dim; i++)
     		{	
     			for(int j = 0; j < test_dim; j++)
     			{
     				result[i][j] = 0; 
     				result_3[i][j] = 0; 
     				if(j > i)
     				{
     					rand_n = ((float)(rand()))/((float)(RAND_MAX)); 
     					if(rand_n < p)
     					{
     						d[j][i] = 1; 
     						d[i][j] = 1; 
     						c[j][i] = 1; 
     						c[i][j] = 1; 
     					}
     					else
     					{
     						d[j][i] = 0; 
     						d[i][j] = 0; 
     						c[j][i] = 0;
     						c[i][j] = 0; 
     					}
     				}
     			}
     		}
     		//run 10 trials
     		count = 0; 
     		for(int i = 0; i < trials; i++)
     		{
     			strassen_mat_mul(d,d,result,test_dim);
     			strassen_mat_mul(result,d,result_3,test_dim);
     			for(int j= 0; j < test_dim; j++)
     			{
     				count += result_3[j][j]; 
     			}
     			zero_buf(result, test_dim); 
     			zero_buf(result_3, test_dim); 
     		}
     		count /= (6 * trials); 
     		cout << "Triangle Count for p = " << p << " is " << count << endl;
     		p += 0.01; 
     	}


     }


    }
    return 0; 
}