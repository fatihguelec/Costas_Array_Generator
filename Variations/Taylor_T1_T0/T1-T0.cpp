//Program written by Fatih GULEC
//to find (q-1) and (q)th order Costas Arrays with Taylor's Variant to Golomb Method
//for q which is a prime or prime power and q ≠ 2^k

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <vector>
#include <iostream>
using std::vector;
#define N 100

struct polynomial
{
    //int *poly = new int[N];
    int poly[N] = {0};
};
struct long_polynomial
{
    int poly[200]= {0};
};

void primefactor(int number,int *p, int *k);
int phi (int i);
int gfprimck(polynomial a, int p);
polynomial gfdeconv(polynomial a, polynomial b, int p);
polynomial poly_mul(polynomial f, polynomial g, int p);
polynomial poly_sub(polynomial a, polynomial b, int p);
int GF(int a, int p);
void costas(polynomial g, int p, int k, int q);
polynomial GF_poly_add(polynomial f, polynomial g, int p);
polynomial GF_poly_sub(polynomial f, polynomial g, int p);
polynomial zeroize(polynomial f);
int costas_test(vector<int> v);

FILE *out, *out2, *out3, *out4, *out5;
int num_costas = 0;//# of Costas arrays of (q-2)th degree
int num_taylor1 = 0;//# of Costas arrays of (q-3)th degree
int num_taylor0 = 0;//# of Costas arrays of (q-4)th degree

int main ()
{
    int q, p, k, num_prim, i,j,tmp; //
    out = fopen("Primitive Elements.txt", "w");
    out2 = fopen("The Elements of GF.txt", "w");
    out3 = fopen("Golomb2.txt", "w");
    out4 = fopen("Taylor1.txt", "w");
    out5 = fopen("Taylor0.txt", "w");
//    printf("Enter a prime or prime power bigger than 2\n");
//    scanf ("%d",&q);
//
//    primefactor(q,&p,&k);//q = p^k

    p = 2;
    while(p == 2)
    {
        printf("Enter a prime or prime power bigger than 2\n");
        scanf ("%d",&q);
        primefactor(q,&p,&k);//q = p^k
        if(p == 2)
        {
            printf("Taylor's Variant to Golomb Construction is not possible for q = 2 or its powers!\n");
        }
    }

    fprintf(out,"q=%d\np=%d\nk=%d\n",q,p,k);

    //Number of primitive elements = phi(p^k - 1)/k - Referans bul!
    num_prim = phi(q-1)/k;
    fprintf(out,"Number of primitive elements: %d\n", num_prim);

    int test_end = 2*pow(p,k) - 1;

    // 'test_dec' is the scalar representation of
    // the polynomial that will be tested each cycle.
    int test_dec = pow(p,k) + 1;

    // Cycle through all possible polynomials.
    fprintf(out,"Primitive Elements:\n");
    while ( test_dec <= test_end )
    {
        // Check that this polynomial is not divisible by X.
        if ( (test_dec % p) != 0 )
        {
            // Expand the scalar value to a polynomial in GF(P).
            tmp = test_dec;
            struct polynomial test_poly;
            //test_poly.poly = new int[N];
            for (i = 0; i< k+1; i++)
            {
                test_poly.poly[i] = tmp % p;
                tmp = floor(tmp/p);
            }

            // Test the polynomial.
            if ( gfprimck(test_poly, p) == 1 )
            {
                fprintf(out,"{");
                for (i=0;i<=k; i++)
                {
                    fprintf(out,"%d ", test_poly.poly[i]);
                }
                fprintf(out,"} = ");
                //Print primitive element as x^0+x^1+...+x^n
                for(j = k; j >= 0; j--)//find the biggest order coefficent index
                {
                    if(test_poly.poly[j] != 0)
                    {
                        break;
                    }
                }

                for(i = 0; i <= k; i++)
                {
                    if(i == j)
                    {
                       fprintf(out,"%d*x^%d", test_poly.poly[i],i);
                    }
                    else if(test_poly.poly[i] != 0)
                    {
                        fprintf(out,"%d*x^%d + ", test_poly.poly[i],i);
                    }
                }
                fprintf(out,"\n");

                // Create Costas Arrays for every primitive element.
                costas(test_poly,p,k,q);

            }

        }
        test_dec++;
    }
    fprintf(out,"Number of (%d)th order Golomb2 Costas Arrays = %d\n", q-2,num_costas);
    fprintf(out,"Number of (%d)th order Taylor1 Costas Arrays = %d\n", q-1,num_taylor1);
    fprintf(out,"Number of (%d)th order Taylor0 Costas Arrays = %d\n", q,num_taylor0);

    return 0;
}

//Takes a prime or prime power q as input and finds q = p^k
void primefactor(int number,int *p, int *k)
{
    int div = 2;
    int i = 2;

    while(number!=0)
    {
        if(number%div!=0)
            div = div + 1;
        else
        {
            number = number / div;
            i++;
            if(number==1)
                break;
        }
    }
    *p = div;
    *k = i-2;
}

int phi (int i)
{
	int res; /* Sonuç */
	int j;

	if (i==1) return 1;

        res=i;

        /* Check for divisibility by every prime number below the square root. */
        /* Start with 2. */
        if (i%2==0)
        {
		res-=res/2;
		do i/=2; while (i%2==0) ;
        }

        /* Since this doesn't use a list of primes, check every odd number. Ideally, skip past composite numbers.*/
	for (j=3; j*j<=i; j+=2)
	{
		if (i%j==0)
		{
			res-=res/j;
			do i/=j; while (i%j==0) ;
		}
	}

        /* If i>1, then it's the last factor at this point. */
	if (i>1) res-=res/i;

        /* Return the result. */
	return res;
}

// CK = -1   A is not an irreducible polynomial;
// CK =  0   A is irreducible but not a primitive polynomial;
// CK =  1   A is a primitive polynomial.
int gfprimck(polynomial a, int p)
{
      struct polynomial r;
      int k,i,n,n_a, temp, y;
      // length of a.poly[] n_a= the index of the greaa.poly[i] = a.poly[i] + p;test element of the polynomial
      //int n_a = sizeof(a.poly)/sizeof(a.poly[0]);
      for(i=N-1; i>=0; i--)
      {
          if(a.poly[i] != 0)
            {
                n_a = i+1;
                break;
            }
      }


      // Allocate space for the result, assume primitive.
      int ck = 1;

      int m = n_a - 1;

      // The polynomial is divisible by x, hence is reducible.
      // The only exception is when the polynomial is x ...
    if (a.poly[0] == 0)
    {
        if (n_a == 1)
            ck = 0;
        else
            ck = -1;
    }
    // This polynomial is actually a constant.
    else if ( m == 0 )
        ck = 1;

        // The typical case.
    else
    {

        // First test if the current polynomial is irreducible.
        n = pow(p,(floor(m/2)+1))-1;
        // 'test_dec' is a vector containing the decimal(scalar) representations of
        // the polynomials that could be divisors of 'a_t'.
        //test_dec = p+1:n;
        int len_t = n-p; //length of test_dec[]
        //struct long_polynomial test_dec;
        int test_dec[len_t];
        for(i = 0;i<len_t; i++)
        {
            //test_dec.poly[i] = p+i+1;
            test_dec[i] = p+i+1;
        }

        // test_dec's that correspond to polynomials divisible by X can be removed.
        //test_dec = test_dec( mod(test_dec,p)~=0 );
        for(i = 0;i <= len_t; i++)
        {
            //temp = test_dec.poly[i]%p;
            temp = test_dec[i]%p;
            if(temp == 0)
            {
                for(k = i;k < len_t; k++)
                {
                    /*test_dec.poly[k] = test_dec.poly[k+1];
                    test_dec.poly[len_t] = 0;*/
                    test_dec[k] = test_dec[k+1];
                    test_dec[len_t] = 0;
                }
                len_t--;i--;
            }
        }
        len_t++;

        struct polynomial test_poly;

        int idx = 0,tmp,idx2;
        // Loop through all polynomials that could be divisors of 'at'.
        while ( idx <= len_t-1 )
        {
            // Expand the scalar value to a polynomial in GF(P).
            //tmp = test_dec.poly[idx];
            tmp = test_dec[idx];
            for (idx2 = 0; idx2<m ; idx2++)
            {
                test_poly.poly[idx2] = tmp%p;
                tmp = floor(tmp/p);
            }



            r = gfdeconv(a,test_poly,p);


            //Find the largest element of r[]
            int maxi = r.poly[0];
            for(i = 0;i < N-1; i++)
            {
                if(r.poly[i+1] > r.poly[i])
                {
                    maxi = r.poly[i+1];
                }
            }

            if ( maxi == 0 )
            {
                ck = -1;
                break;
            }
            idx = idx + 1;
        }

        if ( ck == 1 )
        {
            // If the current polynomial is irreducible then check if it is primitive.
            // To be primitive, the polynomial must not be a factor of another
            // polynomial of the form X^n + 1 for any value of n in the range
            //    m < n <p^m - 1
            // To check for this we check to see if the polynomial divides X^n
            // with a remainder of 1 for all values of n in this range.
            int test_ord = m;
            //int test_poly[m+1];
            //test_poly = [zeros(1,m) 1];
            for(i = 0; i<=N; i++)
            {
               test_poly.poly[i] = 0;
            }
            test_poly.poly[m] = 1;

            y = pow(p,m)-1;
            while ( test_ord < y )
            {

                //[ignored, r] = gfdeconv(test_poly, at, p); //calculate the remainder
                r = gfdeconv(test_poly,a,p);
                int len_r; //length of r[]
                for(i = N-1; i>=0; i--)
                {
                    if(r.poly[i] != 0)
                    {
                        len_r = i+1;
                        break;
                    }
                }
                if ((r.poly[0]==1) && (len_r ==1))
                {
                    // If we find a value of n in this range for which the remainder is
                    // 1, we can then conclude the test and declare that the polynomial
                    // is not primitive.
                    ck = 0;
                    break;
                }
                else
                {
                    // To reduce the computational load, on each successive test we
                    // simply need to test against the remainder of the previous test
                    // multiplied by X (i.e., a shifted version of the previous remainder).

                    for(i = 0;i < m+1; i++)
                    {
                        test_poly.poly[i] = 0;
                    }
                    for (idx = 0; idx<len_r ;idx++)
                    {
                        test_poly.poly[idx+1] = r.poly[idx];
                    }
                    test_ord = test_ord + 1;
                }
            }

        }
    }


    return ck;
}

polynomial gfdeconv(polynomial a, polynomial b, int p) // a(x) = b(x)*q(x) + r(x) in GF(p). It returns the remaining part.
{
    int i,j,t;
    double temp2;
    struct polynomial q;  //quotient;;
    //j = index of the greatest order of b.poly[] = degree of b.poly[]
    for(j =  N-1;j >= 0; j--)
    {
         if(b.poly[j] != 0)
            {
                break;
            }
    }
    // degree of a.poly[]
    for(i =  N-1;i >= 0; i--)
    {
         if(a.poly[i] != 0)
            {
                break;
            }
    }
    // if the degree of a.poly[] is smaller than b.poly[]
    if(i<j)
    {
        return a;
    }

    for(i =  N-1;i > 0; i--)
    {
        if((a.poly[i]!= 0)&&(i>=j))
        {
                struct polynomial temp;

                while(a.poly[i] < b.poly[j])
                {
                    a.poly[i] = a.poly[i] + p;
                }

                temp2 = a.poly[i]/b.poly[j];
                t = floor(temp2);
                q.poly[i-j] = GF(t,p);
                temp.poly[i-j] = q.poly[i-j];
                a = poly_sub(a, poly_mul(temp, b, p), p);

        }

    }

    return a; //remaining r(x)
}
polynomial poly_mul(polynomial f, polynomial g, int p)//h(x)=f(x)*g(x) in GF(p)
{
    int i,j,tmp;
    polynomial h;

    for (i=0;i<= N;i++)
    {
        tmp=0;
        for (j=0;j<=i;j++)
        {
            tmp += f.poly[j]*g.poly[i-j];
            h.poly[i]= GF(tmp,p);
        }
    }

    return(h);
}

polynomial poly_sub(polynomial a, polynomial b, int p)//c(x) = a(x) - b(x) in GF(p)
{
    int i;
    polynomial c;

    for(i=0; i<=N; i++)
    {
        c.poly[i] = GF((a.poly[i]-b.poly[i]),p);
    }

    return c;
}

int GF(int a, int p)
{
    int x=0;
    x = a%p;

    while(x < 0)
    {
        x = x + p;
    }

    return x;
}

void costas(polynomial g, int p, int k, int q)
{
    int i,j,m, n,r, cnt = 0;
    int idx = 0;
    struct polynomial temp, temp3, temp4, temp5;
    vector<int> costas_array; costas_array.resize(q-2);

    //Allocate for the elements of F for alfa
    vector<vector<int> > F;
    F.resize(q);
    for(i = 0; i<q; i++)
    {
        F[i].resize(k);
    }

    //create smaller elements of F[][] smaller than a^k for alfa
    //For example: For k=2, create 0,1,a
    //             For k=4, create 0,1,a,a^2,a^3
    F[0][0] = 0;
    for(i = 1; i<= k; i++)
    {
        F[i][i-1] = 1;
    }


    //create a^k for alfa
    for (j = k-1; j>=0; j--)
    {
        F[k+1][j] = GF((-1)*g.poly[j], p);
    }

    //create the elements bigger than a^k for alfa

    for(i = k+1; i<= q-2; i++)
    {

        for(j = 0; j<k; j++)
        {
            temp.poly[j] = F[i][j];
        }
        struct polynomial alfa; alfa.poly[1] = 1;// define a^1
        temp = poly_mul(temp, alfa, p);
        if(temp.poly[k] != 0)
        {
            struct polynomial temp2;
            temp2.poly[0] = temp.poly[k];
            for(m = 0; m <= N; m++)
            {
               temp3.poly[m] = F[k+1][m];
            }
            temp3 = poly_mul(temp3, temp2, p);
            temp = GF_poly_add(temp, temp3, p);
        }

        for(j = 0; j < k; j++)
        {
            F[i+1][j] = temp.poly[j];
        }

    }

    fprintf(out2, "The elements of F(%d^%d) for primitive element \n{", p, k);
    for (i=0;i<=k; i++)
    {
        fprintf(out2,"%d ", g.poly[i]);
    }
    fprintf(out2,"} = ");
    //Print primitive element as x^0+x^1+...+x^n
    for(j = k; j >= 0; j--)//find the biggest order coefficent index
    {
        if(g.poly[j] != 0)
        {
            break;
        }
    }

    for(i = 0; i <= k; i++)
    {
        if(i == j)
        {
           fprintf(out2,"%d*x^%d", g.poly[i],i);
        }
        else if(g.poly[i] != 0)
        {
            fprintf(out2,"%d*x^%d + ", g.poly[i],i);
        }
    }
    fprintf(out2," and for beta = a\n");

    for(i = 0; i< q; i++)
    {
        for(j = 0; j < k; j++)
        {
           fprintf(out2, "F[%d][%d] = %d \n", i, j, F[i][j] );
        }
        fprintf(out2, "\n");
    }

    //Golomb algorithm
    //Create Costas arrays of degree(q-2) for beta = a
    // a^i + B^j = 1 <==> the coordinates of Costas array  is (i,j)

    for(i = 2; i < q; i++)
    {
        for(m = 0; m < k; m++)
        {
            temp4.poly[m] = F[i][m];
        }
        struct polynomial one;
        one.poly[0] = 1;
        temp5 = GF_poly_sub(one, temp4, p); // 1-a^i
        //compare for the elements of F(p^k) to find B^j = 1-a^i
        for(n = 2; n < q; n++)
        {
            for(j = 0; j < k; j++)
            {
                if(F[n][j] == temp5.poly[j])
                {
                    cnt++;
                }

            }
            if(cnt == k)
            {
                costas_array[idx] = n-1;
                idx++;cnt = 0;
                break;
            }
            cnt = 0;

        }
    }

    //Print Costas array for beta = a
    fprintf(out3,"\nFor primitive element:{");
    for (i=0;i<=k; i++)
    {
        fprintf(out3,"%d ", g.poly[i]);
    }
    fprintf(out3,"} = ");
    //Print primitive element as x^0+x^1+...+x^n
    for(j = k; j >= 0; j--)//find the biggest order coefficent index
    {
        if(g.poly[j] != 0)
        {
            break;
        }
    }

    for(i = 0; i <= k; i++)
    {
        if(i == j)
        {
           fprintf(out3,"%d*x^%d", g.poly[i],i);
        }
        else if(g.poly[i] != 0)
        {
            fprintf(out3,"%d*x^%d + ", g.poly[i],i);
        }
    }
    fprintf(out3," and for beta = a\n(%d)th order Costas Array :\n{", q-2);

    for(j = 0; j < q-2; j++)
    {
        fprintf(out3,"%d ", costas_array[j]);
    }
    fprintf(out3,"}\n");
    num_costas++;

    //Taylor1 Variant Construction by adding a dot to one of the corners and testing
    fprintf(out4,"Taylor1 Variant for Costas Array :{");
    for(j = 0; j < q-2; j++)
    {
        fprintf(out4,"%d ", costas_array[j]);
    }
    fprintf(out4,"}\n");

    vector<int> costas_array1;
    costas_array1.resize(q-1);
    //Corner 1 adding (q-1,q-1)
    costas_array1[q-2] = q-1;
    for(j = 0; j < q-2; j++)
    {
        costas_array1[j] = costas_array[j];
    }
    if(costas_test(costas_array1) == 1)
    {
        fprintf(out4,"{");
        for(j = 0; j < q-1; j++)
        {
            fprintf(out4,"%d ", costas_array1[j]);
        }
        fprintf(out4,"}\n");
        num_taylor1++;
    }
    //Corner 2 adding (q-1,0)
    for(j = 0; j < q-2; j++)
    {
        costas_array1[j] = costas_array[j]+1;
    }
    costas_array1[q-2] = 1;

    if(costas_test(costas_array1) == 1)
    {
        fprintf(out4,"{");
        for(j = 0; j < q-1; j++)
        {
            fprintf(out4,"%d ", costas_array1[j]);
        }
        fprintf(out4,"}\n");
        num_taylor1++;
    }

    //Corner 3 adding (0,0)
    costas_array1[0] = 1;
    for(j = 1; j < q-1; j++)
    {
        costas_array1[j] = costas_array[j-1]+1;
    }
    if(costas_test(costas_array1) == 1)
    {
        fprintf(out4,"{");
        for(j = 0; j < q-1; j++)
        {
            fprintf(out4,"%d ", costas_array1[j]);
        }
        fprintf(out4,"}\n");
        num_taylor1++;
    }

    //Corner 4 adding (0,q-1)
    costas_array1[0] = q-1;
    for(j = 1; j < q-1; j++)
    {
        costas_array1[j] = costas_array[j-1];
    }
    if(costas_test(costas_array1) == 1)
    {
        fprintf(out4,"{");
        for(j = 0; j < q-1; j++)
        {
            fprintf(out4,"%d ", costas_array1[j]);
        }
        fprintf(out4,"}\n");
        num_taylor1++;
    }
    fprintf(out4,"\n");

    //Taylor0 Variant Construction by adding 2 dots to cross corners and testing
    fprintf(out5,"Taylor0 Variant for Costas Array :{");
    for(j = 0; j < q-2; j++)
    {
        fprintf(out5,"%d ", costas_array[j]);
    }
    fprintf(out5,"}\n");

    vector<int> costas_array0;
    costas_array0.resize(q);

    //Cross Corner 1 adding (0,0) and (q-1,q-1)
    costas_array0[0] = 1;
    for(j = 1; j < q-1; j++)
    {
        costas_array0[j] = costas_array[j-1]+1;
    }
    costas_array0[q-1] = q;

    if(costas_test(costas_array0) == 1)
    {
        fprintf(out5,"{");
        for(j = 0; j < q; j++)
        {
            fprintf(out5,"%d ", costas_array0[j]);
        }
        fprintf(out5,"}\n");
        num_taylor0++;
    }
    //fprintf(out5,"\n");

    //Cross Corner 2 adding (0,q-1) and (q-1,0)
    costas_array0[0] = q;
    for(j = 1; j < q-1; j++)
    {
        costas_array0[j] = costas_array[j-1]+1;
    }
    costas_array0[q-1] = 1;
    if(costas_test(costas_array0) == 1)
    {
        fprintf(out5,"{");
        for(j = 0; j < q; j++)
        {
            fprintf(out5,"%d ", costas_array0[j]);
        }
        fprintf(out5,"}\n");
        num_taylor0++;
    }
    fprintf(out5,"\n");

    //Allocate for the elements of Fb for beta
    vector<vector<int> > Fb;
    Fb.resize(q);
    for(i = 0; i<q; i++)
    {
        Fb[i].resize(k);
    }


    //Create element of Fb for beta = a^r
    for(r = 2; r <= q-2; r++)
    {
        idx = 0;
        bool primitive = true;
        //zeroize used polynomials
        temp = zeroize(temp);
        temp3 = zeroize(temp3);
        temp4 = zeroize(temp4);
        temp5 = zeroize(temp5);

        //create first elements of Fb[][]
        for(j = 0; j < k; j++)
        {
            Fb[0][j] = 0;
            Fb[1][j] = 0;
        }
        Fb[1][0] = 1;

        //create beta = a^r
        for (j = 0; j < k; j++)
        {
            Fb[2][j] = F[r+1][j];
        }

        //create the elements bigger than B^k

        for(i = 2; i<= q-2; i++)
        {
            for(j = 0; j<N; j++)
            {
                if(j>=k)
                {
                    temp.poly[j] = 0;
                }
                else
                {
                    temp.poly[j] = Fb[i][j];
                }
            }

            //define and create B^r
            struct polynomial beta;
            for(j = 0; j<N; j++)
            {
                if(j >= k)
                {
                    beta.poly[j] = 0;
                }
                else
                {
                    beta.poly[j] = F[r+1][j];
                }
            }
            temp = poly_mul(temp, beta, p);
            for(m = r+k-1; m >= k; m--)
            {
                if(temp.poly[m] != 0)
                {
                    for(j = 0; j < N; j++)
                    {
                        if(j >= k)
                        {
                            temp3.poly[j] = 0;
                        }
                        else
                        {
                            temp3.poly[j] = GF(F[m+1][j]*temp.poly[m],p);
                        }

                    }
                    temp.poly[m] = 0;
                    temp = GF_poly_add(temp, temp3, p);
                }
            }

            for(j = 0; j < k; j++)
            {
                Fb[i+1][j] = temp.poly[j];
            }

            temp3 = zeroize(temp3);

            for(j = 0; j < k; j++)
            {
               printf("Fb[%d][%d] = %d \n", i+1, j, Fb[i+1][j] );
            }
            printf("\n");
        }


        //Print the element of Fb for beta = a^r
        fprintf(out2, "The elements of Fb(%d^%d) for primitive element \n{", p, k);
        for (i=0;i<=k; i++)
        {
            fprintf(out2,"%d ", g.poly[i]);
        }
        fprintf(out2,"} = ");
        //Print primitive element as x^0+x^1+...+x^n
        for(j = k; j >= 0; j--)//find the biggest order coefficent index
        {
            if(g.poly[j] != 0)
            {
                break;
            }
        }

        for(i = 0; i <= k; i++)
        {
            if(i == j)
            {
               fprintf(out2,"%d*x^%d", g.poly[i],i);
            }
            else if(g.poly[i] != 0)
            {
                fprintf(out2,"%d*x^%d + ", g.poly[i],i);
            }
        }
        fprintf(out2," and for beta = a^%d\n",r);

        for(i = 0; i< q; i++)
        {
            for(j = 0; j < k; j++)
            {
               fprintf(out2, "Fb[%d][%d] = %d \n", i, j, Fb[i][j] );
            }
            fprintf(out2, "\n");
        }

        //Check if the elements of Fb for beta = a^r primitive or not.
        //To do this check if there is a repeating element.
        int cnt2 = 0;
        for(i = k+2; i < q; i++)
        {
            //compare F[k+1]'s(x^(k+1)) elements with all of the elements beyond it.
            //For example: q=9 k=2 ==> start from a^2
            //             q=16 k=4 ==> start from a^4
            for(j = 0; j < k; j++)
            {
                if(Fb[k+1][j] == Fb[i][j])
                {
                    cnt2++;
                }
            }

            if(cnt2 == k)
            {
                primitive = false;
                cnt2 = 0;
                break;
            }
            cnt2 = 0;
        }

        //Golomb algorithm
        //Create Costas arrays of degree(q-2) for beta = a^r
        // a^i + B^j = 1 <==> the coordinates of Costas array  is (i,j)

        if(primitive == true)
        {
            for(i = 2; i < q; i++)
            {


                for(m = 0; m < k; m++)
                {
                    temp4.poly[m] = F[i][m];
                }
                struct polynomial one;
                one.poly[0] = 1;
                temp5 = GF_poly_sub(one, temp4, p); // 1-a^i
                //compare for the elements of F(p^k) to find B^j = 1-a^i
                for(n = 2; n < q; n++)
                {
                    for(j = 0; j < k; j++)
                    {
                        if(Fb[n][j] == temp5.poly[j])
                        {
                            cnt++;
                        }

                    }
                    if(cnt == k)
                    {
                        costas_array[idx] = n-1;
                        idx++;cnt = 0;
                        break;
                    }
                    cnt = 0;

                }
            }

            for(j = 0; j < q-2; j++)
            {
                //fprintf(out3,"%d ", costas_array.poly[j]);
                printf("%d ", costas_array[j]);
            }

            fprintf(out3,"\nFor primitive element:{");
            for (i=0;i<=k; i++)
            {
                fprintf(out3,"%d ", g.poly[i]);
            }
            fprintf(out3,"} = ");
            //Print primitive element as x^0+x^1+...+x^n
            for(j = k; j >= 0; j--)//find the biggest order coefficent index
            {
                if(g.poly[j] != 0)
                {
                    break;
                }
            }

            for(i = 0; i <= k; i++)
            {
                if(i == j)
                {
                   fprintf(out3,"%d*x^%d", g.poly[i],i);
                }
                else if(g.poly[i] != 0)
                {
                    fprintf(out3,"%d*x^%d + ", g.poly[i],i);
                }
            }
            fprintf(out3," and for beta = a^%d\n(%d)th order Costas Array :\n{",r, q-2);

            for(j = 0; j < q-2; j++)
            {
                fprintf(out3,"%d ", costas_array[j]);
            }
            fprintf(out3,"}\n");
            num_costas++;

            //Taylor1 Variant Construction by adding a dot to one of the corners and testing
            fprintf(out4,"Taylor1 Variant for Costas Array :{");
            for(j = 0; j < q-2; j++)
            {
                fprintf(out4,"%d ", costas_array[j]);
            }
            fprintf(out4,"}\n");

            vector<int> costas_array1;
            costas_array1.resize(q-1);
            //Corner 1 adding (q-1,q-1)
            costas_array1[q-2] = q-1;
            for(j = 0; j < q-2; j++)
            {
                costas_array1[j] = costas_array[j];
            }

            if(costas_test(costas_array1) == 1)
            {
                fprintf(out4,"{");
                for(j = 0; j < q-1; j++)
                {
                    fprintf(out4,"%d ", costas_array1[j]);
                }
                fprintf(out4,"}\n");
                num_taylor1++;
            }
            //Corner 2 adding (q-1,0)
            for(j = 0; j < q-2; j++)
            {
                costas_array1[j] = costas_array[j]+1;
            }
            costas_array1[q-2] = 1;

            if(costas_test(costas_array1) == 1)
            {
                fprintf(out4,"{");
                for(j = 0; j < q-1; j++)
                {
                    fprintf(out4,"%d ", costas_array1[j]);
                }
                fprintf(out4,"}\n");
                num_taylor1++;
            }

            //Corner 3 adding (0,0)
            costas_array1[0] = 1;
            for(j = 1; j < q-1; j++)
            {
                costas_array1[j] = costas_array[j-1]+1;
            }
            if(costas_test(costas_array1) == 1)
            {
                fprintf(out4,"{");
                for(j = 0; j < q-1; j++)
                {
                    fprintf(out4,"%d ", costas_array1[j]);
                }
                fprintf(out4,"}\n");
                num_taylor1++;
            }

            //Corner 4 adding (0,q-1)
            costas_array1[0] = q-1;
            for(j = 1; j < q-1; j++)
            {
                costas_array1[j] = costas_array[j-1];
            }
            if(costas_test(costas_array1) == 1)
            {
                fprintf(out4,"{");
                for(j = 0; j < q-1; j++)
                {
                    fprintf(out4,"%d ", costas_array1[j]);
                }
                fprintf(out4,"}\n");
                num_taylor1++;
            }
            fprintf(out4,"\n");
            //Taylor0 Variant Construction by adding 2 dots to cross corners and testing
            fprintf(out5,"Taylor0 Variant for Costas Array :{");
            for(j = 0; j < q-2; j++)
            {
                fprintf(out5,"%d ", costas_array[j]);
            }
            fprintf(out5,"}\n");

            vector<int> costas_array0;
            costas_array0.resize(q);

            //Cross Corner 1 adding (0,0) and (q-1,q-1)
            costas_array0[0] = 1;
            for(j = 1; j < q-1; j++)
            {
                costas_array0[j] = costas_array[j-1]+1;
            }
            costas_array0[q-1] = q;

            if(costas_test(costas_array0) == 1)
            {
                fprintf(out5,"{");
                for(j = 0; j < q; j++)
                {
                    fprintf(out5,"%d ", costas_array0[j]);
                }
                fprintf(out5,"}\n");
                num_taylor0++;
            }
            //fprintf(out5,"\n");

            //Cross Corner 2 adding (0,q-1) and (q-1,0)
            costas_array0[0] = q;
            for(j = 1; j < q-1; j++)
            {
                costas_array0[j] = costas_array[j-1]+1;
            }
            costas_array0[q-1] = 1;
            if(costas_test(costas_array0) == 1)
            {
                fprintf(out5,"{");
                for(j = 0; j < q; j++)
                {
                    fprintf(out5,"%d ", costas_array0[j]);
                }
                fprintf(out5,"}\n");
                num_taylor0++;
            }
            fprintf(out5,"\n");

        }// end of primitive control

    }

}

polynomial GF_poly_add(polynomial f, polynomial g, int p)//h(x) = f(x)+g(x) in GF(p)
{
    int i,tmp;
    struct polynomial h;

    for (i=0; i< N; i++)
    {
        tmp = f.poly[i] + g.poly[i];
        h.poly[i] = GF(tmp, p);
    }

    return(h);
}

polynomial GF_poly_sub(polynomial f, polynomial g, int p)//h(x) = f(x)-g(x) in GF(p)
{
    int i,tmp;
    struct polynomial h;

    for (i=0; i< N; i++)
    {
        tmp = f.poly[i] - g.poly[i];
        h.poly[i] = GF(tmp, p);
    }

    return(h);
}

polynomial zeroize(polynomial f)
{
    int i;
    for(i = 0; i < N; i++)
    {
        f.poly[i] = 0;
    }
    return f;
}

int costas_test(vector<int> v)
{
    int i,j,k,n,m, res = 1;
    int D = v.size();
    vector<int> diff;
    diff.resize(D-1);

    for (j=1;j<D;j++)
    {
        for (i=0;i<(D-1);i++)
        {
            k=i+j;
            if (k<D)
            {
                diff[i] = v[i+j]-v[i];
            }
        }
        //Similarity control in diff array
        for(n = 0; n < D-j; n++)
        {
            for(m = 0; m < D-j; m++)
            {
                if((diff[n] == diff[m])&&(n != m))
                {
                    res = 0;
                    break;
                }
            }
            if(res == 0)
            {
                break;
            }
        }
    }

    return res;
}
