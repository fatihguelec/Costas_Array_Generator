//Program written by Fatih GULEC
//to find (p)th order Costas Arrays by adding a dot to arrays obtained with Welch Method of order (p-1)
//for p which is a prime.

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <vector>
#include <iostream>
using std::vector;

int asal_kontrol(int num);
int *primitive_root(int p);
int phi (int i);
int gcd(int x, int y);
int mypow( int base, int pow, int mod );
int **matrix(int N, int M);
void free_matrix(int **q);
int costas_test(vector<int> v);

//Welch yontemi ile (p-1)'inci
//dereceden bir Costas dizisi olusturur.

FILE *out, *out0;

int main ()
{
    out = fopen("welch1.txt", "w");
    out0 = fopen("welch0.txt", "w");
    int p, num_welch0 = 0;
    int *g;
    int *ptr,i,k;
    int j,tmp;
    int **costas;
    printf ("Please enter a prime number greater than 2.\n");
    scanf ("%d",&p); // p should be a prime

    //p sayisinin asal olup olmadiðinin kontrolü;
    while ((asal_kontrol(p) == 1) || (p <= 2))
    {
        printf("%d is not a prime number greater than 2. Please enter a prime number greater than 2.\n",p);
        scanf ("%d",&p);
    }

    fprintf(out,"%d is a prime number.\n", p);
    //Primitif kokler bulunur.
    //int g[phi(p-1)]; //primitif köklerin tutulduðu g dizisi olusturulur.

    g=new int [phi(p-1)];
    ptr=new int [phi(p-1)];
    costas=matrix(phi(p),p);        // costas[0...phi(p-1)][0..(p-1)] // [0..9][0..22]
    //costas2=matrix(phi(p),p-1);     // costas2[0...phi(p-1)][0..(p-2)] // [0..9][0..21]
    //costas3=new int [p-2];          // costas3[0...p-3] // [0..20]

    ptr = primitive_root(p);

    for (i=0;i<phi(p-1); i++) {
        g[i]= ptr[i];
        fprintf(out, "Primitive root[%d] = %d\n",i, g[i]);
    }

    //Costas dizisinin olusturulmasi.
    //Burada costas matrisinin her satiri ayri bir Costas dizisidir.
    //int costas[i][p-1];

    fprintf(out,"(p-1)=%d'th order Costas array : \n",p-1);
    for(k=0;k<i;k++)
    {
        for(j=0;j<p-1;j++)
        {
            tmp = mypow(g[k],j,p);
            costas[k][j] = tmp;
        }
    }

    for(k=0;k<i;k++)
    {
        fprintf(out,"{");
        for(j=0;j<p-1;j++)
        {
            fprintf(out,"%2d ",costas[k][j]);
        }
        fprintf(out,"}\n");
    }

    //Welch0
    for(k = 0; k < phi(p-1); k++)
    {
        fprintf(out0,"Welch_0 Variant for Costas Array :{");
        for(j = 0; j < p-1; j++)
        {
            fprintf(out0,"%d ", costas[k][j]);
        }
        fprintf(out0,"}\n");

        vector<int> costas_array0;
        costas_array0.resize(p);
        //Corner 1 adding (p-1,p-1)
        costas_array0[p-1] = p;
        for(j = 0; j < p-1; j++)
        {
            costas_array0[j] = costas[k][j];
        }
        if(costas_test(costas_array0) == 1)
        {
            fprintf(out0,"{");
            for(j = 0; j < p; j++)
            {
                fprintf(out0,"%d ", costas_array0[j]);
            }
            fprintf(out0,"}\n");
            num_welch0++;
        }
        //Corner 2 adding (p-1,0)
        for(j = 0; j < p-1; j++)
        {
            costas_array0[j] = costas[k][j]+1;
        }
        costas_array0[p-1] = 1;

        if(costas_test(costas_array0) == 1)
        {
            fprintf(out0,"{");
            for(j = 0; j < p; j++)
            {
                fprintf(out0,"%d ", costas_array0[j]);
            }
            fprintf(out0,"}\n");
            num_welch0++;
        }

        //Corner 3 adding (0,0)
        costas_array0[0] = 1;
        for(j = 1; j < p; j++)
        {
            costas_array0[j] = costas[k][j-1]+1;
        }
        if(costas_test(costas_array0) == 1)
        {
            fprintf(out0,"{");
            for(j = 0; j < p; j++)
            {
                fprintf(out0,"%d ", costas_array0[j]);
            }
            fprintf(out0,"}\n");
            num_welch0++;
        }

        //Corner 4 adding (0,p-1)
        costas_array0[0] = p;
        for(j = 1; j < p; j++)
        {
            costas_array0[j] = costas[k][j-1];
        }
        if(costas_test(costas_array0) == 1)
        {
            fprintf(out0,"{");
            for(j = 0; j < p; j++)
            {
                fprintf(out0,"%d ", costas_array0[j]);
            }
            fprintf(out0,"}\n");
            num_welch0++;
        }
        fprintf(out0,"\n");
    }




    free_matrix(costas);
    free(ptr);
    free(g);
    fprintf(out0,"Number of (%d)th order Welch_0 Costas Arrays = %d\n", p,num_welch0);
    fclose(out);
    return 0;
} // end of int main()

/* Girilen sayinin asal olup olmadiðinin kontrolü
Girilen sayi asal ise fonksiyon 0, asal deðilse 1 dondürür.*/
int asal_kontrol(int num)
{
   int j,flag=0;
   for(j=2;j<=num/2;++j){
        if(num%j==0){
            flag=1;
            break;
        }
   }
   return flag;
}

int *primitive_root(int q) //primitif_kok fonksiyonu p asal sayisinin primitif koklerini bulur.
{
        //Euler phi fonksiyonu ile phi(p-1) adet primitif kok olduðu bulunur.
        //int *g = malloc(phi(q-1)*sizeof(int));
    int *g;
    int GF[q],i,j,m,k,r,n;
    unsigned long temp;
        //En küçük primitif kok bulunur.
        //g = i için {g^1,g'2,...,g^p-1) kümesinin olusturulmasi. Yani GF(p)nin elemanlari.

    g=new int [phi(q-1)];

    for (i=2;i<q; i++) {
        for (j=0;j<(q-1); j++) {
//                temp = pow(i,j);
//                GF[j]=(temp%q);
            GF[j] = mypow(i,j,q);

            fprintf(out,"i=%d, j=%d, q=%d, GF[%d] = %d\n",i,j,q, j, GF[j]); //GF'yi yazdir.

            fprintf(out,"GF[%d] = %d\n",j, GF[j]); //GF'yi yazdir.
        }

            //1. GF'in elemanlarini siraya diz.
        for(m=0;m<(q-1);m++) {
            for(k=m;k<(q-1);k++) {
                if(GF[m] > GF[k]) {
                    temp = GF[m];
                    GF[m] = GF[k];
                    GF[k] = temp;
                }
            }
        }

        fprintf(out,"siralı GF\n");
        //Siraya dizilmis GF'yi yazdir.
        for(n=0;n<(q-1);n++) {
            fprintf(out,"GF[%d] = %d\n",n, GF[n]);
        }

        //2. GF'in elemanlarini sirayla kontrol et.
        int c=0; //sayaç
        for(m=0;m<(q-1);m++) {
            if(GF[m]==m+1) {
                   c++;
            }
        }
        if(c == m) {
            g[0] = i;
            break;
        }

    }

    //fprintf(out,"g[%d]=%d\n",0,g[0]);

    //Diðer primitif kokleri bulma islemi:
    //g([0] mod q'ya gore primitif kok ise OBEB(r,phi(q))=1 için g[0]^r primitif koktür.
    i = 1;
    //for (r = g[0];r < (q-1); r++) {
    for (r = 2;r < (q-1); r++) { // EA
        if (gcd(r,q-1) == 1) {
            g[i] = mypow(g[0],r,q);
            //fprintf(out,"r=%d, q=%d\n",r,q);
            //fprintf(out,"g[%d]=%d\n",i,g[i]);
            i = i+1;
        }

    }

    return g;
} // of int *primitive_root()

//OBEB bulma fonksiyonu
int gcd(int x, int y)
{
    int m,i;

    if(x>y)
         m=y;
    else
         m=x;

    for(i=m;i>=1;i--)
        {
             if((x%i==0)&&(y%i==0))
             {
                 return i;
             }
        }
    return 0;
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

int mypow( int base, int pow, int mod ){
    if( pow == 0 ) return 1;
    if( pow % 2 == 0 )
    {
        int tmp = mypow( base, pow >> 1, mod );
        return tmp * tmp % mod;
    }
    else
    {
        return base * mypow( base, pow - 1, mod ) % mod;
    }
}

int **matrix(int N, int M)
// Allocates two-dimensional array q[0..(N-1)][0..(M-1)]
{
   int **q;
   int j;

   q = new int* [N];
   if(!q) {
     printf("Out of memory\n");
     exit(1);
   }

   for(j=0; j<N; j++) {
      q[j] = new int [M];
      if(!q[j]) {
        printf("Out of memory\n");
        exit(0);
      }
   }
   return(q);
} // of int **matrix()

void free_matrix(int **q)
{
   delete q;
} // of void free_matrix()

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
