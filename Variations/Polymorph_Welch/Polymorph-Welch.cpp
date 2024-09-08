//Program written by Fatih GULEC
//to find polymorphs of (p-1), (p-2) and (p-3)th order Costas Arrays with Welch Method
//for p which is a prime.

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<vector>
#include <iostream>

using std::vector;

int asal_kontrol(int num);
int *primitive_root(int p);
int phi (int i);
int gcd(int x, int y);
int mypow( int base, int pow, int mod );
int **matrix(int N, int M);
void free_matrix(int **q);
void find_polymorphs(int p, vector<int> polymorph);
vector<int> rotate_90(int n, vector<int> polymorph);

//Welch yontemi ile (p-1)'inci
//dereceden bir Costas dizisi olusturur.

FILE *out;
FILE *out2;

int main ()
{
    out = fopen("welch1.txt", "w");
    out2 = fopen("Polymorph-Welch.txt", "w");
    int p;
    int *g;
    int *ptr,i,k;
    int j,tmp;
    int **costas;
    printf ("Please enter a prime numer greater than 2.\n");
    scanf ("%d",&p); // p should be a prime

    //p sayisinin asal olup olmadiðinin kontrolü;
    while ((asal_kontrol(p) == 1) || (p <= 2))
    {
        printf("%d is not a prime number greater than 2. Please enter a prime numer greater than 2.\n",p);
        scanf ("%d",&p);
    }

    fprintf(out,"%d is a prime number.\n", p);

    vector<int> polymorph;
    polymorph.resize(p-1);
    //Primitif kokler bulunur.
    //int g[phi(p-1)]; //primitif köklerin tutulduðu g dizisi olusturulur.


    g=new int [phi(p-1)];
    ptr=new int [phi(p-1)];
    costas=matrix(phi(p),p-1);        // costas[0...phi(p-1)][0..(p-1)] // [0..9][0..22]
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
            polymorph[j] = costas[k][j];
        }
        fprintf(out,"}\n");
        find_polymorphs(p-1, polymorph);
    }


//(p-2)inci dereceden Costas dizisi. (1,1) koşe noktasının bulunduğu satır ve sütunlar silinerek
//(p-2)'inci dereceden Costas dizisi elde edilir.
    polymorph.resize(p-2);
    fprintf(out,"(p-2)=%d'th order Costas array : \n",p-2);
    for(k=0;k<i;k++)
    {
        fprintf(out,"{");
        for(j=1;j<p-1;j++)
        {
            fprintf(out,"%2d ",costas[k][j]-1);
            polymorph[j-1] = costas[k][j]-1;
        }
        fprintf(out,"}\n");
        find_polymorphs(p-2, polymorph);
    }
    fprintf(out,"\n");

    //Eğer 2 primitif kök olarak var ise (p-2)inci dereceden Costas dizisinin (1,2) elemanının olduğu
    //satır ve sütun silinerek (p-3)uncu dereceden Costas dizisi elde edilir.

    if((costas[0][1] == 2)&& (p>3))
    {
        polymorph.resize(p-3);
        fprintf(out,"(p-3)=%d'th order Costas array : \n{",p-3);
        for(j=2;j<p-1;j++)
        {
            fprintf(out,"%2d ",costas[0][j]-2);
            polymorph[j-2] = costas[0][j]-2;
        }
        fprintf(out,"}\n");
        find_polymorphs(p-3, polymorph);
    }
    else
    {
        fprintf(out,"(p-3)'th order Costas array cannot be constructed.\n");
    }

    free_matrix(costas);
    free(ptr);
    free(g);

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

//Finds and prints the polymorphs of a Costas array.
void find_polymorphs(int n, vector<int> polymorph)
{
    int i,j;
    vector<int> temp2; temp2.resize(n);

    for(i = 0; i < n; i++)
    {
        temp2[i] = polymorph[i];
    }

    fprintf(out2,"Polymorphs of the Costas Array = {");
    for(j=0;j<n;j++)
    {
        fprintf(out2,"%2d ",polymorph[j]);
    }
    fprintf(out2,"}\n");

    //90 degrees CW
    polymorph = rotate_90(n, polymorph);

    //180 degrees
    polymorph = rotate_90(n, polymorph);

    //270 degrees CW
    polymorph = rotate_90(n, polymorph);

    //Diagonal Transposition

    for(i = 0; i < n; i++)
    {
        polymorph[temp2[i]-1] = i+1;
    }

    if(temp2 == polymorph)
    {
        fprintf(out2,"Diagonally Symmetric\n");
    }
    else
    {
        fprintf(out2,"{");
        for(i = 0; i < n; i++)
        {
            fprintf(out2,"%d ",polymorph[i]);
        }
        fprintf(out2,"} //Diagonal\n");


        //Diagonal 90 degrees CW
        polymorph = rotate_90(n, polymorph);

        //Diagonal 180 degrees CW
        polymorph = rotate_90(n, polymorph);

        //Diagonal 270 degrees
        polymorph = rotate_90(n, polymorph);

    }

    fprintf(out2,"\n");
}

//Rotates the Costas array 90 degrees and returns the new rotated Costas array(Polymorph)
vector<int> rotate_90(int n, vector<int> polymorph)
{
    int i;
    vector<int> temp; temp.resize(n);
    for(i = 0; i < n; i++)
    {
        temp[i] = polymorph[n-1-i];
    }

    for(i = 0; i < n; i++)
    {
        polymorph[temp[i]-1] = i+1;
    }
    fprintf(out2,"{");
    for(i = 0; i < n; i++)
    {
        fprintf(out2,"%d ",polymorph[i]);
    }
    fprintf(out2,"}\n");

    return polymorph;
}
