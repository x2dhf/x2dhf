/**************************************************************************
*                                                                         *
*   Copyright (C) 2008 Tomasz Dziubak <tomek@fizyka.umk.pl>               *
*                                                                         *     
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License version 2 as     *
*   published by the Free Software Foundation.                            *
*                                                                         *
***************************************************************************/


#include <math.h>
double AM_HK_ZSimpson(double x, double y, double z, double e0, double w, 
		      double V0, double a, double Nt);


// #include <stdio.h>				

double poth3_(double *z, double *s, int *m, double *V0, double *a, double *precis)

// This function returns a value of the model potential -V0/Sqrt(a^2+r^2)
// transformed into the cylindrical coordinates, i.e. the value
//         -V0/Sqrt(a^2+z^2+s^2)
// plus a m-dependent correction resulting from the azimuthal variable
// being factored out.

// PARAMETERS
// z,s - cylindrical coordinates
//   m - magnetic quantum number
//  V0 - potential depth
//   a - potential core width
//	 if a=0, V0=1 -> Coulomb potential
//	 if a>0, V0>0 -> smoothed Coulomb potential 
//         

{
    double popr_cylind;

    if (*s==0.0 && *z==0.0) 
      {
	popr_cylind=(*m * *m)/(2.0* *precis* *precis);
	return -*V0/sqrt(*a * *a+ *precis* *precis)+popr_cylind;
      }
    
    else 
      {
	//	printf ("%12.6e\n",*precis);
	if (*s==0.0) 
	  popr_cylind=(*m * *m)/(2.0* *precis* *precis);
	else
	  popr_cylind=(*m * *m)/(2.0 * *s * *s);
	return -*V0/sqrt(*a * *a+*z * *z+*s * *s)+popr_cylind;
      }
    
}




double potkh_(double *z, double *s, int *m, double *eps, double *w, 
	      double *V0, double *a, double *N, double *precis)


 
// This function returns a value of the Kramers-Henneberger potential
// at a point (z,s) in cylindrical coordinates
// for the smoothed Coulomb potential plus a m-dependent correction
// correction resulting from the azimuthal variable being factored
// out.

// PARAMETERS
//  z,s - cylindrical coordinates
//    m - magnetic quantum number
//  eps - laser field intensity
//    w - laser cycle frequency
//   V0 - original potential depth
//    a - original potential core width (a>0)
//    N - number of intervals in the Simpson quadrature

{
    	double popr_cylind;

        if (*s==0.0) 
	  {
	    popr_cylind=(*m * *m)/(2.0* *precis* *precis);
            return AM_HK_ZSimpson(*z, *precis, 0.0, *eps, *w, *V0, *a, *N)+popr_cylind;
	  }
	else 
	  {
	    popr_cylind=(*m * *m)/(2.0 * *s * *s);
	    return AM_HK_ZSimpson(*z, *s, 0.0, *eps, *w, *V0, *a, *N)+popr_cylind;
	  }
	


}

// AUXILIARY FUNCTIONS
//-------------------------------------------------
// HK in integral form
double AM_HK(double t, double x, double y, double z, double a, double V0, double e0, double w)
{
    double alpha0=e0/(w*w), bx=x+alpha0+alpha0*(cos(w*t)-1);
    return -V0/((2*M_PI/w)*sqrt(a*a+bx*bx+y*y+z*z));
}


// numerical integration by means of the composite Simpson quadrature
double AM_HK_ZSimpson(double x, double y, double z, double e0, double w, 
		      double V0, double a, double Nt)
{
    double T=2*M_PI/w, h=T/Nt;
    double parts[3];
    double t=T-h;
    int i,k;
    parts[0]=parts[1]=0;
    parts[2]=AM_HK(0,x,y,z,a,V0,e0,w)+AM_HK(T,x,y,z,a,V0,e0,w);
    for (i=Nt,k=0; --i; t-=h,k=!k)
    	parts[k] += AM_HK(t,x,y,z,a,V0,e0,w);

    return (parts[2]+2*parts[1]+4*parts[0])*h/3;
}

/*
#include <stdio.h>

int main (int argc, char* argv)
{
	double z, s, eps, w, V0, a, N;
	double zmin, smin, dz, ds;
	FILE *filename;

	int i,j,m;

	zmin=-20;
	smin=0;
	eps=5.0;
	w=1;
	V0=1;
	a=1;
	N=100;
	m=0;
    dz=0.1;
    ds=0.1;

	filename=fopen("pothk.dat","w");
	fprintf(filename,"#z\ts\tHK(z,s)\n");
	for (i=0; i<400; i++)
	{
		z=zmin+i*dz;
		for (j=0; j<100; j++)
		{
			s=smin+j*ds;
			fprintf(filename,"%5.3f\t%5.3f\t%5.3f\n", z,s,potkh_(&z, &s, &m, &eps, &w, &V0, &a, &N));
		}
	}
	return 0; */



/*

   see http://www.ibiblio.org/pub/languages/fortran/ch2-5.html

*/


