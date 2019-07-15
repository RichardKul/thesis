/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on December 27, 2018, 1:15 PM
 */

#include <cstdlib>
#include <cstdio>
//#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
//#include <math.h>
#include <iostream>

using namespace std;

int main(int argc, char **argv) {
 
FILE *ft;
FILE *fp;
FILE *fu;

fu=fopen("deffaven4.txt","w");//gemitteltes D_eff
ft=fopen("geschwn4.txt","w");
fp=fopen("deffn4.txt","w");

//Initialisierung rng
const gsl_rng_type * T;
gsl_rng * r;

gsl_rng_env_setup();

T=gsl_rng_default;
r=gsl_rng_alloc (T);

char *eptr;

//int N=strtol(argv[1], &eptr,10);

int N=200000;//Anzahl Zeitschritte
int O=400;//Anzahl der Kraftschritte
int P=50000;//Zeitschritte, über die gemittelt wird
double df=0.001;
double sigma=1;
//Speicherplatz reservieren
double *v;
v= (double *) malloc((N+1) *100* sizeof(double));
double *x; 
x= (double *) malloc((N+1) *100 * sizeof(double));
double deff[O];
deff[0]=0;
double deff2[O];
deff2[0]=0;
double vav[O];

/*
double dt=strtod(argv[2], &eptr);
double F=strtod(argv[3], &eptr);
double gam=strtod(argv[4], &eptr);
double kT=strtod(argv[5], &eptr);
*/


double dt=0.1;
double F0=0.5;//kleinste Kraft
double kT=0.094;
double gam=0.4;

int h;
int i;
int j;
int l;
int a,b,c,a2,b2,c2;
double F;
double xav[O];
xav[0]=0;
double x2av[O];
x2av[0]=0;
double vavt[N+1],xavt[N+1],x2avt[N+1],defft[N+1];
int s;

//Schleife über versch Kräfte
for(s=0;s<O;s++){

F=F0+s*df;

//Schleife über alle Zeitschritte
for(j=0;j<N;j++){

//Initialisierung von x
for(h=0;h<100;h++){
x[0+h]=0;
}

//Initialisierung von v
for(i=0;i<100;i++){
v[0+i]=0;
}  
//Schleife über alle Teilchen
    for(l=0;l<100;l++){
	v[(j+1)*100+l]=v[j*100+l]-gam*v[j*100+l]*dt+F*dt-sin(x[j*100+l])*dt+sqrt(2*kT*gam*dt)*gsl_ran_gaussian(r,sigma);
	x[(j+1)*100+l]=x[j*100+l]+v[j*100+l]*dt;
    }
    
    double sum4=0;
    for(a2=0;a2<100;a2++){
    sum4=sum4+v[(j+1)*100+a2];
    }
    vavt[j+1]=sum4/100;//mittl Geschwindigkeit
    
    double sum5=0;
    for(b2=0;b2<100;b2++){
    sum5=sum5+x[(j+1)*100+b2];
    }
    xavt[j+1]=sum5/100;//mittl Position
    
    double sum6=0;
    for(c2=0;c2<100;c2++){
    sum6=sum6+x[(j+1)*100+c2]*x[(j+1)*100+c2];
    }
    x2avt[j+1]=sum6/100;//Mittel der Positionsquadrate
    
    defft[j+1]=(x2avt[j+1]-xavt[j+1]*xavt[j+1])/(2*(j+1)*dt);//Diffusionskoeffizient

}

double min=defft[N-P];
int d;
//Finden des Minimums der letzten Werte
for(d=0;d<P;d++){
    if(defft[N-P+d+1]<min){
        min=defft[N-P+d+1];
    }
}

//Finden des Maximums der letzten Werte
double max=defft[N-P];
int e;
for(e=0;e<P;e++){
    if(defft[N-P+e+1]>max){
        max=defft[N-P+e+1];
    }
}

deff2[s]=(min+max)/2;//mittl Diffusionskoeffizient

double sum=0;
    for(a=0;a<100;a++){
    sum=sum+v[N*100+a];
    }
    vav[s]=sum/100;
    
    double sum2=0;
    for(b=0;b<100;b++){
    sum2=sum2+x[N*100+b];
    }
    xav[s]=sum2/100;
    
    double sum3=0;
    for(c=0;c<100;c++){
    sum3=sum3+x[N*100+c]*x[N*100+c];
    }
    x2av[s]=sum3/100;
    
    deff[s]=(x2av[s]-xav[s]*xav[s])/(2*N*dt);//Diffusionskoeffizient aus den letzten Werten
}

//Ausgabe
std::cout <<"wasss kommt raus %f %f %f",vav[5],vav[O-5],deff[O-5];

//int k;
//for(k=0;k<O;k++){
//fprintf(ft,"%f %f\n",F0+k*df,vav[k]);	
//}

//int z;
//for(z=0;z<O;z++){
//fprintf(fp,"%f %f\n",F0+z*df,deff[z]);	
//}
//
//int y;
//for(y=0;y<O;y++){
//fprintf(fu,"%f %f\n",F0+y*df,deff2[y]);	
//}

//Speicher freigeben und Dokumente schließen
//fclose(ft);
//fclose(fu);
gsl_rng_free (r);
//fclose(fp);
free(v);
free(x);
    return 0;
}



