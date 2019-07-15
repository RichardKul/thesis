/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on January 17, 2019, 12:14 AM
 */

#include <cstdlib>
#include <random>
#include <cmath>
#include <fstream>
#include <ios>
#include <iomanip>
#include <fftw3.h>
#include <math.h>
#include <complex>

#define REAL 0
#define IMAG 1


void fft(fftw_complex *in, fftw_complex *out, int N)
{
	// create a DFT plan
	fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	// execute the plan
	fftw_execute(plan);
	// do some cleaning
	fftw_destroy_plan(plan);
	fftw_cleanup();
}

double average(double x[],int n);
void square(long long int x[],int n);
double ninf(double V, double k, double Vh);
void initd(double x[],int n);
int main(int argc, char** argv) {
    
    double d=atof(argv[1]);
    double D=d*0.01;
    double nr=atof(argv[2]);
    double I=-0.1+0.02*nr;
std::ostringstream d_convert;
std::ostringstream nr_convert;
    d_convert << d;
    nr_convert << nr;
std::	 string date="realfast11jjem2st";
std::    string d_str = d_convert.str();
std::    string nr_str = nr_convert.str();
std::    string filenamef="f"+date+d_str+nr_str+".txt";
std::    string filenamed="d"+date+d_str+nr_str+".txt";
std::    string filenameg="g"+date+d_str+nr_str+".txt";
std::    string filenametime="time"+date+d_str+nr_str+".txt";
std::    string filenametimenew="timenew"+date+d_str+nr_str+".txt";
std::    string filenameparam="param"+date+d_str+nr_str+".txt";
std::    string filenameft="ft"+date+d_str+nr_str+".txt";
std::	 string filenamesp="spike"+date+d_str+nr_str+".txt";
std::ofstream myfile;
std::ofstream myfile2;
std::ofstream myfile3;
std::ofstream myfile4;
std::ofstream myfile5;
std::ofstream myfile6;
std::ofstream myfile7;
std::ofstream myfile8;
    myfile.open (filenamed);
    myfile2.open (filenameg);
    myfile3.open (filenamef);
    myfile4.open (filenametime);
    myfile5.open (filenametimenew);
    myfile6.open (filenameparam); 
    myfile7.open (filenameft); 
    myfile8.open (filenamesp);

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::normal_distribution<> n(0,1);
    const double pi = 3.14159265358979323846;   
    std::uniform_real_distribution<> dis(0, 2*pi);
    //double I0=-1; // 0
    //double dI=0.2;
    double gL=0.3; // 8
    double EL=-80; // -80  
    double gNa=1; // 20
    double ENa=60; // 60
    double gK=0.4; // 9
    double EK=-90; // -90
    double gM=5; // 5
    //Na current
    double km=14; // 15 
    double vm=-18; // -20 
    //fast K+ current
    double kn=5; // 5 
    double vn=-25; // -20 
    double tau=3; // 0.152  
    int N=120000000; //12*10⁷
    int Neq=20000000; //2*10⁷
    
    int points=100000;
    int timepoints=10000;
    int j;
    double dt=0.0005;//0.01;
    
    //int spikes=100;
    //double* spike=new double[spikes];
    double th,thn;
    int runs=500; //250
    int repetitions=20; // 20 
    int sampling=((N-Neq)/points)*repetitions;
    
    
    double epsilon=0.01;//0.001
    double omega=0.001;

    int fpoints=1000000;
    int fsampling=((N-Neq)/fpoints)*repetitions;
    double Tf=(N-Neq)*repetitions*dt;
    double dtint=fsampling*dt;
    double deltapeak=1/dtint;
    double *sqout=new double[fpoints];
    
    int spikecount;
    double *spikeout=new double[fpoints];
    initd(sqout,fpoints);
    int q;
    double tf,phi;

//fftw_complex in[points], out[points];
fftw_complex* in = new fftw_complex[fpoints];
fftw_complex* out = new fftw_complex[fpoints];

fftw_complex* inspike = new fftw_complex[fpoints];
fftw_complex* outspike = new fftw_complex[fpoints];
    // simulate equilibrium point
    int Nsim=5000000;
    double vsim,nfsim,nfsims,vsims,neq,veq;
    nfsims=0.4;
    vsims=-70;
    for(int jsim=0;jsim<Nsim;jsim++){ 
      vsim=I*dt+vsims-gL*(vsims-EL)*dt-gNa*ninf(vsims,km,vm)*(vsims-ENa)*dt-gK*nfsims*(vsims-EK)*dt;
      nfsim=nfsims+(ninf(vsims,kn,vn)-nfsims)*dt/tau;
      vsims=vsim;
      nfsims=nfsim;
    }
    veq=vsims;
    neq=nfsims;
    //simulate instable focus
    int Nfoc=20000000;
    double dtfoc,vfoc,nffoc,nffocs,vfocs;
    dtfoc=-dt;
    nffocs=0.5;
    vfocs=-25;
    for(int jfoc=0;jfoc<Nfoc;jfoc++){
      vfoc=I*dtfoc+vfocs-gL*(vfocs-EL)*dtfoc-gNa*ninf(vfocs,km,vm)*(vfocs-ENa)*dtfoc-gK*nffocs*(vfocs-EK)*dtfoc;
      nffoc=nffocs+(ninf(vfocs,kn,vn)-nffocs)*dtfoc/tau;
      vfocs=vfoc;
      nffocs=nffoc;
    }
    if(vfocs<=0 && vfocs >=-100 && nffocs <= 1 && nffocs >=0){
      th=vfocs;
      thn=nffocs;
    }
    else{
      th=-25;
      thn=0.7;
    }
long long int count;
double Nav,Nava,Nava2,Navaalt[points],N2av[points],N2ava,T;
double mx,nx,nx2,v,nf,nfs,vs,Deff,vav,vavalt,sigmav,Fano,Navz[points],sqNavz[points],state[timepoints],state2[timepoints];
int avalue[timepoints],avalue2[timepoints],lasta,tvcount,j2value[timepoints],j2value2[timepoints],lastj2,jrvalue[timepoints],jrvalue2[timepoints],lastjr;
avalue[0]=0;
j2value[0]=0;
jrvalue[0]=0;
avalue2[0]=0;
j2value2[0]=0;
jrvalue2[0]=0;
lasta=0;lastjr=0;lastj2=0;tvcount=0;
initd(Navz,points);
initd(N2av,points);
initd(Navaalt,points);
int bv,bv2,lastbv,lastbv2,sv,svlocked,t,t2,p,qvar,pvar;
Nava=0;Nava2=0;sv=0;t=0;t2=0;nfs=0;vs=0;svlocked=0;
mx=0.0309-0.0047*I;
nx=1.714+0.01*I;
nx2=mx*63;
for(int a=0;a<runs;a++){ 
  for(j=0;j<Neq;j++){
     tf=dt*a*Neq+j*dt; 
    v=I*dt+vs-gL*(vs-EL)*dt-gNa*ninf(vs,km,vm)*(vs-ENa)*dt-gK*nfs*(vs-EK)*dt+sqrt(2*D*dt)*n(gen)+epsilon*cos(2*pi*omega*tf+phi)*dt;
    nf=nfs+(ninf(vs,kn,vn)-nfs)*dt/tau;
    vs=v;
    nfs=nf;
  }
}
if(nfs>mx*vs+nx2){
  bv=0;
  bv2=0;
}
else{
  bv=1;
  bv2=1;
}
state[0]=bv;state2[0]=bv2;lastbv=bv;lastbv2=bv2;
for(int a=0;a<runs;a++){
  count=0;
  spikecount=0;
  q=0;
  qvar=0;
  p=0;
  pvar=0;
  phi=dis(gen);
  for(int jr=0;jr<repetitions;jr++){
  for(int j2=0;j2<N-Neq;j2++){
    tvcount=tvcount+1;
    pvar=pvar+1;
    qvar=qvar+1;
   tf=dt*(N-Neq)*(a*repetitions+jr)+dt*j2;
    v=I*dt+vs-gL*(vs-EL)*dt-gNa*ninf(vs,km,vm)*(vs-ENa)*dt-gK*nfs*(vs-EK)*dt+sqrt(2*D*dt)*n(gen)+epsilon*cos(2*pi*omega*tf+phi)*dt;
    nf=nfs+(ninf(vs,kn,vn)-nfs)*dt/tau;
    if(v > th && vs < th){
      sv=1;
    }
    if(sv==1 && nf > thn && nfs < thn){
      count=count+1;
      sv=0;
      bv2=1;
      svlocked=0;
      lasta=a;
      lastjr=jr;
      lastj2=j2;
    }
    if(bv2==1 && v<veq+5 && nf<neq+0.1){
      if(svlocked==0){
        if(v < veq && vs > veq){
          svlocked=1;
        }}
      if(svlocked==0){
        if(nf < neq && nfs > neq){
          svlocked=-1;
        }
      }
      if(svlocked==1){
        if(nf < neq && nfs > neq){
          bv2=0;
          svlocked=2;
        }
      }
      if(svlocked==-1){
        if(v < veq && vs > veq){
          bv2=0;
          svlocked=2;
        }
      }
    }
    if((bv==1 && nf>mx*v+nx2) || (bv==1 && (((a-lasta)*repetitions+(jr-lastjr))*dt*(N-Neq)+(j2-lastj2)*dt)>20)){
        bv=0;
    }
    if(bv==0 && nf<mx*v+nx){
	    bv=1;  
    }
    if(tvcount==40000){
      tvcount=0;
    if((t<timepoints) && (0.5<bv+lastbv) && (1.5>bv+lastbv)){
    state[t]=bv;
    avalue[t]=a;
    jrvalue[t]=jr;
    j2value[t]=j2;
    t=t+1;
    }
    lastbv=bv;
    }
    if((t2<timepoints) && (0.5<bv2+lastbv2) && (1.5>bv2+lastbv2)){
      state2[t2]=bv2;
      avalue2[t2]=a;
      jrvalue2[t2]=jr;
      j2value2[t2]=j2;
      t2=t2+1;
      lastbv2=bv2;
    }
    vs=v;
    nfs=nf;
    T=2*(jr*dt*(N-Neq)+j2*dt)+dt;
    if(pvar==sampling){
	Navaalt[p]=(a*Navaalt[p]+2*count/T)/(a+1);
	Navz[p]=(a*Navz[p]+count/sqrt(T))/(a+1);
    	N2av[p]=(N2av[p]*a+(count*count)/T)/(a+1);
    	p=p+1;
	pvar=0;
    }
    if(qvar==fsampling){
	    qvar=0;
    in[q][REAL]=v;
    in[q][IMAG]=0;
    if(count>spikecount+0.5){
      inspike[q][REAL]=deltapeak*(count-spikecount);
      inspike[q][IMAG]=0;
      spikecount=count;
    }
    else{
      inspike[q][REAL]=0;
      inspike[q][IMAG]=0;
    }
    q=q+1;
    }
  }
  }
  Nava=(a*Nava+2*count/T)/(a+1);
  Nava2=(a*Nava2+(2*count/T)*(2*count/T))/(a+1);
  fft(in, out,fpoints);  
for(int ind=0;ind<fpoints;ind++){
    sqout[ind]=(a*sqout[ind]+(out[ind][REAL]*out[ind][REAL]+out[ind][IMAG]*out[ind][IMAG])*dtint*dtint)/(a+1);
}
fft(inspike, outspike,fpoints);  
for(int ind=0;ind<fpoints;ind++){
  spikeout[ind]=(a*spikeout[ind]+(outspike[ind][REAL]*outspike[ind][REAL]+outspike[ind][IMAG]*outspike[ind][IMAG])*dtint*dtint)/(a+1);
}
} 
for(int qua=0;qua<points;qua++){
sqNavz[qua]=Navz[qua]*Navz[qua];
}
N2ava=average(N2av,points);
Nav=average(sqNavz,points);
Deff=(N2ava-Nav);
vav=Nava;
vavalt=average(Navaalt,points);
sigmav=(Nava2-Nava*Nava);
Fano=2*Deff/vav;   
//} 
myfile << I << " " << Deff << "\n";

myfile.close(); 

myfile2 << I << " " << vav << " " << vavalt << " "<< sigmav << "\n";

myfile2.close();     

myfile3 << I << " " << Fano << "\n";
myfile3.close(); 

for(int rv=0;rv<timepoints;rv++){
  myfile4 << avalue[rv] << " " << jrvalue[rv] << " " << j2value[rv] <<" " <<state[rv] << "\n";
}
myfile4.close();
for(int rv=0;rv<timepoints;rv++){
  myfile5 << avalue2[rv] << " " << jrvalue2[rv] << " " << j2value2[rv] <<" " <<state2[rv] << "\n";
}
myfile5.close();
myfile6 << "N" << " " << "Neq" << " " << "dt" <<" " << "runs" << " "<< "repetitions" << " " << "neq" << " " << "veq" << " "<< "thn" << " " << "th" << " "<< "epsilon" << " " << "omega" << "\n";
myfile6 << N <<" " << Neq << " "<< dt << " " << runs << " " << repetitions << " " << neq << " " << veq << " " << thn << " " << th << " "<< epsilon << " "<< omega << "\n";
myfile6.close();

for(int rv=0;rv<fpoints;rv++){
  myfile7 << rv/Tf << " " << sqout[rv] <<"\n";
}
myfile7.close();

for(int rv2=0;rv2<fpoints;rv2++){
  myfile8 << rv2/Tf << " " << spikeout[rv2] <<"\n";
}
myfile8.close();

delete [] in;
delete [] out;
delete [] sqout;
delete [] inspike;
delete [] outspike;
delete [] spikeout;

return 0;
}
double ninf(double V, double k, double Vh){
    double f=1/(1+exp((Vh-V)/k));
    return f;
}

double average(double x[],int n){ 
    int i;
    double sum=0;
    int e2=0;
    for(i=0;i<n;i++){
        sum=(sum*e2+x[i])/(e2+1);
	e2=e2+1;
    }
    return sum;
}

void square(long long int x[],int n){
    int j;
    for(j=0;j<n;j++){
        x[j]=x[j]*x[j];
    }
}

void initd(double x[],int n){
    for(int k2=0;k2<n;k2++){
        x[k2]=0;
    }
}
