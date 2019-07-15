/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: richard
 *
 * Created on December 29, 2018, 9:27 PM
 */

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>

using namespace std;

/*
 * 
 */
double average(double x[],int n);
double square(double x[],int n);
double init(double x[],int n);

int main(int argc, char** argv) {
    
    std::random_device rd;
    std::mt19937 gen(rd(10));
    std::normal_distribution<> n(0,1);
    
    /*double filename=atof(argv[1]);
    //filename << std::setprecision(3);
    ostringstream filename_convert;
    filename_convert << filename;
    string filename_str = filename_convert.str();
    filename_str="aha"+filename_str;
    filename_str=filename_str+".txt";
    ofstream myfile;
    myfile.open(filename_str);*/
    
    double ra[10];
    /*double vec[5] = { 1.4, 2.2, 3.9, 4.5, 5};
    int s=sizeof(vec)/sizeof(vec[3]);
    double x=average(vec,s);
    square(vec,s);
    double z=average(vec,s);
    init(vec,s);
    double a=average(vec,s);   
    double testets=atof(argv[2]);
    std::cout << x << " " << z << " " << a << testets << "\n";
    int start=1;
    double* test=new double[start];
    double var[]={1,2,3,4,5};
    int i;
    for(i=0;i<5;i++){
        int length=sizeof(test)/sizeof(test[0]);
    if (i >= length) {
        start= start*2;            // double the previous size
        double* temp= new double[start]; // create new bigger array.
        for (int j=0; j<i; j++) {
            temp[j] = test[j];       // copy values to new array.
        }
        delete [] test;              // free old array memory.
        test = temp;                 // now a points to new array.
    }              
    if(var[i] < 3 && var[i] > 1){
            test[i]=0.1;
    }
    else if(var[i] < 3 && var[i] <= 1){
            test[i]=-0.9;
    }
    else{
        test[i]=1.1;
    }
}
    int k=0;
    while(test[k]!=0){
        k++;
    }
    double testfinal[k]; // create array of right size.
        for (int j=0; j<k; j++) {
            testfinal[j] = test[j];       // copy values to new array.
        }
        delete [] test;              // free old array memory. 
    int k2=sizeof(testfinal)/sizeof(testfinal[0]);    
    double ave=average(testfinal,k);
    double dividend=0;
    double divisor=21766;
    double test123=(0.0+0)/21766+0;*/
    for (int j=0; j<10; j++) {
            ra[j] = n(gen)[j];       // copy values to new array.
        }
    
    
    std::cout << ra[5];
    /*myfile << filename;
    myfile.close();*/
    return 0;
}

double average(double x[],int n){ 
    int i;
    double sum=0;
    for(i=0;i<n;i++){
        sum=sum+x[i];
    }
    double av=sum/n;
    return av;
}

double square(double x[],int n){
    int j;
    for(j=0;j<n;j++){
        x[j]*=x[j];
    }
}

double init(double x[],int n){
    int k;
    for(k=0;k<n;k++){
        x[k]=0;
    }
}
