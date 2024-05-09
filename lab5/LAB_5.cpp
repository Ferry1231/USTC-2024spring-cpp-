#include"./LAB_3.h"
#include"./LAB_2.h"
#include<iostream>
#include<cmath>
#include<cstdlib>
#include<string>
#include<fstream>

using namespace std;

class M_spline{
    public:
    M_spline(){};
    int n;
    double** location;
    double* MList;
    double* dList;
    double* hList;
    double* lambda;
    double* mu;
    double** yCoefs;
    double** matrix;
    void readData(string filename){
        ifstream file(filename);
        string line;
        n=0;
        while(getline(file,line)){
            n++;
        }
        n--;

        location=new double*[n+1];
        for(int i=0;i<n+1;i++){
            location[i]=new double[2];
        }
        hList=new double[n];
        dList=new double[n-1];
        MList=new double[n+1];
        lambda=new double[n-1];
        yCoefs=new double*[n];
        mu=new double[n-1];
        yCoefs=new double*[n];
        for(int i=0;i<n;i++){
            yCoefs[i]=new double[4];
        }

        file.clear();
        file.seekg(0);

        for(int i=0;i<n+1;i++){
            for(int j=0;j<2;j++){
                file>>location[i][j];
            }
        }
        for(int i=0;i<n;i++){
            hList[i]=location[i+1][0]-location[i][0];
        }
        for(int i=0;i<n-1;i++){
            dList[i]=6/(hList[i+1]+hList[i])*((location[i+2][1]-location[i+1][1])/hList[i+1]-(location[i+1][1]-location[i][1])/hList[i]);
        }
        MList[0]=0;
        MList[n]=0;
        for(int i=0;i<n-1;i++){
            lambda[i]=hList[i+1]/(hList[i+1]+hList[i]);
            mu[i]=1-lambda[i];
        }
        for(int i=0;i<n;i++){
            yCoefs[i]=new double[4];
        }
    }
    void makeMatrix(){
        matrix=new double*[n-1];
        for(int i=0;i<n-1;i++){
            matrix[i]=new double[n-1];
        }
        for(int i=0;i<n-1;i++){
            for(int j=0;j<n-1;j++){
                if(i==j){
                    matrix[i][j]=2;
                }
                else if(i==j-1){
                    matrix[i][j]=lambda[i];
                }
                else if(i==j+1){
                    matrix[i][j]=mu[i];
                }
                else{
                    matrix[i][j]=0;
                }
            }
        }
    }
    void solve(string filename){
        readData(filename);
        makeMatrix();
        inversePower solve_M;
        solve_M.A=matrix;
        solve_M.n=n-1;
        double* temp=new double[n-1];
        temp=solve_M.Doolittle(dList);
        for(int i=0;i<n-1;i++){
            MList[i+1]=temp[i];
        }
        delete [] temp;
        for(int i=0;i<n;i++){
            double Mi=MList[i];
            double Mi1=MList[i+1];
            double xi=location[i][0];
            double yi=location[i][1];
            double xi1=location[i+1][0];
            double yi1=location[i+1][1];
            double hi=hList[i];
            
            yCoefs[i][0]=(Mi*pow(xi1,3)-Mi1*pow(xi,3))/(6*hi)+(xi1*yi-xi*yi1)/hi+(xi*Mi1-xi1*Mi)*hi/6;
            yCoefs[i][1]=(Mi1*pow(xi,2)-Mi*pow(xi1,2))/(2*hi)+(yi1-yi)/hi+(Mi*hi-Mi1*hi)/6;
            yCoefs[i][2]=(Mi*xi1-Mi1*xi)/(2*hi);
            yCoefs[i][3]=(Mi1-Mi)/(6*hi);
        }
    }
    void outputData(string filename){
        ofstream file(filename);
        for(int i=0;i<n;i++){
            for(int j=0;j<4;j++){
                file<<yCoefs[i][j];
            }
            file<<"\n";
        }
    }
};