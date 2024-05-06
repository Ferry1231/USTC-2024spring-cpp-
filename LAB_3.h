#ifndef LAB_3_H
#define LAB_3_H

#include"LAB_2.h"
#include<iostream>
#include<cmath>
#include<string>

using namespace std;

class inversePower{
    public:
    double feature_min;
    double* feature_vec;
    int n;
    double** A;
    double** L;
    double** U;
    void readMatrix(string filename){//ok
        ifstream matrix(filename);
        int count=0;
        string line;
        while(getline(matrix,line)){
            ++count;
        }
        n=count;
        A=new double*[n];
        for(int i=0; i<n; i++){
            A[i]=new double[n];
        }
        matrix.clear();
        matrix.seekg(0);
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                matrix>>A[i][j];
            }
        }
        matrix.close();
    }
    double* solveL(double* b){//ok
        double* y=new double[n];
        for(int i=0;i<n;i++){
            y[i]=b[i];
            for(int j=0;j<i;j++){
                y[i]-=y[j]*L[i][j];
            }
        }
        return y;
    }
    double* solveU(double* y){//ok
        double* x=new double[n];
        for(int i=n-1;i>=0;i--){
            x[i]=y[i]/U[i][i];
            for(int j=i+1;j<n;j++){
                x[i]-=x[j]*U[i][j]/U[i][i];
            }
        }
        return x;
    }
    double* Doolittle(double* y_first){//ok
        L=new double*[n];
        U=new double*[n];
        for(int i=0;i<n;i++){
            L[i]=new double[n];
            U[i]=new double[n];
            for(int j=0;j<n;j++){
                L[i][j]=0;
                U[i][j]=0;
            }
        }
        for(int k=0;k<n;k++){
            for(int j=k;j<n;j++){
                U[k][j]=A[k][j];
                for(int i=0;i<k;i++){
                    U[k][j]-=L[k][i]*U[i][j];
                }
                L[j][k]+=A[j][k]/U[k][k];
                for(int i=0;i<k;i++){
                    L[j][k]-=L[j][i]*U[i][k]/U[k][k];
                }
            }
        }
        double* y=solveL(y_first);
        double* x=solveU(y);
        printMatrix(U,n);
        cout<<endl;
        printMatrix(L,n);
        cout<<endl;
        return x;
    }
    void solve(string filename){//ok
        readMatrix(filename);
        double* x_first=new double[n];
        double* x_next=new double[n];
        double* null=new double[n];
        double* y_first=new double[n];
        double lambda_first,lambda_next;
        int k=0;
        feature_vec=new double[n];
        for(int i=0;i<n;i++){
            null[i]=0;
            x_first[i]=1;
        }
        do{
            k+=1;
            if(x_next[0]!=0){
                lambda_first=y_first[0]/x_next[0];
            }
            else{
                lambda_first=0;
            }
            for(int i=0;i<n;i++){
                y_first[i]=x_first[i]/mis(x_first,null,n);
            }
            x_next=Doolittle(y_first);
            for(int i=0;i<n;i++){
                if(x_next[i]!=0){
                    lambda_next=y_first[i]/x_next[i];
                }
            }
            cout<<"Epoch"<<k<<": "<<endl;
            cout<<"Feature_min: "<<lambda_next<<endl;
            cout<<"x-"<<k<<": ";
            printList(x_first,n);
            cout<<"y-"<<k<<": ";
            printList(y_first,n);
            cout<<endl;
            if(fabs(lambda_first-lambda_next)<pow(10,-5) && fabs(lambda_first-lambda_next)>0) break;
            else{
                for(int i=0;i<n;i++){
                    x_first[i]=x_next[i];
                }
                lambda_first=lambda_next;
            }
        }
        while(fabs(lambda_first-lambda_next)>=pow(10,-5) || fabs(lambda_first-lambda_next)==0);
        feature_min=lambda_next;
        for(int i=0;i<n;i++){
            feature_vec[i]=y_first[i];
        }

    }
    void printFeature(){
        cout<<"feature_vec: ";
        printList(feature_vec,n);
        cout<<"feature_min: "<<feature_min<<endl;
    }
};
#endif
