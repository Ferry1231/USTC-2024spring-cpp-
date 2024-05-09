#ifndef LAB_2_H
#define LAB_2_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>

using namespace std;

double** createMatrix(double epi,double a,int n){
    double h=1.00/n;
    //allocate memory for matrix
    double **Matrix=new double*[n];
    for(int i=0;i<n;i++){
        Matrix[i]=new double[n];
    }
    //assign values
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(i==j) Matrix[i][j]=-2*epi-h;
            else if(j==i+1) Matrix[i][j]=epi+h;
            else if(i==j+1) Matrix[i][j]=epi;
            else Matrix[i][j]=0;
        }
    }
    return Matrix;
}

double** matrixSub(double**m1,double v1,int n){
    double** result=new double*[n];
    for(int i=0;i<n;i++){
        result[i]=new double[n];
    }
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            result[i][j]=m1[i][j];
            if(i==j){
                result[i][j]-=v1;
            }
        }
    }
    return result;
}

double** matrixTimes(double**m1,double** m2,int n){
    double** result=new double*[n];
    //allocate memory
    for(int i=0;i<n;i++){
        result[i]=new double[n];
    }
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            result[i][j]=0;
        }
    }
    //compute
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            for(int k=0;k<n;k++){
                result[i][k]+=m1[i][j]*m2[j][k];
            }
        }
    }
    return result;
}

double* matrixTimesList(double** m,double* l,int n){
    double* result=new double[n];
    for(int i=0;i<n;i++){
        result[i]=0;
    }
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            result[i]+=m[i][j]*l[j];
        }
    }
    return result;
}

int findMax(double* list){
    if(fabs(list[0])>fabs(list[1])) return 0;
    else return 1;
}

void printMatrix(double** M,int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cout<<M[i][j]<<"\t";
        }
        cout<<endl;
    }
}

void printMatrix(double** M,int n,int m){
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            cout<<M[i][j]<<"\t";
        }
        cout<<endl;
    }
}

void printList(double* list,int n){
    for(int i=0;i<n;i++){
        cout<<list[i]<<",";
    }
    cout<<endl;
}

double mis(double* l1,double* l2,int n){
        double* temp=new double[n];
        for(int i=0;i<n;i++){
            temp[i]=fabs(l1[i]-l2[i]);
        }
        double max=0;
        for(int i=0;i<n-1;i++){
            if(temp[i]>temp[i+1]){
                double temp_num;
                temp_num=temp[i];
                temp[i]=temp[i+1];
                temp[i+1]=temp_num;
            }
        }
        max=temp[n-1];
        return max;
    }

class GaussColumn{
    public:
    GaussColumn(){}
    double epi,a;
    int n;
    double** Matrix;
    double* result;
    double* real;
    double** columnFirst(int i_col){
        double* temp_list=new double[2];
        for(int i=i_col;i<n;i++){
            temp_list[i-i_col]=Matrix[i][i_col];
        }
        int temp_max=findMax(temp_list);
        delete [] temp_list;
        if(temp_max+i_col==i_col){
            return Matrix;
        }
        else{
            double** M_temp=new double*[n];
            for(int i=0;i<n;i++){
                M_temp[i]=new double[n];
            }
            for(int i=0;i<n;i++){
                if(i==i_col){
                    M_temp[i][temp_max+i_col]=1;
                    M_temp[i][i]=0;
                }
                else if(i==temp_max+i_col){
                    M_temp[temp_max+i_col][i_col]=1;
                    M_temp[i][i]=0;
                }
                else{
                    M_temp[i][i]=1;
                }
            }
            Matrix=matrixTimes(M_temp,Matrix,n);
            return Matrix;
        }
    }
    double real_solution(int i,double h){
        double x=i*h;
        return (1-a)*(1-exp(-x/epi))/(1-exp(-1/epi))+a*x;
    }
    void solve(){
        double h=1.00/n;
        //Matrix=createMatrix(epi,a,n);
        result=new double[n];
        real=new double[n];
        for(int i=1;i<=n;i++){
            result[i-1]=0;
            double x=i*h;
            real[i-1]=(1-a)*(1-exp(-x/epi))/(1-exp(-1/epi))+a*x;
        }
        for(int i=0;i<n-1;i++){
            Matrix=columnFirst(i);
            double temp=Matrix[i+1][i];
            for(int j=i;j<n;j++){
                Matrix[i+1][j]-=Matrix[i][j]*temp/Matrix[i][i];
            }
        }
        for(int i=n-1;i>=0;i--){
            int j=i;
            if(j+1<n){
                for(int k=j+1;k<n;k++){
                    result[j]-=result[k]*Matrix[j][k]/Matrix[j][j];
                }
                result[j]+=a*pow(h,2)/Matrix[j][j];
            }
            else if(j+1==n){
                result[j]=real_solution(n,h);
            }
        }
    }
    void print(){
        for(int i=0;i<n;i++){
            cout<<"y_"<<i+1<<" : "<<result[i]<<"  "<<real[i]<<"  ";
            if(i%5==0){
                cout<<endl;
            }
        }
        
    }
};

class GaussSeidel{
    public:
    GaussSeidel(){}
    double epi,a;
    int n;
    double** Matrix;
    double* b;
    double* result;
    double* real;
    double** reverse_DL(double** Matrix){
        double h=1.00/n;
        double** reverse=new double*[n];
        for(int i=0;i<n;i++){
            reverse[i]=new double[n];
        }
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                if(i<j){
                    reverse[i][j]=0;
                }
                else{
                    reverse[i][j]=pow((epi)/(2*epi+h),i-j);
                    reverse[i][j]/=(2*epi+h);
                }
            }
        }
        return reverse;
    }
    double real_solution(int i,double h){
        double x=i*h;
        return (1-a)*(1-exp(-x/epi))/(1-exp(-1/epi))+a*x;
    }
    void solve(){
        double* x_first=new double[n];
        double* x_next=new double[n];
        double** DL=new double*[n];
        double** U=new double*[n];
        real=new double[n];
        Matrix=createMatrix(epi,a,n);
        result=new double[n];
        b=new double[n];
        double h=1.00/n;
        for(int i=0;i<n;i++){
            x_first[i]=(i+1)/10.00;
            b[i]=a*pow(h,2);
            if(i==n-1){
                x_first[i]=1;
                b[i]=epi*real_solution(n-1,h)-(2*epi+h)*real_solution(n,h);
            } 
            DL[i]=new double[n];
            U[i]=new double[n];
            real[i]=real_solution(i+1,h);
        }
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                if(i>=j){
                    DL[i][j]=Matrix[i][j];
                    U[i][j]=0;
                }
                else{
                    DL[i][j]=0;
                    U[i][j]=Matrix[i][j];
                }
            }
        }
        double** S=matrixTimes(reverse_DL(DL),U,n);
        double* f=matrixTimesList(reverse_DL(DL),b,n);
        for(int i=0;i<n;i++){
            f[i]=-f[i];
        }
        for(int i=0;i<n;i++){
            x_next[i]=0;
        }
        while(mis(x_first,x_next,n)>pow(10,-5) || mis(x_first,x_next,n)==0){
            double* x_temp=matrixTimesList(S,x_first,n);
            for(int i=0;i<n;i++){
                x_next[i]=x_temp[i]+f[i];
            }
            if(mis(x_first,x_next,n)<=pow(10,-5) && mis(x_first,x_next,n)>0){
                for(int i=0;i<n;i++){
                    result[i]=x_next[i];
                }
                return;
            }
            for(int i=0;i<n;i++){
                x_first[i]=x_next[i];
            }
            delete [] x_temp;
        }
    }
        void print(){
        for(int i=0;i<n;i++){
            cout<<"y_"<<i+1<<" : "<<setprecision(4)<<result[i]<<"  "<<real[i]<<"  ";
            if(i%5==0){
                cout<<endl;
            }
        }
        
    }
};  

#endif