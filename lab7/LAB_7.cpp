#include<iostream>
#include<cmath>
#include<string>
#include<cstdlib>
#include"../lab2/LAB_2.h"

using namespace std;

double ax(double t){
    return sin(t)/(pow(t,0.5)+1);
}

double ay(double t){
    return log(t+1)/(t+1);
}

double dx(double t,double c){
    return (c-t)*sin(t)/(pow(t,0.5)+1);
}

double dy(double t,double c){
    return (c-t)*log(t+1)/(t+1);
}

int near(double t){
    if(fabs(t-int(t))>0.5){
        return int(t)+1;
    }
    else return int(t);
}

class Romberg_locomotion{
    public:
    Romberg_locomotion(){}
    double* Ax;
    double* Ay;
    double* Vx;
    double* Vy;
    double* Dx;
    double* Dy;
    double* points;
    double error;// to be initialized
    int n_points;// to be initialized
    int M;// to be initialized
    typedef double(*func)(double);
    func f1;
    func f2;
    typedef double (*func_d)(double,double);
    func_d f1_d;
    func_d f2_d;

    double romberg(func f,double* X_points,int n_count){
        double** R=new double*[M];
        for(int i=0;i<M;i++){
            R[i]=new double[M];
        }
        double h=X_points[n_count-1]-X_points[0];
        R[0][0]=(f(X_points[0])+f(X_points[n_count-1]))*h/2;
        if(n_count==1){
            return f(X_points[0])*X_points[0];
        }
        for(int k=1;k<M;k++){
            double f_counts=0;
            for(int i=0;i<pow(2,k-1);i++){
                f_counts+=f((2*i+1)*h*pow(0.5,k));
            }
            R[k][0]=(R[k-1][0]+h*pow(0.5,k-1)*f_counts)/2;
            for(int j=1;j<=k;j++){
                R[k][j]=R[k][j-1]+(R[k][j-1]-R[k-1][j-1])/(pow(4.0,j)-1);
                if(fabs(R[k][k]-R[k-1][k-1])<error || (k==M-1 && j==M-1)){
                    return R[k][k];
                }
            }
        }
        
        cout<<endl;
    }

    double romberg(func_d f,double* X_points,int n_count){
        double** R=new double*[M];
        for(int i=0;i<M;i++){
            R[i]=new double[M];
        }
        double h=X_points[n_count-1]-X_points[0];
        R[0][0]=(f(X_points[0],X_points[n_count-1])+f(X_points[n_count-1],X_points[n_count-1]))*h/2;
        if(n_count==1){
            return f(X_points[0],X_points[n_count-1])*X_points[0];
        }
        for(int k=1;k<M;k++){
            double f_counts=0;
            for(int i=0;i<pow(2,k-1);i++){
                f_counts+=f((2*i+1)*h*pow(0.5,k),X_points[n_count-1]);
            }
            R[k][0]=(R[k-1][0]+h*pow(0.5,k-1)*f_counts)/2;
            for(int j=1;j<=k;j++){
                R[k][j]=R[k][j-1]+(R[k][j-1]-R[k-1][j-1])/(pow(4.0,j)-1);
                if(fabs(R[k][k]-R[k-1][k-1])<error || (k==M-1 && j==M-1)){
                    return R[k][k];
                }
            }
        }
    }
    void init(){
        points=new double[n_points];
        Ax=new double[n_points];
        Ay=new double[n_points];
        Vx=new double[n_points];
        Vy=new double[n_points];
        Dx=new double[n_points];
        Dy=new double[n_points];
        for(int i=0;i<n_points;i++){
            points[i]=0.1*(i+1);
            Ax[i]=f1(points[i]);
            Ay[i]=f2(points[i]);
        }
        return;
    }

    void V_romberg(){
        for(int i=0;i<n_points;i++){
            Vx[i]=romberg(ax,points,i+1);
            Vy[i]=romberg(ay,points,i+1);
        }
        return;
    }

    void D_romberg(){
        for(int i=0;i<n_points;i++){
            Dx[i]=romberg(dx,points,i+1);
            Dy[i]=romberg(dy,points,i+1);
        }
        return;
    }

    void solve(){
        init();
        V_romberg();
        D_romberg();
        cout<<"Vx: ";
        printList(Vx,n_points);
        cout<<"Vy: ";
        printList(Vy,n_points);
        cout<<"Dx: ";
        printList(Dx,n_points);
        cout<<"Dy: ";
        printList(Dy,n_points);
        return;
    }
};

int main(){
    Romberg_locomotion test;
    test.n_points=100;
    test.M=8;
    test.error=pow(10,-6);
    test.f1=ax;
    test.f2=ay;
    test.f1_d=dx;
    test.f2_d=dy;
    test.solve();
    return 0;
}