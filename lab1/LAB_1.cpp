#include <iostream>
#include <cmath>

using namespace std;

struct Coefficients_14{
    double c4,c3,c2,c1,c0;
};

struct Coordinates{
    double x,y;
};

double innerProduct2d(const Coordinates& P,const Coordinates& Q){
    double result=P.x*Q.x+P.y*Q.y;
    return result; 
}

double f(Coordinates coef,Coordinates& P,Coordinates& Q){
    return 2*innerProduct2d(P,Q)*coef.x*coef.y+2*innerProduct2d(Q,Q)*pow(coef.y,2)-coef.y-1+coef.x;
}
double g(Coordinates coef,Coordinates& P,Coordinates& Q){
    return innerProduct2d(P,P)*pow(coef.x,2)+2*innerProduct2d(P,Q)*coef.x*coef.y+innerProduct2d(Q,Q)*pow(coef.y,2)-1;
}
double f_x(Coordinates coef,Coordinates& P,Coordinates& Q){
    return 2*innerProduct2d(P,Q)*coef.y+1;
}
double f_y(Coordinates coef,Coordinates& P,Coordinates& Q){
    return 2*innerProduct2d(P,Q)*coef.x+4*innerProduct2d(Q,Q)*coef.y-1;
}
double g_x(Coordinates coef,Coordinates& P,Coordinates& Q){
    return 2*innerProduct2d(P,P)*coef.x+2*innerProduct2d(P,Q)*coef.y;
}
double g_y(Coordinates coef,Coordinates& P,Coordinates& Q){
    return 2*innerProduct2d(Q,Q)*coef.y+2*innerProduct2d(P,Q)*coef.x;
}

Coordinates solve_matrix(double a1,double b1,double a2,double b2,double c1,double c2){//ok
    Coordinates Solution;
    Solution.x=(c2*b1-c1*b2)/(b1*a2-b2*a1);
    Solution.y=(c1*a2-c2*a1)/(b1*a2-b2*a1);
    return Solution;
}

Coordinates solveQuadraticEquation_NewtonMethod(Coordinates coef,Coordinates& P,Coordinates& Q){//ok
    Coordinates coef_delta;
    int i=0;
    do{
        coef_delta=solve_matrix(f_x(coef,P,Q),f_y(coef,P,Q),g_x(coef,P,Q),g_y(coef,P,Q),-f(coef,P,Q),-g(coef,P,Q));
        coef.x+=coef_delta.x; 
        coef.y+=coef_delta.y;
        i++;
    }
    while((coef_delta.x>=0.001 || coef_delta.y>=0.001) && i<10);
    return coef;
}

Coordinates reflectedPoint(Coordinates& T,const Coordinates& P,const Coordinates& Q,Coordinates& R,const double& x,const double& y){
    T.x=x*P.x+y*Q.x;
    T.y=x*P.y+y*Q.y;
    R=solve_matrix(T.x,T.y,T.y,-T.x,2-Q.x*T.x-Q.y*T.y,Q.x*T.y-Q.y*T.x);
    return R;
}

int main(){
    /*Coordinates P,Q,R,T,Solution;
    P.x=-1024;
    P.y=0;
    Q.x=-8;
    Q.y=4;
    Solution.x=1;
    Solution.y=0;
    Solution=solveQuadraticEquation_NewtonMethod(Solution,P,Q);
    cout<<"x:"<<Solution.x<<"  y:"<<Solution.y<<endl;
    reflectedPoint(T,P,Q,R,Solution.x,Solution.y);
    cout<<"T:("<<T.x<<","<<T.y<<")"<<",";
    cout<<"R:("<<R.x<<","<<R.y<<")";
    cout<<endl;*/
    Coordinates temp=solve_matrix(0.002,87.13,4.453,-7.26,87.15,37.27);
    cout<<temp.x<<" "<<temp.y<<endl;
    return 0;
}