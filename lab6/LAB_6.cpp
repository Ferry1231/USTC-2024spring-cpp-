#include<iostream>
#include<complex>
#include<cmath>
#include<typeinfo>
#include<string>
#include<cstdlib>

#define PI acos(-1)

using namespace std;


double f1(double t){
    return 0.7 * sin(4 * PI * t.get_d()) + sin(10 * PI * t.get_d());
}

double f2(double t){
    return 0.7 * cos(4 * PI * t.get_d()) + sin(10 * PI * t.get_d()) + 0.3 * (rand() % 2);
}

class FFT{
public:
    FFT(){}
    int n;
    typedef double(*func)(double);
    func f;
    complex<double> two;
    complex<double>* fn;
    complex<double>* fn_1;
    complex<double>* gn;

    void init(){
        two = complex<double>(2.0, 0.0);
        fn = new complex<double>[n];
        gn = new complex<double>[n];
        fn_1 = new complex<double>[n];
        for(int i = 0; i < n; i++){
            fn[i] = complex<double>(f(double(i) / n), 0.0);
        }
    }

    complex<double>* FT(complex<double>* Fn, int num){
        if(num == 1){
            return Fn;
        }
        complex<double> wn(cos(2 * PI / num), sin(-2 * PI / num));
        complex<double> w(1.0, 0.0);
        complex<double>* f0 = new complex<double>[num / 2];
        complex<double>* f1 = new complex<double>[num / 2];
        complex<double>* g = new complex<double>[num];
        for(int i = 0; i < num / 2; i++){
            f0[i] = Fn[2 * i];
            f1[i] = Fn[2 * i + 1];
        }
        complex<double>* g0 = FT(f0, num / 2);
        complex<double>* g1 = FT(f1, num / 2);
        for(int k = 0; k < num / 2; k++){
            g[k] = (g0[k] + w * g1[k]) / two;
            g[k + num / 2] = (g0[k] - w * g1[k]) / two;
            w = w * wn;
        }
        if(num == n){
            for(int i = 0; i < num; i++){
                gn[i] = g[i];
            }
        }
        return g;
    }

    complex<double>* IFT(complex<double>* Gn, int num){
        if(num == 1){
            return Gn;
        }
        complex<double> wn(cos(2 * PI / num), sin(-2 * PI / num));
        complex<double> w(1.0, 0.0);
        complex<double>* f0 = new complex<double>[num / 2];
        complex<double>* f1 = new complex<double>[num / 2];
        complex<double>* g = new complex<double>[num];
        for(int i = 0; i < num / 2; i++){
            f0[i] = Gn[2 * i];
            f1[i] = Gn[2 * i + 1];
        }
        complex<double>* g0 = IFT(f0, num / 2);
        complex<double>* g1 = IFT(f1, num / 2);
        for(int k = 0; k < num / 2; k++){
            g[k] = (g0[k] + w * g1[k]);
            g[k + num / 2] = (g0[k] - w * g1[k]);
            w = w * wn;
        }
        if(num == n){
            for(int i = 0; i < num; i++){
                fn_1[i] = g[i];
            }
        }
        return g;
    }

    void print(){
        cout << "F: ";
        for(int i = 0; i < n; i++){
            cout << abs(fn[i]) << " ";
        }
        cout << endl;
        cout << "G: ";
        for(int i = 0; i < n; i++){
            cout << abs(gn[i]) << " ";
        }
        cout << endl;
        cout << "F_1: ";
        for(int i = 0; i < n; i++){
            cout << abs(fn_1[i]) << " ";
        }
    }
};

int main(){
    FFT test;
    test.f = f1;
    test.n = int(pow(2, 4));
    test.init();
    test.FT(test.fn, test.n);
    test.IFT(test.gn, test.n);
    test.print();
    return 0;
}
