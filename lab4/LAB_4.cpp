#include"../lab3/LAB_3.h"
#include"../lab2/LAB_2.h"
#include<iostream>
#include<cmath>
#include<cstdlib>
#include<string>
#include <fstream>

using namespace std;

// n->column; m->row

double** createRandomM(int n,int m){//ok
    double** matrix=new double*[n];
    for(int i=0;i<n;i++){
        matrix[i]=new double[m];
    }
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            matrix[i][j]=(rand()%100)/100.00;
        }
    }
    return matrix;
}

double** readMatrix(string filename, int rows, int cols) {//ok
    double** matrix = new double*[rows];
    for (int i = 0; i < rows; ++i) {
        matrix[i] = new double[cols];
    }

    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Failed to open file." << endl;
        return nullptr;
    }

    for (int i = 0; i < rows; ++i) {
        string line;
        if (!getline(file, line)) {
            cerr << "Failed to read line " << i << " from file." << endl;
            return nullptr;
        }

        size_t pos = 0;
        for (int j = 0; j < cols - 1; ++j) {
            size_t nextPos = line.find(',', pos);
            if (nextPos == string::npos) {
                cerr << "Failed to find delimiter in line " << i << "." << endl;
                return nullptr;
            }

            matrix[i][j] = stod(line.substr(pos, nextPos - pos));
            pos = nextPos + 1;
        }

        // Read the last value without looking for a comma
        matrix[i][cols - 1] = stod(line.substr(pos));
    }

    file.close();
    return matrix;
}

double** AAT(double** matrix,int n,int m,bool reverse){//ok
    double** result;
    if(reverse==true){
        result=new double*[m];
        for(int i=0;i<m;i++){
            result[i]=new double[m];
        }
        for(int i=0;i<m;i++){
            for(int j=0;j<m;j++){
                for(int k=0;k<n;k++){
                    result[i][j]+=matrix[k][i]*matrix[k][j];
                }
            }
        }
    }
    else if(reverse==false){
        result=new double*[n];
        for(int i=0;i<n;i++){
            result[i]=new double[n];
        }
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                for(int k=0;k<m;k++){
                    result[i][j]+=matrix[i][k]*matrix[j][k];
                }
            }
        }
    }
    return result;
}

double** transposition(double** matrix,int n,int m){//ok
    double** matrix_new=new double*[m];
    for(int i=0;i<m;i++){
        matrix_new[i]=new double[n];
    }
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            matrix_new[i][j]=matrix[j][i];
        }
    }
    return matrix_new;
}

bool compare(double a,double b){//ok
    return(a>b);
}

void swap(double& a,double& b){//ok
    double temp=a;
    a=b;
    b=temp;
}

void swap(int& a,int& b){
    int temp=a;
    a=b;
    b=temp;
}

double countSum(double** A,int n,int m){
    double sum=0;
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            if(i!=j) sum+=pow(A[i][j],2);
        }
    }
    return sum;
}

class JacobiDecom{//ok
    public:
    JacobiDecom(){}
    double** matrix;
    double** result;
    double** eigenvec;
    int n;
    int* findMax(double** A,int n){//ok
        int* res=new int[2];
        double maxValue=0;
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                if(fabs(A[i][j])>fabs(maxValue) && i!=j){
                    maxValue=A[i][j];
                    res[0]=i;
                    res[1]=j;
                }
            }
        }
        return res;
    }
    bool isEnd(double**A,int n){//ok
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                if(i!=j && fabs(A[i][j])>=pow(10,-6)) return false;
            }      
        }
        return true;
    }
    double** solve(){//ok
        result=new double*[n];
        eigenvec=new double*[n];
        for(int i=0;i<n;i++){
            result[i]=new double[n];
            eigenvec[i]=new double[n];
        }
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                result[i][j]=matrix[i][j];
                if(i==j) eigenvec[i][j]=1;
                else eigenvec[i][j]=0;
            }
        }
        while(!isEnd(result,n)){
            int i=findMax(result,n)[0];
            int j=findMax(result,n)[1];
            double** G=new double*[n];
            for(int f=0;f<n;f++){
                G[f]=new double[n];
            }
            double s=(result[j][j]-result[i][i])/(2*(result[i][j]));
            double t;
            if(s>0){
                t=pow(pow(s,2)+1,0.5)-s;
            }
            else if(s<0){
                t=-pow(pow(s,2)+1,0.5)-s;
            }
            else if(s==0){
                t=1;
            }
            for(int e=0;e<n;e++){
                for(int f=0;f<n;f++){
                    if(e==f) G[e][f]=1;
                    else G[e][f]=0;
                }
            }
            G[i][i]=1/pow(pow(t,2)+1,0.5);
            G[j][j]=1/pow(pow(t,2)+1,0.5);
            G[i][j]=t/pow(pow(t,2)+1,0.5);
            G[j][i]=-t/pow(pow(t,2)+1,0.5);
            result=matrixTimes(transposition(G,n,n),result,n);
            result=matrixTimes(result,G,n);
            result[i][j]=0;
            result[j][i]=0;
            eigenvec=matrixTimes(eigenvec,G,n);
            cout<<"sum: "<<countSum(result,n,n)<<endl;
        }
        eigenvec=transposition(eigenvec,n,n);
        for(int i=0;i<n-1;i++){
            for(int t=0;t<n-i-1;t++){
                if(result[t][t]<result[t+1][t+1]){
                    swap(result[t][t],result[t+1][t+1]);
                    for(int j=0;j<n;j++){
                    swap(eigenvec[t][j],eigenvec[t+1][j]);
                    }
                }
            }   
        }
        for(int i=0;i<n;i++){
            double sum=0;
            for(int j=0;j<n;j++){
                sum+=pow(eigenvec[i][j],2);
            }
            sum=pow(sum,0.5);
            for(int j=0;j<n;j++){
                eigenvec[i][j]/=sum;
            }
        }
        eigenvec=transposition(eigenvec,n,n);
        return result;
    }
    void print(){//ok
        for(int i=0;i<n;i++){
            cout<<result[i][i]<<endl;
        }
    }
};

class SVD{//ok
    public:
    SVD(){}
    double* eigenvalues;
    int n;
    int m;
    double** A;
    double** U;
    double** V;
    double** E;//E->Sigma
    double* evalues;
    double** solve_E(){
        E=new double*[n];
        for(int i=0;i<n;i++){
            E[i]=new double[n];
        }
        double** AAt=AAT(A,n,m,0);
        eigenvalues=new double[n];
        JacobiDecom temp_E;
        temp_E.matrix=AAt;
        temp_E.n=n;
        temp_E.solve();
        for(int i=0;i<n;i++){
            eigenvalues[i]=temp_E.result[i][i];
        }
        for(int i=0;i<n;i++){
            E[i][i]=pow(eigenvalues[i],0.5);
            if(i==m-1 && m<n) break;
        }
        return E;
    }
    double** solve_U(){
        JacobiDecom AAt;
        AAt.matrix=AAT(A,n,m,0);
        AAt.n=n;
        AAt.solve();
        U=AAt.eigenvec;
        return U;
    }
    double** solve_V(){
        JacobiDecom AtA;
        AtA.matrix=AAT(A,n,m,1);
        AtA.n=m;
        AtA.solve();
        V=AtA.eigenvec;
        return V;
    }
    void solve(){
        solve_E();
        solve_U();
        solve_V();
    }
    void print(){
        cout<<"matrix U: "<<endl;
        printMatrix(U,n);
        cout<<"matrix E: "<<endl;
        printMatrix(E,n,m);
        cout<<"matrix Vt: "<<endl;
        printMatrix(transposition(V,m,m),m);
    }
};

class PCA{
    public:
    PCA(){}
    double** data;
    double* eigenvec_1;
    double* eigenvec_2;
    double eigenvalue_1;
    double eigenvalue_2;
    int n;
    int m;
    double** covMatrix;
    double** decentra(){//ok
        double* avg=new double[n];
        for(int i=0;i<n;i++){
            for(int j=0;j<m;j++){
                avg[i]+=data[i][j];
            }
        }
        for(int i=0;i<n;i++){
            avg[i]/=m;
        }
        for(int i=0;i<n;i++){
            for(int j=0;j<m;j++){
                data[i][j]-=avg[i];
            }
        }
        return data;
    }
    void solve(){
        data=decentra();
        covMatrix=AAT(data,n,m,0);
        printMatrix(covMatrix,n);
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                covMatrix[i][j]/=m;
            }
        }
        JacobiDecom PCA_solve;
        PCA_solve.matrix=covMatrix;
        PCA_solve.n=n;
        PCA_solve.solve();
        double** values=PCA_solve.result;
        double* eigenvalues=new double[n];
        for(int i=0; i<n;i++){
            eigenvalues[i]=values[i][i];
        }
        eigenvalue_1=eigenvalues[0];
        eigenvalue_2=eigenvalues[1];
        double** eivec=PCA_solve.eigenvec;
        eivec=transposition(eivec,n,n);
        eigenvec_1=eivec[0];
        eigenvec_2=eivec[1];
    }
    void print(){
        cout<<eigenvalue_1<<": ";
        printList(eigenvec_1,n);
        cout<<endl;
        cout<<eigenvalue_2<<": ";
        printList(eigenvec_2,n);
        cout<<endl;
    }
};

int main(){
    // Q(b)
    /*JacobiDecom Qb;
    double** A=createRandomM(4,3);
    Qb.matrix=AAT(A,4,3,0);
    Qb.n=4;
    printMatrix(A,4,3);
    printMatrix(Qb.matrix,4);
    cout<<endl;
    Qb.solve();
    Qb.print();
    for(int i=0;i<4;i++){
        printMatrix(matrixSub(Qb.matrix,Qb.result[i][i],4),4);
        cout<<endl;
    }*/
    /*SVD Qb_2;
    Qb_2.A=A;
    Qb_2.n=4;
    Qb_2.m=3;
    Qb_2.solve();
    Qb_2.print();*/
    //Q(c)
    PCA test;
    test.data=transposition(readMatrix("iris.txt",150,5),150,4);
    test.n=4;
    test.m=150;
    test.solve();
    printMatrix(test.covMatrix,4);
    test.print();
    return 0;
}
