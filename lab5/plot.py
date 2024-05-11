import matplotlib.pyplot as plt
import numpy as np
import math


Xplain=[i-9 for i in range(21)]
Yplain=[]

X=[]

with open('Computational Methods/lab5/point2.txt',"r") as f1:
    for line in f1:
        values = list(map(float, line.split()))
        X.append(values[0])
        Yplain.append(values[1])
        if X[-1]==11:
            X.pop()
            break
        for i in range(1,1000):
            X.append(values[0]+i/1000)


F=[]

with open('Computational Methods/lab5/result2.txt',"r") as f2:
    for line in f2:
        values = list(map(float, line.split()))
        F.append(values)

n=len(F)


Y=[(F[math.floor(x)+9][0]+F[math.floor(x)+9][1]*x+F[math.floor(x)+9][2]*x**2+F[math.floor(x)+9][3]*x**3) for x in X]

co=[math.floor(x)+9 for x in X]


X=np.array(X)
Y=np.array(Y)
Xplain=np.array(Xplain)
Yplain=np.array(Yplain)

plt.scatter(X,Y,s=1,color='b')
plt.scatter(Xplain,Yplain,s=10,color='r')
plt.xlabel('X Axis')
plt.ylabel('Y Axis')
#plt.ylim(-2,3)
plt.title("$M spline$")
plt.grid(True)

plt.show()
