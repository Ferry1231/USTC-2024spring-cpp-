import matplotlib.pyplot as plt
import numpy as np

with open("Computational Methods/iris.txt","r") as f:
    content=f.read()

lines=content.split("\n")
data=[]
for line in lines:
    values=line.split(',')
    row=[float(value) for value in values]
    data.append(row)

data=np.array(data)
E1=np.array([0.36159,-0.0822689,0.856572,0.358844])
E2=np.array([0.65654,0.729712,-0.175767,-0.0747065])

data=data.T
avg=[]
for i in range(4):
    avg.append(sum(data[i])/150)
    #print(sum(data[i])-data[i][4])
for i in range(4):
    for j in range(150):
        data[i][j]-=avg[i]

data=data.T
print(data)

dict={0:"red",1:"green",2:"blue"}

for i in range(150):
    x1=np.dot(data[i][:4],E1)
    x2=np.dot(data[i][:4],E2)
    plt.scatter(x1,x2,color=dict[data[i][4]])
plt.xlabel("$x_1$")
plt.ylabel("$x_2$")
plt.show()
