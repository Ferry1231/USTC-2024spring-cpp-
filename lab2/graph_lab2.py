import matplotlib.pyplot as plt

def txtToList(filename)->list:
    result=[]
    with open(filename,"r") as f:
        for line in f:
            row=[float(x) for x in line.split()]
            result.append(row)  
    return result

GC=txtToList("output_gc.txt")
GS=txtToList("output_gs.txt")
RE=txtToList("output_re.txt")

SUB_GC=[]
SUB_GS=[]

count=[i+1 for i in range(len(RE[0]))]

for i in range(4):
    sub_gc=[GC[i][j]-RE[i][j] for j in range(len(RE[i]))]
    sub_gs=[GS[i][j]-RE[i][j] for j in range(len(RE[i]))]
    SUB_GC.append(sub_gc)
    SUB_GS.append(sub_gs)

plt.plot(count,SUB_GS[0],label='$\epsilon=1$')
plt.plot(count,SUB_GS[1],label='$\epsilon=0.1$')
plt.plot(count,SUB_GS[2],label='$\epsilon=0.01$')
plt.plot(count,SUB_GS[3],label='$\epsilon=0.0001$')
plt.title("Gauss Seidel Iteration Errors")
plt.xlabel("count")
plt.ylabel("errors")
plt.xticks(range(0,101,5))
#plt.yticks(range(0,0.5,0.001))
plt.legend()
plt.show()