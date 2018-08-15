from vertex import *
def Reading(fileName):
    with open(fileName) as file:
        data = file.readlines()
        x = []
        y = []
        for i in data:
            x.append(float(i.split('\t')[0]))
            y.append(float(i.split('\t')[1]))
    return [x,y]
#get vertex list (X)
    # print('enter file name:   ')
    # fileName = input()
fileName = 'NACA4415.txt'
data = Reading(fileName)
X = [Vertex]
for i in range(0,len(data[0])):
    X.append(Vertex(data[0][i], data[1][i]))
n = len(X)-1
i = range(n-1)
# computing u
u = []
u.append(0)
sum1 = 0.0
sum2 = 0.0
e = 0.5
for k in i:
    for j in range(0,i[-1],1):
        asd = j
        sum1 = sum1 + X[j+1].length(X[j])**e
    for j in range(0,n):
        sum2 = sum2 + X[j+1].length(X[j])**e
    u.append(sum1/sum2)
print(len(u))

