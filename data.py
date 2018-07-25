from numpy import genfromtxt, savetxt
data = genfromtxt('nature16965-s2.txt')
savetxt('transpose.txt',data.T)