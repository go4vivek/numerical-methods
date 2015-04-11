n=3
m=n-1

a = 2./3*n**3 - 1./2*n**2 + 5./6*n - 1
b = 2./3*m**3 + 3./2*m**2 + 11./6*m
c = (2.*m**2+5.*m+3.)*n - 4./3*m**3 - 11./2*m**2 - 37./6*m - 3

print a, b, c

Sum = 0.
for i in range(n-m,n-1+1):
	Sum += 2*(n-i)**2 + (n-i+1)
print Sum  