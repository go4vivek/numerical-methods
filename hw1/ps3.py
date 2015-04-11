import numpy as np 
from scipy.linalg import lu

A = np.matrix([[2,-1,1], [-2,1,3], [4,0,-1]])
print lu(A)