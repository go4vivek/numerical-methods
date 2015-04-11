import numpy as np 
from scipy.linalg import lu

A = np.matrix([[1,0,0,0,1],[-1,1,0,0,1],[-1,-1,1,0,1],[-1,-1,-1,1,1],[-1,-1,-1,-1,1]])
print lu(A)
