from scipy.linalg import lu
import numpy as np 
A = np.matrix([[1,2,-7,-1.5,2],[4,4,0,-6,-2],[-2,-1,2,6,2.5],[0,-2,2,-4,1],[2,1,13,-5,3.5]])
print lu(A)