import condition_number as cn
import matrix_helpers as mh
import numpy as np

# Creating a 2X2 matrix 
matrix = np.array([[4, 2], [3, 1]]) 
inversed = np.linalg.inv(matrix)
vector = [2, 3, 5]
  
print("Original matrix:") 
print(matrix) 
print("Inversed matrix:")
print(inversed)
  
# Output 
result =  np.linalg.cond(matrix, 1) 
print("Condition number of the matrix:") 
print(result) 

print("Condition number of the matrix by our method:", mh.get_matrix_cond(matrix))

# Good
#print("Norm by numpy (vec):", np.linalg.norm(vector))
#print("Norm by my (vec):", mh.get_vector_norm_euc(vector))
# ---

print("Norm by numpy (matrix):", np.linalg.norm(matrix, 1))
print("Norm by my (matrix):", mh.get_matrix_norm(matrix))

# Good
#print("Norm vector =", np.linalg.norm([1, 2, 3]))
#print("my Norm vector =", mh.get_vector_norm_euc([1, 2, 3]))
# ---

matrix_2 = [[0.1, -0.4, 0],
            [0.2, 0, -0.3],
            [0, 0.1, 0.3 ]]
print("Matrix norm (numpy)=", np.linalg.norm(matrix_2))
print("my matrix norm=", mh.get_matrix_norm(matrix_2))
