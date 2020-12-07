import matrix_helpers as mh
import numpy as np

class ConditionNumber:
    def __init__(self, matrixAB=None, matrixA=None, vectorB=None):
        self.whole_matrix = matrixAB

        if matrixAB is None:
            self.matrixA = matrixA
            self.vectorB = vectorB
        else:
            self.matrixA = mh.get_matrixA(matrixAB)
            self.vectorB = mh.get_vectorB_unpacked(matrixAB)
        
        # Get the condition number of the matrix
        if matrixAB is None:
            self.cond = mh.get_matrix_cond(matrixA)
        else:
            self.cond = mh.get_matrix_cond(mh.get_matrixA(matrixAB))
        
    def experimental(self, vectorX, delta_vectorX, printall=False):
        print("\n--- Experimental condition number check")
        # Get the norms of the vectorB and delta_vectorB
        vectorB_norm = mh.get_vector_norm_euc(self.vectorB)
        delta_vectorB = mh.get_deltaB_from_vectorB(self.vectorB)
        deltaB_norm = mh.get_vector_norm_euc(delta_vectorB)

        # Get the norms of the vectorX and delta_vectorX
        vectorX_norm = mh.get_vector_norm_euc(vectorX)
        deltaX_norm = mh.get_vector_norm_euc(delta_vectorX)

        if printall:
            print(f'\n --- Norms:\n{"||B||":<12} {"=":>4} {vectorB_norm}')
            print(f'{"||deltaB||":<12} {"=":>4} {deltaB_norm}')
            print(f'{"||X||":<12} {"=":>4} {vectorX_norm}')
            print(f'{"||deltaX||":<12} {"=":>4} {deltaX_norm}')
            print("\nMatrix determinant =", mh.get_matrix_determinant(self.matrixA))
            print("\nChecking Cond(A) for matrix (1):")
            print("Cond(A) >= ( ||deltaX|| / ||X|| ) * ( ||B|| / ||deltaB|| )")
            print(self.cond, ">= (", deltaX_norm, "/", vectorX_norm, ") * (", vectorB_norm, "/", deltaB_norm, ")")
            print(self.cond, ">=", (deltaX_norm / vectorX_norm) * (vectorB_norm / deltaB_norm))

            print("\nChecking Cond(A) for matrix (2):")
            print("||deltaX|| / ||X|| <= cond(A) * ( ||deltaB|| / ||B||)")
            print(deltaX_norm, "/", vectorX_norm, "<=", self.cond, "* (", deltaB_norm, "/", vectorB_norm, ")")
            print(deltaX_norm / vectorX_norm, "<=", self.cond * (deltaB_norm / vectorB_norm))

        condition_1 = self.cond >= (deltaX_norm / vectorX_norm) * (vectorB_norm / deltaB_norm)
        condition_2 = (deltaX_norm / vectorX_norm) <= self.cond * (deltaB_norm / vectorB_norm)

        if printall:
            print("\n --- Condition returned (1)", condition_1, "and (2)", condition_2, "\n")
            print("Condition number of matrix Cond(A) =", self.cond)
            print("Condition number from np.linalg.cond(..., 1) =", np.linalg.cond(self.matrixA, 1))
            if self.cond < 10:
                print("Cond(A) of the matrix falls into 1-10 range. This means that the results of this matrix calculation are adequate.")
            elif self.cond < 1000:
                print("Cond(A) of the matrix falls into 10-1000 range. This means that the results of this matrix calculation CAN BE INADEQUATE.")
            elif self.cond > 1000:
                print("Cond(A) of the matrix falls into 1000-infinity range. This means that the results of this matrix calculation ARE INADEQUATE.")
        
        return condition_1, condition_2
        
