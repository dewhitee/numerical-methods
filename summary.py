import matrix_helpers as mh
from numpy import array
import lab1_gauss_elimination.gauss_elimination as ge
import lab1_gauss_seidel.gauss_seidel as gs
import condition_number as cn

class Summary:
    def __init__(self, matrixAB: list, adjusted_matrix: list, vars: list, X: list, print_only_results: bool = False, matrix_name: str = ""):
        self.matrixAB = matrixAB
        self.adjusted_matrix = adjusted_matrix
        self.vars = vars
        self.X = X
        self.print_only_results = print_only_results
        self.matrix_name = matrix_name
        self.AX = mh.get_matrixAX(matrixAB, X)
        self.B = mh.get_vectorB_unpacked(matrixAB)

        if not print_only_results:
            mh.show_A_B(matrixAB)
            mh.show_AX_B(matrixAB, X)

        self.condition_number = cn.ConditionNumber(matrixAB)
        self.condA = self.condition_number.cond
        print("\nCondition number Cond(A) of this matrix is =", self.condA)
        
        self.deltaB = mh.get_deltaB_from_whole_matrix(matrixAB)
        print("\nDeltaB vector:", self.deltaB)

        adjusted_matrix = mh.append_vectorB_to_whole_matrix(
            matrixAB, array(self.B) + array(self.deltaB))

        if not print_only_results:
            mh.full_print(adjusted_matrix, vars, 'Modified (adjusted) matrix')
            
    def compare_with_adjusted_gauss_elimination(self, printall):
        print("\n--- (Gauss Elimination) Comparing initial matrix with adjusted ---")
        deltaX = ge.GaussElimination(
            matrixAB=self.adjusted_matrix, 
            vars=self.vars, 
            print_only_results=self.print_only_results, 
            matrix_name="(Adjusted) " + self.matrix_name + " with B + deltaB").solution_vectorX

        if not self.print_only_results:
            print("DeltaX vector (Adjusted matrix Gauss Elimination result) = ", deltaX)

        self.condition_number.experimental(self.X, deltaX, printall)

    def compare_with_adjusted_gauss_seidel(self, equations, e, use_matrix, printall):
        print("\n--- (Gauss Seidel) Comparing initial matrix with adjusted ---")
        deltaX = gs.GaussSeidel(
            matrixAB=self.adjusted_matrix if use_matrix else None,
            equations=equations, 
            vars=self.vars, 
            e=e, 
            print_only_results=self.print_only_results,
            matrix_name="(Adjusted) " + self.matrix_name + " with B + deltaB").solution_vectorX

        if not self.print_only_results:
            print("DeltaX vector (Adjusted matrix Gauss Seidel result) = ", deltaX)
        
        self.condition_number.experimental(self.X, deltaX, printall)
