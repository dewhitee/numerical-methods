import matrix_norm as norm
import matrix_inverse as inverse

#
#
#
#

def matrix_cond(matrix: list) -> float:
    """

    """
    return norm.matrix_norm(matrix) * norm.matrix_norm(inverse.matrix_inverse(matrix))
