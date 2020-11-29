#include <vector>
#include <map>

using namespace std;

template <typename T>
using Matrix = vector<vector<T>>;

template <typename T>
class Gauss
{
    Matrix<T>       _matrix;
    vector<T>       _answers;

public:

    Gauss(Matrix<T> matrix)
        : _matrix(matrix)
    {
        Calculate();
    }

    vector<T> GetAnswers() { return _answers; }

private:

    void Calculate()
    {
        size_t initSize = _matrix.size();
        Matrix<T> clone = _matrix;

        size_t n = initSize;
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n + 1; j++)
            {
                clone[i][j] = _matrix[i][j];
            }
        }

        StraightStep(clone, n);
        ReversedStep(clone, n);
        SeparateAnswers(clone, n);
    }

    void StraightStep(Matrix<T>& clone, size_t n)
    {
        // Straight step (Nulling the bottom-left corner)
        for (size_t k = 0; k < n; k++) // k-th row number
        {
            for (size_t i = 0; i < n + 1; i++) // i-th column number
            {
                clone[k][i] /= _matrix[k][k]; // Dividing the k-th row by the first term != 0 for transforming it into 1
            }

            for (size_t i = k + 1; i < n; i++) // i-th number of the next row after k
            {
                T K = clone[i][k] / clone[k][k]; // Coefficient

                for (size_t j = 0; j < n + 1; j++)
                {
                    clone[i][j] -= clone[k][j] * K; // Nulling of matrix elements under the first term transformed into 1
                }

                for (size_t i = 0; i < n; i++) // Updating the initial matrix
                {
                    for (size_t j = 0; j < n + 1; j++)
                    {
                        _matrix[i][j] = clone[i][j];
                    }
                }
            }
        }
    }

    void ReversedStep(Matrix<T>& clone, size_t n)
    {
        // Reversed step (Nulling upper-left corner)
        for (size_t k = n - 1; k > -1; k--)
        {
            for (size_t i = n; i > -1; i--)
            {
                clone[k][i] /= _matrix[k][k];
            }

            for (size_t i = k - 1; i > -1; i--)
            {
                T K = clone[i][k] / clone[k][k];
                for (size_t j = n; j > -1; j--)
                {
                    clone[i][j] -= clone[k][j] * K;
                }
            }
        }
    }

    void SeparateAnswers(Matrix<T> &clone, size_t n)
    {
        for (size_t i = 0; i < n; i++)
        {
            _answers.push_back(clone[i][n]);
        }
    }
};