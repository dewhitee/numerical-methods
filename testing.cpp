#include "lab1_gauss_elimination/gauss_elimination.h"
#include <iostream>

using namespace std;

int main()
{
    auto gauss = Gauss<double>(Matrix<double>(
        {{1.3, 1.33}, {2.3, 5.3}, {0.5, 3.003}}
    ));
    for (auto &&a : gauss.GetAnswers())
    {
        cout << "Answer = " << a << endl;
    }

    return 0;
}