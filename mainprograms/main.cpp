

#include <iostream>
#include <math.h>
#include <TMatrix.h>
#include <IntRule1d.h>

using std::cout;
using std::endl;
using std::cin;

int main ()
{
    int intnum;
    cout << "Polinomio:\n";
    cout << "     f(x) = 2 + x -x^2\n\n";
    cout << "Integrando de 0 a 0.8\n\n";
    VecDouble f = {2, 1, -1};
    IntRule1d polinomio(intnum);
    polinomio.Print(std::&cout);
    return 0;
}