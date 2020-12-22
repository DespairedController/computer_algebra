#include <iostream>
#include <tuple>
#include "test.h"

using namespace Givaro;

int main() {
    std::cout << "q,d,time" << std::endl;
    std::vector<std::tuple<long, long, long>> ans;
    CantorZassenhaus::GaloisField Z3(7, 1);  // integers modulo 13
    Poly1Dom<GFqDom<int>, Dense> DP29(Z3, Indeter("X"));
}
