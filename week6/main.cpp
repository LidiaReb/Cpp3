#include <iostream>
#include <cassert>
#include "grid.h"

int main() {
    Grid<float,3> const g3(2, 3, 4, 1.0f);
    g3.print();
    assert(1.0f == g3(1, 1, 1));

    Grid<float,2> g2(2, 5, 2.0f);
    assert(2.0f == g2(1, 1));

    g2 = g3[1];
    assert(1.0f == g2(1, 1));
    
    std::cout << "Done\n";
}