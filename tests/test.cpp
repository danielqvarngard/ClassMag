#include <iostream>

int main(){
    for (auto ii = 0u; ii < 2; ++ii)
        std::cout << 2.0*static_cast<double>(ii % 2) - 1.0 << " ";
    std::cout << "\n";
}