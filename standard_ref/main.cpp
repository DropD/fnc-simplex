
#define VERBOSE
//~ #include "simple.hpp"
#include "general.hpp"

#include <iostream>
#include <string>

using namespace std;

int main(int argc, char ** argv) {

    string fname;
    
    if(argc != 2) {
        cout << "Requires input file" << endl;
        return 0;
    }
    else
        fname = argv[1];


    Simplex<double> s;
    s.load(fname);
    s.print();

    s.solve();

    vector<double> sol = s.solutions();
    cout << "Optimal value: " << sol[0] << endl;
    cout << "Variables:";
    for(int i = 1; i < sol.size(); ++i)
        cout << "  " << sol[i];
    cout << endl;


    return 0;

}
