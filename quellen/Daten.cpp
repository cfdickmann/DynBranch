#include "AmericanOption.h"

using namespace std;

void AmericanOption::Daten() {
    int Example =3;

    X0 = DoubleFeld(100); //genug Platz
    sigma = DoubleFeld(100);

    if (Example == 3) { //Glasserman Example MaxCall
        PfadModell = ITO;
        option = MAX_CALL;
        delta = 0.1;
        D = 3;//Achtung!
        for (int j = 0; j < D; ++j) {
            X0[j] = 90.;
            sigma[j] = 0.2;
        }
        Strike = 100.;
        r = 0.05;
        T = 3;
        N = 10;
        number_of_Exercise_Dates = 10;
    }
}
