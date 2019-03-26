/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "IntRule1d.h"
#include "tpanic.h"

#define PI 3.1415926535898

IntRule1d::IntRule1d() {
    fOrder = -1;
    fPoints.Resize(1, 1);
    fWeights.resize(0);
}

IntRule1d::IntRule1d(int order) {
    SetOrder(order);
}

void IntRule1d::SetOrder(int order) {
    if (order < 0) {
        DebugStop();
    }

    fOrder = order;
    int nPoints = order + 1;
    fPoints.Resize(nPoints, 1);
    fWeights.resize(nPoints);

    switch(order) {

        case 0:
            fPoints(0, 0) =  0.00000000000000; fWeights[0] = 2.00000000000000;
            break;

        case 1:
            fPoints(0, 0) = -0.57735026918963; fWeights[0] = 1.00000000000000;
            fPoints(1, 0) =  0.57735026918963; fWeights[1] = 1.00000000000000;
            break;

        case 2:
            fPoints(0, 0) = -0.77459666924148; fWeights[0] = 0.55555555555555;
            fPoints(1, 0) =  0.00000000000000; fWeights[1] = 0.88888888888889;
            fPoints(2, 0) =  0.77459666924148; fWeights[2] = 0.55555555555555;
            break;

        case 3:
            fPoints(0, 0) = -0.86113631159405; fWeights[0] = 0.34785484513745;
            fPoints(1, 0) = -0.33998104358486; fWeights[1] = 0.65214515486255;
            fPoints(2, 0) =  0.33998104358486; fWeights[2] = 0.65214515486255;
            fPoints(3, 0) =  0.86113631159405; fWeights[3] = 0.34785484513745;
            break;

        case 4:
            fPoints(0, 0) = -0.90617984593866; fWeights[0] = 0.23692688505618;
            fPoints(1, 0) = -0.53846931010568; fWeights[1] = 0.47862867049937;
            fPoints(2, 0) =  0.00000000000000; fWeights[2] = 0.56888888888889;
            fPoints(3, 0) =  0.53846931010568; fWeights[3] = 0.47862867049937;
            fPoints(4, 0) =  0.90617984593866; fWeights[4] = 0.23692688505618;
            break;

        case 5:
            fPoints(0, 0) = -0.93246951420315; fWeights[0] = 0.17132449237916;
            fPoints(1, 0) = -0.66120938646626; fWeights[1] = 0.36076157304814;
            fPoints(2, 0) = -0.23861918608320; fWeights[2] = 0.46791393457269;
            fPoints(3, 0) =  0.23861918608320; fWeights[3] = 0.46791393457269;
            fPoints(4, 0) =  0.66120938646626; fWeights[4] = 0.36076157304814;
            fPoints(5, 0) =  0.93246951420315; fWeights[5] = 0.17132449237916;
            break;

        case 6:
            fPoints(0, 0) = -0.94910791234276; fWeights[0] = 0.12948496616886;
            fPoints(1, 0) = -0.74153118559939; fWeights[1] = 0.27970539148928;
            fPoints(2, 0) = -0.40584515137740; fWeights[2] = 0.38183005050512;
            fPoints(3, 0) =  0.00000000000000; fWeights[3] = 0.41795918367347;
            fPoints(4, 0) =  0.40584515137740; fWeights[4] = 0.38183005050512;
            fPoints(5, 0) =  0.74153118559939; fWeights[5] = 0.27970539148928;
            fPoints(6, 0) =  0.94910791234276; fWeights[6] = 0.12948496616886;
            break;

        case 7:
            fPoints(0, 0) = -0.96028985649754; fWeights[0] = 0.10122853629037;
            fPoints(1, 0) = -0.79666647741363; fWeights[1] = 0.22238103445337;
            fPoints(2, 0) = -0.52553240991633; fWeights[2] = 0.31370664587789;
            fPoints(3, 0) = -0.18343464249565; fWeights[3] = 0.36268378337836;
            fPoints(4, 0) =  0.18343464249565; fWeights[4] = 0.36268378337836;
            fPoints(5, 0) =  0.52553240991633; fWeights[5] = 0.31370664587789;
            fPoints(6, 0) =  0.79666647741363; fWeights[6] = 0.22238103445337;
            fPoints(7, 0) =  0.96028985649754; fWeights[7] = 0.10122853629037;
            break;

        case 8:
            fPoints(0, 0) = -0.96816023950763; fWeights[0] = 0.08127438836157;
            fPoints(1, 0) = -0.83603110732664; fWeights[1] = 0.18064816069486;
            fPoints(2, 0) = -0.61337143270059; fWeights[2] = 0.26061069640294;
            fPoints(3, 0) = -0.32425342340381; fWeights[3] = 0.31234707704000;
            fPoints(4, 0) =  0.00000000000000; fWeights[4] = 0.33023935500126;
            fPoints(5, 0) =  0.32425342340381; fWeights[5] = 0.31234707704000;
            fPoints(6, 0) =  0.61337143270059; fWeights[6] = 0.26061069640294;
            fPoints(7, 0) =  0.83603110732664; fWeights[7] = 0.18064816069486;
            fPoints(8, 0) =  0.96816023950763; fWeights[8] = 0.08127438836157;
            break;

        case 9:
            fPoints(0, 0) = -0.97390652851717; fWeights[0] = 0.06667134430868;
            fPoints(1, 0) = -0.86506336668898; fWeights[1] = 0.14945134915058;
            fPoints(2, 0) = -0.67940956829902; fWeights[2] = 0.21908636251598;
            fPoints(3, 0) = -0.43339539412925; fWeights[3] = 0.26926671930999;
            fPoints(4, 0) = -0.14887433898163; fWeights[4] = 0.29552422471475;
            fPoints(5, 0) =  0.14887433898163; fWeights[5] = 0.29552422471475;
            fPoints(6, 0) =  0.43339539412925; fWeights[6] = 0.26926671930999;
            fPoints(7, 0) =  0.67940956829902; fWeights[7] = 0.21908636251598;
            fPoints(8, 0) =  0.86506336668898; fWeights[8] = 0.14945134915058;
            fPoints(9, 0) =  0.97390652851717; fWeights[9] = 0.06667134430868;
            break;

        case 10:
            fPoints(0, 0) = -0.97822865814606; fWeights[0] = 0.05566856711617;
            fPoints(1, 0) = -0.88706259976810; fWeights[1] = 0.12558036946490;
            fPoints(2, 0) = -0.73015200557405; fWeights[2] = 0.18629021092773;
            fPoints(3, 0) = -0.51909612920681; fWeights[3] = 0.23319376459199;
            fPoints(4, 0) = -0.26954315595234; fWeights[4] = 0.26280454451025;
            fPoints(5, 0) =  0.00000000000000; fWeights[5] = 0.27292508677790;
            fPoints(6, 0) =  0.26954315595234; fWeights[6] = 0.26280454451025;
            fPoints(7, 0) =  0.51909612920681; fWeights[7] = 0.23319376459199;
            fPoints(8, 0) =  0.73015200557405; fWeights[8] = 0.18629021092773;
            fPoints(9, 0) =  0.88706259976810; fWeights[9] = 0.12558036946490;
            fPoints(10, 0) =  0.97822865814606; fWeights[10] = 0.05566856711617;
            break;

        default:
            gauleg(-1, 1, fPoints, fWeights);
    }
}

/* Reference: Press et al. Numerical Recipes: The Art of Scientific Computing. (3 ed.). Page: 184
 * */
void IntRule1d::gauleg(const double x1, const double x2, Matrix &co, VecDouble &w) {
    auto nPoints = co.Rows();
    double m = (nPoints + 1) / 2.;
    double xm = (x1 + x2) / 2.;
    double xl = (x2 - x1) / 2.;

    double z, z1, pp, p1, p2, p3;

    for (auto i = 0; i < m; i++) {

        // Initial root approximation
        z = cos(PI * (i + 0.75) / (nPoints + 0.5));

        do {
            p1 = 1.0;
            p2 = 0.0;

            for (auto j = 0; j < nPoints; j++) {
                p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1);
            }

            pp = nPoints * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        }
        while (std::abs(z - z1) > EPS);

        co(i, 0) = xm - xl * z;
        co(nPoints - 1 - i, 0) = xm + xl * z;
        w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
        w[nPoints - 1 - i] = w[i];
    }
}
