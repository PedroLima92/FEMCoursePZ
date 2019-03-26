/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "IntRule.h"
#include "tpanic.h"
#include <iomanip>

IntRule::IntRule() : fPoints(1, 1), fWeights(0), fOrder(-1) {
}

IntRule::IntRule(int order) {
    DebugStop();
}

IntRule::~IntRule() {
}

IntRule::IntRule(const IntRule& copy) {
    fOrder = copy.fOrder;
    fWeights = copy.fWeights;
    fPoints = copy.fPoints;
}

IntRule &IntRule::operator=(const IntRule &cp) {
    DebugStop();
}

int IntRule::NPoints() const {
    DebugStop();
}

void IntRule::Print(std::ostream &out) {

    out << ":: Integration Rule Info" << std::endl;
    out << " Order: " << fOrder << std::endl;

    for (auto i = 0; i < fPoints.Rows();  i++) {

        out << " Point: " << std::right << std::setw(2) << i
            << std::fixed << std::setprecision(printPrecision)
            << ",    x: " << std::right << std::setw(printPrecision + 3) << fPoints(i, 0)
            << ",    w: " << std::right << std::setw(printPrecision + 3) << fWeights[i]
            << std::endl;
    }
}

void IntRule::Point(int p, VecDouble& co, double& w) const {
    DebugStop();
}
