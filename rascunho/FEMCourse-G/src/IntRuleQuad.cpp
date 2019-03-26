/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "IntRuleQuad.h"
#include "tpanic.h"

IntRuleQuad::IntRuleQuad() {
    fOrder = -1;
    fPoints.Resize(1, 1);
    fWeights.resize(0);
}

IntRuleQuad::IntRuleQuad(int order) {
    SetOrder(order);
}

void IntRuleQuad::SetOrder(int order) {
    DebugStop();  
}

void IntRuleQuad::gaulegQuad(const double x1, const double x2, VecDouble &co, VecDouble &w) {
    DebugStop();   
}
