/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeoElementSide.h"
#include "tpanic.h"

GeoElementSide::GeoElementSide() {
    DebugStop();
}

GeoElementSide::GeoElementSide(const GeoElementSide &copy) {
    DebugStop();
}

GeoElementSide &GeoElementSide::operator=(const GeoElementSide &copy) {
    DebugStop();
}

GeoElementSide GeoElementSide::Neighbour() const {
    DebugStop();
}

void GeoElementSide::SetNeighbour(const GeoElementSide &neighbour) {
    DebugStop();
}

bool GeoElementSide::IsNeighbour(const GeoElementSide &candidate) {
    DebugStop();
}

void GeoElementSide::IsertConnectivity(GeoElementSide &candidate) {
    DebugStop();
}

void GeoElementSide::AllNeighbours(std::vector<GeoElementSide> &allneigh) {
    DebugStop();
}

void GeoElementSide::ComputeNeighbours(std::vector<GeoElementSide> &compneigh) {
    DebugStop();
}

bool GeoElementSide::DataConsistency(GeoElementSide &candidate) {
    DebugStop();
}
