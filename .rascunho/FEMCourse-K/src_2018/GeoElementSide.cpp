//
//  GeoElementSide.cpp
//  FemSC
//
//  Created by Karolinne Oliveira Coelho on 4/28/18.
//
//

#include <stdio.h>
#include "GeoElementSide.h"
#include "GeoElement.h"
#include "GeoMesh.h"
#include <algorithm>

GeoElementSide::GeoElementSide() : fElement(0), fSide(0){
    
}

GeoElementSide::GeoElementSide(const GeoElementSide &copy) {
    fSide = copy.fSide;
    fElement = copy.fElement;
}

GeoElementSide &GeoElementSide::operator=(const GeoElementSide &copy) {
    fElement = copy.fElement;
    fSide = copy.fSide;
    return *this;
}

GeoElementSide GeoElementSide::Neighbour() const {
    return fElement ? fElement->Neighbour(fSide) : GeoElementSide();
}

void GeoElementSide::SetNeighbour(const GeoElementSide &neighbour) {
    fElement->SetNeighbour(fSide,neighbour);
}

bool GeoElementSide::IsNeighbour(const GeoElementSide &candidate) {
    if(*this == candidate) return true;
    GeoElementSide neigh = Neighbour();
    while((neigh == *this) == false)
    {
        if(neigh == candidate) return true;
        neigh = neigh.Neighbour();
    }
    return false;
}

void GeoElementSide::IsertConnectivity(GeoElementSide &candidate) {
    GeoElementSide myneigh = Neighbour();
    GeoElementSide neighneigh = candidate.Neighbour();
    
    if (IsNeighbour(candidate)) {
        return;
    }
    candidate.SetNeighbour(myneigh);
    SetNeighbour(neighneigh);
    
    //    if ((candidate.Neighbour() == myneigh) == false) {
    //        neighbour.SetNeighbour(myneigh);
    //    }
    //
    //    if ((myneigh.Neighbour() == candidate) == false) {
    //        myneigh.SetNeighbour(candidate);
    //    }
    
}

void GeoElementSide::AllNeighbours(std::vector<GeoElementSide> &compneigh){
    GeoElementSide neigh = Neighbour();
    int64_t pos = 0;
    while(!(neigh == *this))
    {
        compneigh[pos] = neigh;
        neigh = neigh.Neighbour();
        pos++;
    };
    return;
}

void GeoElementSide::ComputeNeighbours(std::vector<GeoElementSide> &neighbour) {
    
    std::vector<GeoElementSide> GeoElSideSet(100);
    VecInt GeoElSet(100,0);
    VecInt nodes;
    VecInt nodesresult;
    
    if(fSide < fElement->NCornerNodes())
    {
        AllNeighbours(neighbour);
    }
    
    // Vetor para guardar todos os elementos associados aos sets
    int nsnodes = fElement->NSideNodes(fSide);
    VecInt nodesside(nsnodes);
    
    int64_t nel = 0;
    for(int in=0; in<nsnodes; in++)
    {
        int locnod = fElement->SideNodeIndex(fSide,in);
        GeoElementSide locside(fElement,locnod);
        locside.AllNeighbours(GeoElSideSet);
        
        while (GeoElSideSet[nel].fElement != NULL) {
            int64_t index = GeoElSideSet[nel].Element()->GetIndex();
            GeoElSet[nel] = index;
            nel++;
        }
    }
    GeoElSet.erase(GeoElSet.begin()+nel, GeoElSet.end());
    
    // Interseção dos elementos
    std::sort(GeoElSet.begin(), GeoElSet.end());
    auto last = std::unique(GeoElSet.begin(), GeoElSet.end());
    GeoElSet.erase(last, GeoElSet.end());
    
    // Preenchendo o vetor de vizinhos
    GeoMesh * gmesh = fElement->GetMesh();
    
    fElement->GetNodes(nodes);
    for (int ins=0; ins<nsnodes; ins++) {
        nodesside[ins] = nodes[fElement->SideNodeIndex(fSide, ins)];
    }
    std::sort(nodesside.begin(), nodesside.end()); // ordenando de forma crescente o vetor de índices de nós do side
    
    VecInt nodessideresult(nel);
    for(int iel=0; iel<nel; iel++) {
        GeoElement * gelresult = gmesh->Element(GeoElSet[iel]);
        gelresult->GetNodes(nodesresult);
        
        for (int ins=0; ins<gelresult->NSides(); ins++) {
            int nslocal = gelresult->NSideNodes(ins);
            nodessideresult.resize(nslocal, 0);
            
            for (int ilocal=0; ilocal<nslocal; ilocal++) {
                nodessideresult[ilocal] = nodesresult[gelresult->SideNodeIndex(ins, ilocal)];
            }
            std::sort(nodessideresult.begin(), nodessideresult.end());
            
            if (nodessideresult==nodesside) {
                neighbour.push_back(GeoElementSide(gelresult,ins));
            }
        }
    }
    
    GeoElSet.clear();
    GeoElSideSet.clear();
    
    nodes.clear();
    nodesside.clear();
    
    nodesresult.clear();
    nodessideresult.clear();
}
