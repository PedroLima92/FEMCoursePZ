//
//  GeoMesh.cpp
//  FemSC
//
//  Created by Karolinne Oliveira Coelho on 4/28/18.
//
//

#include <stdio.h>
#include <fstream>

#include "GeoMesh.h"
#include "GeoElementSide.h"


GeoMesh::GeoMesh(){
    Nodes.resize(0);
    Elements.resize(0);
    Reference = NULL;
    fDim = 0;
};

GeoMesh::GeoMesh(const GeoMesh &cp){
    this->operator =(cp);
}

GeoMesh &GeoMesh::operator=(const GeoMesh &cp){
    Nodes = cp.Nodes;
    Elements = cp.Elements;
    Reference = cp.Reference;
    fDim = cp.fDim;
    return *this;
}

void GeoMesh::SetNumNodes(int nnodes){
    Nodes.resize(nnodes);
}

void GeoMesh::SetNumElements(int numelements){
    Elements.resize(numelements);
}

int GeoMesh::NumNodes(){
    return Nodes.size();
}

int GeoMesh::NumElements(){
    return Elements.size();
}

GeoNode &GeoMesh::Node(int node){
    return Nodes[node];
}

void GeoMesh::SetElement(int elindex, GeoElement *gel){
    if (elindex+1 >= Elements.size()) {
        Elements.resize(elindex+1);
    }
    Elements[elindex] = gel;
}

GeoElement *GeoMesh::Element(int elindex){
    return Elements[elindex];
}

void GeoMesh::BuildConnectivity(){
    
    int64_t nelem = NumElements();
    
    std::vector<GeoElement*> NeighNode(NumNodes(),0);
    std::vector<GeoElementSide> neighbours(0);
//    VecInt SideNum(NumNodes(),-1);
    
    int64_t SideNum[NumNodes()];
    for (int64_t inodes=0; inodes<NumNodes(); inodes++) {
        SideNum[inodes] = -1;
    }
    
    for(int64_t iel=0; iel<nelem; iel++)
    {
        GeoElement *gel = Elements[iel];
        if(!gel) continue;
        
        int ncor = gel->NCornerNodes();
        
        for(int in=0; in<ncor; in++)
        {
            
            int64_t nod = gel->NodeIndex(in);
            if(SideNum[nod] == -1) // Se o side não tá definido
            {
                SideNum[nod] = in;
                NeighNode[nod] = gel;
            }
            else // Caso o side já esteja definido
            {
                GeoElementSide neigh(NeighNode[nod], SideNum[nod]); // Cria um ElementSide com base na definição (vetores)
                GeoElementSide gelside(gel,in); // ElementSide criado para teste com o corner node do for

                GeoElementSide neighbour = gelside.Neighbour();
                if (neighbour.Element() == 0) DebugStop();
                
                if (gelside.IsNeighbour(neigh) == false)
                {
                    gelside.IsertConnectivity(neigh);
                }
            }
        }
    }
    
    int64_t SideNum0[NumElements()];
    for (int64_t iel=0; iel<NumElements(); iel++) {
        SideNum0[iel] = -1;
    }
    for(int64_t iel=0; iel<nelem; iel++)
    {
        GeoElement *gel = Elements[iel];
        if(!gel) continue;
        
        int ncor = gel->NCornerNodes();
        int nsides = gel->NSides();
        
        for(int is=ncor; is<nsides; is++)
        {
            SideNum0[iel] = is;
            
            int64_t SideNum1[NumElements()];
            for (int64_t iel=0; iel<NumElements(); iel++) {
                SideNum1[iel] = -1;
            }
            
            GeoElementSide gelside(gel,is);
            gelside.ComputeNeighbours(neighbours);
            
            int64_t nneigh = neighbours.size();
            
            for(int64_t in=0; in<nneigh; in++) {
                if(neighbours[in].Side() == -1)
                {
                    std::cout << "GeoMesh::BuildConnectivity : Inconsistent mesh detected analysing element/side:" << std::endl;
                    continue;
                }
                
                if(SideNum1[iel] == -1) // ERRO
                {
                    SideNum1[iel] = in;
                }
                gelside.IsertConnectivity(neighbours[in]);
            }
            neighbours.clear();
        }
    }
    NeighNode.clear();
    std::cout << "Geometric Mesh: BuildConnectivity done!"<< std::endl;
    return;
}

void GeoMesh::Print(std::ostream &out){
    out << "\n\t\t GEOMETRIC GeoMesh INFORMATIONS:\n\n";
    out << "number of nodes               = " << Nodes.size() << "\n";
    out << "number of elements            = " << Elements.size() << "\n";
    
    out << "\n\tGeometric Node Information:\n\n";
    int64_t i;
    int64_t nnodes = Nodes.size();
    for(i=0; i<nnodes; i++)
    {
        out << "Index: " << i << " ";
        Nodes[i].Print(out);
    }
    out << "\n\tGeometric Element Information:\n\n";
    int64_t nelem = Elements.size();
    for(i=0; i<nelem; i++)
    {
        if(Elements[i]) Elements[i]->Print(out);
        out << "\n";
    }
}
