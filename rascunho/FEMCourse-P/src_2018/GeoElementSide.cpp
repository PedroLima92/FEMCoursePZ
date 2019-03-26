//
//  GeoElementSide.cpp
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#include "GeoElementSide.h"
#include "GeoMesh.h"

    
    GeoElementSide::GeoElementSide() : fElement(0), fSide(-1) {

    }

    
    GeoElementSide::GeoElementSide(const GeoElementSide &copy){
        fElement=copy.fElement;
        fSide = copy.fSide;
    }
    
    GeoElementSide &GeoElementSide::operator=(const GeoElementSide &copy){
        fElement=copy.fElement;
        fSide = copy.fSide;
        return *this;
    }

    GeoElementSide GeoElementSide::Neighbour() const{
        return fElement ? fElement->Neighbour(fSide) : GeoElementSide();
    }


    /** @brief Fill in the data structure for the neighbouring information*/
    void GeoElementSide::SetNeighbour(const GeoElementSide &neighbour){
        fElement->SetNeighbour(fSide,neighbour);
    }

    bool GeoElementSide::DataConsistency(GeoElementSide &candidate){
        int a = IsNeighbour(candidate);
        int b = candidate.IsNeighbour(*this);
        if((a && !b) || (!a && b))
        {
            std::cout << "Wrong data structure" <<std::endl;
            return false;
        }else{
            return true;
        }
    }


    bool GeoElementSide::IsNeighbour(const GeoElementSide &candidate){
        if(candidate == *this) return 1;
        GeoElementSide neighbour = Neighbour();
        if(!(neighbour.Element()!=0 && neighbour.Side()>-1)) return 0;
        while((neighbour == *this)==false) {
            if(candidate == neighbour) return 1;
            neighbour = neighbour.Neighbour();
        }
        return 0;
    }

    void GeoElementSide::IsertConnectivity(GeoElementSide &candidate){
#ifdef PADEBUG
        if(!DataConsistency(candidate)) DebugStop();
#endif
        if(IsNeighbour(candidate)) return;
        GeoElementSide neighneigh, currentneigh;
        neighneigh = candidate.Neighbour();
        currentneigh = Neighbour();
        SetNeighbour(neighneigh);
        candidate.SetNeighbour(currentneigh);
        
    }

    void GeoElementSide::AllNeighbours(std::vector<GeoElementSide> &allneigh) {
        GeoElementSide neigh = Neighbour();

        while((neigh == *this)==false)
        {
            allneigh.push_back(neigh);
            neigh = neigh.Neighbour();
        }
    }


    void GeoElementSide::ComputeNeighbours(std::vector<GeoElementSide> &neighbour){
        if(fSide < fElement->NCornerNodes())
        {
            AllNeighbours(neighbour);
            return;
        }
        
        int nsidenodes = fElement->NSideNodes(fSide);
        std::vector<GeoElementSide> GeoElSideSet;
        std::vector<int> Set[27];
        VecInt nodeindexes(nsidenodes,0.);
        for (int in=0; in<nsidenodes; in++) {
            
            nodeindexes[in] = fElement->NodeIndex(fElement->SideNodeIndex(fSide, in));
            int locnod = fElement->SideNodeIndex(fSide, in);
            GeoElSideSet.resize(0);
            GeoElementSide locside(fElement,locnod);
            locside.AllNeighbours(GeoElSideSet);
            int nel = GeoElSideSet.size();
            for (int el=0; el<nel; el++) {
                Set[in].push_back(GeoElSideSet[el].Element()->GetIndex());
            }
            std::sort(Set[in].begin(),Set[in].end());
        }
        std::vector<int> result;
        switch (nsidenodes) {
            case 1:
            {
                result = Set[0];
            }
                break;
            case 2:
            {
                std::set_intersection(Set[0].begin(), Set[0].end(), Set[1].begin(), Set[1].end(), std::back_inserter(result));
            }
                break;
            case 3:
            {
                std::vector<int> inter1;
                std::set_intersection(Set[0].begin(), Set[0].end(), Set[1].begin(), Set[1].end(), std::back_inserter(inter1));
                std::set_intersection(inter1.begin(), inter1.end(), Set[2].begin(), Set[2].end(), std::back_inserter(result));
            }
                break;
            case 4:
            {
                std::vector<int> inter1, inter2;
                std::set_intersection(Set[0].begin(), Set[0].end(), Set[2].begin(), Set[2].end(), std::back_inserter(inter1));
                if(inter1.size()==0) break;
                std::set_intersection(Set[1].begin(), Set[1].end(), Set[3].begin(), Set[3].end(), std::back_inserter(inter2));
                if(inter2.size()==0) break;
                std::set_intersection(inter1.begin(), inter1.end(), inter2.begin(), inter2.end(), std::back_inserter(result));
            }
                break;
            default:
            {
                std::vector<int> inter1, inter2;
                inter1 = Set[0];
                for(int in=0; in<nsidenodes-1; in++) {
                    inter2.resize(0);
                    std::set_intersection(inter1.begin(), inter1.end(), Set[in+1].begin(), Set[in+1].end(), std::back_inserter(inter2));
                    if(inter2.size() == 0) break;
                    inter1 = inter2;
                }
                result = inter2;
            }
                
        }
        int el,nel = result.size();
        GeoMesh * geoMesh = fElement->GetMesh();
        for(el=0; el<nel; el++) {
            GeoElement * gelResult = geoMesh->Element(result[el]);
            int whichSd = gelResult->WhichSide(nodeindexes);
            if(whichSd > 0)
            {
                neighbour.push_back(GeoElementSide(gelResult, whichSd));
            }
        }
        
    }


