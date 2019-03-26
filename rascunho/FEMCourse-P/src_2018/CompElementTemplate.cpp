//
//  CompElementTemplate.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "CompElementTemplate.h"
#include "CompElement.h"
#include "CompMesh.h"
#include "GeoMesh.h"
#include "GeoElementTemplate.h"
#include "GeoElement.h"
#include "Shape1d.h"
#include "ShapeQuad.h"
#include "ShapeTetrahedron.h"
#include "ShapeTriangle.h"
#include "CompMesh.h"
#include "tpanic.h"
#include "DataTypes.h"
#include "MathStatement.h"
#include "DOF.h"
#include "GeoNode.h"


    template<class Shape>
    CompElementTemplate<Shape>::CompElementTemplate() : CompElement(), dofindexes(0),intrule(0){
        
    }

    template<class Shape>
    CompElementTemplate<Shape>::CompElementTemplate(int64_t ind, CompMesh *cmesh,  GeoElement *geo) : CompElement(ind,cmesh,geo){
        int Nelem = cmesh->GetElementVec().size();
        cmesh->SetNumberElement(Nelem);
        //Nelem+=1;
        //ind = Nelem-1;
        cmesh->SetElement(ind, this);
        geo->SetReference(this);
        int nsidesw = geo->NSides();
        this->SetIndex(ind);
        int order = cmesh->GetDefaultOrder();
        intrule.SetOrder(order*2);
        SetIntRule(&intrule);
        MathStatement *mat = cmesh->GetMath(ind);
        this->SetStatement(mat);
        int nsides = Shape::nSides;
        SetNDOF(nsides);
        for (int is=0; is<nsides; is++) {
            GeoElementSide gelside(GetGeoElement(),is);
            GeoElementSide neighbour = gelside.Neighbour();
            while(neighbour != gelside)
            {
                if (neighbour.Element()->GetReference()) {
                    break;
                }
                neighbour = neighbour.Neighbour();
            }
            if(neighbour != gelside)
            {
                int indexneigh = neighbour.Element()->GetIndex();
                if (indexneigh<ind) {
                    CompElement *cel = neighbour.Element()->GetReference();
                    this->SetDOFIndex(is, cel->GetDOFIndex(neighbour.Side()));
                }else{
                    
                    int order = cmesh->GetDefaultOrder();
                    int nshape = Shape::NShapeFunctions(is,order);
                    int nstate = mat->NState();
                    int64_t ndof = cmesh->GetNumberDOF();
                    cmesh->SetNumberDOF(ndof+1);
                    DOF dof;
                    dof.SetNShapeStateOrder(nshape,nstate,order);
                    cmesh->SetDOF(ndof, dof);
                    this->SetDOFIndex(is, ndof);
                }

            }
            else
            {
                int order = cmesh->GetDefaultOrder();
                int nshape = Shape::NShapeFunctions(is,order);
                int nstate = mat->NState();
                int64_t ndof = cmesh->GetNumberDOF();
                cmesh->SetNumberDOF(ndof+1);
                DOF dof;
                dof.SetNShapeStateOrder(nshape,nstate,order);
                cmesh->SetDOF(ndof, dof);
                this->SetDOFIndex(is, ndof);
            }
        }
    }

    template<class Shape>
    CompElementTemplate<Shape>::CompElementTemplate(const CompElementTemplate &copy) : CompElement(copy){
        dofindexes=copy.dofindexes;
        intrule=copy.intrule;
    }

    template<class Shape>
    CompElementTemplate<Shape> &CompElementTemplate<Shape>::operator=(const CompElementTemplate &copy){
        dofindexes=copy.dofindexes;
        intrule=copy.intrule;
        return *this;
    }

    template<class Shape>
    CompElementTemplate<Shape>::~CompElementTemplate(){
        
    }

    template<class Shape>
    void CompElementTemplate<Shape>::SetNDOF(int64_t ndof){
        dofindexes.resize(ndof);
    }

    template<class Shape>
    void CompElementTemplate<Shape>::SetDOFIndex(int i, int64_t dofindex){
        dofindexes[i]=dofindex;
    }

    template<class Shape>
    int64_t CompElementTemplate<Shape>::GetDOFIndex(int i){
        return dofindexes[i];
    }

    template<class Shape>
    CompElement * CompElementTemplate<Shape>::Clone() const{
        //CompElementTemplate<Shape>* result = new CompElementTemplate<Shape>(*this);
        //return result;
    }

    template<class Shape>
    void  CompElementTemplate<Shape>::ShapeFunctions(const VecDouble &intpoint, VecDouble &phi, Matrix &dphi) const{
        VecInt orders(NDOF());
        int ndof = NDOF();
        CompMesh *cmesh =this->GetCompMesh();
        for (int ic=0; ic<ndof; ic++) {
            orders[ic]=cmesh->GetDOF(ic).GetOrder();
        }
        Shape::Shape(intpoint, orders, phi, dphi);

    }
    template<class Shape>
    void CompElementTemplate<Shape>::GetMultiplyingCoeficients(VecDouble &coefs) const{
        
        CompMesh *cmesh = this->GetCompMesh();
        int neq = cmesh->Solution().size();
        VecInt iGlob;
        coefs.resize(0);
        
        int ndofel = NDOF();
        int indexdof = 0;
        for (int idof =0; idof<ndofel; idof++) {
            int idcon = dofindexes[idof];
            DOF dof = cmesh->GetDOF(idcon);
            int nshape = dof.GetNShape();
            int nstat = dof.GetNState();
            for(int i=0; i<nshape*nstat; i++) {
                iGlob.resize(indexdof+1);
                coefs.resize(indexdof+1);
                iGlob[indexdof] = dof.GetFirstEquation()+i;
                coefs[indexdof] = cmesh->Solution()[iGlob[indexdof]];
                indexdof++;
            }
        }
           
    }

    template<class Shape>
    int CompElementTemplate<Shape>::NShapeFunctions() const{
        //CompMesh cmesh = *this->GetCompMesh();
        int ndof = NDOF();
        int nshape = 0;
        for (int ic=0; ic<ndof; ic++) {
            nshape +=NShapeFunctions(ic);
        }
        return nshape;
    }

    template<class Shape>
    int CompElementTemplate<Shape>::NDOF() const{
        return dofindexes.size();
    }
    
    /// returns the number of shape functions stored in the DOF data structure
    template<class Shape>
    int CompElementTemplate<Shape>::NShapeFunctions(int doflocindex) const{
        CompMesh cmesh = *this->GetCompMesh();
        DOF dofex = cmesh.GetDOF(doflocindex);
        return cmesh.GetDOF(doflocindex).GetNShape();
    }
    
    /// uses the Shape template class to compute the number of shape functions
    template<class Shape>
    int CompElementTemplate<Shape>::ComputeNShapeFunctions(int doflocindex, int order){
        dofindexes.resize(doflocindex+1);
        dofindexes[doflocindex]=doflocindex;
        return Shape::NShapeFunctions(doflocindex,order);
    }

    template<class Shape>
    void CompElementTemplate<Shape>::Print(std::ostream &out){
        out << __PRETTY_FUNCTION__ << std::endl;

        out << "\nOutput for a computable element index: " << GetIndex();
        //  if(geoel)
        //  {
        //     out << "\nCenter coordinate: ";
        //  }
        if(this->GetStatement())
        {
            out << "\nMaterial id " << this->GetStatement()->GetMatID() << "\n";
        }
        else {
            out << "\nNo material\n";
        }

        out << "Number of connects = " << this->NDOF();
        out<< "\nConnect indexes : ";
        int nod;
        for(nod=0; nod< NDOF(); nod++)
        {
            out << GetDOFIndex(nod) <<  ' ' ;
        }
        out << std::endl;

        out << "Index = " << GetIndex();
        out << " - Center coordinate:\n";

        for (int i=0;i< Shape::nCorners;i++)
        {
            VecDouble center(3,0.);
            GeoElement *geo = GetGeoElement();

            for (int j=0;j<3;j++) center[j] = geo->GetMesh()->Node(i).Coord(j);
            out << "[" <<  i << "] ";
            for (int j=0;j<3;j++){
                out << center[j] << ", ";
            }
        }

        out << std::endl;
        int nconects=this->NDOF();
        out << "Number of connects = " << this->NDOF() << " Connect indexes/NShape : ";
        
        for(nod=0; nod< nconects; nod++)
        {
            DOF &c = GetCompMesh()->GetDOF(nod);
            out << GetDOFIndex(nod) <<  '/' << c.GetNShape() << ' ' ;
        }
        out << std::endl;

        CompMesh *Pcmesh = GetCompMesh();


        MathStatement * material = this->GetStatement();
        if(!material)
        {
            out << " no material " << std::endl;
        }

        if (material) {
            out << "material id " << material->GetMatID() << std::endl;
        }
        IntRule *intr = GetIntRule();
        out << "Integration orders : \t";
        int introrder = Pcmesh->GetDefaultOrder()*2;
        out << "material id " << introrder << std::endl;
        intr->Print(out);
        
        out << std::endl;
    }




template class CompElementTemplate<Shape1d>;
template class CompElementTemplate<ShapeQuad>;
template class CompElementTemplate<ShapeTriangle>;
template class CompElementTemplate<ShapeTetrahedron>;
