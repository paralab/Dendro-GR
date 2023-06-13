//
// Created by milinda on 11/22/18.
//

#ifndef DENDRO_5_0_HEATVEC_H
#define DENDRO_5_0_HEATVEC_H

#include "oda.h"
#include "feVector.h"

namespace HeatEq
{
    class HeatVec : public feVector<HeatVec>{

    private:

        double * imV1;
        double * imV2;

    public:
        HeatVec(ot::DA* da,unsigned int dof=1);
        ~HeatVec();

        /**@biref elemental compute vec for rhs*/
        virtual void elementalComputVec(const VECType* in,VECType* out, double*coords=NULL,double scale=1.0);


        bool preComputeVec(const VECType* in,VECType* out, double scale=1.0);

        bool postComputeVec(const VECType* in,VECType* out, double scale=1.0);

        /**@brief octree grid x to domin x*/
        double gridX_to_X(double x);
        /**@brief octree grid y to domin y*/
        double gridY_to_Y(double y);
        /**@brief octree grid z to domin z*/
        double gridZ_to_Z(double z);



    };
}



#endif //DENDRO_5_0_HEATVEC_H
