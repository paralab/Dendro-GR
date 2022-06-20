//
// Created by milinda on 7/25/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains a basic structure to store the black hole parameters.
*/
//

#ifndef SFCSORTBENCH_BH_H
#define SFCSORTBENCH_BH_H

#include "point.h"
namespace bssn
{

    struct BH
    {
        private:
            /**@brief mass of the black hole*/
            double m_uiMass;
            /**@brief coordinate of the black hole*/
            Point m_uiCoord;
            /**@brief */
            Point m_uiV;
            /**@brief spin parameter of the black hole. */
            double m_uiSpin;
            /**@brief spin parameter of the black hole. */
            double m_uiSpinTheta;
            /**@brief spin parameter of the black hole. */
            double m_uiSpinPhi;

        public:
            /**@brief default constructor*/
            BH()
            {
                m_uiMass=0.0;
                m_uiCoord=Point(0.0,0.0,0.0);
                m_uiV=Point(0.0,0.0,0.0);
                m_uiSpin=0;
                m_uiSpinTheta=0;
                m_uiSpinPhi=0;

            }

            /**@brief constructor to BH structure*/
            BH(double pMass,double pCx,double pCy,double pCz,double pVx,double pVy,double pVz,double pSpin,double pSpinTheta, double pSpinPhi)
            {
                m_uiMass=pMass;
                m_uiCoord=Point(pCx,pCy,pCz);
                m_uiV=Point(pVx,pVy,pVz);
                m_uiSpin=pSpin;
                m_uiSpinTheta=pSpinTheta;
                m_uiSpinPhi=pSpinPhi;
            }

            /**@brief returns mass*/
            inline double getBHMass() const {return m_uiMass;}
            /** @brief returns BH coordinates*/
            inline Point getBHCoord() const { return m_uiCoord;}
            /**@brief returns BH x coordinate*/
            inline double getBHCoordX() const {return m_uiCoord.x();}
            /**@brief returns BH y coordinate*/
            inline double getBHCoordY() const {return m_uiCoord.y();}
            /**@brief returns BH z coordinate*/
            inline double getBHCoordZ() const {return m_uiCoord.z();}
            /**@brief returns m_uiV*/
            inline Point getV() const {return m_uiV;}
            /**@brief returns m_uiV x coordinate*/
            inline double getVx() const {return m_uiV.x();}
            /**@brief returns m_uiV y coordinate*/
            inline double getVy() const {return m_uiV.y();}
            /**@brief returns m_uiV z coordinate*/
            inline double getVz() const {return m_uiV.z();}
            /**@brief returns spin of BH*/
            inline double getBHSpin() const {return m_uiSpin;}
            /**@brief returns spin(theta) of BH */
            inline double getBHSpinTheta() const {return m_uiSpinTheta;}
            /**@brief returns spin(phi) of BH */
            inline double getBHSpinPhi() const {return m_uiSpinPhi;}






    };



}// end of namespace bssn


#endif //SFCSORTBENCH_BH_H
