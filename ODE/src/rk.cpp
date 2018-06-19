//
// Created by milinda on 5/25/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//

#include "rk.h"

namespace ode
{
    namespace solver
    {
        RK::RK(ot::Mesh *pMesh, double pTBegin, double pTEnd,double pTh)
        {
            m_uiMesh=pMesh;
            m_uiComm=m_uiMesh->getMPIGlobalCommunicator();
            m_uiOrder=pMesh->getElementOrder();
            m_uiTimeBegin=pTBegin;
            m_uiTimeEnd=pTEnd;
            m_uiT_h=pTh;
            m_uiT_h_prev=m_uiT_h;

            m_uiCurrentStep=0;
            m_uiCurrentTime=m_uiTimeBegin;
            m_uiNrp=m_uiOrder+1;



        }

        RK::~RK()
        {

        }




    }


}


