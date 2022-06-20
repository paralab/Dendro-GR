#include "enutsOp.h"
namespace ts
{
    
    ENUTSOp::ENUTSOp(ETSType type)
    {
        m_uiType = type;

        if(type == ETSType::RK3)
        {
            m_uiNumStages = 3;
            const unsigned int d1 = m_uiNumStages;
            const unsigned int d2 = (m_uiNumStages*m_uiNumStages);


            const DendroScalar ETS_U[] =        { 0.0     ,0.0      , 0.0, 
                                                  1.0     ,0.0      , 0.0,
                                                  1.0/4.0 ,1.0/4.0  , 0.0 };

            // coefficient matrix, (picewise linear approx) 
            const DendroScalar ETS_C[] =        { 1.0     , 0.0      , 0.0, 
                                                  1.0     , 1.0      , 0.0,
                                                  1.0     , 1.0/2.0  , 1.0/4.0 };

            const DendroScalar ETS_InvC[] =     { 1.0     , 0.0      , 0.0, 
                                                  -1.0     , 1.0      , 0.0,
                                                  -2.0     , -2.0      , 4.0 };


            
            m_uiAij.resize(d2);
            m_uiBij.resize(d2);
            m_uiCij.resize(d2);
            m_uiInvCij.resize(d2);
            m_uiPi.resize(d1);
            m_uiInvPi.resize(d1);
            m_uiMij.resize(d2);
            m_uiVin.resize(d1);
            m_uiVout.resize(d1);

            std::memcpy(m_uiAij.data(), ETS_U, sizeof(DendroScalar)*d2);
            std::memcpy(m_uiCij.data(), ETS_C, sizeof(DendroScalar)*d2);
            std::memcpy(m_uiInvCij.data(), ETS_InvC, sizeof(DendroScalar)*d2);
            

            
        
        }else if( type == ETSType::RK4)
        {
            m_uiNumStages = 4; 
            const unsigned int d1 = m_uiNumStages;
            const unsigned int d2 = (m_uiNumStages*m_uiNumStages);
            
            const DendroScalar ETS_U[] =        { 0.0   , 0.0       , 0.0   , 0.0,
                                                1.0/2.0 , 0.0       , 0.0   , 0.0,
                                                0.0     , 1.0/2.0   , 0.0   , 0.0,
                                                0.0     , 0.0       , 1.0   , 0.0 };

            // coefficient matrix, (picewise linear approx) 
            const DendroScalar ETS_C[] =        { 1.0     , 0.0          , 0.0,      0.0,
                                                  1.0     , 1.0/2.0      , 0.0,      0.0,
                                                  1.0     , 1.0/2.0      , 1.0/4.0,  0.0,
                                                  1.0     , 1.0          , 1.0/2.0,  1.0/4.0};


            const DendroScalar ETS_InvC[] =     { 1.0     , 0.0          , 0.0,      0.0,
                                                  -2.0    , 2.0          , 0.0,      0.0,
                                                  0.0     , -4.0         , 4.0,      0.0,
                                                  4.0     , 0.0          , -8.0,     4.0};


            
            m_uiAij.resize(d2);
            m_uiBij.resize(d2);
            m_uiCij.resize(d2);
            m_uiInvCij.resize(d2);
            m_uiPi.resize(d1);
            m_uiInvPi.resize(d1);
            m_uiMij.resize(d2);
            m_uiVin.resize(d1);
            m_uiVout.resize(d1);

            std::memcpy(m_uiAij.data(), ETS_U, sizeof(DendroScalar)*d2);
            std::memcpy(m_uiCij.data(), ETS_C, sizeof(DendroScalar)*d2);
            std::memcpy(m_uiInvCij.data(), ETS_InvC, sizeof(DendroScalar)*d2);

        }else if ( type == ETSType::RK5)
        {
            // these coefficients looks wrong. - Milinda. (need to fix those and enable the rk5 method. )
            return;
            /*m_uiNumStages = 5;
            const unsigned int d1 = m_uiNumStages;
            const unsigned int d2 = (m_uiNumStages*m_uiNumStages);
            
            const DendroScalar ETS_U[] =        {   0.0       ,  0.0        , 0.0       , 0.0        , 0.0,
                                                    1.0/8.0   ,  1.0/8.0    , 0.0       , 0.0        , 0.0,
                                                    0.0       , -1.0/2.0    , 1.0       , 0.0        , 0.0,
                                                    3.0/16.0  ,  0.0        , 0.0       , 9.0/16.0   , 0.0,
                                                    -3.0/7.0  ,  2.0/7.0    , 12.0/7.0  , -12.0/7.0  , 8.0/7.0};
            
            m_uiAij.resize(d2);
            m_uiBij.resize(d2);
            m_uiCij.resize(d2);
            m_uiPi.resize(d1);

            std::memcpy(m_uiAij.data(), ETS_U, sizeof(DendroScalar)*(d2));*/

        }

        //std::cout<<" m_uinum stages : "<<m_uiNumStages<<std::endl;
        return;

    }

    void ENUTSOp::Pdt(DendroScalar dt, unsigned int rk_s)
    {
        
        m_uiPi[0]    = 1.0;
        m_uiInvPi[0] = 1.0;
        for(unsigned int s=2; s <= rk_s; s++)
        {
            m_uiPi[s-1] = m_uiPi[s-2]*dt;
            m_uiInvPi[s-1] = 1.0/m_uiPi[s-1];
        }
            

        
    }



    void ENUTSOp::Bdt(DendroScalar dt, unsigned int rk_s)
    {
        
        if(m_uiType == ETSType::RK3)
        {

            /*const DendroScalar ETS_B[] =        { 1.0     , 1.0      , 1.0/2.0, 
                                                  0.0     , 1.0      , 1.0,
                                                  0.0     , 0.0      , 1.0  };*/

            m_uiBij[0 * m_uiNumStages + 0 ] = 1.0;
            m_uiBij[0 * m_uiNumStages + 1 ] = 1.0 * dt;
            m_uiBij[0 * m_uiNumStages + 2 ] = (1.0/2.0) *dt*dt;

            m_uiBij[1 * m_uiNumStages + 0 ] = 0.0;
            m_uiBij[1 * m_uiNumStages + 1 ] = 1.0;
            m_uiBij[1 * m_uiNumStages + 2 ] = 1.0 * dt;


            m_uiBij[2 * m_uiNumStages + 0 ] = 0.0;
            m_uiBij[2 * m_uiNumStages + 1 ] = 0.0;
            m_uiBij[2 * m_uiNumStages + 2 ] = 1.0;




        }else if(m_uiType == ETSType::RK4)
        {

            /*const DendroScalar ETS_B[] =       { 1.0     , 1.0      , 1.0/2.0 , 1.0/6.0,
                                                  0.0     , 1.0      , 1.0     , 1.0/2.0,    
                                                  0.0     , 0.0      , 1.0     , 1.0    ,
                                                  0.0     , 0.0      , 0.0     , 1.0    };*/
            
            m_uiBij[0 * m_uiNumStages + 0 ] = 1.0;
            m_uiBij[0 * m_uiNumStages + 1 ] = 1.0 * dt;
            m_uiBij[0 * m_uiNumStages + 2 ] = (1.0/2.0) * dt*dt;
            m_uiBij[0 * m_uiNumStages + 3 ] = (1.0/6.0) * dt*dt*dt;

            m_uiBij[1 * m_uiNumStages + 0 ] = 0.0;
            m_uiBij[1 * m_uiNumStages + 1 ] = 1.0;
            m_uiBij[1 * m_uiNumStages + 2 ] = 1.0 *dt;
            m_uiBij[1 * m_uiNumStages + 3 ] = (1.0/2.0)*dt*dt;


            m_uiBij[2 * m_uiNumStages + 0 ] = 0.0;
            m_uiBij[2 * m_uiNumStages + 1 ] = 0.0;
            m_uiBij[2 * m_uiNumStages + 2 ] = 1.0;
            m_uiBij[2 * m_uiNumStages + 3 ] = 1.0 * dt;

            m_uiBij[3 * m_uiNumStages + 0 ] = 0.0;
            m_uiBij[3 * m_uiNumStages + 1 ] = 0.0;
            m_uiBij[3 * m_uiNumStages + 2 ] = 0.0;
            m_uiBij[3 * m_uiNumStages + 3 ] = 1.0;



        }
        
    }


    void ENUTSOp::Cfc(DendroScalar** out, const DendroScalar** in, const unsigned int* sz, unsigned int rk_s, DendroScalar dt_c, DendroScalar dt_f, DendroScalar dt, unsigned int dof)
    {

        const unsigned int nx = sz[0];
        const unsigned int ny = sz[1];
        const unsigned int nz = sz[2];
        const unsigned int NN = nx*ny*nz;
        const unsigned int d = m_uiNumStages;
        
        Bdt(dt,rk_s);
        
        const double inv_dtf = (1.0/dt_f);
        const double dtc_s = pow(dt_c,rk_s-1);

        // B x P^{-1} x C^{-1}
        for(unsigned int i=0; i < d; i++)
            for(unsigned int j=0; j < d; j++)
            {
                m_uiMij[i*d + j] = 0.0;
                double tmp = 1.0;
                for(unsigned int k=0; k < d; k++)
                {
                    m_uiMij[i*d + j] += m_uiBij[i*d+k] * m_uiInvCij[k*d + j]*tmp;
                    tmp = tmp*inv_dtf; 
                }
                    
            }
                
        std::swap(m_uiMij,m_uiBij);

        // C x P x B x P^{-1} x C^{-1}
        for(unsigned int i=0; i < d; i++)
            for(unsigned int j=0; j < d; j++)
            {
                m_uiMij[i*d + j] = 0.0;
                double tmp=1.0;
                for(unsigned int k=0; k < d; k++)
                {
                    m_uiMij[i*d + j] += m_uiCij[i*d+k] * m_uiBij[k*d + j]*tmp;
                    tmp = tmp * dt_c; 
                }
                    
            }
        
        // std::cout<<"M f2c : "<<std::endl; 
        // for(unsigned int i=0; i < d; i++)
        // {
        //     for(unsigned int j=0; j < d; j++)
        //     {
        //         std::cout<<m_uiMij[i*d + j]<<",";
        //     }
        //     std::cout<<std::endl;
        // }
            

        

        // now M = C x B x C^{-1}
       

        for(unsigned int v=0; v < dof; v++)
        {
        
            for(unsigned int n=0; n < NN; n++)
            {
                for(unsigned int s=0; s <rk_s; s++)
                    m_uiVin[s] = in[v*rk_s + (s)][n];
                   
                m_uiVout[(rk_s-1)]=0;
                for(unsigned int k=0; k < rk_s; k++)
                    m_uiVout[rk_s-1] += m_uiMij[(rk_s-1)*d + k ] * m_uiVin[k];

                //m_uiVout[rk_s-1] *=dtc_s;
                //out[v*rk_s + rk_s-1][n] = in[v*rk_s + (rk_s-1)][n]; //m_uiVout[rk_s-1];
                // if(m_uiVout[rk_s-1]>3)
                // {  
                //     std::cout<<"vout c2f: "<<m_uiVout[rk_s-1]<<std::endl;;
                //     for(unsigned int s=0; s <rk_s; s++)
                //      std::cout<<" vin : "<<m_uiVin[s]<<std::endl;

                //     std::cout<<"dt_s: "<<dtc_s<<std::endl;
                //     std::cout<<"dt_f: "<<dt_f<<" dt : "<<dt <<" dt_c: "<<dt_c<<std::endl;
                //     exit(0);

                // }
                out[v*rk_s + rk_s-1][n] = m_uiVout[rk_s-1];
                //printf("out s %d n %d  val %f \n", rk_s-1,n, out[v*rk_s + rk_s-1][n]);

            }

        }

    }


    void ENUTSOp::Ccf(DendroScalar** out, const DendroScalar** in, const unsigned int* sz, unsigned int rk_s, DendroScalar dt_c, DendroScalar dt_f, DendroScalar dt, unsigned int dof)
    {

        const unsigned int nx = sz[0];
        const unsigned int ny = sz[1];
        const unsigned int nz = sz[2];

        const unsigned int NN = nx*ny*nz;
        const unsigned int d = m_uiNumStages;
        
        Bdt(dt,rk_s);
        // std::cout<<"Bij : "<<dt<<std::endl; 
        // for(unsigned int i=0; i < d; i++)
        // {
        //     for(unsigned int j=0; j < d; j++)
        //     {
        //         std::cout<<m_uiBij[i*d + j]<<",";
        //     }
        //     std::cout<<std::endl;
        // }

        const double inv_dtc = (1.0/dt_c);
        const double dtf_s = pow(dt_f,rk_s-1);

        // B x P^{-1} x C^{-1}
        for(unsigned int i=0; i < d; i++)
            for(unsigned int j=0; j < d; j++)
            {
                m_uiMij[i*d + j] = 0.0;
                double tmp = 1.0;
                for(unsigned int k=0; k < d; k++)
                {
                    m_uiMij[i*d + j] += m_uiBij[i*d+k] * m_uiInvCij[k*d + j]*tmp;
                    tmp =tmp *inv_dtc; 
                }
                    
            }
                
        std::swap(m_uiMij,m_uiBij);

        // C x P x B x P^{-1} x C^{-1}
        for(unsigned int i=0; i < d; i++)
            for(unsigned int j=0; j < d; j++)
            {
                m_uiMij[i*d + j] = 0.0;
                double tmp=1.0;
                for(unsigned int k=0; k < d; k++)
                {
                    m_uiMij[i*d + j] += m_uiCij[i*d+k] * m_uiBij[k*d + j]*tmp; 
                    tmp=tmp*dt_f;
                }
                    
            }
        
        // now M = C x B x C^{-1}
        


        // std::cout<<"M c2f : "<<dt<<std::endl; 
        // for(unsigned int i=0; i < d; i++)
        // {
        //     for(unsigned int j=0; j < d; j++)
        //     {
        //         std::cout<<m_uiMij[i*d + j]<<",";
        //     }
        //     std::cout<<std::endl;
        // }

        for(unsigned int v=0; v < dof; v++)
        {
        
            for(unsigned int n=0; n < NN; n++)
            {
                for(unsigned int s=0; s <rk_s; s++)
                    m_uiVin[s] = in[v*rk_s + s][n];
                   
                m_uiVout[(rk_s-1)]=0;
                for(unsigned int k=0; k < rk_s; k++)
                    m_uiVout[rk_s-1] += m_uiMij[(rk_s-1)*d + k ] * m_uiVin[k];

                // if(m_uiVout[rk_s-1]>3)
                // {  
                //     std::cout<<"vout c2f: "<<m_uiVout[rk_s-1]<<std::endl;;
                //     for(unsigned int s=0; s <rk_s; s++)
                //      std::cout<<" vin : "<<m_uiVin[s]<<std::endl;

                //     std::cout<<"dt_s: "<<dtf_s<<std::endl;
                //     std::cout<<"dt_f: "<<dt_f<<" dt : "<<dt <<" dt_c: "<<dt_c<<std::endl;
                //     exit(0);

                // }
                    
                //out[v*rk_s + rk_s-1][n] = in[v*rk_s + (rk_s-1)][n]; // m_uiVout[rk_s-1];
                out[v*rk_s + rk_s-1][n] = m_uiVout[rk_s-1];

            }

        }

    }


    void ENUTSOp::coarser_finer_ut_correction(DendroScalar** out, const DendroScalar** in, const unsigned int* sz, DendroScalar dt_c, DendroScalar dt_f, DendroScalar dt, unsigned int dof)
    {
        const unsigned int nx = sz[0];
        const unsigned int ny = sz[1];
        const unsigned int nz = sz[2];
        const unsigned int NN = nx*ny*nz;
        const unsigned int d = m_uiNumStages;
        const unsigned int totalVars = d + 1;
        
        
        // note: when dt=0, then the B = I;
        Bdt(dt,m_uiNumStages);

        const double inv_dtc = (1.0/dt_c);
        
        // B x P^{-1} x C^{-1}
        for(unsigned int i=0; i < d; i++)
            for(unsigned int j=0; j < d; j++)
            {
                m_uiMij[i*d + j] = 0.0;
                double tmp = 1.0;
                for(unsigned int k=0; k < d; k++)
                {
                    m_uiMij[i*d + j] += m_uiBij[i*d+k] * m_uiInvCij[k*d + j]*tmp;
                    tmp = tmp*inv_dtc; 
                }
                    
            }


        // M_ij has the correct operator. 
        // Vin  : 0, 1, 2, .. m_uiNumStages-1 these are RK stages, m_uiNumStages contains the U_f or U_c. 
        // Vout : 0, 1, 2, .. m_uiNumStages-1 contains the correct stages to compute correct U form U_f or U_c
        for(unsigned int v=0; v < dof; v++)
        {
        
            for(unsigned int n=0; n < NN; n++)
            {
                for(unsigned int s=0; s <d; s++)
                {
                    m_uiVin[s]  = in[v*(d+1) + (s)][n];
                    m_uiVout[s] = 0;
                }

                for(unsigned int si=0; si < m_uiNumStages; si++)
                 for(unsigned int sj=0; sj < m_uiNumStages; sj++)
                    m_uiVout[si] += m_uiMij[si*d + sj ] * m_uiVin[sj];
                
                double dt_by_fac_p = 1.0;
                
                out[v*(d+1) + m_uiNumStages][n] = in[v*(d+1) + m_uiNumStages][n];
            
                for(unsigned int si=0; si < m_uiNumStages; si++)
                {
                    dt_by_fac_p *= (dt_f / ((double)(si+1)));
                    out[v*(d+1) + m_uiNumStages][n] -= (dt_by_fac_p*m_uiVout[si]);  
                }
                
                //printf("out s %d n %d  val %f \n", rk_s-1,n, out[v*rk_s + rk_s-1][n]);
                //printf("coarset to finer: in value: %f out value: %f \n",in[v*(d+1) + m_uiNumStages][n],out[v*(d+1) + m_uiNumStages][n]);
            }
        }


        return;

    }

    void ENUTSOp::finer_coarser_ut_correction(DendroScalar** out, const DendroScalar** in, const unsigned int* sz, DendroScalar dt_c, DendroScalar dt_f, DendroScalar dt, unsigned int dof)
    {
        //std::cout<<"caling finer_coarser"<<std::endl;
        const unsigned int nx = sz[0];
        const unsigned int ny = sz[1];
        const unsigned int nz = sz[2];
        const unsigned int NN = nx*ny*nz;
        const unsigned int d = m_uiNumStages;
        const unsigned int totalVars = d + 1;
        
        
        // note: when dt=0, then the B = I;
        Bdt(dt,m_uiNumStages);

        const double inv_dtf = (1.0/dt_f);
        
        // B x P^{-1} x C^{-1}
        for(unsigned int i=0; i < d; i++)
            for(unsigned int j=0; j < d; j++)
            {
                m_uiMij[i*d + j] = 0.0;
                double tmp = 1.0;
                for(unsigned int k=0; k < d; k++)
                {
                    m_uiMij[i*d + j] += m_uiBij[i*d+k] * m_uiInvCij[k*d + j]*tmp;
                    tmp = tmp*inv_dtf; 
                }
                    
            }


        // M_ij has the correct operator. 
        // Vin  : 0, 1, 2, .. m_uiNumStages-1 these are RK stages, m_uiNumStages contains the U_f or U_c. 
        // Vout : 0, 1, 2, .. m_uiNumStages-1 contains the correct stages to compute correct U form U_f or U_c
        for(unsigned int v=0; v < dof; v++)
        {
        
            for(unsigned int n=0; n < NN; n++)
            {
                for(unsigned int s=0; s <d; s++)
                {
                    m_uiVin[s]  = in[v*(d+1) + (s)][n];
                    m_uiVout[s] = 0;
                }

                for(unsigned int si=0; si < m_uiNumStages; si++)
                 for(unsigned int sj=0; sj < m_uiNumStages; sj++)
                    m_uiVout[si] += m_uiMij[si*d + sj ] * m_uiVin[sj];
                
                double dt_by_fac_p = 1.0;
                
                out[v*(d+1) + m_uiNumStages][n] = in[v*(d+1) + m_uiNumStages][n];
            
                for(unsigned int si=0; si < m_uiNumStages; si++)
                {
                    dt_by_fac_p *= (dt_f / ((double)(si+1)));
                    out[v*(d+1) + m_uiNumStages][n] += (dt_by_fac_p*m_uiVout[si]);  
                }
                //printf("out s %d n %d  val %f \n", rk_s-1,n, out[v*rk_s + rk_s-1][n]);
                //printf("finer to coarser: in value: %f out value: %f \n",in[v*(d+1) + m_uiNumStages][n],out[v*(d+1) + m_uiNumStages][n]);

            }

        }

        return ;
        
    }

}