#include "fluidMath.h"
#include "root.h"
#include "geom.h"
#include "math.h"

using namespace fluid;
namespace fluid_math
{


  Error prim_to_con(const ot::Mesh * pMesh,double**v,double **u)
  {
    
    double gd[3][3], gu[3][3], detg,sdetg;
    double pos[3];

    unsigned int ownerID, ii_x, jj_y, kk_z;
    const ot::TreeNode* allElements = &(*(pMesh->getAllElements().begin()));
    const unsigned int* cg_to_dg = &(*(pMesh->getCG2DGMap().begin()));
    const unsigned int elementOrder = pMesh->getElementOrder();
    
    unsigned int lookUp;
    Error status;
    unsigned int failCount=0;
    unsigned int sz;
    for(unsigned int node=pMesh->getNodeLocalBegin();node<pMesh->getNodeLocalEnd();node++)
    {

      lookUp=cg_to_dg[node];
      pMesh->dg2eijk(lookUp,ownerID,ii_x,jj_y,kk_z);

      sz=(1u<<(m_uiMaxDepth-allElements[ownerID].getLevel()));

      pos[0]=GRIDX_TO_X(allElements[ownerID].minX() + ii_x*(sz/((double)elementOrder)));
      pos[1]=GRIDY_TO_Y(allElements[ownerID].minY() + jj_y*(sz/((double)elementOrder)));
      pos[2]=GRIDZ_TO_Z(allElements[ownerID].minZ() + kk_z*(sz/((double)elementOrder)));

      metric_vars(gu, &detg, gd, pos);
      sdetg = sqrt(detg);

      double u_pp[FLUID_NUM_EVOL_VARS];
      double v_pp[FLUID_NUM_PRIM_VARS]={v[PVAR::V_RHO][node], v[PVAR::V_VX][node],v[PVAR::V_VY][node],v[PVAR::V_VZ][node],v[PVAR::V_P][node],v[PVAR::V_U1][node],v[PVAR::V_U2][node],v[PVAR::V_U3][node]};

      status=prim_to_con_pt(v_pp,u_pp,gd,sdetg);

      for(unsigned int var=0;var<FLUID_NUM_EVOL_VARS;var++)
        u[var][node]=u_pp[var];


      if(status != Error::FLUID_SOLVER_SUCCESS)
        failCount++;

    }

    if(failCount==0)
      return Error::FLUID_SOLVER_SUCCESS;
    else
    {
      std::cout<<"[Error]: primitive to conserved variable solver failed. "<<std::endl;
      print_stacktrace();  
      return Error::FLUID_PRIM_TO_CONS_FAILED;
    }

      
    

  }
 
  Error prim_to_con_pt(double v[], double u[], double gd[3][3], double sdetg){
    // This probably occurs after an interpolation, so let's go ahead and floor
    // the density and pressure.
    // FIXME: Take a look at the floor.
    if(v[PVAR::V_RHO] < FLUID_VACUUM_RESET){
      v[PVAR::V_RHO] = FLUID_VACUUM_RESET;
      v[PVAR::V_VX] = 0.0;
      v[PVAR::V_VY] = 0.0;
      v[PVAR::V_VZ] = 0.0;
    }
    v[PVAR::V_P] = fmax(FLUID_VACUUM_TAU_RESET, v[PVAR::V_P]);

    double rho = v[PVAR::V_RHO];
    double P   = v[PVAR::V_P];

    /*if(rho < 0.0){
      printf("prim_to_con_pt: rho < 0. rho=%g\n",rho);
      rho=FLUID_VACUUM_RESET;
      return Error::FLUID_NEGATIVE_RHO;
    }
    if(P < 0.0){
      printf("prim_to_con_pt: P < 0. P=%g\n",P);
      P = FLUID_VACUUM_TAU_RESET;
      return Error::FLUID_NEGATIVE_P;

    }*/

    double vu[3], vd[3];
    /*double v4[3]={v[PVAR::V_U1],v[PVAR::V_U2],v[PVAR::V_U3]};
    
    double v3[3];
    double alpha=1.0;
    double Beta[3]={0.0,0.0,0.0};
    Vvec4_to_Vvec3(v4,v3,gd,alpha,Beta);

    vu[0] = v3[0];//v[PVAR::V_VX];
    vu[1] = v3[1];//v[PVAR::V_VY];
    vu[2] = v3[2];//v[PVAR::V_VZ];*/

    vu[0] = v[PVAR::V_VX];
    vu[1] = v[PVAR::V_VY];
    vu[2] = v[PVAR::V_VZ];


    double vsq = square_vector(vu, gd);
    vector_lower_index(vd, vu, gd);

    double Wsq = 1.0/(1.0-vsq); // Squared Lorentz factor
    //double h = rho + FLUID_GAMMA*P / (FLUID_GAMMA - 1.0); // total enthalpy
    double h = rho*(1.0 + FLUID_GAMMA/(FLUID_GAMMA - 1.0)*P/rho);
    double hWsq = h*Wsq;

    if(Wsq < 1.0){
      printf("prim_to_con_pt: Wsq < 1. Wsq=%f vsq=%f. v=(%f, %f, %f)\n",Wsq,vsq, vu[0], vu[1], vu[2]);
      printf("prim_to_con_pt: metric: %f, %f, %f, %f, %f, %f\n",gd[0][0],gd[0][1],gd[0][2],gd[1][1],gd[1][2],gd[2][2]);
      print_stacktrace();
      return Error::FLUID_W_SQUARE_LT_1;
    }
    else{
      u[VAR::U_D  ] = sdetg*(rho *sqrt(Wsq));
      u[VAR::U_SX ] = sdetg*(hWsq * vd[0]);
      u[VAR::U_SY ] = sdetg*(hWsq * vd[1]);
      u[VAR::U_SZ ] = sdetg*(hWsq * vd[2]);
      u[VAR::U_TAU] = sdetg*(hWsq - P) - u[VAR::U_D];
    }

    return Error::FLUID_SOLVER_SUCCESS;

  }


  Error con_to_prim(const ot::Mesh * pMesh, double**u,double **v, std::set<unsigned int>& badBounds)  
  {
    
    double gd[3][3], gu[3][3], detg,sdetg;
    double pos[3];

    unsigned int ownerID, ii_x, jj_y, kk_z;
    const ot::TreeNode* allElements = &(*(pMesh->getAllElements().begin()));
    const unsigned int* cg_to_dg = &(*(pMesh->getCG2DGMap().begin()));
    const unsigned int elementOrder = pMesh->getElementOrder();
    
    unsigned int lookUp;
    Error status;
    unsigned int failCount=0;
    unsigned int sz;
    std::set<unsigned int> badPoints;
    for(unsigned int node=pMesh->getNodeLocalBegin();node<pMesh->getNodeLocalEnd();node++)
    {

      lookUp=cg_to_dg[node];
      pMesh->dg2eijk(lookUp,ownerID,ii_x,jj_y,kk_z);

      sz=(1u<<(m_uiMaxDepth-allElements[ownerID].getLevel()));

      pos[0]=allElements[ownerID].minX() + ii_x*(sz/((double)elementOrder));
      pos[1]=allElements[ownerID].minY() + jj_y*(sz/((double)elementOrder));
      pos[2]=allElements[ownerID].minZ() + kk_z*(sz/((double)elementOrder));

      metric_vars(gu, &detg, gd, pos);
      sdetg = sqrt(detg);

      double u_pp[FLUID_NUM_EVOL_VARS]={u[VAR::U_D][node], u[VAR::U_SX][node], u[VAR::U_SY][node], u[VAR::U_SZ][node], u[VAR::U_TAU][node]};
      double v_pp[FLUID_NUM_PRIM_VARS];
      
      for(unsigned int var=0; var< FLUID_NUM_PRIM_VARS;var++)
        v_pp[var]=v[var][node];

      
      /*if(u_pp[VAR::U_TAU] < 0){
        std::cout << "Tau is negative!\n";
        std::cout << "  tau = " << u_pp[VAR::U_TAU] << "\n";
        std::cout << "  point = (" << ii_x << ", " << jj_y << ", " << kk_z << ")\n";
        std::cout << "  pos = (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")\n";
      }*/
      status=con_to_prim_pt(u_pp,v_pp,gd,gu,sdetg,pos);

      // Sometimes unphysical variables are thrown into the primitive solver.
      // These variables are usually rescaled to be physical, but this rescaling
      // is a bit arbitrary and can cause other issues. So, under those circumstances
      // we add points that have been rescaled to a set. These points are then
      // interpolated using linear interpolation with up to four nearest neighbors.
      // If those points are unavailable, the code leaves the points alone. Here's
      // the basic overview for implementation:
      // 1. Flag unphysical points and add them to the bad point set.
      // 2. For each point, do the following:
      //    - If it's on a processor boundary, leave it alone.
      //    - If it's an interior point...
      //      - Check both x neighbors. If neither is in the bad point set, add them
      //        to a vector of averaging points.
      //      - Check both y neighbors. If neither is in the bad point set, add them
      //        to a vector of averaging points.
      //      - If the averaging vector is greater than zero, use compensated
      //        summation to average each the values of each point for each primitve
      //        variable. Add an option for rho to be averaged with log(rho), and 
      //        remember to average the four-velocity instead of the three-velocity.
      // 3. Bad points on processor boundaries should be averaged after the variables
      //    are unzipped so that they have access to an updated padding region.

      /*if(status!=Error::FLUID_SOLVER_SUCCESS)
      {
        std::cout<<"rank: "<<(pMesh->getMPIRank())<<"node : "<<node<<" con_to prim failed : "<<status<<std::endl;
        for(unsigned int var=0;var<FLUID_NUM_EVOL_VARS;var++)
          std::cout<<"cons var: "<<FLUID_EVOL_VAR_NAMES[var]<<" value: "<<u_pp[var]<<std::endl;
      }*/
      

      for(unsigned int var=0;var<FLUID_NUM_PRIM_VARS;var++)
        v[var][node]=v_pp[var];
      // U can change inside con_to_prim because of the floor.
      for(unsigned int var=0;var<FLUID_NUM_EVOL_VARS;var++){
        u[var][node]=u_pp[var];
      }
      
      /*if(status != Error::FLUID_SOLVER_SUCCESS)
        failCount++;*/

    }

    if(failCount==0)
      return Error::FLUID_SOLVER_SUCCESS;
    else
    {
      // need to interpolate for the failed points. 
      //@todo: Milinda write this part. 
      std::cout<<"[Error]: conserved to primitive variable solver failed. "<<std::endl;
      print_stacktrace();  
      return Error::FLUID_CONS_TO_PRIM_FAILED;
    }
      
    

  }



Error con_to_prim_pt(double *dupt, double *vpt,double gd[3][3], double gu[3][3], double sdetg,double pos[]){

  int num = FLUID_NUM_EVOL_VARS;
  double upt[num];
  // FIXME: Something to experiment with -- this
  // factor of 1.0e-15 actually removes energy from
  // the system, albeit at a very small rate. Over
  // long periods of time, this seems to have an
  // effect, particularly on simulations that
  // develop instabilities, even for flat Cartesian
  // spaces. For alternative coordinate systems and
  // non-flat spaces, this will certainly have a
  // larger effect.
  double isdetg = 1.0/(sdetg + 1.0e-15);
  /*double isdetg = 1.0/sdetg;
  if(sdetg == 0){
    isdetg = 1.0e-16;
  }*/

  double alpha = 1.0;
  double Beta[3] = {0.0, 0.0, 0.0};
  bool rescaled = false;

  // Undensitize the variables
  upt[VAR::U_D  ] = isdetg * dupt[VAR::U_D  ];
  upt[VAR::U_SX ] = isdetg * dupt[VAR::U_SX ];
  upt[VAR::U_SY ] = isdetg * dupt[VAR::U_SY ];
  upt[VAR::U_SZ ] = isdetg * dupt[VAR::U_SZ ];
  upt[VAR::U_TAU] = isdetg * dupt[VAR::U_TAU];
  
  // Apply floor -- floor should be applied to the undensitized vars
  // FIXME  The flooring given below depends on D and tau being below vacuum
  // values independently.  This should probably be imposed in a more
  // intelligent way that considers more of the possibilities.

  // FIXME: Eric -- look at consistency of the floor.
  if(upt[VAR::U_D] < FLUID_VACUUM){
    upt[VAR::U_D  ] = FLUID_VACUUM_RESET;
    upt[VAR::U_SX ] = 0.0;
    upt[VAR::U_SY ] = 0.0;
    upt[VAR::U_SZ ] = 0.0;
    upt[VAR::U_TAU] = fmax(upt[U_TAU], FLUID_VACUUM_TAU_RESET);
    rescaled = true;
  }
  else{
    if (upt[VAR::U_TAU] < FLUID_VACUUM_TAU){
      upt[VAR::U_TAU] = FLUID_VACUUM_TAU_RESET;
      rescaled = true;
    }
    double Sd[3];
    Sd[0] = upt[VAR::U_SX];
    Sd[1] = upt[VAR::U_SY];
    Sd[2] = upt[VAR::U_SZ];
    // FIXME: This scaling is arbitrary. Consider other
    // possibilities.
    double Ssq = square_form(Sd, gu);
    double Ssq_max = (2.0*upt[VAR::U_D] + upt[VAR::U_TAU])*upt[VAR::U_TAU];
    if(Ssq_max < Ssq){
      double t_min = VELOCITY_FACTOR * sqrt(Ssq_max/Ssq);
      upt[VAR::U_SX] *= t_min;
      upt[VAR::U_SY] *= t_min;
      upt[VAR::U_SZ] *= t_min;
      rescaled = true;
    }
  }

  // Call primary EOS solver
  double fpar[NFPAR];
  c2p_work_vars(fpar, upt, vpt, gd, gu, FLUID_GAMMA, pos);

  int eos_solver = 0;
  int rc = mhd_solver_4(upt, vpt, gd, gu, fpar);
  
  if( rc < 0){
    eos_solver = 1;
    rc = mhd_solver_isentropic(upt,vpt,pos,gd,gu,fpar);
  }
  
  // computes the 4 velocity vector
  double v3[]={vpt[PVAR::V_VX],vpt[PVAR::V_VY],vpt[PVAR::V_VZ]};
  double v4[3];

  Vvec3_to_Vvec4(v3,v4,gd,alpha,Beta);

  vpt[PVAR::V_U1]=v4[0];
  vpt[PVAR::V_U2]=v4[1];
  vpt[PVAR::V_U3]=v4[2];

  Vvec4_to_Vvec3(v4,v3,gd,alpha,Beta);
  vpt[PVAR::V_VX]=v3[0];
  vpt[PVAR::V_VY]=v3[1];
  vpt[PVAR::V_VZ]=v3[2];
  
  // Densitize the conserved variables
  dupt[VAR::U_D  ] = sdetg*upt[VAR::U_D  ];
  dupt[VAR::U_SX ] = sdetg*upt[VAR::U_SX ];
  dupt[VAR::U_SY ] = sdetg*upt[VAR::U_SY ];
  dupt[VAR::U_SZ ] = sdetg*upt[VAR::U_SZ ];
  dupt[VAR::U_TAU] = sdetg*upt[VAR::U_TAU];

  if(vpt[PVAR::V_P] < FLUID_VACUUM_TAU){
    vpt[PVAR::V_P] = FLUID_VACUUM_TAU_RESET;
    prim_to_con_pt(dupt, vpt, gd, sdetg);
    rescaled = true;
  }
  
  // left boundary.  
  if(FLUID_COORDS == 1 &&  (pos[0]==0.0)){
    // Reflection symmetries about the azimuthal axis
    const double coef[8] = {1.0, -1.0, 1.0, 1.0, 1.0,-1.0,1.0,1.0};

    for(unsigned int var=0;var<FLUID_NUM_PRIM_VARS;var++)
      vpt[var]=coef[var]*vpt[var];
    
    for(unsigned int var=0;var<FLUID_NUM_EVOL_VARS;var++)
      dupt[var]=coef[var]*dupt[var];
  
  }
    
  if(rc<0)
    return Error::FLUID_CONS_TO_PRIM_FAILED;
    
  if(rescaled){
    return Error::FLUID_CONS_TO_PRIM_RESCALED;
  }

  return Error::FLUID_SOLVER_SUCCESS;
    
  
}

Error Vvec3_to_Vvec4(double* v3, double* v4, double gd[3][3],double alpha, double Beta[3])
{
    const double vSq=square_vector(v3,gd);
    if(vSq>1.0)
      return Error::FLUID_V_SQUARE_GT_1;

    const double W = sqrt(1.0/(1-vSq));

    v4[0]=W*(v3[0]-Beta[0]/alpha);
    v4[1]=W*(v3[1]-Beta[1]/alpha);
    v4[2]=W*(v3[2]-Beta[2]/alpha);

    return Error::FLUID_SOLVER_SUCCESS;


}


Error Vvec4_to_Vvec3(double* v4, double* v3, double gd[3][3], double alpha, double Beta[3])
{

    double Betad[3];
    vector_lower_index(Betad, Beta, gd);

    const double bu = Betad[0]*v4[0] + Betad[1]*v4[1] + Betad[2]*v4[2];
    const double usq = square_vector(v4, gd);
    const double bsq = square_vector(Beta, gd);

    const double inv_alpha=1.0/alpha;
    double inv_W = (alpha - bsq*inv_alpha) / (bu + sqrt(bu*bu + (alpha*alpha - bsq)*(1.0 + usq)));

    v3[0]=v4[0]*inv_W + Beta[0]*inv_alpha;
    v3[1]=v4[1]*inv_W + Beta[1]*inv_alpha;
    v3[2]=v4[2]*inv_W + Beta[2]*inv_alpha;

   /* debugging check */
    const double vSq=square_vector(v3,gd);
    if(vSq>1.0) {
      printf("Vsq: %f v=(%f,%f,%f)\n",vSq,v3[0],v3[1],v3[2]);
      print_stacktrace();
      return Error::FLUID_V_SQUARE_GT_1;
    }
    else {
      return Error::FLUID_SOLVER_SUCCESS;
    }

}


// c2p_work_vars {{{
void c2p_work_vars(double *fpar, double *u, double *v, double gd[3][3],
                   double gu[3][3], double gamma, double *pos){
  double kappa;
  double Sd[3], vu[3];
  double vsq, Ssq;
  double P, rho_0;

  // Primitive variables before the inversion. These are used below for some
  // auxiliary variables.
  rho_0 = v[PVAR::V_RHO];
  vu[0] = v[PVAR::V_VX ];
  vu[1] = v[PVAR::V_VY ];
  vu[2] = v[PVAR::V_VZ ];
  P     = v[PVAR::V_P  ];

  // Some of the conserved variables before the inversion. These are used 
  // below for some auxiliary variables.
  Sd[0] = u[VAR::U_SX];
  Sd[1] = u[VAR::U_SY];
  Sd[2] = u[VAR::U_SZ];

  Ssq = square_form(Sd, gu);
  vsq = square_vector(vu, gd);
  if (vsq >= 1.0){
    printf("  c2p_work_vars:  The velocity is superluminal.  To correct \n");
    printf("  c2p_work_vars:  this, we simply set it to 0.99999999 \n");
    printf("  c2p_work_vars:  which seems sort of stupid. \n");
    printf("  c2p_work_vars:       vsq = %25.20e \n", vsq);
    //FIXME It would be well to come up with something more intelligent
    vsq = 0.99999999;
  }

  // Define some temporary and auxiliary variables.
  fpar[FP_SSQ]   = Ssq;
  fpar[FP_D]     = u[VAR::U_D];
  fpar[FP_TAU]   = u[VAR::U_TAU];
  fpar[FP_GAMMA] = gamma;
  fpar[FP_VSQ]   = vsq;
  fpar[FP_COORD_X] = pos[0];
  fpar[FP_COORD_Y] = pos[1];
  fpar[FP_COORD_Z] = pos[2];
  kappa            = P/pow(rho_0,gamma);
  fpar[FP_KAPPA]   = (kappa < 1000.0) ? kappa : 1000.0;

  // TRACE is reserved for internal use in the specific solvers */
  fpar[FP_TRACE]   = 0.0;
  fpar[FP_C2PWARN] = FLUID_CONTOPRIMWARN;
  fpar[FP_VACUUM]  = FLUID_VACUUM;
}
// }}}

// mhd_solver_4 {{{
// This is the most robust and accurate primitive solver without a magnetic field.
int mhd_solver_4(double *u, double *v, double gd[3][3], double gu[3][3], double *fpar){
  const int ltrace = 0;

  double gamma = FLUID_GAMMA;

  double Ssq = fpar[FP_SSQ];
  double D   = u[VAR::U_D];
  double tau = u[VAR::U_TAU];
  double Sd[3];

  Sd[0] = u[VAR::U_SX];
  Sd[1] = u[VAR::U_SY];
  Sd[2] = u[VAR::U_SZ];
  
  int j;
  double qf_1, qf_2, qg_2, qf_max;
  double ql, qh, f_at_qh, f_at_qh_v2;
  double f_at_ql = 0.0;
  double rescale_factor1;
  double a = tau + D;
  double t1 = 0.5 * gamma;
  double beta_sq = Ssq/(a*a);
  double t2 = t1*t1 - (gamma-1.0)*beta_sq;
  // t2 shows up under a square root, so we check to see if it is < 0.
  if(t2 < 0.0){
    printf("mhd_solver_4:  problem with upper bound. t2 < 0. t2 = %g\n",t2);
    exit(2);
  }
  
  double sqrt_t2 = sqrt(t2);
  qf_1 = t1 - sqrt_t2;  // the smaller root of f
  qf_2 = t1 + sqrt_t2;   // the larger root of f
  qf_max = t1;          // the location of the max of f
  //qg_2 = sqrt(beta_sq); // the larger root of g
  qg_2 = 1.0;

  // If the velocity vector (or S^i) is the zero vector, the conserved to
  // primitive inversion is analytic.  This translates to beta_sq being zero.
  // We do this case separately and return.
  if (beta_sq < 1.e-16){
    v[PVAR::V_RHO] = D;
    v[PVAR::V_P  ] = (gamma - 1.0) * tau;
    v[PVAR::V_VX ] = 0.0;
    v[PVAR::V_VY ] = 0.0;
    v[PVAR::V_VZ ] = 0.0;

    return 1;
  }

  if(fabs((qg_2 - qf_1) / (qg_2 + qf_1 + 1.0e-15)) < 1.e-6){
    if (ltrace){
      printf("  mhd_solver_4:  Roots are colliding near q=0. \n");
      printf("  mhd_solver_4:  Adding a 1.e-5 piece to the larger \n");
      printf("  mhd_solver_4:  hoping this resolves the problem. \n");
    }
    qg_2 += 1.e-5;
  }

  // Set the lower bound to qg_2 unless qg_2 is smaller than the location of
  // the maximum of f.
  ql = fmax(qg_2, qf_max);

  // This temporarily initializes the upper bound to a possible value.
  // The following code checks this to make sure i works.
  qh = qf_2;

  /* ... documentation on the cases that need to be considered ... {{{
  We consider various possibilities for the relative location of qg_2 and
  qf_2 as they should be the lower and upper bounds, respectively, for our
  root, q, of F(q).

  Case 1 is that qg_2 < qf_2 and a physical inequality is obeyed.  All is well
  in this case and the Newton solve should proceed without issue.

  Case 2 considers the violation of another physical inequality and adjusts
  (rescales) the components of S in order to enforce the inequality.

  Case 3 considers the numerical "vanishing" of D/a relative to beta_sq which
  can also numerically violate physical inequalities.

  Case 4 considers the (seemingly remote but nonetheless realized) possibility
  that qg_2 and qf_2 are exactly the same (to machine precision).

  No claim is made that these cases are exhaustive.  They were simply ones
  that arose in debugging the code and trying to make this primitive solver as
  robust as possible.  There is a sort of "ad hoc-ness" to some of these.
  Could they be improved?  Probably.  But for now they seem to work.
  }}} */
  // CASE 1
  if (qg_2 < qf_2 && beta_sq < 1.0){
    qh = qf_2;
  }
  // CASE 2
  else if (beta_sq + (D/a)*(D/a) - 1.0 >= 0.0){
    if(ltrace){
      // {{{
      printf("  mhd_solver_4:  Our bounds are colliding near q=1.  \n");
      printf("  mhd_solver_4:  In particular, we have an inequality \n");
      printf("  mhd_solver_4:  violation.  In particular, Ssq/a^2 + D^2/a^2\n");
      printf("  mhd_solver_4:  should be less than 1.  Instead it is \n") ;
      printf("  mhd_solver_4:   Ssq/a^2+D^2/a^2=%35.30e\n",beta_sq+(D/a)*(D/a));
      printf("  mhd_solver_4:  where:     \n" );
      printf("  mhd_solver_4:      Ssq = %25.20e \n", Ssq );
      printf("  mhd_solver_4:      a^2 = %25.20e \n", a*a );
      printf("  mhd_solver_4:      D^2 = %25.20e \n", D*D );
      printf("  mhd_solver_4:      tau = %25.20e \n", tau );
      printf("  mhd_solver_4:       D  = %25.20e \n", D );
      printf("     \n");
      printf("  mhd_solver_4:  We will rescale S_i and Ssq in hopes this \n");
      printf("  mhd_solver_4:  resolves the problem. \n");
      // }}}
    }

    
    /*
    We rescale the components of Sd (and hence Ssq) in order that
    the inequality is enforced.

    FIXME:  Note that this scaling factor is a bit arbitrary.  It would
    FIXME:  be well to experiment with the value added in.
    */
    //double  eps0 = 1.e-14 ;
    //rescale_factor1 = 1.0 / (sqrt(beta_sq + (D/a)*(D/a)) + eps0) ;
    rescale_factor1 = 1.0 / (sqrt(beta_sq + (D/a)*(D/a)) + 1.e-14) ;

    // Now multiply Sd by the scaling factor.
    for ( j=0; j<3; j++ )
    { Sd[j] *= rescale_factor1 ; }

    // Recalculate Ssq.
    Ssq = square_form(Sd, gu);

    // Now we have to recalculate the roots of f and g with this new Ssq.
    beta_sq = Ssq/(a*a);
    t2 = t1*t1 - (gamma-1.0)*beta_sq ;
    if (t2 < 0.0) { printf("mhd_solver_4: t2 < 0, i.e. t2=%25.20e\n",t2); }
    qf_2 = t1 + sqrt(t2) ; // new root of f; upper bound on q
    qg_2 = sqrt(beta_sq) ; // new root of g; lower bound on q

    // Check that these new roots are indeed good upper and lower bounds.
    if ( qg_2 < qf_2 )
    {
      ql = qg_2 ;
      qh = qf_2 ;
    }
    else
    {
      if (ltrace) {
//{{{
        printf("  mhd_solver_4:  Even after rescaling, we still have a \n");
        printf("  mhd_solver_4:  problem with the relative values of the \n");
        printf("  mhd_solver_4:  bounds qg_2 and qf_2.  In particular, we n");
        printf("  mhd_solver_4:  should have that  qg_2 < qf_2. Instead, we\n");
        printf("  mhd_solver_4:  have:  \n");
        printf("  mhd_solver_4:      qg_2 = %25.20e \n", qg_2 );
        printf("  mhd_solver_4:      qf_2 = %25.20e \n", qf_2 );
//}}}
      }
    }

    // Because we made a rescaling of the components of S and of Ssq, we need
    // to save these.
    u[VAR::U_SX] = Sd[0];
    u[VAR::U_SY] = Sd[1];
    u[VAR::U_SZ] = Sd[2];
    fpar[FP_SSQ] = Ssq;

    //{{{ Having rescaled, check again that the inequality that failed before is now satisfied.
    if ( beta_sq + (D/a)*(D/a) >= 1.0 )
    {
      if (ltrace) {
        printf("  mhd_solver_4:  Even after rescaling S_i the physical \n");
        printf("  mhd_solver_4:  inequality that should be satisfied is not\n");
        printf("  mhd_solver_4:  In particular, \n");
        printf("  mhd_solver_4:   (Ssq+D^2)/a^2=%35.30e\n",beta_sq+(D/a)*(D/a));
      }
    }
    else
    {
      if (ltrace) {
        printf("  mhd_solver_4:  After what appears to be a successful \n");
        printf("  mhd_solver_4:  rescaling, we have  \n");
        printf("  mhd_solver_4:   (Ssq+D^2)/a^2=%35.30e\n",beta_sq+(D/a)*(D/a));
      }
    }
  // }}}
  }
  
  // CASE 3
  else if (fabs((beta_sq + (D/a)*(D/a)) - beta_sq) <= 6.0e-16){
    // {{{
    printf("  mhd_solver_4:  SURPRISE! I deleted the treatment of this \n");
    printf("  mhd_solver_4:  possibility as I didn't think it ever arose. \n");
    printf("  mhd_solver_4:  I guess it does after all. \n");
/*
    printf("  mhd_solver_4:  D/a is so small as to be lost in roundoff.\n");
    printf("  mhd_solver_4:  We will adjust Ssq in hopes this \n");
    printf("  mhd_solver_4:  resolves the problem. \n");

    //Instead of a true (multiplicative) rescaling, in this case, we have
    //found that simply subtracting a small amount from each of the
    //components of S can bring Ssq back into the physical regime.
    //  FIXME:  Again, the value 1.e-15 seems a bit arbitrary.
    //  FIXME:  It might be well to play with this some.
    //  FIXME:  Actually, this could be a real problem.  If Sd is already negative, we could be making the magnitude,
 Ssq, in fact larger.  this is definitely a potential problem and should be addressed ...
    for ( j=0; j<3; j++ ) { Sd[j] -= 1.e-15 * a ; }

    // Recalculate Ssq.
    Ssq = square_form(Sd, gu);

    // Now we have to recalculate the roots of f and g with this new Ssq.
    beta_sq = Ssq/(a*a);
    t2 = 1.0 - (gamma-1.0)*beta_sq/(t1*t1) ;
    if (t2 < 0.0) { printf("mhd_solver_4: t2 < 0, i.e. t2=%25.20e\n",t2); }
    qf_2 = t1*( 1.0 + sqrt(t2)) ; // new root of f; upper bound on q
    qg_2 = sqrt(beta_sq) ;        // new root of g; lower bound on q

    // Check that these new roots are indeed good upper and lower bounds.
    if ( qg_2 < qf_2 )
    {
      ql = qg_2 ;
      qh = qf_2 ;
    }
    else
    {
      printf("  mhd_solver_4:  Even after adjusting (via subtracting) \n") ;
      printf("  mhd_solver_4:  the components of S_i, we still have a \n");
      printf("  mhd_solver_4:  problem with the relative values of the \n");
      printf("  mhd_solver_4:  bounds  qg_2  and  qf_2.  In particular,\n");
      printf("  mhd_solver_4:  we should have  qg_2 < qf_2.  Instead, we \n");
      printf("  mhd_solver_4:  have \n");
      printf("  mhd_solver_4:      qg_2 = %25.20e \n", qg_2 );
      printf("  mhd_solver_4:      qf_2 = %25.20e \n", qf_2 );
    }

    // Because we made a change to the S components and Ssq, we need to save
    // them.
    u[U_SX] = Sd[0] ;
    u[U_SY] = Sd[1] ;
    u[U_SZ] = Sd[2] ;
    fpar[FP_SSQ] = Ssq;

    // Having modified the conserved variables to force the inequality
    // to be satisfied, let's check that that actually happened.
    if ( fabs( (beta_sq + (D/a)*(D/a)) - beta_sq ) <= 6.0e-16 )
    {
      printf("  mhd_solver_4:  On adjusting componets of S, we still \n");
      printf("  mhd_solver_4:  have a problem with D being too small: \n");
      printf("  mhd_solver_4:      fabs( (beta_sq + (D/a)*(D/a)) - beta_sq ) = %25.20e \n", fabs((beta_sq+(D/a)*(D/a)
) - beta_sq) );
      printf("  mhd_solver_4:      Ssq/a^2 + D^2/a^2 = %35.30e \n",beta_sq+(D/a)*(D/a));
    }*/
    // }}}
  }
  // CASE 4
  else if (qg_2 == qf_2){
    qg_2 -= 1.e-16;
    qf_2 += 1.e-16;

    ql = qg_2;
    qh = qf_2;
  }
  else{
    // {{{
    // This else is effectively a garbage dump for all the other possibilities
    // not delineated in the if/else above.  It is a sort of "none of the
    // above". As a result, if we do get here, we essentially don't know why...
    printf("  mhd_solver_4:   something unexpected (1) \n");
    printf("        qg_2 = %25.20e \n", qg_2 );
    printf("        qf_2 = %25.20e \n", qf_2 );
    printf("     beta_sq = %25.20e \n", beta_sq );
    printf("        D    = %25.20e \n", D );
    printf("        tau  = %25.20e \n", tau );
    printf("        Ssq  = %25.20e \n", Ssq );
    printf("    beta_sq + (D/a)^2 = %25.20e \n", beta_sq+(D/a)*(D/a) );
    printf("   dump fpar  \n");
    printf("      fpar[FP_SSQ]     = %25.20e \n", fpar[FP_SSQ] );
    printf("      fpar[FP_D  ]     = %25.20e \n", fpar[FP_D  ] );
    printf("      fpar[FP_TAU]     = %25.20e \n", fpar[FP_TAU] );
    printf("      fpar[FP_GAMMA]   = %25.20e \n", fpar[FP_GAMMA] );
    printf("      fpar[FP_VSQ]     = %25.20e \n", fpar[FP_VSQ] );
    printf("      fpar[FP_COORD_X] = %25.20e \n", fpar[FP_COORD_X] );
    printf("      fpar[FP_COORD_Y] = %25.20e \n", fpar[FP_COORD_Y] );
    printf("      fpar[FP_COORD_Z] = %25.20e \n", fpar[FP_COORD_Z] );
    printf("      fpar[FP_INDEX_X] = %25.20e \n", fpar[FP_INDEX_X] );
    // }}}
  }

  // At this point, we should have the two values, namely qg_2 and qf_2,
  // that shoud form a legitimate bracket for finding the root of F(q).
  // Let's check.
  //
  // We want to check that at ql the function F(q) := f(q)-(const)*sqrt(g(q))
  // is positive.  First just calculate F(ql).
  int rc=0;
  rc = func4_p( &f_at_ql, ql, fpar );
  if (rc < 0)
  {
    printf(" mhd_solver_4:  Something's amiss with ql \n");
    printf(" mhd_solver_4:  We couldn't calculate f(ql). \n");
  }

  // Now check that at qh the function F(q) := f(q)-(const)*sqrt(g(q))
  // is negative. For now, just calculate F(qh).
  rc = 0;
  rc = func4_p( &f_at_qh, qh, fpar );
  if (rc < 0)
  {
    printf(" mhd_solver_4:  Something's amiss with qh \n");
    printf(" mhd_solver_4:  We couldn't calculate f(qh). \n");
  }

  // If, indeed, F(ql) > 0 and F(qh) < 0, then we should have a 
  // bracket for the root.  We will check for the failure of this and do
  // what we can.
  // FIXME Actually, if there is a problem, we don't do anything here.
  if ( f_at_ql < 0.0 || f_at_qh > 0.0 )
  {
    if (ltrace) {
//{{{
      printf("  mhd_solver_4:  There is a problem with the bracketing of \n");
      printf("  mhd_solver_4:  the root of F(q).    \n");
      printf("  mhd_solver_4:         ql = %25.20e \n",ql);
      printf("  mhd_solver_4:    f_at_ql = %25.20e should be > 0.0\n",f_at_ql);
      printf("  mhd_solver_4:         qh = %25.20e \n",qh);
      printf("  mhd_solver_4:    f_at_qh = %25.20e should be < 0.0\n",f_at_qh);
      printf("  mhd_solver_4:  We are sort of at a loss at this point.  We \n");
      printf("  mhd_solver_4: will go on, but without much hope of success.\n");
//}}}
    }
  }

  if ( f_at_ql > 0.0  &&  f_at_qh > 0.0 )
  {
    if (ltrace) {
      printf("  mhd_solver_4:  At the lower and upper bounds, the \n");
      printf("  mhd_solver_4:  function is positive so we fail to \n");
      printf("  mhd_solver_4:  get a bracket at all.  Needless to \n");
      printf("  mhd_solver_4:  say, this is bad and may mean no   \n");
      printf("  mhd_solver_4:  solution exists.                   \n");
      printf("  mhd_solver_4:  Optimistically, it may meant that  \n");
      printf("  mhd_solver_4:  the upper bound is not large enough. \n");
      printf("  mhd_solver_4:  We will try to shift the current value \n");
      printf("  mhd_solver_4:  of qh (which is %25.20e)up by 6.e-16,i.e.\n",qh);
    }

    qh += 6.e-16 ;

    if (ltrace) {
      printf("  mhd_solver_4:      qh = %25.20e \n",qh );
      printf("  mhd_solver_4:  Checking to see if that worked ... \n");
    }

    rc = func4_p( &f_at_qh_v2, qh, fpar ) ;

    if (ltrace) {
      printf("  mhd_solver_4:  ... well, among other things, the new \n");
      printf("  mhd_solver_4:  value of f (which should now be negative) \n");
      printf("  mhd_solver_4:  at the new value of qh is: \n");
      printf("  mhd_solver_4:    old      qh = %25.20e \n",qh-6.e-16);
      printf("  mhd_solver_4:    old f_at_qh = %25.20e \n",f_at_qh);
      printf("  mhd_solver_4:    new      qh = %25.20e \n",qh);
      printf("  mhd_solver_4:    new f_at_qh = %25.20e (should be neg) \n",
                f_at_qh_v2);
    }
  }

  /*
  At this point, we should have a legitimate bracket for our root, namely
  [ql, qh] := [qg_2, qf_2].  We should now be able to proceed with the
  Newton solve.

  As our initial guess for the root, take the midpoint of the bracket.
  */
  //double qi = 0.5*(ql + qh);
  // Actually, we'll use false position here.
  double qi = (ql*f_at_qh - qh*f_at_ql)/(f_at_qh - f_at_ql);
  double q = qi;

  //  Incorporate a little redundancy in the event things go bad.
  double ql0 = ql;
  double qh0 = qh;
  double q0 = q;

  //  Set trace parameter for rtsafe to 0 (no tracing), set the tolerance for
  //  the Newton solve, and call rtsafe to get the root, q.
  int   xtrace = 0;
  double  xtol = 1.0e-14;
  //int xrc = froot::rtsafe(func4_p_d, &q, ql, qh, xtol, fpar, xtrace);
  double fq = 1.0;
  double dfq = 1.0;
  double ddfq = 1.0;
  unsigned int count = 0;
  int held = 0;
  int maxIterations = 50;
  while(fabs(fq) > xtol && count < maxIterations){
    double oldq = q;
    // Newton's method
    //func4_p_d(&fq, &dfq, q, fpar);
    //q = q - fq/dfq;
    //Halley's method
    func4_p_dd(&fq, &dfq, &ddfq, q, fpar);
    q = q - 2.0*fq*dfq/(2.0*dfq*dfq - fq*ddfq);
    // Check for convergence by ensuring that the new q is
    // still inside the bounds. If not, apply the Illinois
    // variant of the false position method.
    if(q > qh || q < ql){
      // Calculate the new root.
      if(fq*f_at_qh > 0){
        double m;
        if(held == -1){
          m = 0.5;
        }
        else{
          m = 1.0;
        }
        q = (f_at_qh*ql - m*f_at_ql*qh)/(f_at_qh - m*f_at_ql);
      }
      else{
        double m;
        if(held == 1){
          m = 0.5;
        }
        else{
          m = 1.0;
        }
        q = (m*f_at_qh*ql - f_at_ql*qh)/(m*f_at_qh - f_at_ql);
      }
    }
    if(fq < 0){
      held = -1;
      qh = oldq;
      f_at_qh = fq;
    }
    else{
      held = 1;
      ql = oldq;
      f_at_ql = fq;
    }
    count++;
  }

  // In the event of failure, a postmortem ...
  //if (xrc < 0)
  if(fabs(fq > xtol))
  {
    printf("  mhd_solver_4:  The Newton solve (rtsafe) failed.\n");
    printf("  mhd_solver_4:  Current values include  \n");
    printf("    ql = %25.20e \n",ql);
    printf("    qh = %25.20e \n",qh);
    printf("    qi = %25.20e \n",qi);
    printf("    q   =%25.20e \n",q);
    printf("    ql0 = %25.20e \n",ql0);
    printf("    qh0 = %25.20e \n",qh0);
    printf("    q0   =%25.20e \n",q0);
    printf("  mhd_solver_4:  Some conserved quantities and  \n");
    printf("  mhd_solver_4:  inequalities,  \n");
    printf("    D               = %25.20e \n", D);
    printf("    tau+D           = %25.20e \n", a);
    printf("    D^2/a^2         = %25.20e \n", (D*D)/(a*a));
    printf("    Ssq/(a*a)       = %25.20e \n", Ssq/(a*a));
    printf("    (Ssq+D^2)/(a*a) = %25.20e \n", (Ssq+D*D)/(a*a));
    //printf("  mhd_solver_4:  After rtsafe failure ... redoing ... \n");

    // Retrying rtsafe but with a different trace parameter so as to get
    // some more info out of rtsafe.
    /*ql = ql0;
    qh = qh0;
    q = q0;

    xtrace = 1;
    xrc = froot::rtsafe(func4_p_d, &q, ql, qh, xtol, fpar, xtrace);*/
    exit(2);
    return -4;
  }

  //  If we get to this point, we may actually have a solution.
  // We have to rescale back by a = tau+D.
  q *= a;
  
  // Knowing the value for q (:= hW^2, i.e. related to the enthalpy), we
  // can now get the remaining primitive variables.  Start with the velocity
  // squared:  v^2.
  double vsq = Ssq / (q*q);
  // We have the "down" components of S; now get the "up" components
  double Su[3];
  form_raise_index(Su, Sd, gu);
  
  // Now calculate rho, pressure, and the "up" components of v^i.
  v[PVAR::V_RHO] = D*sqrt(1.0 - vsq);
  v[PVAR::V_P  ] = (gamma-1.0)/gamma*(q*(1.0-vsq) - v[PVAR::V_RHO]);
  //v[PVAR::V_P  ] = q - a;
  v[PVAR::V_VX ] = Su[0] / q;
  v[PVAR::V_VY ] = Su[1] / q;
  v[PVAR::V_VZ ] = Su[2] / q;
  
  if(v[PVAR::V_P] == 0){
    if(ltrace){
      printf("  mhd_solver_4:  P = 0 after solve.\n");
      printf("  mhd_solver_4:  The velocity is probably really\n");
      printf("  mhd_solver_4:  small, so we'll call it zero and \n");
      printf("  mhd_solver_4:  find P analytically.\n");
    }
    v[PVAR::V_P] = (gamma - 1.0)*tau;
  }

  // Check some primitive values
  if(v[PVAR::V_RHO] < 0.0 || v[PVAR::V_P] < 0.0 || vsq < 0.0 || vsq > 1.0 ){
    if(ltrace){
      // {{{
      printf(" mhd_solver_4:  Problematic values after the primitive. \n");
      printf(" mhd_solver_4:  solve: \n");
      printf(" mhd_solver_4:     rho = %15.10e \n", v[V_RHO] );
      printf(" mhd_solver_4:      P  = %15.10e \n", v[V_P] );
      printf(" mhd_solver_4:     v^2 = %15.10e \n", vsq );
      // }}}
    }
  }
  if(v[PVAR::V_RHO] < FLUID_VACUUM){
    if(ltrace){
      printf("mhd_solver_4:  rho is negative:  rho = %15.10e \n", v[PVAR::V_RHO]);
      printf("mhd_solver_4:  we reset it to the vacuum.  \n");
    }
    v[PVAR::V_RHO] = FLUID_VACUUM_RESET;
    v[PVAR::V_VX] = 0.0;
    v[PVAR::V_VY] = 0.0;
    v[PVAR::V_VZ] = 0.0;
    u[U_D] = FLUID_VACUUM_RESET;
    u[U_SX] = 0.0;
    u[U_SY] = 0.0;
    u[U_SZ] = 0.0;
  }
  if(v[PVAR::V_P] < FLUID_VACUUM_TAU){
    if(ltrace){
      printf("mhd_solver_4:  P is negative:  P = %15.10e \n", v[PVAR::V_P]);
      printf("mhd_solver_4:  we reset it to the vacuum.  \n");
    }
    v[PVAR::V_P] = FLUID_VACUUM_TAU_RESET;
    prim_to_con_pt(v, u, gd, 1.0);
  }
  assert(D > 0.0);
  
  return 0;
}
// }}}

// mhd_solver_isentropic {{{
int mhd_solver_isentropic(double *u, double *v, double *pos, double gd[3][3],
                          double gu[3][3], double *fpar){
  int rc;
  double y, yi, vsq, F, hW2;
  double xi;

  double Sd[3], Su[3];
  
  double gamma = FLUID_GAMMA;
  
  double Ssq = fpar[FP_SSQ];
  double kappa = fpar[FP_KAPPA];
  
  double D     = u[VAR::U_D];

  const int ltrace2 = 0;
  int fail;

  if(ltrace2){
    printf("mhd_solver_isentropic begin...\n");
  }

  Sd[0] = u[VAR::U_SX];
  Sd[1] = u[VAR::U_SY];
  Sd[2] = u[VAR::U_SZ];
  
  form_raise_index(Su, Sd, gu);

  // An initial guess or rho0
  vsq = fpar[FP_VSQ];
  yi = D*sqrt(1.0-vsq);
  xi = yi;

  rc = SolvePrimVarsIsentropic(&xi, u, fpar);
  
  fail = 0;
  if (rc < 0){
    if(FLUID_CONTOPRIMWARN == 1){
      printf("Isentropic Solver failed. Die here.\n");
    }
    fail = 1;
    return -4;
  }

  y = xi;
    
  if(y <= 0.0){
    if(FLUID_CONTOPRIMWARN == 1){
      printf("### Isentropic Solver has unphysical root, y=%20.14e\n", y);
    }
  }
  
  if (y <= 0.0) {
    if ( FLUID_CONTOPRIMWARN == 1 ) {
      printf("### y is unphysical, y= %15.12f\n",y);
      printf("### SSq                    = %15.12f\n",Ssq);
      rc = func_p_isentropic(&F, y, fpar);
      printf("###  f                     = %e\n",F);
    }
    fail = 1;
  }

  if ( fail ) {
    //if (DIE_ON_C2P_ERRORS) {
      MathDump(fpar);
      exit(-1);
    //}
    return -4;
  }
  else{
    if(y < FLUID_VACUUM){
      y = FLUID_VACUUM_RESET;
      v[PVAR::V_P] = FLUID_VACUUM_TAU_RESET;
    }
    else{
      v[PVAR::V_RHO] = y;
      v[PVAR::V_P]   = kappa*pow(y,gamma);
    }

    hW2 = (y + gamma/(gamma - 1.0)*v[PVAR::V_P])*D*D/(y*y);

    // Note:  The momensum S is stored with index down.  The velocity
    //        components are index up.
    v[PVAR::V_VX] = Su[0] / hW2;
    v[PVAR::V_VX] = Su[1] / hW2;
    v[PVAR::V_VX] = Su[2] / hW2;

    if (ltrace2) {
      printf(">>> Isentropic solution >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
      printf(">>>    rho          = %g\n",v[PVAR::V_RHO]);
      printf(">>>    vx           = %g\n",v[PVAR::V_VX]);
      printf(">>>    vy           = %g\n",v[PVAR::V_VY]);
      printf(">>>    vz           = %g\n",v[PVAR::V_VZ]);
      printf(">>>    P            = %g\n",v[PVAR::V_P]);
    }
  }

  return 0;
}
// }}}

// SolvePrimVarsIsentropic {{{
int SolvePrimVarsIsentropic(double *qa, double *upt, double *fpar){
  int rc;
  double q = *qa;
  const double tol = 1.0e-10;
  const int ltrace = 1;

  double D = fpar[FP_D];

  // Bracket the root
  int trace = 0;
  double rl = 0.1*D;
  double rh = 2.0*(fpar[FP_D] + fpar[FP_TAU]);
  if(ltrace){
    fprintf(stdout,"@@ SolvePrimVarsIsentropic: ... rl=%g, rh=%g\n",rl,rh);
  }
  rc = froot::rbrac(func_p_isentropic, &rl, &rh, fpar, trace);
  if (rc < 0){
    if (FLUID_CONTOPRIMWARN == 1){
      fprintf(stdout,"SolvePrimVarsIsentropic: Could not bracket the root for hW^2\n");
      fprintf(stdout,"      rl = %g, rh = %g\n",qa[1],qa[2]);
    }
    rl = 0.1*D;
    rh = 2.0*(fpar[FP_D] + fpar[FP_TAU]);
    trace = 1;
    rc = froot::rbrac(func_p_isentropic, &rl, &rh, fpar, trace);
    return -1;
  }

  // All set up...should now be able to get the root.
  rc = froot::rtsafe(func_p_d_isentropic, &q, rl, rh, tol, fpar, trace);
  if(rc < 0){
    if(FLUID_CONTOPRIMWARN == 1){
      fprintf(stdout,"SolvePrimVarsIsentropic:  Error in rtsafe\n");
    }
    return -1;
  }
  
  *qa = q;
  return 0;
}
// }}}

// func4_p {{{
int func4_p(double *f, double q, double *fpar){
  double Ssq   = fpar[FP_SSQ];
  double D     = fpar[FP_D];
  double tau   = fpar[FP_TAU];
  double gamma = fpar[FP_GAMMA];

  double dis;
  
  double a = tau+D;
  double betasq = Ssq/(a*a);
  dis = q*q - betasq;
  if(dis < 1.e-18){
    dis = 1.e-18;
  }

  double qoff = q - 0.5*gamma;

  *f = -qoff*qoff + gamma*gamma/4.0
     - (gamma-1.0)*betasq - (gamma-1.0)*D/a*sqrt(dis);

  return 1;
}
// }}}

// func4_p_d {{{
int func4_p_d(double *f, double *df, double q, double *fpar){
  double Ssq   = fpar[FP_SSQ];
  double D     = fpar[FP_D];
  double tau   = fpar[FP_TAU];
  double gamma = fpar[FP_GAMMA];
  
  double dis;
  
  double a = tau+D;
  dis = q*q - Ssq/(a*a);
  if ( dis < 1.e-18){
    dis = 1.e-18;
  }

  *f = - (q - 0.5*gamma)*(q - 0.5*gamma) + gamma*gamma/4.0
       - (gamma-1.0)*Ssq/(a*a) - (gamma-1.0)*D/a*sqrt(dis);

  *df = -2.0*(q - 0.5*gamma) - (gamma-1.0)*D/a*q/sqrt(dis);

  return 1;
}
// }}}

// func4_p_dd {{{
int func4_p_dd(double *f, double *df, double *ddf, double q, double *fpar){
  double Ssq   = fpar[FP_SSQ];
  double D     = fpar[FP_D];
  double tau   = fpar[FP_TAU];
  double gamma = fpar[FP_GAMMA];
  
  double dis;
  
  double a = tau+D;
  double betasq = Ssq/(a*a);
  dis = q*q - betasq;
  if ( dis < 1.e-18){
    dis = 1.e-18;
  }
  double sdis = sqrt(dis);
  double qoff = q - 0.5*gamma;
  double coeff = (gamma-1.0)*D/a;

  *f = -qoff*qoff + gamma*gamma/4.0 - (gamma-1.0)*betasq - coeff*sdis;

  *df = -2.0*qoff - coeff*q/sdis;

  *ddf = -2.0 + coeff*betasq/(dis*sdis);

  return 1;
}
// }}}

// MathDump {{{
void MathDump(double *fpar){
  int ii = (int) fpar[FP_INDEX_X];
  int jj = (int) fpar[FP_INDEX_Y];
  int kk = (int) fpar[FP_INDEX_Z];

  double tau   = fpar[FP_TAU];
  double d     = fpar[FP_D];
  double Ssq   = fpar[FP_SSQ];
  double gamma = fpar[FP_GAMMA];

  fprintf(stdout,"INDEX:  (i,j) = (%d, %d, %d)\n",ii,jj,kk);
  fprintf(stdout,"COORDS: (x,y,z) = (%g, %g, %g)\n",
          fpar[FP_COORD_X],fpar[FP_COORD_Y],fpar[FP_COORD_Z]);
  fprintf(stdout,"---------- MATHEMATICA TRACE ----------\n");
  fprintf(stdout,"tau = %20.16f\n",tau);
  fprintf(stdout,"d = %20.16f\n",d);
  fprintf(stdout,"g = %20.16f\n",gamma);
  fprintf(stdout,"Ssq = %20.16f\n",Ssq);
}
// }}}

// func_p_isentropic {{{
int func_p_isentropic(double *f, double rho0, double *fpar)
{

  double Ssq   = fpar[FP_SSQ];
  double D     = fpar[FP_D];
  double Gamma = fpar[FP_GAMMA];
  double kappa = fpar[FP_KAPPA];

  double h = rho0 + Gamma/(Gamma - 1.0)*kappa*pow(rho0,Gamma);
  double W = D/rho0;

  *f  = (h*h*W*W)*(W*W - 1.0) - Ssq;

  return 1;

}
// }}}

// func_p_d_isentropic {{{
int func_p_d_isentropic(double *f, double *df, double rho0, double *fpar)
{

  double Ssq   = fpar[FP_SSQ];
  double D     = fpar[FP_D];
  double g     = fpar[FP_GAMMA];
  double kappa = fpar[FP_KAPPA];

  double P = kappa*pow(rho0,g);
  double h = rho0 + g/(g - 1.0)*P;
  double W = D/rho0;

  /* Maple vars */
  double t2, t5, t7, t8, t10, t11, t12, t17, t21, t22, t23, t24;
  double t51, t52, t53;

  *f  = (h*h*W*W)*(W*W - 1.0) - Ssq;

      t2 = 1.0/(g-1.0);
      t5 = P;
      t7 = rho0+g*t2*t5;
      t8 = D*D;
      t10 = rho0*rho0;
      t11 = 1.0/t10;
      t12 = g*g;
      t17 = 1.0+t12*t2*t5/rho0;
      t21 = t7*t7;
      t22 = t21*t8;
      t23 = t10*rho0;
      t24 = 1.0/t23;
      t51 = t8*t8;
      t52 = 1.0/t51;
      t53 = t10*t10;
     *df = (2.0*t7*t8*t11*t17-2.0*t22*t24)*(t8*t11-0.1E1)-2.0*(t22*t11)*t8*t24;


  return 1;
}
// }}}


} // end of namespace fluid math

