//
// Created by milinda on 10/5/17.
/**
 *@author Milinda Fernando
 *School of Computing, University of Utah
 *@brief check all teh bssn derivs for user specified function.
 */
//

#ifndef SFCSORTBENCH_CHECKDERIVATIVES_H
#define SFCSORTBENCH_CHECKDERIVATIVES_H

#include <functional>

#include "TreeNode.h"
#include "derivs.h"
#include "iostream"
#include "mesh.h"
#include "oct2vtk.h"
#include "octUtils.h"

#define NUM_VARS 10

#define FX    0
#define Dx_FX 1
#define Dy_FX 2
#define Dz_FX 3

#define DxDx_FX 4
#define DxDy_FX 5
#define DxDz_FX 6
#define DyDy_FX 7
#define DyDz_FX 8
#define DzDz_FX 9

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    if (argc < 4) {
        if (!rank)
            std::cout << "Usage: " << argv[0]
                      << " maxDepth wavelet_tol partition_tol eleOrder"
                      << std::endl;
        return 0;
    }

    m_uiMaxDepth         = atoi(argv[1]);
    double wavelet_tol   = atof(argv[2]);
    double partition_tol = atof(argv[3]);
    unsigned int eOrder  = atoi(argv[4]);

    if (!rank) {
        std::cout << YLW << "maxDepth: " << m_uiMaxDepth << NRM << std::endl;
        std::cout << YLW << "wavelet_tol: " << wavelet_tol << NRM << std::endl;
        std::cout << YLW << "partition_tol: " << partition_tol << NRM
                  << std::endl;
        std::cout << YLW << "eleOrder: " << eOrder << NRM << std::endl;
    }

    _InitializeHcurve(m_uiDim);

    const double d_min = -0.5;
    const double d_max = 0.5;

    std::function<double(double, double, double)> fx[NUM_VARS];

    /*fx[FX]=[d_min,d_max](const double x,const double y,const double z){ return
    sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z);};
    // 1st order derivatives
    fx[Dx_FX]=[d_min,d_max](const double x,const double y,const double z){
    return (2*M_PI*cos(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z));};
    fx[Dy_FX]=[d_min,d_max](const double x,const double y,const double z){
    return (2*M_PI*sin(2*M_PI*x)*cos(2*M_PI*y)*sin(2*M_PI*z));};
    fx[Dz_FX]=[d_min,d_max](const double x,const double y,const double z){
    return (2*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y)*cos(2*M_PI*z));};

    // 2nd order derivatives
    fx[DxDx_FX]=[d_min,d_max](const double x,const double y,const double z){
    return (-4*M_PI*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z));};
    fx[DxDy_FX]=[d_min,d_max](const double x,const double y,const double z){
    return ( 4*M_PI*M_PI*cos(2*M_PI*x)*cos(2*M_PI*y)*sin(2*M_PI*z));};
    fx[DxDz_FX]=[d_min,d_max](const double x,const double y,const double z){
    return ( 4*M_PI*M_PI*cos(2*M_PI*x)*sin(2*M_PI*y)*cos(2*M_PI*z));};
    fx[DyDy_FX]=[d_min,d_max](const double x,const double y,const double z){
    return (-4*M_PI*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z));};
    fx[DxDz_FX]=[d_min,d_max](const double x,const double y,const double z){
    return ( 4*M_PI*M_PI*sin(2*M_PI*x)*cos(2*M_PI*y)*cos(2*M_PI*z));};
    fx[DzDz_FX]=[d_min,d_max](const double x,const double y,const double z){
    return (-4*M_PI*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z));};*/

    fx[FX] = [d_min, d_max](const double x, const double y, const double z) {
        return sin(2 * M_PI *
                   ((x / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)) *
               sin(2 * M_PI *
                   ((y / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)) *
               sin(2 * M_PI *
                   ((z / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min));
    };
    // 1st order derivatives
    fx[Dx_FX] = [d_min, d_max](const double x, const double y, const double z) {
        return (2 * M_PI *
                cos(2 * M_PI *
                    ((x / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)) *
                sin(2 * M_PI *
                    ((y / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)) *
                sin(2 * M_PI *
                    ((z / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)));
    };
    fx[Dy_FX] = [d_min, d_max](const double x, const double y, const double z) {
        return (2 * M_PI *
                sin(2 * M_PI *
                    ((x / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)) *
                cos(2 * M_PI *
                    ((y / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)) *
                sin(2 * M_PI *
                    ((z / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)));
    };
    fx[Dz_FX] = [d_min, d_max](const double x, const double y, const double z) {
        return (2 * M_PI *
                sin(2 * M_PI *
                    ((x / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)) *
                sin(2 * M_PI *
                    ((y / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)) *
                cos(2 * M_PI *
                    ((z / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)));
    };

    // 2nd order derivatives
    fx[DxDx_FX] = [d_min, d_max](const double x, const double y,
                                 const double z) {
        return (-4 * M_PI * M_PI *
                sin(2 * M_PI *
                    ((x / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)) *
                sin(2 * M_PI *
                    ((y / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)) *
                sin(2 * M_PI *
                    ((z / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)));
    };
    fx[DxDy_FX] = [d_min, d_max](const double x, const double y,
                                 const double z) {
        return (4 * M_PI * M_PI *
                cos(2 * M_PI *
                    ((x / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)) *
                cos(2 * M_PI *
                    ((y / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)) *
                sin(2 * M_PI *
                    ((z / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)));
    };
    fx[DxDz_FX] = [d_min, d_max](const double x, const double y,
                                 const double z) {
        return (4 * M_PI * M_PI *
                cos(2 * M_PI *
                    ((x / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)) *
                sin(2 * M_PI *
                    ((y / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)) *
                cos(2 * M_PI *
                    ((z / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)));
    };
    fx[DyDy_FX] = [d_min, d_max](const double x, const double y,
                                 const double z) {
        return (-4 * M_PI * M_PI *
                sin(2 * M_PI *
                    ((x / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)) *
                sin(2 * M_PI *
                    ((y / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)) *
                sin(2 * M_PI *
                    ((z / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)));
    };
    fx[DyDz_FX] = [d_min, d_max](const double x, const double y,
                                 const double z) {
        return (4 * M_PI * M_PI *
                sin(2 * M_PI *
                    ((x / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)) *
                cos(2 * M_PI *
                    ((y / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)) *
                cos(2 * M_PI *
                    ((z / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)));
    };
    fx[DzDz_FX] = [d_min, d_max](const double x, const double y,
                                 const double z) {
        return (-4 * M_PI * M_PI *
                sin(2 * M_PI *
                    ((x / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)) *
                sin(2 * M_PI *
                    ((y / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)) *
                sin(2 * M_PI *
                    ((z / (1u << m_uiMaxDepth)) * (d_max - d_min) + d_min)));
    };

    /* std::function<double(double,double,double)> fx1[NUM_VARS];
     fx1[FX]=[d_min,d_max](const double x,const double y,const double z){ return
     (sin(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min)));};
     fx1[Dx_FX]=[d_min,d_max](const double x,const double y,const double z){
     return
     (2*M_PI*(1.0/(1u<<m_uiMaxDepth)*(d_max-d_min)))*(cos(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min)));};
     fx1[Dy_FX]=[d_min,d_max](const double x,const double y,const double z){
     return
     (2*M_PI*(1.0/(1u<<m_uiMaxDepth)*(d_max-d_min)))*(sin(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*cos(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min)));};
     fx1[Dz_FX]=[d_min,d_max](const double x,const double y,const double z){
     return
     (2*M_PI*(1.0/(1u<<m_uiMaxDepth)*(d_max-d_min)))*(sin(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*cos(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min)));};
     fx1[DxDx_FX]=[d_min,d_max](const double x,const double y,const double z){
     return
     (-pow(2*M_PI*((1.0/(1u<<m_uiMaxDepth))*(d_max-d_min)),2)*(sin(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))));};
     fx1[DxDy_FX]=[d_min,d_max](const double x,const double y,const double z){
     return (
     pow(2*M_PI*((1.0/(1u<<m_uiMaxDepth))*(d_max-d_min)),2)*(cos(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*cos(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))));};
     fx1[DxDz_FX]=[d_min,d_max](const double x,const double y,const double z){
     return (
     pow(2*M_PI*((1.0/(1u<<m_uiMaxDepth))*(d_max-d_min)),2)*(cos(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*cos(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))));};
     fx1[DyDy_FX]=[d_min,d_max](const double x,const double y,const double z){
     return
     (-pow(2*M_PI*((1.0/(1u<<m_uiMaxDepth))*(d_max-d_min)),2)*(sin(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))));};
     fx1[DyDz_FX]=[d_min,d_max](const double x,const double y,const double z){
     return (
     pow(2*M_PI*((1.0/(1u<<m_uiMaxDepth))*(d_max-d_min)),2)*(sin(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*cos(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*cos(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))));};
     fx1[DzDz_FX]=[d_min,d_max](const double x,const double y,const double z){
     return
     (-pow(2*M_PI*((1.0/(1u<<m_uiMaxDepth))*(d_max-d_min)),2)*(sin(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))));};*/

    ot::TreeNode root(m_uiDim, m_uiMaxDepth);

    std::vector<ot::TreeNode> grid;
    function2Octree(fx[FX], grid, m_uiMaxDepth, wavelet_tol, eOrder, comm);

    std::vector<ot::TreeNode> tmpOcts;
    SFC::parSort::SFC_treeSort(grid, tmpOcts, tmpOcts, tmpOcts, partition_tol,
                               m_uiMaxDepth, root, ROOT_ROTATION, 1,
                               TS_REMOVE_DUPLICATES, 2, comm);
    std::swap(tmpOcts, grid);
    tmpOcts.clear();

    SFC::parSort::SFC_treeSort(grid, tmpOcts, tmpOcts, tmpOcts, partition_tol,
                               m_uiMaxDepth, root, ROOT_ROTATION, 1,
                               TS_CONSTRUCT_OCTREE, 2, comm);
    std::swap(tmpOcts, grid);
    tmpOcts.clear();

    SFC::parSort::SFC_treeSort(grid, tmpOcts, tmpOcts, tmpOcts, partition_tol,
                               m_uiMaxDepth, root, ROOT_ROTATION, 1,
                               TS_BALANCE_OCTREE, 2, comm);
    std::swap(tmpOcts, grid);
    tmpOcts.clear();

    DendroIntL localSz;
    DendroIntL globalSz;

    localSz = grid.size();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) std::cout << "Num balance octs: " << globalSz << std::endl;

    ot::Mesh pMesh(grid, 1, eOrder, comm);

    localSz = pMesh.getNumLocalMeshNodes();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) std::cout << "Num total vertices: " << globalSz << std::endl;

    double** varFx       = new double*[NUM_VARS];  // original vars
    double** varSFx      = new double*[NUM_VARS];  // stencil computed.
    double** varDiff     = new double*[NUM_VARS];  // stencil computed.
    double** unzipVarSFx = new double*[NUM_VARS];

    for (unsigned int i = 0; i < NUM_VARS; i++) {
        varFx[i]       = pMesh.createVector<double>(fx[i]);
        varSFx[i]      = pMesh.createVector<double>(0);
        varDiff[i]     = pMesh.createVector<double>(0);
        unzipVarSFx[i] = pMesh.createUnZippedVector<double>();
    }

    for (unsigned int i = 0; i < NUM_VARS; i++) {
        pMesh.performGhostExchange(varFx[i]);
    }

    pMesh.unzip(varFx[FX], unzipVarSFx[FX]);

    const std::vector<ot::Block>& blockList = pMesh.getLocalBlockList();

    unsigned int offset;
    unsigned int sz[3];
    double h[3];

    const Point p_min(d_min, d_min, d_min);
    const Point p_max(d_max, d_max, d_max);

    unsigned int bflag;

    for (unsigned int blk = 0; blk < blockList.size(); blk++) {
        offset = blockList[blk].getOffset();
        bflag  = blockList[blk].getBlkNodeFlag();

        sz[0]  = blockList[blk].getAllocationSzX();
        sz[1]  = blockList[blk].getAllocationSzY();
        sz[2]  = blockList[blk].getAllocationSzZ();

        h[0]   = blockList[blk].computeDx(p_min, p_max);
        h[1]   = blockList[blk].computeDy(p_min, p_max);
        h[2]   = blockList[blk].computeDz(p_min, p_max);

        /*h[0]=blockList[blk].computeGridDx();
        h[1]=blockList[blk].computeGridDy();
        h[2]=blockList[blk].computeGridDz();*/

        deriv42_x(unzipVarSFx[Dx_FX] + offset, unzipVarSFx[FX] + offset, h[0],
                  sz, bflag);
        deriv42_y(unzipVarSFx[Dy_FX] + offset, unzipVarSFx[FX] + offset, h[1],
                  sz, bflag);
        deriv42_z(unzipVarSFx[Dz_FX] + offset, unzipVarSFx[FX] + offset, h[2],
                  sz, bflag);

        deriv42_xx(unzipVarSFx[DxDx_FX] + offset, unzipVarSFx[FX] + offset,
                   h[0], sz, bflag);
        deriv42_yy(unzipVarSFx[DyDy_FX] + offset, unzipVarSFx[FX] + offset,
                   h[1], sz, bflag);
        deriv42_zz(unzipVarSFx[DzDz_FX] + offset, unzipVarSFx[FX] + offset,
                   h[2], sz, bflag);
    }

    for (unsigned int blk = 0; blk < blockList.size(); blk++) {
        offset = blockList[blk].getOffset();
        bflag  = blockList[blk].getBlkNodeFlag();

        sz[0]  = blockList[blk].getAllocationSzX();
        sz[1]  = blockList[blk].getAllocationSzY();
        sz[2]  = blockList[blk].getAllocationSzZ();

        h[0]   = blockList[blk].computeDx(p_min, p_max);
        h[1]   = blockList[blk].computeDy(p_min, p_max);
        h[2]   = blockList[blk].computeDz(p_min, p_max);

        /*h[0]=blockList[blk].computeGridDx();
        h[1]=blockList[blk].computeGridDy();
        h[2]=blockList[blk].computeGridDz();*/

        deriv42_y(unzipVarSFx[DxDy_FX] + offset, unzipVarSFx[Dx_FX] + offset,
                  h[1], sz, bflag);
        deriv42_z(unzipVarSFx[DxDz_FX] + offset, unzipVarSFx[Dx_FX] + offset,
                  h[2], sz, bflag);
        deriv42_z(unzipVarSFx[DyDz_FX] + offset, unzipVarSFx[Dy_FX] + offset,
                  h[2], sz, bflag);
    }

    double norm_inf[NUM_VARS];

    if (pMesh.isActive())
        for (unsigned int i = 0; i < NUM_VARS; i++) {
            pMesh.zip(unzipVarSFx[i], varSFx[i]);
            norm_inf[i] = normLInfty(
                varFx[i] + pMesh.getNodeLocalBegin(),
                varSFx[i] + pMesh.getNodeLocalBegin(),
                (pMesh.getNodeLocalEnd() - pMesh.getNodeLocalBegin()),
                pMesh.getMPICommunicator());

            if (!rank) {
                std::cout << "var: " << i << " l_inf: " << norm_inf[i]
                          << std::endl;
            }

            for (unsigned int node = pMesh.getNodeLocalBegin();
                 node < pMesh.getNodeLocalEnd(); node++) {
                varDiff[i][node] = varFx[i][node] - varSFx[i][node];
            }
        }

    for (unsigned int i = 0; i < NUM_VARS; i++) {
        pMesh.performGhostExchange(varFx[i]);
        pMesh.performGhostExchange(varSFx[i]);
        pMesh.performGhostExchange(varDiff[i]);
    }

    const char* pDataNames[] = {"fx",      "dx_fx",   "dy_fx",   "dz_fx",
                                "dxdx_fx", "dxdy_fx", "dxdz_fx", "dydy_fx",
                                "dydz_fx", "dzdz_fx"};
    io::vtk::mesh2vtuFine(&pMesh, "fx", 0, NULL, NULL, NUM_VARS, pDataNames,
                          (const double**)varFx);
    io::vtk::mesh2vtuFine(&pMesh, "fx_s", 0, NULL, NULL, NUM_VARS, pDataNames,
                          (const double**)varSFx);
    io::vtk::mesh2vtuFine(&pMesh, "fx_d", 0, NULL, NULL, NUM_VARS, pDataNames,
                          (const double**)varDiff);

    for (unsigned int i = 0; i < NUM_VARS; i++) {
        delete[] varFx[i];
        delete[] varSFx[i];
        delete[] varDiff[i];
        delete[] unzipVarSFx[i];
    }

    delete[] varFx;
    delete[] varSFx;
    delete[] varDiff;
    delete[] unzipVarSFx;

    MPI_Finalize();
    return 0;
}

#endif  // SFCSORTBENCH_CHECKDERIVATIVES_H
