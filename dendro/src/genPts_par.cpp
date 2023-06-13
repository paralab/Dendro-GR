/*
 *
 * @author: Milinda Fernando
 * School of Computing, University of Utah
 * @date: 11/10/2015
 *
 *
 * Contains the Normal(Gaussian) and Logarithmic normal random number (octants) generator based on new c+11 rand engines.
 *
 *
 * */

#include "genPts_par.h"


void genGauss(const double& sd, const DendroIntL numPts, int dim, char * filePrefix,MPI_Comm comm)
{


    int rank=0,size=1;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&size);

    char ptsFileName[256];

//    MPI_Comm_size(comm,&size);
//    MPI_Comm_rank(comm,&rank);


        std::random_device rd;
        std::mt19937 gen(rd());

    sprintf(ptsFileName, "%s%d_%d.pts", filePrefix, rank, size);

    std::normal_distribution<> d(0.5,sd);
    long int ptsSize=numPts;
    double* xyz = new double[ptsSize*dim];
    double temp=0;
    for(DendroIntL i=0;i<dim*numPts;i++)
    {
        temp=(double)d(gen);

        if(temp<0) temp=-1*temp;
        if(temp>=1) temp=1.0/(2*temp);

        xyz[i]=temp;

    }

    FILE* outfile;
    outfile = fopen(ptsFileName,"wb");
    fwrite(&numPts, sizeof(unsigned int), 1, outfile);
    fwrite(xyz, sizeof(double), (dim*numPts), outfile);
    fclose(outfile);

    delete [] xyz;
    xyz=NULL;


}


void genGauss(const double& sd, const DendroIntL numPts, int dim,double *xyz)
{

    std::random_device rd;
    std::mt19937_64 gen(rd());
    //std::mt19937 gen;


    std::normal_distribution<double> d(0.5,sd);
    //long int ptsSize=numPts;
    double temp=0;
    DendroIntL totPts=numPts*dim;
    for(DendroIntL i=0;i<totPts;i++)
    {
        temp=(double)d(gen);

        if(temp<0) temp=-1*temp;
        if(temp>=1) temp=1.0/(2*temp);

        xyz[i]=temp;

    }
     return ;

}


void genUniformRealDis(const double& sd, const DendroIntL numPts, int dim,double *xyz)
{

    std::random_device rd;
    std::mt19937_64 gen(rd());
    //std::mt19937 gen;


    //std::normal_distribution<>d(0.5,sd);
    std::uniform_real_distribution<> d(0,1);
    //std::uniform_real_distribution<> d1(0.6,0.7);

    double temp=0;
    DendroIntL totPts=numPts*dim;
    for(DendroIntL i=0;i<totPts;i++)
    {
        temp=(double)d(gen);
        if(i>(totPts/2))
        {
            temp=(double)d(gen);
        }

        if(temp<0) temp=-1*temp;
        if(temp>=1) temp=1.0/(2*temp);

        xyz[i]=temp;

    }
    return ;

}


void genLogarithmicGauss(const double& sd, const DendroIntL numPts, int dim,double *xyz)
{

    std::random_device rd;
    std::mt19937_64 gen(rd());


    std::lognormal_distribution<> d(0.5,sd);

    //long int ptsSize=numPts;
    double temp=0;
    DendroIntL totPts=numPts*dim;
    for(DendroIntL i=0;i<totPts;i++)
    {
        temp=(double)d(gen);

        if(temp<0) temp=-1*temp;
        if(temp>=1) temp=1.0/(2*temp);

        xyz[i]=temp;

    }
    return ;

}




void genLogarithmicGauss(const double& sd, const DendroIntL numPts, int dim, char * filePrefix)
{


    int rank=0,size=1;

    char ptsFileName[256];

//    MPI_Comm_size(comm,&size);
//    MPI_Comm_rank(comm,&rank);


    std::random_device rd;
    std::mt19937 gen(rd());

    sprintf(ptsFileName, "%s%d_%d.pts", filePrefix, rank, size);

    std::lognormal_distribution<>d(0.5,sd);

    
    double* xyz = new double[dim*numPts];
    double temp=0;
    for(DendroIntL i=0;i<(long long)dim*numPts;i++)
    {
        temp=(double)d(gen);

        if(temp<0) temp=-1*temp;
        if(temp>=1) temp=1.0/(2*temp);

        xyz[i]=temp;

    }

    FILE* outfile;
    outfile = fopen(ptsFileName,"wb");
    fwrite(&numPts, sizeof(unsigned int), 1, outfile);
    fwrite(xyz, sizeof(double), (dim*numPts), outfile);
    fclose(outfile);

    delete [] xyz;



}



void pts2Octants(std::vector<ot::TreeNode> & pNodes,double * pts, DendroIntL totPts, unsigned int dim ,unsigned int maxDepth)
{
    pNodes.clear();
#ifdef DIM_2

    for (DendroIntL i = 0; i < totPts; i += 2) {
        if ((pts[i] > 0.0) &&
            (pts[i + 1] > 0.0) &&
            (((unsigned int) (pts[i] * ((double) (1u << maxDepth)))) < (1u << maxDepth)) &&
            (((unsigned int) (pts[i + 1] * ((double) (1u << maxDepth)))) < (1u << maxDepth)) ) {

            pNodes.push_back(ot::TreeNode((unsigned int) (pts[i] * (double) (1u << maxDepth)),
                                            (unsigned int) (pts[i + 1] * (double) (1u << maxDepth)),
                                            0,maxDepth, dim, maxDepth));
        }
    }

#else
    for (DendroIntL i = 0; i < totPts; i += 3) {
        if ((pts[i] > 0.0) &&
            (pts[i + 1] > 0.0)
            && (pts[i + 2] > 0.0) &&
            (((unsigned int) (pts[i] * ((double) (1u << maxDepth)))) < (1u << maxDepth)) &&
            (((unsigned int) (pts[i + 1] * ((double) (1u << maxDepth)))) < (1u << maxDepth)) &&
            (((unsigned int) (pts[i + 2] * ((double) (1u << maxDepth)))) < (1u << maxDepth))) {

            pNodes.push_back(ot::TreeNode((unsigned int) (pts[i] * (double) (1u << maxDepth)),
                                          (unsigned int) (pts[i + 1] * (double) (1u << maxDepth)),
                                          (unsigned int) (pts[i + 2] * (double) (1u << maxDepth)),
                                          maxDepth, dim, maxDepth));
        }
    }
#endif





}