#include <zoltan_hilbert.h>
#include <assert.h>
#include <parUtils.h>
#include "../include/zoltan_hilbert.h"

#define NUM_GID_SIZE 1
#define NUM_LID_SIZE 1




void user_sizes_node(void *data,int num_gid_entries,int num_lid_entries, int nids, ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int* nbytes, int *ierr)
{
    // MESH_DATA * mesh =(MESH_DATA*) data;
    *ierr = ZOLTAN_OK;
    for (int i=0; i<nids; ++i) nbytes[i] = sizeof(Node_Type);
}


void user_pack_node(void *data,int num_gid_entries,int num_lid_entries, int nids, ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int* dest, int* size, int* index, char *buf, int *ierr)
{
    *ierr = ZOLTAN_OK;
    MESH_DATA * mesh =(MESH_DATA*) data;
    Node_Type *node= (Node_Type* )buf;

    ZOLTAN_ID_TYPE lid;

    for (int i=0; i<nids; ++i) {
        lid = local_id[i];
        node[i].x=mesh->x[lid];
        node[i].y=mesh->y[lid];
        node[i].z=mesh->z[lid];

        node[i].globalID=mesh->globalIds[lid];


    }

    return;




}

void user_unpack_node(void *data,int gidSize,int num_ids, ZOLTAN_ID_PTR global_id, int * size,int * idx , char *buf, int *ierr)
{
    *ierr = ZOLTAN_OK;
    MESH_DATA * mesh =(MESH_DATA*) data;

    Node_Type * node =(Node_Type *)(buf);

    mesh->globalIds = (ZOLTAN_ID_TYPE *) realloc(mesh->globalIds,sizeof(ZOLTAN_ID_TYPE)*num_ids);
    mesh->localIds = (ZOLTAN_ID_TYPE *) realloc(mesh->localIds,sizeof(ZOLTAN_ID_TYPE)*num_ids);

    mesh->x=(ZOLTAN_ID_TYPE *) realloc(mesh->x,sizeof(ZOLTAN_ID_TYPE)*num_ids);
    mesh->y=(ZOLTAN_ID_TYPE *) realloc(mesh->x,sizeof(ZOLTAN_ID_TYPE)*num_ids);
    mesh->z=(ZOLTAN_ID_TYPE *) realloc(mesh->x,sizeof(ZOLTAN_ID_TYPE)*num_ids);


    mesh->numLocalPts=num_ids;


    for(int i=0;i<num_ids;i++)
    {
        mesh->x[i]=node[i].x;
        mesh->y[i]=node[i].y;
        mesh->z[i]=node[i].z;
        mesh->localIds[i]=i;
        mesh->globalIds[i]=node[i].globalID;

   }

}




void ZoltanHibertLBandSort(DendroIntL grainSz,unsigned int dim ,unsigned int maxDepth,double tol, std::vector<double> &stats,unsigned int distribution,MPI_Comm comm)
{

    int rank,npes;
    int rc ;// error code
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);




    struct Zoltan_Struct *zz;
    ZOLTAN_VOID_FN *fn = NULL;
    int *fndata = NULL;

    DendroIntL numGlobalPts;
    MPI_Allreduce(&grainSz, &numGlobalPts, 1, MPI_LONG_LONG, MPI_SUM, comm);

    //double tol= ((double)grainSz/(double)tolOCt);
    MESH_DATA mesh;

    /******************************************************************
     ** Generate mesh points.
     ******************************************************************/

    DendroIntL totalgrainSz=grainSz*dim;

    double *pts =new double[totalgrainSz];

    if(distribution==0)
        genGauss(0.5,grainSz,dim,pts);
    else if(distribution==1)
    {
        genLogarithmicGauss(0.5,grainSz,dim,pts);

    }else if(distribution==2)
    {
        genUniformRealDis(0.5,grainSz,dim,pts);
    }else
    {
        genGauss(0.5,grainSz,dim,pts);
    }

    mesh.numLocalPts =grainSz;
    mesh.numGlobalPoints=numGlobalPts;

    mesh.globalIds =new ZOLTAN_ID_TYPE[grainSz];
    mesh.localIds =new ZOLTAN_ID_TYPE[grainSz];
    mesh.x=new ZOLTAN_ID_TYPE[grainSz];
    mesh.y=new ZOLTAN_ID_TYPE[grainSz];
    mesh.z=new ZOLTAN_ID_TYPE[grainSz];


    for (DendroIntL i = 0; i < grainSz; i ++) {
        if ((pts[3*i] > 0.0) &&
            (pts[3*i + 1] > 0.0)
            && (pts[3*i + 2] > 0.0) &&
            (((unsigned int) (pts[3*i] * ((double) (1u << maxDepth)))) < (1u << maxDepth)) &&
            (((unsigned int) (pts[3*i + 1] * ((double) (1u << maxDepth)))) < (1u << maxDepth)) &&
            (((unsigned int) (pts[3*i + 2] * ((double) (1u << maxDepth)))) < (1u << maxDepth))) {


            mesh.x[i] = (unsigned int) (pts[3 * i] * (double) (1u << maxDepth));
            mesh.y[i] = (unsigned int) (pts[3 * i + 1] * (double) (1u << maxDepth));
            mesh.z[i] = (unsigned int) (pts[3 * i + 2] * (double) (1u << maxDepth));


            mesh.globalIds[i]= rank * i;
            mesh.localIds[i]=i;

        }

    }

    delete[](pts);
    auto zz_start=std::chrono::system_clock::now();//MPI_Wtime();





       zz = Zoltan_Create(comm);
       /* General parameters */

       Zoltan_Set_Param(zz, "DEBUG_LEVEL", "1");
       Zoltan_Set_Param(zz, "LB_METHOD", "HSFC");
       Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", std::to_string(NUM_GID_SIZE).c_str());
       Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", std::to_string(NUM_LID_SIZE).c_str());
       Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
       Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");
       Zoltan_Set_Param(zz, "IMBALANCE_TOL", std::to_string(tol).c_str());
       Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", std::to_string(npes).c_str());
       Zoltan_Set_Param(zz, "NUM_LOCAL_PARTS", "1");
       Zoltan_Set_Param(zz, "ORDER_METHOD", "LOCAL_HSFC");


       Zoltan_Set_Param(zz, "AUTO_MIGRATE", "FALSE");


       /******************************************************************
       ** Create a Zoltan library structure for this instance of load
       ** balancing.  Set the parameters and query functions that will
       ** govern the library's calculation.  See the Zoltan User's
       ** Guide for the definition of these and many other parameters.
       ******************************************************************/

       //Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS","8");
       // HSFC Parameters
       // Defaults: KEEP_CUTS = 0
       //           REDUCE_DIMENSIONS = 0
       //           DEGENERATE_RATIO = 10

       Zoltan_Set_Param(zz, "KEEP_CUTS", "1"); // KEEP_CUTS = 0 means don't keep cuts, KEEP_CUTS = 1 means keep cuts,

       // REDUCE_DIMENSIONS: When a 3 dimensional geometry is almost flat, it may make more sense to treat it as a 2 dimensional geometry when applying the HSFC algorithm.
       // (Coordinate values in the omitted direction are ignored for the purposes of partitioning.) If this parameter
       // is set to 1, a 3 dimensional geometry will be treated as 2 dimensional if is very flat, or 1 dimensional if it very thin.
       // And a 2 dimensional geometry will be treated as 1 dimensional if it is very thin. Turning this parameter on
       // removes the possibility that disconnected parts will appear on the surface of a flat 3 dimensional object.
       //zz->Set_Param( "REDUCE_DIMENSIONS", "0");

       // DEGENERATE_RATIO:
       // If the REDUCE_DIMENSIONS parameter is set, then this parameter determines when a geometry is considered to be flat.
       // A bounding box which is oriented to the geometry is constructed, and the lengths of its sides are tested against a ratio of 1
       //zz->Set_Param( "DEGENERATE_RATIO", "10"); // default value

       /* Query functions, to provide geometry to Zoltan */

       Zoltan_Set_Num_Obj_Fn(zz, get_number_of_objects, &mesh);
       Zoltan_Set_Obj_List_Fn(zz, get_object_list, &mesh);
       Zoltan_Set_Num_Geom_Fn(zz, get_num_geometry, &mesh);
       Zoltan_Set_Geom_Multi_Fn(zz, get_geometry_list, &mesh);


       int changes = 0, numGidEntries = 1, numLidEntries = 1, numImport = 0, numExport = 0;
       ZOLTAN_ID_PTR importGlobalGids = NULL, importLocalGids = NULL, exportGlobalGids = NULL, exportLocalGids = NULL;
       int *importProcs = NULL, *importToPart = NULL, *exportProcs = NULL, *exportToPart = NULL;


    if(npes<=1) {


        auto zz_start=std::chrono::system_clock::now();//MPI_Wtime();

        ZOLTAN_ID_TYPE * permutedIds=new ZOLTAN_ID_TYPE[grainSz];
        Zoltan_Order(zz,NUM_GID_SIZE,mesh.numLocalPts,mesh.globalIds,permutedIds);

        auto zz_end=std::chrono::system_clock::now();//MPI_Wtime();

        double zz_time= (std::chrono::duration_cast<std::chrono::milliseconds>((zz_end-zz_start)).count())/(MILLISECOND_CONVERSION);


        stats.push_back(0);
        stats.push_back(0);
        stats.push_back(0);

        stats.push_back(0);
        stats.push_back(0);
        stats.push_back(0);

        stats.push_back(0);
        stats.push_back(0);
        stats.push_back(0);

        stats.push_back(zz_time);
        stats.push_back(zz_time);
        stats.push_back(zz_time);






        return ;
     }

/*    *//******************************************************************
    ** Zoltan can now partition the vertices in the simple mesh.
    ** In this simple example, we assume the number of partitions is
    ** equal to the number of processes.  Process rank 0 will own
    ** partition 0, process rank 1 will own partition 1, and so on.
    ******************************************************************/



        auto zz_splitter_start=std::chrono::system_clock::now();//MPI_Wtime();

        rc = Zoltan_LB_Partition(zz,  /*input (all remaining fields are output)*/
                                 &changes,        /* 1 if partitioning was changed, 0 otherwise */
                                 &numGidEntries,   /*Number of integers used for a global ID */
                                 &numLidEntries,   /*Number of integers used for a local ID */
                                 &numImport,       /*Number of vertices to be sent to me */
                                 &importGlobalGids,  /* Global IDs of vertices to be sent to me*/
                                 &importLocalGids,    /*Local IDs of vertices to be sent to me*/
                                 &importProcs,     /*Process rank for source of each incoming vertex*/
                                 &importToPart,    /*New partition for each incoming vertex*/
                                 &numExport,       /*Number of vertices I must send to other processes*/
                                 &exportGlobalGids,   /*Global IDs of the vertices I must send*/
                                 &exportLocalGids,    /*Local IDs of the vertices I must send*/
                                 &exportProcs,     /*Process to which I send each of the vertices*/
                                 &exportToPart);   /*Partition to which each vertex will belong*/

        auto zz_splitter_end=std::chrono::system_clock::now();//MPI_Wtime();


        /*if(rc==ZOLTAN_OK)
        {
            if(!rank) std::cout<<"=================================================================="<<std::endl;
            if(!rank) std::cout<<"                     Zoltan LB PARTIION FINISHED                  "<<std::endl;
            if(!rank) std::cout<<"=================================================================="<<std::endl;
        }*/

        //MPI_Barrier(comm);


        //std::cout << "rank : " << rank << " numPts: " << mesh.numLocalPts << " sendCnt : " << numExport << " recvCnt: " << numImport << std::endl;


        auto zz_all2all_start=std::chrono::system_clock::now();//MPI_Wtime();

        ZOLTAN_ID_TYPE *sendCnts = new ZOLTAN_ID_TYPE[npes];
        ZOLTAN_ID_TYPE *recvCnts = new ZOLTAN_ID_TYPE[npes];

        ZOLTAN_ID_TYPE *Cnts = new ZOLTAN_ID_TYPE[npes];

        for (unsigned int i = 0; i < npes; i++) {
            Cnts[i] = 0;
            sendCnts[i] = 0;
            recvCnts[i] = 0;
        }


        ZOLTAN_ID_TYPE totalSendCnt = 0;

        for (unsigned int i = 0; i < numExport; i++) {
//        if(!rank) std::cout <<" Exp Proc: "<<exportProcs[i]<<std::endl;
            sendCnts[exportProcs[i]]++;


        }

        for (unsigned int i = 0; i < numImport; i++)
            recvCnts[importProcs[i]]++;


        ZOLTAN_ID_TYPE *senCntsOffset = new ZOLTAN_ID_TYPE[npes];
        ZOLTAN_ID_TYPE *recvCntsOffset = new ZOLTAN_ID_TYPE[npes];

        senCntsOffset[0] = 0;
        recvCntsOffset[0] = 0;

        totalSendCnt = sendCnts[0];

        for (int i = 1; i < npes; i++) {
            senCntsOffset[i] = senCntsOffset[i - 1] + sendCnts[i - 1];
            recvCntsOffset[i] = recvCntsOffset[i - 1] + recvCntsOffset[i - 1];
            totalSendCnt += sendCnts[i];
        }


       /* MPI_Comm newComm;
        MPI_Comm commOld=comm;
        par::splitComm2way(!numExport,&newComm,comm);


        comm=newComm;*/


        //assert(!numExport);

        std::vector<Node_Type> sendData;
        std::vector<Node_Type> recvData;

//        std::vector<Node_Type> sendData(numExport);
//        std::vector<Node_Type> recvData(numImport);


        if(numExport) sendData.resize(numExport);
        if(numImport) recvData.resize(numImport);

      /*  MPI_Datatype MPI_NODETYPE;
        MPI_Type_contiguous(4, MPI_UNSIGNED, &MPI_NODETYPE);
        MPI_Type_commit(&MPI_NODETYPE);*/


        Node_Type temp;

        for (unsigned int i = 0; i < numExport; i++) {

            //        std::cout<<" Index : "<<(senCntsOffset[exportProcs[i]]+Cnts[exportProcs[i]])<<std::endl;
            temp.x = mesh.x[exportLocalGids[i]];
            temp.y = mesh.y[exportLocalGids[i]];
            temp.z = mesh.z[exportLocalGids[i]];
            temp.globalID = mesh.globalIds[exportLocalGids[i]];
            sendData[senCntsOffset[exportProcs[i]] + Cnts[exportProcs[i]]] = temp;
            Cnts[exportProcs[i]]++;

        }



      /*  MPI_Alltoallv(&sendData.front(), (int *) sendCnts, (int *) senCntsOffset, MPI_NODETYPE, &recvData.front(),
                      (int *) recvCnts, (int *) recvCntsOffset, MPI_NODETYPE, comm);*/

        //par::Mpi_Alltoallv(&sendData.front(),(int *) sendCnts,(int *) senCntsOffset,&recvData.front(), (int *) recvCnts, (int *) recvCntsOffset,comm);

        par::Mpi_Alltoallv_Kway(&sendData.front(),(int *) sendCnts,(int *) senCntsOffset,&recvData.front(), (int *) recvCnts, (int *) recvCntsOffset,comm);

        auto zz_all2all_end=std::chrono::system_clock::now();//MPI_Wtime();

        ZOLTAN_ID_TYPE localSz = recvData.size();
        ZOLTAN_ID_TYPE glbalSz = 0;
        MPI_Allreduce(&localSz, &glbalSz, 1, MPI_UNSIGNED, MPI_SUM, comm);

        //if (mesh.numLocalPts > 0){

        if (mesh.globalIds != NULL) delete[] (mesh.globalIds);
        if (mesh.localIds != NULL) delete[] (mesh.localIds);
        if (mesh.x != NULL) delete[] (mesh.x);
        if (mesh.y != NULL) delete[] (mesh.y);
        if (mesh.z != NULL) delete[] (mesh.z);

        //}

        delete[](sendCnts);
        delete[](recvCnts);
        delete[](senCntsOffset);
        delete[](recvCntsOffset);

        sendData.clear();
        recvData.clear();


        mesh.numGlobalPoints = glbalSz;
        mesh.numLocalPts = localSz;

        mesh.globalIds = new ZOLTAN_ID_TYPE[localSz];
        mesh.localIds = new ZOLTAN_ID_TYPE[localSz];

        mesh.x = new ZOLTAN_ID_TYPE[localSz];
        mesh.y = new ZOLTAN_ID_TYPE[localSz];
        mesh.z = new ZOLTAN_ID_TYPE[localSz];

        for (unsigned int i = 0; i < recvData.size(); i++) {

            mesh.x[i] = recvData[i].x;
            mesh.y[i] = recvData[i].y;
            mesh.z[i] = recvData[i].z;
            mesh.globalIds[i] = recvData[i].globalID;
            mesh.localIds[i] = i;

        }

    //}



     /*   Zoltan_Set_Obj_Size_Multi_Fn(zz, user_sizes_node, &mesh);
        Zoltan_Set_Pack_Obj_Multi_Fn(zz, user_pack_node, &mesh);
        Zoltan_Set_Unpack_Obj_Multi_Fn(zz, user_unpack_node, &mesh);


        Zoltan_Migrate(zz, numImport, importGlobalGids, importLocalGids, importProcs, importToPart, numExport,
                       exportGlobalGids, exportLocalGids, exportProcs, exportToPart);*/


    //}

    auto zz_local_start=std::chrono::system_clock::now();//MPI_Wtime();

    ZOLTAN_ID_TYPE localPts=mesh.numLocalPts;

    ZOLTAN_ID_TYPE * permutedIds=new ZOLTAN_ID_TYPE[localPts];
    Zoltan_Order(zz,NUM_GID_SIZE,mesh.numLocalPts,mesh.globalIds,permutedIds);

    auto zz_local_end=std::chrono::system_clock::now();//MPI_Wtime();

    /*
    MESH_DATA  sortedMesh;
    sortedMesh.numGlobalPoints = mesh.numGlobalPoints;
    sortedMesh.numLocalPts=localPts;

    sortedMesh.localIds=new ZOLTAN_ID_TYPE[localPts];
    sortedMesh.globalIds=new ZOLTAN_ID_TYPE[localPts];
    sortedMesh.x=new ZOLTAN_ID_TYPE[localPts];
    sortedMesh.y=new ZOLTAN_ID_TYPE[localPts];
    sortedMesh.z=new ZOLTAN_ID_TYPE[localPts];

    int count=0;
    ZOLTAN_ID_TYPE minRank=permutedIds[0];
    for(unsigned int i=1;i<localPts;i++) {
        if (minRank > permutedIds[i]) minRank = permutedIds[i];
    }

    for(unsigned int i=0;i<localPts;i++)
    {

        //if(!rank)

        if( (permutedIds[i] - minRank) > localPts){
            std::cout<<"rank "<<rank<<" permute id: "<<permutedIds[i]-minRank<<std::endl;
            count++;
            continue;
        }

        sortedMesh.x[permutedIds[i] - minRank]=mesh.x[i];
        sortedMesh.y[permutedIds[i] - minRank]=mesh.y[i];
        sortedMesh.z[permutedIds[i] - minRank]=mesh.z[i];
        sortedMesh.localIds[i]=i;
        sortedMesh.globalIds[permutedIds[i] - minRank]=mesh.globalIds[i];

    }

    std::cout<<" rank "<<rank<<" Missed pts count: "<<count<<std::endl;*/



    auto zz_end=std::chrono::system_clock::now();//MPI_Wtime();

    double zz_time= (std::chrono::duration_cast<std::chrono::milliseconds>((zz_end-zz_start)).count())/(MILLISECOND_CONVERSION);
    double zz_splitter= (std::chrono::duration_cast<std::chrono::milliseconds>((zz_splitter_end-zz_splitter_start)).count())/(MILLISECOND_CONVERSION);
    double zz_all2all =(std::chrono::duration_cast<std::chrono::milliseconds>((zz_all2all_end-zz_all2all_start)).count())/(MILLISECOND_CONVERSION);
    double zz_local =(std::chrono::duration_cast<std::chrono::milliseconds>((zz_local_end-zz_local_start)).count())/(MILLISECOND_CONVERSION);


    double zz_time_g[3];
    double zz_splitter_g[3];
    double zz_all2all_g[3];
    double zz_local_g[3];

    MPI_Reduce(&zz_time,zz_time_g,1,MPI_DOUBLE,MPI_MIN,0,comm);
    MPI_Reduce(&zz_time,zz_time_g+1,1,MPI_DOUBLE,MPI_SUM,0,comm);
    MPI_Reduce(&zz_time,zz_time_g+2,1,MPI_DOUBLE,MPI_MAX,0,comm);


    MPI_Reduce(&zz_splitter,zz_splitter_g,1,MPI_DOUBLE,MPI_MIN,0,comm);
    MPI_Reduce(&zz_splitter,zz_splitter_g+1,1,MPI_DOUBLE,MPI_SUM,0,comm);
    MPI_Reduce(&zz_splitter,zz_splitter_g+2,1,MPI_DOUBLE,MPI_MAX,0,comm);

    MPI_Reduce(&zz_all2all,zz_all2all_g,1,MPI_DOUBLE,MPI_MIN,0,comm);
    MPI_Reduce(&zz_all2all,zz_all2all_g+1,1,MPI_DOUBLE,MPI_SUM,0,comm);
    MPI_Reduce(&zz_all2all,zz_all2all_g+2,1,MPI_DOUBLE,MPI_MAX,0,comm);

    MPI_Reduce(&zz_local,zz_local_g,1,MPI_DOUBLE,MPI_MIN,0,comm);
    MPI_Reduce(&zz_local,zz_local_g+1,1,MPI_DOUBLE,MPI_SUM,0,comm);
    MPI_Reduce(&zz_local,zz_local_g+2,1,MPI_DOUBLE,MPI_MAX,0,comm);


    zz_time_g[1]=zz_time_g[1]/(double)(npes);
    zz_splitter_g[1]=zz_splitter_g[1]/(double)(npes);
    zz_all2all_g[1]=zz_all2all_g[1]/(double)(npes);
    zz_local_g[1]=zz_local_g[1]/(double)(npes);


    if(!rank) {
     /*   std::cout << "zz_min\tzz_mean\tzz_max"<<std::endl;
        std::cout << zz_time_g[0]<<"\t"<<zz_time_g[1]<<"\t"<<zz_time_g[2]<<std::endl;*/


        stats.push_back(zz_splitter_g[0]);
        stats.push_back(zz_splitter_g[1]);
        stats.push_back(zz_splitter_g[2]);

        stats.push_back(zz_all2all_g[0]);
        stats.push_back(zz_all2all_g[1]);
        stats.push_back(zz_all2all_g[2]);

        stats.push_back(zz_local_g[0]);
        stats.push_back(zz_local_g[1]);
        stats.push_back(zz_local_g[2]);


        stats.push_back(zz_time_g[0]);
        stats.push_back(zz_time_g[1]);
        stats.push_back(zz_time_g[2]);




    }


    if(mesh.globalIds!=NULL) delete[] mesh.globalIds;
    if(mesh.localIds!=NULL) delete[] mesh.localIds;
    if(mesh.x!=NULL) delete[] mesh.x;
    if(mesh.y!=NULL) delete[] mesh.y;
    if(mesh.z!=NULL) delete[] mesh.z;



    Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids,
                        &importProcs, &importToPart);
    Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids,
                        &exportProcs, &exportToPart);

    Zoltan_Destroy(&zz);







}



int main(int argc, char *argv[])
{

    int rc;
    int rank,npes;

    float ver;




    /******************************************************************
    ** Initialize MPI and Zoltan
    ******************************************************************/

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    rc = Zoltan_Initialize(argc, argv, &ver);

    if (rc != ZOLTAN_OK){
        printf("Zoltan Initialization Failed...\n");
        MPI_Finalize();
        exit(0);
    }


    if (argc <5){
        if (rank == 0)
            std::cout << "Usage: mpirun -np x ./zoltan numPts dim maxDepth tol dis (0- gauss 1-logNormal 2-Uniform)" << std::endl;
        MPI_Finalize();
        exit(0);

    }

    DendroIntL grainSize = atoll(argv[1]);
    unsigned int dim= atoi(argv[2]);
    unsigned int maxDepth= atoi(argv[3]);
    double tol=atof(argv[4]);

    unsigned int distribution=atoi(argv[5]);

    std::vector<double> warm_up;

    if (!rank) std::cout << "======================================================================" << std::endl;
    if (!rank) std::cout << "     Warm Up run  (par)                                " << std::endl;
    if (!rank) std::cout << "======================================================================" << std::endl;
    ZoltanHibertLBandSort(grainSize,dim,maxDepth,tol,warm_up,distribution,MPI_COMM_WORLD);

    std::vector<double> stats_zz;

    if(!rank)
    {
        stats_zz.push_back(grainSize);
        stats_zz.push_back(npes);
        stats_zz.push_back(dim);
        stats_zz.push_back(maxDepth);
        stats_zz.push_back(tol);
    }




    if (!rank) std::cout << "======================================================================" << std::endl;
    if (!rank) std::cout << "     Zoltan Hilbert Partitioning Run  (par)                                " << std::endl;
    if (!rank) std::cout << "======================================================================" << std::endl;


    ZoltanHibertLBandSort(grainSize,dim,maxDepth,tol,stats_zz,distribution,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

#ifdef RUN_WEAK_SCALING
    int proc_group = 0;
    int min_np = 1;
    for (int i = npes; rank < i && i >= min_np; i = i >> 1) proc_group++;
    MPI_Comm comm;
    MPI_Comm_split(MPI_COMM_WORLD, proc_group, rank, &comm);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);



    //std::cout <<"Rank "<<rank<<"from "<<npes<<std::endl;


  /*  if(!rank) {
        std::cout << "numPts: " << numPts << std::endl;
        std::cout << "npes : "<<npes<<std::endl;
        std::cout << "dim: " << dim << std::endl;
        std::cout << "maxDepth: " << maxDepth << std::endl;
        std::cout << "tolerance: " << tolerance << std::endl;


    }*/


    if (!rank) {


        stats_zz.push_back(grainSize);
        stats_zz.push_back(npes);
        stats_zz.push_back(dim);
        stats_zz.push_back(maxDepth);
        stats_zz.push_back(tol);

    }

    ZoltanHibertLBandSort(grainSize,dim,maxDepth,tol,stats_zz,distribution,comm);


    //}
#endif


    int splitStatus=stats_zz.empty(); // this is same as stas_ss.empty
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    MPI_Comm statComm;

    MPI_Comm_split(MPI_COMM_WORLD,splitStatus,rank,&statComm);

    if(splitStatus==0) {


        MPI_Comm_size(statComm, &npes);
        MPI_Comm_rank(statComm, &rank);
        int localSz = stats_zz.size();
        // @milinda change to stack alloc.
        int *localSz_g = new int[npes];



        //std::cout<<"Local Size: "<<localSz<<std::endl;
        MPI_Allgather(&localSz, 1, MPI_INT, localSz_g, 1, MPI_INT, statComm);

        int *gatherCountScan = new int[npes];
        int *gatherCountOfffset = new int[npes];
        gatherCountScan[0] = localSz_g[0];
        gatherCountOfffset[0] = 0;
        for (int i = 1; i < npes; i++) {
            gatherCountScan[i] = gatherCountScan[i - 1] + localSz_g[i];
            gatherCountOfffset[i] = gatherCountOfffset[i - 1] + localSz_g[i - 1];
        }


        std::vector<double> stat_g_zz(gatherCountScan[npes - 1]);


        MPI_Gatherv(&stats_zz[0], localSz, MPI_DOUBLE, &stat_g_zz[0], localSz_g, gatherCountOfffset, MPI_DOUBLE, 0,
                    statComm);


        if(!rank)
        {


            std::cout<<"======================================================================"<<std::endl;
            std::cout<<"                        TABLE FORMAT                                  "<<std::endl;
            std::cout<<"======================================================================"<<std::endl;

            if(distribution==0)
                std::cout<<"Using Normal distribution"<<std::endl;
            else if(distribution==1)
                std::cout<<"Using Uniform distribution"<<std::endl;
            else if(distribution==2)
                std::cout<<"Using Log Normal distribution"<<std::endl;
            else
                std::cout<<"Using Normal distribution"<<std::endl;

            std::cout<<"npes\tzz_splitter_min\tzz_splitter_mean\tzz_splitter_max\tzz_all2all_min\tzz_all2all_mean\tzz_all2all_max\tzz_local_min\tzz_local_mean\tzz_local_max\tzz_min\tzz_mean\tz_max"<<std::endl;
            std::cout<<"Weak Scaling Results for : "<<stat_g_zz[0]<<std::endl;

            /*for(int i=0;i<stat_g_zz.size();i++)
                std::cout<<"stat : "<<i<<" "<<stat_g_zz[i]<<std::endl;*/



            int seqCount=0;
            for(int i=17;i<stat_g_zz.size();i=i+17) {

                if(stat_g_zz[i+1]==1) {
                    seqCount++;
                }
                if(seqCount==1)
                    continue;

                std::cout<<stat_g_zz[i+1]<<"\t"<<stat_g_zz[i+5]<<"\t"<<stat_g_zz[i+6]<<"\t"<<stat_g_zz[i+7]<<"\t"<<stat_g_zz[i+8]<<"\t"<<stat_g_zz[i+9]<<"\t"<<stat_g_zz[i+10]<<"\t"<<stat_g_zz[i+11]<<"\t"<<stat_g_zz[i+12]<<"\t"<<stat_g_zz[i+13]<<"\t"<<stat_g_zz[i+14]<<"\t"<<stat_g_zz[i+15]<<"\t"<<stat_g_zz[i+16]<<std::endl;

            }

            for(int i=0;i<17;i=i+17) {

                std::cout<<stat_g_zz[i+1]<<"\t"<<stat_g_zz[i+5]<<"\t"<<stat_g_zz[i+6]<<"\t"<<stat_g_zz[i+7]<<"\t"<<stat_g_zz[i+8]<<"\t"<<stat_g_zz[i+9]<<"\t"<<stat_g_zz[i+10]<<"\t"<<stat_g_zz[i+11]<<"\t"<<stat_g_zz[i+12]<<"\t"<<stat_g_zz[i+13]<<"\t"<<stat_g_zz[i+14]<<"\t"<<stat_g_zz[i+15]<<"\t"<<stat_g_zz[i+16]<<std::endl;

            }








        }




    }




    MPI_Finalize();



    return 0;
}

/* Application defined query functions */

 int get_number_of_objects(void *data, int *ierr)
{
    MESH_DATA *mesh= (MESH_DATA *)data;
    *ierr = ZOLTAN_OK;
    return mesh->numLocalPts;
}

 void get_object_list(void *data, int sizeGID, int sizeLID,
                            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                            int wgt_dim, float *obj_wgts, int *ierr)
{
    int i;
    MESH_DATA *mesh= (MESH_DATA *)data;
    *ierr = ZOLTAN_OK;

    /* In this example, return the IDs of our objects, but no weights.
     * Zoltan will assume equally weighted objects.
     */

    for (i=0; i<mesh->numLocalPts; i++){
        globalID[i] = mesh->globalIds[i];
        localID[i] = mesh->localIds[i];
    }
}

 int get_num_geometry(void *data, int *ierr)
{
    *ierr = ZOLTAN_OK;
    return 3;
}

void get_geometry_list(void *data, int sizeGID, int sizeLID,
                              int num_obj,
                              ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                              int num_dim, double *geom_vec, int *ierr)
{
    int i;

    MESH_DATA *mesh= (MESH_DATA *)data;

    if ( (sizeGID != 1) || (sizeLID != 1) || (num_dim != 3)){
        *ierr = ZOLTAN_FATAL;
        return;
    }

    *ierr = ZOLTAN_OK;

    for (i=0;  i < num_obj ; i++){
        geom_vec[3*i] = (double)mesh->x[i];
        geom_vec[3*i + 1] = (double)mesh->y[i];
        geom_vec[3*i + 2] = (double)mesh->z[i];
    }

    return;
}

unsigned int getParticleDataSize(int dim)
{
    unsigned int dataSz= sizeof(float)*dim + 2*sizeof(ZOLTAN_ID_TYPE);
    return dataSz;
}