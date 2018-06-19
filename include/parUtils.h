
/**
  @file parUtils.h
  @brief A set of parallel utilities.	
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  @author Hari Sundar, hsundar@gmail.com
  @author Shravan Veerapaneni, shravan@seas.upenn.edu
  @author Santi Swaroop Adavani, santis@gmail.com
  */ 

#ifndef __PAR_UTILS_H_
#define __PAR_UTILS_H_

#define KEEP_HIGH 100
#define KEEP_LOW  101

#ifdef __DEBUG__
#ifndef __DEBUG_PAR__
#define __DEBUG_PAR__
#endif
#endif

#include "mpi.h"
#include <vector>
#include "dendro.h"

#ifndef KWAY
#define KWAY 128
#endif



#ifdef PETSC_USE_LOG

#include "petscsys.h"

namespace par {
  extern int sortEvent;
  extern int concatEvent;
  extern int remdupEvent;
  extern int partwEvent;
  extern int searchEvent;
  extern int parScatterEvent;
  extern int a2avWaitEvent;
  extern int all2AllvSparseEvent;
  extern int all2AllvDenseEvent;
  extern int allGatherEvent;
  extern int reduceEvent;
  extern int sendRecvEvent;
  extern int allReduceEvent;
  extern int all2AllEvent;
  extern int allGathervEvent;
  extern int gatherEvent;
  extern int scanEvent;
  extern int bcastEvent;
  extern int splitComm2wayEvent;
  extern int splitCommEvent;
}

#define PROF_A2AV_WAIT_BEGIN PetscLogEventBegin(a2avWaitEvent,0,0,0,0);
#define PROF_A2AV_WAIT_END  PetscLogEventEnd(a2avWaitEvent,0,0,0,0); 

#define PROF_SPLIT_COMM_2WAY_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(splitComm2wayEvent,0,0,0,0);
#define PROF_SPLIT_COMM_2WAY_END         \
  PetscLogEventEnd(splitComm2wayEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_SPLIT_COMM_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(splitCommEvent,0,0,0,0);
#define PROF_SPLIT_COMM_END         \
  PetscLogEventEnd(splitCommEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_SEARCH_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(searchEvent,0,0,0,0);
#define PROF_SEARCH_END         \
  PetscLogEventEnd(searchEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PAR_CONCAT_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(concatEvent,0,0,0,0);
#define PROF_PAR_CONCAT_END         \
  PetscLogEventEnd(concatEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PAR_SCATTER_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(parScatterEvent,0,0,0,0);
#define PROF_PAR_SCATTER_END         \
  PetscLogEventEnd(parScatterEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PAR_SENDRECV_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(sendRecvEvent,0,0,0,0);
#define PROF_PAR_SENDRECV_END \
  PetscLogEventEnd(sendRecvEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PAR_BCAST_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(bcastEvent,0,0,0,0);
#define PROF_PAR_BCAST_END \
  PetscLogEventEnd(bcastEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PAR_SCAN_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(scanEvent,0,0,0,0);
#define PROF_PAR_SCAN_END \
  PetscLogEventEnd(scanEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PAR_GATHER_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(gatherEvent,0,0,0,0);
#define PROF_PAR_GATHER_END \
  PetscLogEventEnd(gatherEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PAR_REDUCE_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(reduceEvent,0,0,0,0);
#define PROF_PAR_REDUCE_END \
  PetscLogEventEnd(reduceEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PAR_ALLREDUCE_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(allReduceEvent,0,0,0,0);
#define PROF_PAR_ALLREDUCE_END \
  PetscLogEventEnd(allReduceEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PAR_ALL2ALL_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(all2AllEvent,0,0,0,0);
#define PROF_PAR_ALL2ALL_END \
  PetscLogEventEnd(all2AllEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PAR_ALLGATHER_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(allGatherEvent,0,0,0,0);
#define PROF_PAR_ALLGATHER_END \
  PetscLogEventEnd(allGatherEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PAR_ALLGATHERV_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(allGathervEvent,0,0,0,0);
#define PROF_PAR_ALLGATHERV_END \
  PetscLogEventEnd(allGathervEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PAR_ALL2ALLV_SPARSE_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(all2AllvSparseEvent,0,0,0,0);
#define PROF_PAR_ALL2ALLV_SPARSE_END         \
  PetscLogEventEnd(all2AllvSparseEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PAR_ALL2ALLV_DENSE_BEGIN \
  PetscFunctionBegin; \
PetscLogEventBegin(all2AllvDenseEvent,0,0,0,0);
#define PROF_PAR_ALL2ALLV_DENSE_END \
  PetscLogEventEnd(all2AllvDenseEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_SORT_BEGIN 	\
  PetscFunctionBegin; \
PetscLogEventBegin(sortEvent,0,0,0,0);
#define PROF_SORT_END         \
  PetscLogEventEnd(sortEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_REMDUP_BEGIN 	\
  PetscFunctionBegin; \
PetscLogEventBegin(remdupEvent,0,0,0,0);
#define PROF_REMDUP_END         \
  PetscLogEventEnd(remdupEvent,0,0,0,0); \
PetscFunctionReturn(0);

#define PROF_PARTW_BEGIN 	\
  PetscFunctionBegin; \
PetscLogEventBegin(partwEvent,0,0,0,0);
#define PROF_PARTW_END         \
  PetscLogEventEnd(partwEvent,0,0,0,0); \
PetscFunctionReturn(0);

#else

#define PROF_A2AV_WAIT_BEGIN 
#define PROF_SPLIT_COMM_2WAY_BEGIN
#define PROF_SPLIT_COMM_BEGIN
#define PROF_SEARCH_BEGIN
#define PROF_PAR_SCATTER_BEGIN
#define PROF_PAR_SENDRECV_BEGIN 
#define PROF_PAR_BCAST_BEGIN 
#define PROF_PAR_GATHER_BEGIN 
#define PROF_PAR_SCAN_BEGIN 
#define PROF_PAR_REDUCE_BEGIN
#define PROF_PAR_ALLREDUCE_BEGIN
#define PROF_PAR_ALL2ALL_BEGIN
#define PROF_PAR_ALLGATHERV_BEGIN
#define PROF_PAR_ALLGATHER_BEGIN
#define PROF_PAR_ALL2ALLV_SPARSE_BEGIN
#define PROF_PAR_ALL2ALLV_DENSE_BEGIN
#define PROF_PAR_CONCAT_BEGIN
#define PROF_SORT_BEGIN
#define PROF_REMDUP_BEGIN
#define PROF_PARTW_BEGIN

#define PROF_A2AV_WAIT_END  
#define PROF_SPLIT_COMM_2WAY_END return 1; 
#define PROF_SPLIT_COMM_END return 1; 
#define PROF_SEARCH_END return 1; 
#define PROF_PAR_SCATTER_END return 1; 
#define PROF_PAR_SENDRECV_END return 1; 
#define PROF_PAR_BCAST_END return 1; 
#define PROF_PAR_GATHER_END return 1; 
#define PROF_PAR_SCAN_END return 1; 
#define PROF_PAR_REDUCE_END return 1; 
#define PROF_PAR_ALLREDUCE_END return 1; 
#define PROF_PAR_ALL2ALL_END return 1; 
#define PROF_PAR_ALLGATHERV_END return 1; 
#define PROF_PAR_ALLGATHER_END return 1; 
#define PROF_PAR_ALL2ALLV_SPARSE_END return 1; 
#define PROF_PAR_ALL2ALLV_DENSE_END return 1; 
#define PROF_PAR_CONCAT_END return 1; 
#define PROF_SORT_END return 1; 
#define PROF_REMDUP_END return 1;
#define PROF_PARTW_END return 1; 

#endif

/**
  @namespace par
  @author Rahul Sampath
  @author Hari Sundar
  @brief Collection of Generic Parallel Functions: Sorting, Partitioning, Searching,...
  */
namespace par {

  template <typename T>
    int Mpi_Isend(T* buf, int count, int dest, int tag, MPI_Comm comm, MPI_Request* request);

  template <typename T>
    int Mpi_Issend(T* buf, int count, int dest, int tag, MPI_Comm comm, MPI_Request* request);

  template <typename T>
    int Mpi_Recv(T* buf, int count, int source, int tag, MPI_Comm comm, MPI_Status* status);

  template <typename T>
    int Mpi_Irecv(T* buf, int count, int source, int tag, MPI_Comm comm, MPI_Request* request);

  /**
   * @author Rahul S. Sampath
   */
  template <typename T>
    int Mpi_Gather( T* sendBuffer, T* recvBuffer, int count, int root, MPI_Comm comm);

  /**
   * @author Rahul S. Sampath
   */
  template <typename T, typename S>
    int Mpi_Sendrecv( T* sendBuf, int sendCount, int dest, int sendTag,
        S* recvBuf, int recvCount, int source, int recvTag,
        MPI_Comm comm, MPI_Status* status);

  /**
   * @author Rahul S. Sampath
   */
  template <typename T>
    int Mpi_Bcast( T* buffer, int count, int root, MPI_Comm comm);

  /**
   * @author Rahul S. Sampath
   */
  template <typename T>
    int Mpi_Scan( T* sendbuf, T* recvbuf, int count, MPI_Op op, MPI_Comm comm);

  /**
   * @author Rahul S. Sampath
   */
  template <typename T>
    int Mpi_Reduce( T* sendbuf, T* recvbuf, int count, MPI_Op op, int root, MPI_Comm comm);

  /**
   * @author Rahul S. Sampath
   */
  template <typename T> 
    int Mpi_Allreduce( T* sendbuf, T* recvbuf, int count, MPI_Op op, MPI_Comm comm);

  /**
   * @author Rahul S. Sampath
   */
  template <typename T>
    int Mpi_Alltoall(T* sendbuf, T* recvbuf, int count, MPI_Comm comm); 

  /**
   * @author Rahul S. Sampath
   */
  template <typename T>
    int Mpi_Allgatherv(T* sendbuf, int sendcount, T* recvbuf,
        int* recvcounts, int* displs, MPI_Comm comm);

  /**
    @author Rahul S. Sampath 
    */
  template <typename T>
    int Mpi_Allgather(T* sendbuf, T* recvbuf, int count, MPI_Comm comm);

  /**
    @author Rahul S. Sampath 
    */
  template <typename T>
    int Mpi_Alltoallv_sparse(T* sendbuf, int* sendcnts, int* sdispls, 
        T* recvbuf, int* recvcnts, int* rdispls, MPI_Comm comm);

  /**
    @author Rahul S. Sampath 
    */
  template <typename T>
    int Mpi_Alltoallv_dense(T* sendbuf, int* sendcnts, int* sdispls, 
        T* recvbuf, int* recvcnts, int* rdispls, MPI_Comm comm);


    template <typename T>
    int Mpi_Alltoallv_Kway(T* sbuff_, int* s_cnt_, int* sdisp_,
                           T* rbuff_, int* r_cnt_, int* rdisp_, MPI_Comm c);



  /**
    @brief Re-distributes a STL vector, preserving the relative ordering of the
    elements. 
    @author Rahul S. Sampath
    @param in The input vector
    @param out The output vector. Memory for this vector will be allocated within the
    function.
    @param outSz The local size of the output vector on the calling processor.
    @param comm The communicator
    @return error flag
    */
  template <typename T> 
    int scatterValues(std::vector<T> & in, std::vector<T> & out, 
        DendroIntL outSz, MPI_Comm comm );

  /**
   *
   @brief A parallel search function.
   @author Rahul Sampath
   * @param keys locally sorted unique list of keys
   * @param searchList globally sorted unique list. No processor 
   *                   must call this function with an empty list.
   * @param results maximum lower bound in searchList for the 
   *                corresponding key
   * @param comm MPI communicator
   * 
   * @return int errorcode
   */
  template <typename T>
    int maxLowerBound(const std::vector<T> & keys, const std::vector<T> & searchList,
        std::vector<T> & results, MPI_Comm comm);

  template<typename T>
    unsigned int defaultWeight(const T *a);

  /**
    @brief A parallel weighted partitioning function. In our implementation, we do not pose any 
    restriction on the input or the number of processors. This function can be used with an odd number of processors as well.
    Some processors can pass an empty vector as input. The relative ordering of the elements is preserved.
    @author Hari Sundar
    @author Rahul Sampath
    @param vec the input vector
    @param getWeight function pointer to compute the weight of each element. If you pass NULL, 
    then every element will get a weight equal to 1.    
    @param comm the communicator
    */
  template<typename T>
    int partitionW(std::vector<T>& vec,
        unsigned int (*getWeight)(const T *), MPI_Comm comm);

  /**
    @brief A parallel concatenation function. listB is appended (globally)
    to listA and the result is stored in listA. An useful  application
    of this function is when listA and listB are individually sorted (globally)
    and the smallest element in listB is greater than the largest element
    in listA and we want to create a merged list that is sorted.
    @param listA a distributed vector, the result is stored in listA
    @param listB another distributed vector that is appended to listA
    listA must not be empty on any of the calling processors. listB can be empty on some of the calling processors. listB will be cleared within the function.
    @param comm the communicator
    @author Rahul Sampath
    */
  template<typename T>
    int concatenate(std::vector<T> & listA,
        std::vector<T> & listB, MPI_Comm comm);

  /**
    @brief A parallel sample sort implementation. In our implementation, we do not pose any 
    restriction on the input or the number of processors. This function can be used with an odd number of processors as well.
    Some processors can pass an empty vector as input. If the total number of elements in the vector (globally) is fewer 
    than 10*p^2, where p is the number of processors, then we will use bitonic sort instead of sample sort to sort the vector.
    We use a paralle bitonic sort to sort the samples in the sample sort algorithm. Hence, the complexity of the algorithm
    is O(n/p log n/p) + O(p log p). Here, n is the global length of the vector and p is the number of processors.
    @author Hari Sundar
    @author Rahul Sampath
    @author Santi Swaroop Adavani
    @author Shravan Veerapaneni
    @param in the input vector
    @param out the output vector
    @param comm the communicator
    */
  template<typename T>
    int sampleSort(std::vector<T>& in, std::vector<T> & out, std::vector<double>& stats,MPI_Comm comm);

  /**
    @brief Removes duplicates in parallel. If the input is not sorted, sample sort will be called 
    within the function to sort the vector and then duplicates will be removed.
    @param nodes the input vector.
    @param isSorted pass 'true' if the vector is globally sorted.
    @param comm The communicator
    @author Rahul Sampath
    */
  /*template<typename T>
    int removeDuplicates(std::vector<T>& nodes, bool isSorted, MPI_Comm comm);*/

  /**
    @brief Splits a communication group into two, one containing processors that passed a value of 'false' 
    for the parameter 'iAmEmpty' and the another containing processors that passed a value of 'true' 
    for the parameter. Both the groups are sorted in the ascending order of their ranks in the old comm.
    @author Rahul Sampath
    @param iAmEmpty     Some flag to determine which group the calling processor will be combined into. 	
    @param orig_comm    The comm group that needs to be split.
    @param new_comm     The new comm group.
    */
  int splitComm2way(bool iAmEmpty, MPI_Comm* new_comm, MPI_Comm orig_comm);

  /**
    @brief Splits a communication group into two depending on the values in isEmptyList.
    Both the groups are sorted in the ascending order of their ranks in the old comm.
    All processors must call this function with the same 'isEmptyList' array.
    @author Rahul Sampath
    @param isEmptyList  flags (of length equal to the number of processors) to determine whether each processor is active or not. 	
    @param orig_comm    The comm group that needs to be split.
    @param new_comm     The new comm group.
    */
  int splitComm2way(const bool* isEmptyList, MPI_Comm* new_comm, MPI_Comm orig_comm);

  /*
     @author Rahul Sampath
     @brief Splits a communication group into two, processors with 
     ranks less than splittingRank form one group and the other 
     processors form the second group. Both the groups are sorted in 
     the ascending order of their ranks in the old comm.
     @param splittingRank The rank used for splitting the communicator
     @param orig_comm    The comm group that needs to be split.
     @param new_comm     The new comm group.
     */
  int splitCommUsingSplittingRank(int splittingRank, MPI_Comm* new_comm, MPI_Comm orig_comm);

  /** 
   * @brief Splits a communication group into two, the first having a power of 2
   * number of processors and the other having the remainder. The first group
   * is sorted in the ascending order of their ranks in the old comm and the second group
   * is sorted in the descending order of their ranks in the old comm
   * @author Hari Sundar
   * @param orig_comm    The comm group that needs to be split.
   * @param new_comm     The new comm group.
   */
  unsigned int splitCommBinary( MPI_Comm orig_comm, MPI_Comm* new_comm);


  /** 
   * @brief Splits a communication group into two, the first having a power of 2
   * number of processors and the other having the remainder. Both the groups
   * are sorted in the ascending order of their ranks in the old comm.
   * @author Hari Sundar
   * @param orig_comm    The comm group that needs to be split.
   * @param new_comm     The new comm group.
   */
  unsigned int splitCommBinaryNoFlip( MPI_Comm orig_comm, MPI_Comm* new_comm);

  /** 
    @author Hari Sundar
   * @brief Merges lists A, and B, retaining either the low or the High in list A.      
   * 
   * @param listA   	Input list, and where the output is stored.
   * @param listB   	Second input list.
   * @param KEEP_WHAT 	determines whether to retain the High or the low values
   * 			from A and B. One of KEEP_HIGH or KEEP_LOW.
   *
   * Merging the two lists when their sizes are not the same is a bit involved.
   * The major condition that needs to be used is that all elements that are less
   * than max(min(A), min(B)) are retained by the KEEP_LOW processor, and
   * similarly all elements that are larger larger than min(max(A), max(B)) are
   * retained by the KEEP_HIGH processor.
   *
   * The reason for this is that, on the Keep_Low side,
   *
   *   max(min(A), min(B)) > min(A) > max(A-)
   *
   * and similarly on the Keep_high side,
   *  
   *   min(max(A), max(B)) < max(A) < min(A+)
   *
   * which guarantees that the merged lists remain bitonic.   
   */
  template <typename T>
    void MergeLists( std::vector<T> &listA, std::vector<T> &listB, int KEEP_WHAT) ;

  /**
    @author Hari Sundar
    @brief The main operation in the parallel bitonic sort algorithm. This implements the compare-split operation.
   * @param which_keys is one of KEEP_HIGH or KEEP_LOW
   * @param partner    is the processor with which to Merge and Split.
   @param local_list the input vector
   @param comm the communicator
   */
  template <typename T>
    void MergeSplit( std::vector<T> &local_list, int which_keys, int partner, MPI_Comm  comm);

  /**
    @author Hari Sundar
    */
  template <typename T>
    void Par_bitonic_sort_incr( std::vector<T> &local_list, int proc_set_size, MPI_Comm  comm );

  /**
    @author Hari Sundar
    */
  template <typename T>
    void Par_bitonic_sort_decr( std::vector<T> &local_list, int proc_set_size, MPI_Comm  comm);

  /**
    @author Hari Sundar
    */
  template <typename T>
    void Par_bitonic_merge_incr( std::vector<T> &local_list, int proc_set_size, MPI_Comm  comm );

  /**
    @brief An implementation of parallel bitonic sort that expects the number of processors to be a power of 2.
    However, unlike most implementations, we do not expect the length of the vector 
    (neither locally nor globally) to be a power of 2 or even. Moreover, each processor can call 
    this with a different number of elements. However, we do expect that 'in' atleast has 1
    element on each processor.
    @param in the vector to be sorted		
    @author Hari Sundar  	
    */
  template <typename T>
    void bitonicSort_binary(std::vector<T> & in, MPI_Comm comm) ;

  /**
    @brief An implementation of parallel bitonic sort that does not expect the number of processors to
    be a power of 2. In fact, the number of processors can even be odd.
    Moreover, we do not even expect the length of the vector 
    (neither locally nor globally) to be a power of 2 or even. Moreover, each processor can call 
    this with a different number of elements. However, we do expect that 'in' atleast has 1
    element on each processor. This recursively calls the function bitonicSort_binary, followed by a
    special parallel merge.
    @param in  the vector to be sorted
    @author Hari Sundar
    @see bitonicSort_binary
    */
  template <typename T>
    void bitonicSort(std::vector<T> & in, MPI_Comm comm) ;

}//end namespace

#ifdef USE_OLD_SORT
#include "parUtils_old.tcc"
#else
#include "parUtils.tcc"
#endif

#endif















