/**
  @file parUtils.txx
  @brief Definitions of the templated functions in the par module.
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  @author Hari Sundar, hsundar@gmail.com
  @author Shravan Veerapaneni, shravan@seas.upenn.edu
  @author Santi Swaroop Adavani, santis@gmail.com
 */

#include "binUtils.h"
#include "seqUtils.h"
#include "dtypes.h"
#include <cassert>
#include <iostream>
#include <algorithm>
#include "testUtils.h"
#include "dendro.h"
#include "ompUtils.h"
#include "TreeNode.h"
#include <chrono>
#include <assert.h>

#include "treenode2vtk.h"

#ifdef __DEBUG__
#ifndef __DEBUG_PAR__
#define __DEBUG_PAR__
#endif
#endif

#define MILLISECOND_CONVERSION 1e3

namespace par
{

template <typename T>
inline int Mpi_Isend(T *buf, int count, int dest, int tag, MPI_Comm comm, MPI_Request *request)
{
  int error_code = MPI_Isend(buf, count, par::Mpi_datatype<T>::value(), dest, tag, comm, request);
  __MPI_CHECK_ERROR__(error_code,comm);
  return error_code;
}

template <typename T>
inline int Mpi_Issend(T *buf, int count, int dest, int tag, MPI_Comm comm, MPI_Request *request)
{

  int error_code = MPI_Issend(buf, count, par::Mpi_datatype<T>::value(),dest, tag, comm, request);
  __MPI_CHECK_ERROR__(error_code,comm);
  return error_code;

}

template <typename T>
inline int Mpi_Recv(T *buf, int count, int source, int tag, MPI_Comm comm, MPI_Status *status)
{

  int error_code = MPI_Recv(buf, count, par::Mpi_datatype<T>::value(),source, tag, comm, status);
  __MPI_CHECK_ERROR__(error_code,comm);
  return error_code;

}

template <typename T>
inline int Mpi_Irecv(T *buf, int count, int source, int tag, MPI_Comm comm, MPI_Request *request)
{

  int error_code = MPI_Irecv(buf, count, par::Mpi_datatype<T>::value(),source, tag, comm, request);
  __MPI_CHECK_ERROR__(error_code,comm);
  return error_code;

}

template <typename T, typename S>
inline int Mpi_Sendrecv(T *sendBuf, int sendCount, int dest, int sendTag, S *recvBuf, int recvCount, int source, int recvTag, MPI_Comm comm, MPI_Status *status)
{
  PROF_PAR_SENDRECV_BEGIN

  error_code = MPI_Sendrecv(sendBuf, sendCount, par::Mpi_datatype<T>::value(), dest, sendTag,recvBuf, recvCount, par::Mpi_datatype<S>::value(), source, recvTag, comm, status);

  PROF_PAR_SENDRECV_END
  

}

template <typename T>
inline int Mpi_Scan(T *sendbuf, T *recvbuf, int count, MPI_Op op, MPI_Comm comm)
{
  #ifdef __PROFILE_WITH_BARRIER__
    MPI_Barrier(comm);
  #endif
  PROF_PAR_SCAN_BEGIN

  error_code = MPI_Scan(sendbuf, recvbuf, count, par::Mpi_datatype<T>::value(), op, comm);

  PROF_PAR_SCAN_END

}

template <typename T>
inline int Mpi_Allreduce(T *sendbuf, T *recvbuf, int count, MPI_Op op, MPI_Comm comm)
{
  #ifdef __PROFILE_WITH_BARRIER__
    MPI_Barrier(comm);
  #endif
  PROF_PAR_ALLREDUCE_BEGIN

  error_code = MPI_Allreduce(sendbuf, recvbuf, count, par::Mpi_datatype<T>::value(), op, comm);

  PROF_PAR_ALLREDUCE_END
  

}

template <typename T>
inline int Mpi_Alltoall(T *sendbuf, T *recvbuf, int count, MPI_Comm comm)
{
  #ifdef __PROFILE_WITH_BARRIER__
    MPI_Barrier(comm);
  #endif
  PROF_PAR_ALL2ALL_BEGIN

  error_code = MPI_Alltoall(sendbuf, count, par::Mpi_datatype<T>::value(),recvbuf, count, par::Mpi_datatype<T>::value(), comm);

  PROF_PAR_ALL2ALL_END
  
}

template <typename T>
inline int Mpi_Alltoallv(T *sendbuf, int *sendcnts, int *sdispls, T *recvbuf, int *recvcnts, int *rdispls, MPI_Comm comm)
{
  #ifdef __PROFILE_WITH_BARRIER__
    MPI_Barrier(comm);
  #endif

  int error_code = MPI_Alltoallv(sendbuf, sendcnts, sdispls, par::Mpi_datatype<T>::value(), recvbuf, recvcnts, rdispls, par::Mpi_datatype<T>::value(), comm);
  
  __MPI_CHECK_ERROR__(error_code,comm);
  return error_code;
  
}

template <typename T>
inline int Mpi_Gather(T *sendBuffer, T *recvBuffer, int count, int root, MPI_Comm comm)
{
  #ifdef __PROFILE_WITH_BARRIER__
    MPI_Barrier(comm);
  #endif
  PROF_PAR_GATHER_BEGIN

  error_code = MPI_Gather(sendBuffer, count, par::Mpi_datatype<T>::value(), recvBuffer, count, par::Mpi_datatype<T>::value(), root, comm);

  PROF_PAR_GATHER_END

  

}

template <typename T>
inline int Mpi_Gatherv(T *sendBuffer, int sendcnts, T *recvBuffer, int* recvcnts, int* rdispls, int root, MPI_Comm comm)
{
  #ifdef __PROFILE_WITH_BARRIER__
    MPI_Barrier(comm);
  #endif
  PROF_PAR_GATHER_BEGIN

  error_code = MPI_Gatherv(sendBuffer, sendcnts, par::Mpi_datatype<T>::value(), recvBuffer, recvcnts, rdispls, par::Mpi_datatype<T>::value(), root, comm);

  PROF_PAR_GATHER_END

}

template <typename T>
inline int Mpi_Bcast(T *buffer, int count, int root, MPI_Comm comm)
{
  #ifdef __PROFILE_WITH_BARRIER__
    MPI_Barrier(comm);
  #endif
  PROF_PAR_BCAST_BEGIN

  error_code = MPI_Bcast(buffer, count, par::Mpi_datatype<T>::value(), root, comm);

  PROF_PAR_BCAST_END
  

}

template <typename T>
inline int Mpi_Reduce(T *sendbuf, T *recvbuf, int count, MPI_Op op, int root, MPI_Comm comm)
{
  #ifdef __PROFILE_WITH_BARRIER__
    MPI_Barrier(comm);
  #endif
  PROF_PAR_REDUCE_BEGIN

  error_code = MPI_Reduce(sendbuf, recvbuf, count, par::Mpi_datatype<T>::value(), op, root, comm);

  PROF_PAR_REDUCE_END
  
}

template <typename T>
int Mpi_Allgatherv(T *sendBuf, int sendCount, T *recvBuf, int *recvCounts, int *displs, MPI_Comm comm)
{
  #ifdef __PROFILE_WITH_BARRIER__
    MPI_Barrier(comm);
  #endif
  PROF_PAR_ALLGATHERV_BEGIN
  
  #ifdef __USE_A2A_FOR_MPI_ALLGATHER__

  int maxSendCount;
  int npes, rank;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &rank);

  par::Mpi_Allreduce<int>(&sendCount, &maxSendCount, 1, MPI_MAX, comm);

  T *dummySendBuf = new T[maxSendCount * npes];
  assert(dummySendBuf);

  #pragma omp parallel for
  for (int i = 0; i < npes; i++)
  {
    for (int j = 0; j < sendCount; j++)
    {
      dummySendBuf[(i * maxSendCount) + j] = sendBuf[j];
    }
  }

  T *dummyRecvBuf = new T[maxSendCount * npes];
  assert(dummyRecvBuf);

  par::Mpi_Alltoall<T>(dummySendBuf, dummyRecvBuf, maxSendCount, comm);

  #pragma omp parallel for
  for (int i = 0; i < npes; i++)
  {
    for (int j = 0; j < recvCounts[i]; j++)
    {
      recvBuf[displs[i] + j] = dummyRecvBuf[(i * maxSendCount) + j];
    }
  }

  delete[] dummySendBuf;
  delete[] dummyRecvBuf;

  #else

  error_code = MPI_Allgatherv(sendBuf, sendCount, par::Mpi_datatype<T>::value(),recvBuf, recvCounts, displs, par::Mpi_datatype<T>::value(), comm);

  #endif
  PROF_PAR_ALLGATHERV_END
  

}

template <typename T>
int Mpi_Allgather(T *sendBuf, T *recvBuf, int count, MPI_Comm comm)
{
  #ifdef __PROFILE_WITH_BARRIER__
    MPI_Barrier(comm);
  #endif
  PROF_PAR_ALLGATHER_BEGIN

  #ifdef __USE_A2A_FOR_MPI_ALLGATHER__

  int npes;
  MPI_Comm_size(comm, &npes);
  T *dummySendBuf = new T[count * npes];
  assert(dummySendBuf);
  #pragma omp parallel for
  for (int i = 0; i < npes; i++)
  {
    for (int j = 0; j < count; j++)
    {
      dummySendBuf[(i * count) + j] = sendBuf[j];
    }
  }
  par::Mpi_Alltoall<T>(dummySendBuf, recvBuf, count, comm);
  delete[] dummySendBuf;

  #else

  error_code = MPI_Allgather(sendBuf, count, par::Mpi_datatype<T>::value(), recvBuf, count, par::Mpi_datatype<T>::value(), comm);

  #endif

  PROF_PAR_ALLGATHER_END
  

}

template <typename T>
int Mpi_Alltoallv_sparse(T *sendbuf, int *sendcnts, int *sdispls, T *recvbuf, int *recvcnts, int *rdispls, MPI_Comm comm)
{
  #ifdef __PROFILE_WITH_BARRIER__
    MPI_Barrier(comm);
  #endif
  PROF_PAR_ALL2ALLV_SPARSE_BEGIN

  #ifndef ALLTOALLV_FIX
    error_code = Mpi_Alltoallv(sendbuf, sendcnts, sdispls, recvbuf, recvcnts, rdispls, comm);
  #else

  int npes, rank;
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &rank);

  int commCnt = 0;

  #pragma omp parallel for reduction(+ \
                                    : commCnt)
    for (int i = 0; i < rank; i++)
    {
      if (sendcnts[i] > 0)
      {
        commCnt++;
      }
      if (recvcnts[i] > 0)
      {
        commCnt++;
      }
    }

  #pragma omp parallel for reduction(+ \
                                    : commCnt)
    for (int i = (rank + 1); i < npes; i++)
    {
      if (sendcnts[i] > 0)
      {
        commCnt++;
      }
      if (recvcnts[i] > 0)
      {
        commCnt++;
      }
    }

    MPI_Request *requests = new MPI_Request[commCnt];
    assert(requests);

    MPI_Status *statuses = new MPI_Status[commCnt];
    assert(statuses);

    commCnt = 0;

    //First place all recv requests. Do not recv from self.
    for (int i = 0; i < rank; i++)
    {
      if (recvcnts[i] > 0)
      {
        par::Mpi_Irecv<T>(&(recvbuf[rdispls[i]]), recvcnts[i], i, 1,
                          comm, &(requests[commCnt]));
        commCnt++;
      }
    }

    for (int i = (rank + 1); i < npes; i++)
    {
      if (recvcnts[i] > 0)
      {
        par::Mpi_Irecv<T>(&(recvbuf[rdispls[i]]), recvcnts[i], i, 1,
                          comm, &(requests[commCnt]));
        commCnt++;
      }
    }

    //Next send the messages. Do not send to self.
    for (int i = 0; i < rank; i++)
    {
      if (sendcnts[i] > 0)
      {
        par::Mpi_Issend<T>(&(sendbuf[sdispls[i]]), sendcnts[i], i, 1,
                          comm, &(requests[commCnt]));
        commCnt++;
      }
    }

    for (int i = (rank + 1); i < npes; i++)
    {
      if (sendcnts[i] > 0)
      {
        par::Mpi_Issend<T>(&(sendbuf[sdispls[i]]), sendcnts[i],
                          i, 1, comm, &(requests[commCnt]));
        commCnt++;
      }
    }

    //Now copy local portion.
  #ifdef __DEBUG_PAR__
    assert(sendcnts[rank] == recvcnts[rank]);
  #endif

  #pragma omp parallel for
    for (int i = 0; i < sendcnts[rank]; i++)
    {
      recvbuf[rdispls[rank] + i] = sendbuf[sdispls[rank] + i];
    }

    
    error_code = MPI_Waitall(commCnt, requests, statuses);
    // 27/05/20(Milinda) : A2A_sparse was crashing on impi19 in frontera. I put a barrier here since all2all v should be a blocking call, otherwise repeated calls might cause problems. 
    //MPI_Barrier(comm); 

    delete[] requests;
    delete[] statuses;
  #endif
  PROF_PAR_ALL2ALLV_SPARSE_END
  

}

template <typename T>
int Mpi_Alltoallv_dense(T *sendbuf, int *sendcnts, int *sdispls, T *recvbuf, int *recvcnts, int *rdispls, MPI_Comm comm)
{
  #ifdef __PROFILE_WITH_BARRIER__
    MPI_Barrier(comm);
  #endif
  PROF_PAR_ALL2ALLV_DENSE_BEGIN

  #ifndef ALLTOALLV_FIX
    error_code = Mpi_Alltoallv(sendbuf, sendcnts, sdispls, recvbuf, recvcnts, rdispls, comm);
  #else
  int npes, rank;
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &rank);

  //Processors may send a lot of information to themselves and a lesser
  //amount to others. If so, we don't want to waste communication by
  //including the local copy size in the max message size.
  int maxNumElemSend = 0;
  for (int i = 0; i < rank; i++)
  {
    if (sendcnts[i] > maxNumElemSend)
    {
      maxNumElemSend = sendcnts[i];
    }
  }

  for (int i = (rank + 1); i < npes; i++)
  {
    if (sendcnts[i] > maxNumElemSend)
    {
      maxNumElemSend = sendcnts[i];
    }
  }

  int allToAllCount;
  par::Mpi_Allreduce<int>(&maxNumElemSend, &allToAllCount, 1, MPI_MAX, comm);

  T *tmpSendBuf = new T[allToAllCount * npes];
  assert(tmpSendBuf);

  T *tmpRecvBuf = new T[allToAllCount * npes];
  assert(tmpRecvBuf);

  for (int i = 0; i < rank; i++)
  {
    for (int j = 0; j < sendcnts[i]; j++)
    {
      tmpSendBuf[(allToAllCount * i) + j] = sendbuf[sdispls[i] + j];
    }
  }

  for (int i = (rank + 1); i < npes; i++)
  {
    for (int j = 0; j < sendcnts[i]; j++)
    {
      tmpSendBuf[(allToAllCount * i) + j] = sendbuf[sdispls[i] + j];
    }
  }

  error_code = par::Mpi_Alltoall<T>(tmpSendBuf, tmpRecvBuf, allToAllCount, comm);

  for (int i = 0; i < rank; i++)
  {
    for (int j = 0; j < recvcnts[i]; j++)
    {
      recvbuf[rdispls[i] + j] = tmpRecvBuf[(allToAllCount * i) + j];
    }
  }

  //Now copy local portion.
  #ifdef __DEBUG_PAR__
  assert(sendcnts[rank] == recvcnts[rank]);
  #endif

  for (int j = 0; j < recvcnts[rank]; j++)
  {
    recvbuf[rdispls[rank] + j] = sendbuf[sdispls[rank] + j];
  }

  for (int i = (rank + 1); i < npes; i++)
  {
    for (int j = 0; j < recvcnts[i]; j++)
    {
      recvbuf[rdispls[i] + j] = tmpRecvBuf[(allToAllCount * i) + j];
    }
  }

  delete[] tmpSendBuf;
  delete[] tmpRecvBuf;
  
  #endif

  PROF_PAR_ALL2ALLV_DENSE_END
}

template <typename T>
int Mpi_Alltoallv_Kway(T *sbuff_, int *s_cnt_, int *sdisp_, T *rbuff_, int *r_cnt_, int *rdisp_, MPI_Comm comm)
{

  #ifdef __PROFILE_WITH_BARRIER__
    MPI_Barrier(comm);
  #endif
  PROF_PAR_ALL2ALLV_DENSE_BEGIN
  #ifndef ALLTOALLV_FIX 
    error_code = Mpi_Alltoallv(sbuff_, s_cnt_, sdisp_, rbuff_, r_cnt_, rdisp_, comm);
  #else

  int kway = KWAY;
  int np, pid;
  MPI_Comm_size(comm, &np);
  MPI_Comm_rank(comm, &pid);
  //std::cout<<" Kway: "<<kway<<std::endl;
  int range[2] = {0, np};
  int split_id, partner;

  std::vector<int> s_cnt(np);
  #pragma omp parallel for
  for (int i = 0; i < np; i++)
  {
    s_cnt[i] = s_cnt_[i] * sizeof(T) + 2 * sizeof(int);
  }
  std::vector<int> sdisp(np);
  sdisp[0] = 0;
  omp_par::scan(&s_cnt[0], &sdisp[0], np);

  char *sbuff = new char[sdisp[np - 1] + s_cnt[np - 1]];
  #pragma omp parallel for
  for (int i = 0; i < np; i++)
  {
    ((int *)&sbuff[sdisp[i]])[0] = s_cnt[i];
    ((int *)&sbuff[sdisp[i]])[1] = pid;
    memcpy(&sbuff[sdisp[i]] + 2 * sizeof(int), &sbuff_[sdisp_[i]], s_cnt[i] - 2 * sizeof(int));
  }

  //int t_indx=0;
  int iter_cnt = 0;
  while (range[1] - range[0] > 1)
  {
    iter_cnt++;
    if (kway > range[1] - range[0])
      kway = range[1] - range[0];

    std::vector<int> new_range(kway + 1);
    for (int i = 0; i <= kway; i++)
      new_range[i] = (range[0] * (kway - i) + range[1] * i) / kway;
    int p_class = (std::upper_bound(&new_range[0], &new_range[kway], pid) - &new_range[0] - 1);
    int new_np = new_range[p_class + 1] - new_range[p_class];
    int new_pid = pid - new_range[p_class];

    //Communication.
    {
      std::vector<int> r_cnt(new_np * kway, 0);
      std::vector<int> r_cnt_ext(new_np * kway, 0);
      //Exchange send sizes.
      for (int i = 0; i < kway; i++)
      {
        MPI_Status status;
        int cmp_np = new_range[i + 1] - new_range[i];
        int partner = (new_pid < cmp_np ? new_range[i] + new_pid : new_range[i + 1] - 1);
        assert((new_pid < cmp_np ? true : new_range[i] + new_pid == new_range[i + 1])); //Remove this.
        MPI_Sendrecv(&s_cnt[new_range[i] - new_range[0]], cmp_np, MPI_INT, partner, 0,
                     &r_cnt[new_np * i], new_np, MPI_INT, partner, 0, comm, &status);

        //Handle extra communication.
        if (new_pid == new_np - 1 && cmp_np > new_np)
        {
          int partner = new_range[i + 1] - 1;
          std::vector<int> s_cnt_ext(cmp_np, 0);
          MPI_Sendrecv(&s_cnt_ext[0], cmp_np, MPI_INT, partner, 0,
                       &r_cnt_ext[new_np * i], new_np, MPI_INT, partner, 0, comm, &status);
        }
      }

      //Allocate receive buffer.
      std::vector<int> rdisp(new_np * kway, 0);
      std::vector<int> rdisp_ext(new_np * kway, 0);
      int rbuff_size, rbuff_size_ext;
      char *rbuff, *rbuff_ext;
      {
        omp_par::scan(&r_cnt[0], &rdisp[0], new_np * kway);
        omp_par::scan(&r_cnt_ext[0], &rdisp_ext[0], new_np * kway);
        rbuff_size = rdisp[new_np * kway - 1] + r_cnt[new_np * kway - 1];
        rbuff_size_ext = rdisp_ext[new_np * kway - 1] + r_cnt_ext[new_np * kway - 1];
        rbuff = new char[rbuff_size];
        rbuff_ext = new char[rbuff_size_ext];
      }

      //Sendrecv data.
      //*
      int my_block = kway;
      while (pid < new_range[my_block])
        my_block--;
      //      MPI_Barrier(comm);
      for (int i_ = 0; i_ <= kway / 2; i_++)
      {
        int i1 = (my_block + i_) % kway;
        int i2 = (my_block + kway - i_) % kway;

        for (int j = 0; j < (i_ == 0 || i_ == kway / 2 ? 1 : 2); j++)
        {
          int i = (i_ == 0 ? i1 : ((j + my_block / i_) % 2 ? i1 : i2));
          MPI_Status status;
          int cmp_np = new_range[i + 1] - new_range[i];
          int partner = (new_pid < cmp_np ? new_range[i] + new_pid : new_range[i + 1] - 1);

          int send_dsp = sdisp[new_range[i] - new_range[0]];
          int send_dsp_last = sdisp[new_range[i + 1] - new_range[0] - 1];
          int send_cnt = s_cnt[new_range[i + 1] - new_range[0] - 1] + send_dsp_last - send_dsp;

          
          MPI_Sendrecv(&sbuff[send_dsp], send_cnt, MPI_BYTE, partner, 0,
                       &rbuff[rdisp[new_np * i]], r_cnt[new_np * (i + 1) - 1] + rdisp[new_np * (i + 1) - 1] - rdisp[new_np * i], MPI_BYTE, partner, 0, comm, &status);
          
          //Handle extra communication.
          if (pid == new_np - 1 && cmp_np > new_np)
          {
            int partner = new_range[i + 1] - 1;
            std::vector<int> s_cnt_ext(cmp_np, 0);
            MPI_Sendrecv(NULL, 0, MPI_BYTE, partner, 0,
                         &rbuff[rdisp_ext[new_np * i]], r_cnt_ext[new_np * (i + 1) - 1] + rdisp_ext[new_np * (i + 1) - 1] - rdisp_ext[new_np * i], MPI_BYTE, partner, 0, comm, &status);
          }
        }
      }
      

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      //Rearrange received data.
      {
        if (sbuff != NULL)
          delete[] sbuff;
        sbuff = new char[rbuff_size + rbuff_size_ext];

        std::vector<int> cnt_new(2 * new_np * kway, 0);
        std::vector<int> disp_new(2 * new_np * kway, 0);
        for (int i = 0; i < new_np; i++)
          for (int j = 0; j < kway; j++)
          {
            cnt_new[(i * 2) * kway + j] = r_cnt[j * new_np + i];
            cnt_new[(i * 2 + 1) * kway + j] = r_cnt_ext[j * new_np + i];
          }
        omp_par::scan(&cnt_new[0], &disp_new[0], 2 * new_np * kway);

        #pragma omp parallel for
        for (int i = 0; i < new_np; i++)
          for (int j = 0; j < kway; j++)
          {
            memcpy(&sbuff[disp_new[(i * 2) * kway + j]], &rbuff[rdisp[j * new_np + i]], r_cnt[j * new_np + i]);
            memcpy(&sbuff[disp_new[(i * 2 + 1) * kway + j]], &rbuff_ext[rdisp_ext[j * new_np + i]], r_cnt_ext[j * new_np + i]);
          }

        //Free memory.
        if (rbuff != NULL)
          delete[] rbuff;
        if (rbuff_ext != NULL)
          delete[] rbuff_ext;

        s_cnt.clear();
        s_cnt.resize(new_np, 0);
        sdisp.resize(new_np);
        for (int i = 0; i < new_np; i++)
        {
          for (int j = 0; j < 2 * kway; j++)
            s_cnt[i] += cnt_new[i * 2 * kway + j];
          sdisp[i] = disp_new[i * 2 * kway];
        }
      }

      
    }

    range[0] = new_range[p_class];
    range[1] = new_range[p_class + 1];
   
  }

  //Copy data to rbuff_.
  std::vector<char *> buff_ptr(np);
  char *tmp_ptr = sbuff;
  for (int i = 0; i < np; i++)
  {
    int &blk_size = ((int *)tmp_ptr)[0];
    buff_ptr[i] = tmp_ptr;
    tmp_ptr += blk_size;
  }
  
  #pragma omp parallel for
  for (int i = 0; i < np; i++)
  {
    int &blk_size = ((int *)buff_ptr[i])[0];
    int &src_pid = ((int *)buff_ptr[i])[1];
    assert(blk_size - 2 * sizeof(int) <= r_cnt_[src_pid] * sizeof(T));
    memcpy(&rbuff_[rdisp_[src_pid]], buff_ptr[i] + 2 * sizeof(int), blk_size - 2 * sizeof(int));
  }
  
  //Free memory.
  if (sbuff != NULL)
    delete[] sbuff;
  #endif

  PROF_PAR_ALL2ALLV_DENSE_END
}

template <typename T>
unsigned int defaultWeight(const T *a)
{
  return 1;
}

template <typename T>
int scatterValues(std::vector<T> &in, std::vector<T> &out, DendroIntL outSz, MPI_Comm comm)
{
  #ifdef __PROFILE_WITH_BARRIER__
    MPI_Barrier(comm);
  #endif
  PROF_PAR_SCATTER_BEGIN

  int rank, npes;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &rank);

  #if 0
    // Milinda(26/05/2020) : simple scatter values code that wrote, because scatterValues was crashing on Frontera
    // for larger number of cores. May be it was due to not setting I_MPI_ALLREDUCE=1. There are some all_reduce (intel mpi) errrors in 
    // Frontera.  
    std::vector<int> sc;
    std::vector<int> sfst;
    std::vector<int> rc;
    std::vector<int> rfst;

    sc.resize(npes,0);
    sfst.resize(npes,0);
    rc.resize(npes,0);
    rfst.resize(npes,0);

    int npes_act = in.size()/DENDRO_DEFAULT_GRAIN_SZ;
    
    if (npes_act > npes)
      npes_act = npes;

    if(npes_act==0)
      npes_act=2;

    const DendroIntL inSz = in.size();

    for(unsigned int i=0; i < npes_act; i++)
      sc[i] = ((i+1)*inSz)/npes_act - ((i)*inSz)/npes_act;

    par::Mpi_Alltoall(sc.data(),rc.data(),1,comm);

    sfst[0] = 0;
    rfst[0] = 0;

    omp_par::scan(sc.data(),sfst.data(),npes);
    omp_par::scan(rc.data(),rfst.data(),npes);

    out.resize(rfst[npes-1] + rc[npes-1]);
    par::Mpi_Alltoallv(in.data(), sc.data(), sfst.data(),out.data(),rc.data(),rfst.data(),comm);
  #endif
  
  #if 1

    MPI_Request request;
    MPI_Status status;

    DendroIntL inSz = in.size();
    out.resize(outSz);

    DendroIntL off1 = 0, off2 = 0;
    DendroIntL *scnIn = NULL;
    if (inSz)
    {
      scnIn = new DendroIntL[inSz];
      assert(scnIn);
    }

    // perform a local scan first ...
    DendroIntL zero = 0;
    if (inSz)
    {
      scnIn[0] = 1;
      for (DendroIntL i = 1; i < inSz; i++)
      {
        scnIn[i] = scnIn[i - 1] + 1;
      } //end for
      // now scan with the final members of
      par::Mpi_Scan<DendroIntL>(scnIn + inSz - 1, &off1, 1, MPI_SUM, comm);
    }
    else
    {
      par::Mpi_Scan<DendroIntL>(&zero, &off1, 1, MPI_SUM, comm);
    }

    // communicate the offsets ...
    if (rank < (npes - 1))
    {
      par::Mpi_Issend<DendroIntL>(&off1, 1, (rank + 1), 0, comm, &request);
    }
    if (rank)
    {
      par::Mpi_Recv<DendroIntL>(&off2, 1, (rank - 1), 0, comm, &status);
    }
    else
    {
      off2 = 0;
    }

    // add offset to local array
    for (DendroIntL i = 0; i < inSz; i++)
    {
      scnIn[i] = scnIn[i] + off2; // This has the global scan results now ...
    }                             //end for

    //Gather Scan of outCnts
    DendroIntL *outCnts;
    outCnts = new DendroIntL[npes];
    assert(outCnts);

    if (rank < (npes - 1))
    {
      MPI_Status statusWait;
      MPI_Wait(&request, &statusWait);
    }

    if (outSz)
    {
      par::Mpi_Scan<DendroIntL>(&outSz, &off1, 1, MPI_SUM, comm);
    }
    else
    {
      par::Mpi_Scan<DendroIntL>(&zero, &off1, 1, MPI_SUM, comm);
    }

    par::Mpi_Allgather<DendroIntL>(&off1, outCnts, 1, comm);

    int *sendSz = new int[npes];
    assert(sendSz);

    int *recvSz = new int[npes];
    assert(recvSz);

    int *sendOff = new int[npes];
    assert(sendOff);

    int *recvOff = new int[npes];
    assert(recvOff);

    // compute the partition offsets and sizes so that All2Allv can be performed.
    // initialize ...
    for (int i = 0; i < npes; i++)
    {
      sendSz[i] = 0;
    }

    //The Heart of the algorithm....
    //scnIn and outCnts are both sorted
    DendroIntL inCnt = 0;
    int pCnt = 0;
    while ((inCnt < inSz) && (pCnt < npes))
    {
      if (scnIn[inCnt] <= outCnts[pCnt])
      {
        sendSz[pCnt]++;
        inCnt++;
      }
      else
      {
        pCnt++;
      }
    }

    // communicate with other procs how many you shall be sending and get how
    // many to recieve from whom.
    par::Mpi_Alltoall<int>(sendSz, recvSz, 1, comm);

    int nn = 0; // new value of nlSize, ie the local nodes.
    for (int i = 0; i < npes; i++)
    {
      nn += recvSz[i];
    }

    // compute offsets ...
    sendOff[0] = 0;
    recvOff[0] = 0;
    for (int i = 1; i < npes; i++)
    {
      sendOff[i] = sendOff[i - 1] + sendSz[i - 1];
      recvOff[i] = recvOff[i - 1] + recvSz[i - 1];
    }

    assert(static_cast<unsigned int>(nn) == outSz);
    // perform All2All  ...
    T *inPtr = NULL;
    T *outPtr = NULL;
    if (!in.empty())
    {
      inPtr = &(*(in.begin()));
    }
    if (!out.empty())
    {
      outPtr = &(*(out.begin()));
    }
    par::Mpi_Alltoallv_sparse<T>(inPtr, sendSz, sendOff,
                                outPtr, recvSz, recvOff, comm);

    // clean up...
    if (scnIn)
    {
      delete[] scnIn;
      scnIn = NULL;
    }

    delete[] outCnts;
    outCnts = NULL;

    delete[] sendSz;
    sendSz = NULL;

    delete[] sendOff;
    sendOff = NULL;

    delete[] recvSz;
    recvSz = NULL;

    delete[] recvOff;
    recvOff = NULL;
  #endif

  PROF_PAR_SCATTER_END


  

}

template <typename T>
int concatenate(std::vector<T> &listA, std::vector<T> &listB, MPI_Comm comm)
{
  #ifdef __PROFILE_WITH_BARRIER__
    MPI_Barrier(comm);
  #endif
  PROF_PAR_CONCAT_BEGIN

  int rank;
  int npes;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &npes);

  assert(!(listA.empty()));

  //1. First perform Allreduce to get total listA size
  //and total listB size;

  DendroIntL locAsz_locBsz[2];
  DendroIntL globAsz_globBsz[2];

  locAsz_locBsz[0] = listA.size();
  locAsz_locBsz[1] = listB.size();
  globAsz_globBsz[0] = 0;
  globAsz_globBsz[1] = 0;

  par::Mpi_Allreduce<DendroIntL>(locAsz_locBsz, globAsz_globBsz, 2, MPI_SUM, comm);

  //2. Re-distribute A and B independently so that
  //B is distributed only on the high rank processors
  //and A is distribute only on the low rank processors.

  DendroIntL avgTotalSize = ((globAsz_globBsz[0] + globAsz_globBsz[1]) / npes);

  //since listA is not empty on any of the active procs,
  //globASz > npes so avgTotalSize >= 1

  DendroIntL remTotalSize = ((globAsz_globBsz[0] + globAsz_globBsz[1]) % npes);

  int numSmallProcs = (npes - remTotalSize);

  //In the final merged list, there will be exactly remTotalSize number
  //of processors each having (avgTotalSize + 1) elements and there will
  //be exactly numSmallProcs number of processors each having
  //avgTotalSize elements.
  //Also, len(A) + len(B) = (numSmallProcs*avg) + (remTotalSize*(avg+1))

  std::vector<T> tmpA;
  std::vector<T> tmpB;

  int numAhighProcs;
  int numAlowProcs;
  int numBothProcs;
  int numBhighProcs;
  int numBlowProcs;
  DendroIntL aSizeForBoth;
  DendroIntL bSizeForBoth;

  if (globAsz_globBsz[1] <= (numSmallProcs * avgTotalSize))
  {
    numBhighProcs = 0;
    numBlowProcs = ((globAsz_globBsz[1]) / avgTotalSize);
    bSizeForBoth = ((globAsz_globBsz[1]) % avgTotalSize);

    assert(numBlowProcs <= numSmallProcs);

    //remBsize is < avgTotalSize. So it will fit on one proc.
    if (bSizeForBoth)
    {
      numBothProcs = 1;
      if (numBlowProcs < numSmallProcs)
      {
        //We don't know if remTotalSize is 0 or not.
        //So, let the common proc be a low proc.
        aSizeForBoth = (avgTotalSize - bSizeForBoth);
        numAhighProcs = remTotalSize;
        numAlowProcs = (numSmallProcs - (1 + numBlowProcs));
      }
      else
      {
        //No more room for small procs. The common has to be a high proc.
        aSizeForBoth = ((avgTotalSize + 1) - bSizeForBoth);
        numAhighProcs = (remTotalSize - 1);
        numAlowProcs = 0;
      }
    }
    else
    {
      numBothProcs = 0;
      aSizeForBoth = 0;
      numAhighProcs = remTotalSize;
      numAlowProcs = (numSmallProcs - numBlowProcs);
    }
  }
  else
  {
    //Some B procs will have (avgTotalSize+1) elements
    DendroIntL numBusingAvgPlus1 = ((globAsz_globBsz[1]) / (avgTotalSize + 1));
    DendroIntL remBusingAvgPlus1 = ((globAsz_globBsz[1]) % (avgTotalSize + 1));
    if (numBusingAvgPlus1 <= remTotalSize)
    {
      //Each block can use (avg+1) elements each, since there will be some
      //remaining for A
      numBhighProcs = numBusingAvgPlus1;
      numBlowProcs = 0;
      bSizeForBoth = remBusingAvgPlus1;
      if (bSizeForBoth)
      {
        numBothProcs = 1;
        if (numBhighProcs < remTotalSize)
        {
          //We don't know if numSmallProcs is 0 or not.
          //So, let the common proc be a high proc
          aSizeForBoth = ((avgTotalSize + 1) - bSizeForBoth);
          numAhighProcs = (remTotalSize - (numBhighProcs + 1));
          numAlowProcs = numSmallProcs;
        }
        else
        {
          //No more room for high procs. The common has to be a low proc.
          aSizeForBoth = (avgTotalSize - bSizeForBoth);
          numAhighProcs = 0;
          numAlowProcs = (numSmallProcs - 1);
        }
      }
      else
      {
        numBothProcs = 0;
        aSizeForBoth = 0;
        numAhighProcs = (remTotalSize - numBhighProcs);
        numAlowProcs = numSmallProcs;
      }
    }
    else
    {
      //Since numBusingAvgPlus1 > remTotalSize*(avg+1)
      //=> len(B) > remTotalSize*(avg+1)
      //=> len(A) < numSmallProcs*avg
      //This is identical to the first case (except for
      //the equality), with A and B swapped.

      assert(globAsz_globBsz[0] < (numSmallProcs * avgTotalSize));

      numAhighProcs = 0;
      numAlowProcs = ((globAsz_globBsz[0]) / avgTotalSize);
      aSizeForBoth = ((globAsz_globBsz[0]) % avgTotalSize);

      assert(numAlowProcs < numSmallProcs);

      //remAsize is < avgTotalSize. So it will fit on one proc.
      if (aSizeForBoth)
      {
        numBothProcs = 1;
        //We don't know if remTotalSize is 0 or not.
        //So, let the common proc be a low proc.
        bSizeForBoth = (avgTotalSize - aSizeForBoth);
        numBhighProcs = remTotalSize;
        numBlowProcs = (numSmallProcs - (1 + numAlowProcs));
      }
      else
      {
        numBothProcs = 0;
        bSizeForBoth = 0;
        numBhighProcs = remTotalSize;
        numBlowProcs = (numSmallProcs - numAlowProcs);
      }
    }
  }

  assert((numAhighProcs + numAlowProcs + numBothProcs + numBhighProcs + numBlowProcs) == npes);

  assert((aSizeForBoth + bSizeForBoth) <= (avgTotalSize + 1));

  if (numBothProcs)
  {
    assert((aSizeForBoth + bSizeForBoth) >= avgTotalSize);
  }
  else
  {
    assert(aSizeForBoth == 0);
    assert(bSizeForBoth == 0);
  }

  if ((aSizeForBoth + bSizeForBoth) == (avgTotalSize + 1))
  {
    assert((numAhighProcs + numBothProcs + numBhighProcs) == remTotalSize);
    assert((numAlowProcs + numBlowProcs) == numSmallProcs);
  }
  else
  {
    assert((numAhighProcs + numBhighProcs) == remTotalSize);
    assert((numAlowProcs + numBothProcs + numBlowProcs) == numSmallProcs);
  }

  //The partition is as follow:
  //1. numAhighProcs with (avg+1) elements each exclusively from A,
  //2. numAlowProcs with avg elements each exclusively from A
  //3. numBothProcs with aSizeForBoth elements from A and
  // bSizeForBoth elements from B
  //4. numBhighProcs with (avg+1) elements each exclusively from B.
  //5. numBlowProcs with avg elements each exclusively from B.

  if (rank < numAhighProcs)
  {
    par::scatterValues<T>(listA, tmpA, (avgTotalSize + 1), comm);
    par::scatterValues<T>(listB, tmpB, 0, comm);
  }
  else if (rank < (numAhighProcs + numAlowProcs))
  {
    par::scatterValues<T>(listA, tmpA, avgTotalSize, comm);
    par::scatterValues<T>(listB, tmpB, 0, comm);
  }
  else if (rank < (numAhighProcs + numAlowProcs + numBothProcs))
  {
    par::scatterValues<T>(listA, tmpA, aSizeForBoth, comm);
    par::scatterValues<T>(listB, tmpB, bSizeForBoth, comm);
  }
  else if (rank <
           (numAhighProcs + numAlowProcs + numBothProcs + numBhighProcs))
  {
    par::scatterValues<T>(listA, tmpA, 0, comm);
    par::scatterValues<T>(listB, tmpB, (avgTotalSize + 1), comm);
  }
  else
  {
    par::scatterValues<T>(listA, tmpA, 0, comm);
    par::scatterValues<T>(listB, tmpB, avgTotalSize, comm);
  }

  listA = tmpA;
  listB = tmpB;
  tmpA.clear();
  tmpB.clear();

  //3. Finally do a simple concatenation A = A + B. If the previous step
  //was performed correctly, there will be atmost 1 processor, which has both
  //non-empty A and non-empty B. On other processors one of the two lists
  //will be empty
  if (listA.empty())
  {
    listA = listB;
  }
  else
  {
    if (!(listB.empty()))
    {
      listA.insert(listA.end(), listB.begin(), listB.end());
    }
  }

  listB.clear();

  PROF_PAR_CONCAT_END
}

template <typename T>
int maxLowerBound(const std::vector<T> &keys, const std::vector<T> &searchList, std::vector<T> &results, MPI_Comm comm)
{
  PROF_SEARCH_BEGIN

  int rank, npes;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &rank);

  // allocate memory for the mins array
  std::vector<T> mins(npes);
  assert(!searchList.empty());

  T *searchListPtr = NULL;
  T *minsPtr = NULL;
  if (!searchList.empty())
  {
    searchListPtr = &(*(searchList.begin()));
  }
  if (!mins.empty())
  {
    minsPtr = &(*(mins.begin()));
  }
  par::Mpi_Allgather<T>(searchListPtr, minsPtr, 1, comm);

  //For each key decide which processor to send to
  unsigned int *part = NULL;

  if (keys.size())
  {
    part = new unsigned int[keys.size()];
    assert(part);
  }

  for (unsigned int i = 0; i < keys.size(); i++)
  {
    //maxLB returns the smallest index in a sorted array such
    //that a[ind] <= key and  a[index +1] > key
    bool found = par::maxLowerBound<T>(mins, keys[i], part + i, NULL, NULL);
    if (!found)
    {
      //This key is smaller than the mins from every processor.
      //No point in searching.
      part[i] = rank;
    }
  }

  mins.clear();

  int *numKeysSend = new int[npes];
  assert(numKeysSend);

  int *numKeysRecv = new int[npes];
  assert(numKeysRecv);

  for (int i = 0; i < npes; i++)
  {
    numKeysSend[i] = 0;
  }

  // calculate the number of keys to send ...
  for (unsigned int i = 0; i < keys.size(); i++)
  {
    numKeysSend[part[i]]++;
  }

  // Now do an All2All to get numKeysRecv
  par::Mpi_Alltoall<int>(numKeysSend, numKeysRecv, 1, comm);

  unsigned int totalKeys = 0; // total number of local keys ...
  for (int i = 0; i < npes; i++)
  {
    totalKeys += numKeysRecv[i];
  }

  // create the send and recv buffers ...
  std::vector<T> sendK(keys.size());
  std::vector<T> recvK(totalKeys);

  // the mapping ..
  unsigned int *comm_map = NULL;

  if (keys.size())
  {
    comm_map = new unsigned int[keys.size()];
    assert(comm_map);
  }

  // Now create sendK
  int *sendOffsets = new int[npes];
  assert(sendOffsets);
  sendOffsets[0] = 0;

  int *recvOffsets = new int[npes];
  assert(recvOffsets);
  recvOffsets[0] = 0;

  int *numKeysTmp = new int[npes];
  assert(numKeysTmp);
  numKeysTmp[0] = 0;

  // compute offsets ...
  for (int i = 1; i < npes; i++)
  {
    sendOffsets[i] = sendOffsets[i - 1] + numKeysSend[i - 1];
    recvOffsets[i] = recvOffsets[i - 1] + numKeysRecv[i - 1];
    numKeysTmp[i] = 0;
  }

  for (unsigned int i = 0; i < keys.size(); i++)
  {
    unsigned int ni = numKeysTmp[part[i]];
    numKeysTmp[part[i]]++;
    // set entry ...
    sendK[sendOffsets[part[i]] + ni] = keys[i];
    // save mapping .. will need it later ...
    comm_map[i] = sendOffsets[part[i]] + ni;
  }

  if (part)
  {
    delete[] part;
  }

  assert(numKeysTmp);
  delete[] numKeysTmp;
  numKeysTmp = NULL;

  T *sendKptr = NULL;
  T *recvKptr = NULL;
  if (!sendK.empty())
  {
    sendKptr = &(*(sendK.begin()));
  }
  if (!recvK.empty())
  {
    recvKptr = &(*(recvK.begin()));
  }

  par::Mpi_Alltoallv_sparse<T>(sendKptr, numKeysSend, sendOffsets,
                               recvKptr, numKeysRecv, recvOffsets, comm);

  std::vector<T> resSend(totalKeys);
  std::vector<T> resRecv(keys.size());

  //Final local search.
  for (unsigned int i = 0; i < totalKeys; i++)
  {
    unsigned int idx;
    bool found = par::maxLowerBound<T>(searchList, recvK[i], &idx, NULL, NULL);
    if (found)
    {
      resSend[i] = searchList[idx];
    }
  } //end for i

  //Exchange Results
  //Return what you received in the earlier communication.
  T *resSendPtr = NULL;
  T *resRecvPtr = NULL;
  if (!resSend.empty())
  {
    resSendPtr = &(*(resSend.begin()));
  }
  if (!resRecv.empty())
  {
    resRecvPtr = &(*(resRecv.begin()));
  }
  par::Mpi_Alltoallv_sparse<T>(resSendPtr, numKeysRecv, recvOffsets,
                               resRecvPtr, numKeysSend, sendOffsets, comm);

  assert(sendOffsets);
  delete[] sendOffsets;
  sendOffsets = NULL;

  assert(recvOffsets);
  delete[] recvOffsets;
  recvOffsets = NULL;

  assert(numKeysSend);
  delete[] numKeysSend;
  numKeysSend = NULL;

  assert(numKeysRecv);
  delete[] numKeysRecv;
  numKeysRecv = NULL;

  for (unsigned int i = 0; i < keys.size(); i++)
  {
    results[i] = resRecv[comm_map[i]];
  } //end for

  // Clean up ...
  if (comm_map)
  {
    delete[] comm_map;
  }

  PROF_SEARCH_END
}

template <typename T>
int partitionW(std::vector<T> &nodeList, unsigned int (*getWeight)(const T *), MPI_Comm comm)
{
  #ifdef __PROFILE_WITH_BARRIER__
    MPI_Barrier(comm);
  #endif
  PROF_PARTW_BEGIN

  int npes;

  MPI_Comm_size(comm, &npes);

  if (npes == 1)
  {
    PROF_PARTW_END
  }

  if (getWeight == NULL)
  {
    getWeight = par::defaultWeight<T>;
  }

  int rank;

  MPI_Comm_rank(comm, &rank);

  MPI_Request request;
  MPI_Status status;
  const bool nEmpty = nodeList.empty();

  DendroIntL off1 = 0, off2 = 0, localWt = 0, totalWt = 0;

  DendroIntL *wts = NULL;
  DendroIntL *lscn = NULL;
  DendroIntL nlSize = nodeList.size();
  if (nlSize)
  {
    wts = new DendroIntL[nlSize];
    assert(wts);

    lscn = new DendroIntL[nlSize];
    assert(lscn);
  }

  // First construct arrays of id and wts.
  #pragma omp parallel for reduction(+ \
                                    : localWt)
    for (DendroIntL i = 0; i < nlSize; i++)
    {
      wts[i] = (*getWeight)(&(nodeList[i]));
      localWt += wts[i];
    }

  #ifdef __DEBUG_PAR__
    MPI_Barrier(comm);
    if (!rank)
    {
      std::cout << "Partition: Stage-1 passed." << std::endl;
    }
    MPI_Barrier(comm);
  #endif

    // compute the total weight of the problem ...
    par::Mpi_Allreduce<DendroIntL>(&localWt, &totalWt, 1, MPI_SUM, comm);

    // perform a local scan on the weights first ...
    DendroIntL zero = 0;
    if (!nEmpty)
    {
      lscn[0] = wts[0];
      //        for (DendroIntL i = 1; i < nlSize; i++) {
      //          lscn[i] = wts[i] + lscn[i-1];
      //        }//end for
      omp_par::scan(&wts[1], lscn, nlSize);
      // now scan with the final members of
      par::Mpi_Scan<DendroIntL>(lscn + nlSize - 1, &off1, 1, MPI_SUM, comm);
    }
    else
    {
      par::Mpi_Scan<DendroIntL>(&zero, &off1, 1, MPI_SUM, comm);
    }

    // communicate the offsets ...
    if (rank < (npes - 1))
    {
      par::Mpi_Issend<DendroIntL>(&off1, 1, rank + 1, 0, comm, &request);
    }
    if (rank)
    {
      par::Mpi_Recv<DendroIntL>(&off2, 1, rank - 1, 0, comm, &status);
    }
    else
    {
      off2 = 0;
    }

    // add offset to local array
  #pragma omp parallel for
    for (DendroIntL i = 0; i < nlSize; i++)
    {
      lscn[i] = lscn[i] + off2; // This has the global scan results now ...
    }                           //end for

  #ifdef __DEBUG_PAR__
    MPI_Barrier(comm);
    if (!rank)
    {
      std::cout << "Partition: Stage-2 passed." << std::endl;
    }
    MPI_Barrier(comm);
  #endif

    int *sendSz = new int[npes];
    assert(sendSz);

    int *recvSz = new int[npes];
    assert(recvSz);

    int *sendOff = new int[npes];
    assert(sendOff);
    sendOff[0] = 0;

    int *recvOff = new int[npes];
    assert(recvOff);
    recvOff[0] = 0;

    // compute the partition offsets and sizes so that All2Allv can be performed.
    // initialize ...

  #pragma omp parallel for
    for (int i = 0; i < npes; i++)
    {
      sendSz[i] = 0;
    }

    // Now determine the average load ...
    DendroIntL npesLong = npes;
    DendroIntL avgLoad = (totalWt / npesLong);

    DendroIntL extra = (totalWt % npesLong);

    //The Heart of the algorithm....
    if (avgLoad > 0)
    {
      for (DendroIntL i = 0; i < nlSize; i++)
      {
        if (lscn[i] == 0)
        {
          sendSz[0]++;
        }
        else
        {
          int ind = 0;
          if (lscn[i] <= (extra * (avgLoad + 1)))
          {
            ind = ((lscn[i] - 1) / (avgLoad + 1));
          }
          else
          {
            ind = ((lscn[i] - (1 + extra)) / avgLoad);
          }
          assert(ind < npes);
          sendSz[ind]++;
        } //end if-else
      }   //end for */


    #ifdef __DEBUG_PAR__
      int tmp_sum = 0;
      for (int i = 0; i < npes; i++)
        tmp_sum += sendSz[i];
      assert(tmp_sum == nlSize);
    #endif
    }
    else
    {
      sendSz[0] += nlSize;
    } //end if-else

    #ifdef __DEBUG_PAR__
    MPI_Barrier(comm);
    if (!rank)
    {
      std::cout << "Partition: Stage-3 passed." << std::endl;
    }
    MPI_Barrier(comm);
    #endif

    if (rank < (npes - 1))
    {
      MPI_Status statusWait;
      MPI_Wait(&request, &statusWait);
    }

    // communicate with other procs how many you shall be sending and get how
    // many to recieve from whom.
    par::Mpi_Alltoall<int>(sendSz, recvSz, 1, comm);

    #ifdef __DEBUG_PAR__
    DendroIntL totSendToOthers = 0;
    DendroIntL totRecvFromOthers = 0;
    for (int i = 0; i < npes; i++)
    {
      if (rank != i)
      {
        totSendToOthers += sendSz[i];
        totRecvFromOthers += recvSz[i];
      }
    }
    #endif

    DendroIntL nn = 0; // new value of nlSize, ie the local nodes.
    #pragma omp parallel for reduction(+ \
                                      : nn)
    for (int i = 0; i < npes; i++)
    {
      nn += recvSz[i];
    }

    // compute offsets ...
    //      for (int i = 1; i < npes; i++) {
    //        sendOff[i] = sendOff[i-1] + sendSz[i-1];
    //        recvOff[i] = recvOff[i-1] + recvSz[i-1];
    //      }
    omp_par::scan(sendSz, sendOff, npes);
    omp_par::scan(recvSz, recvOff, npes);

    #ifdef __DEBUG_PAR__
    MPI_Barrier(comm);
    if (!rank)
    {
      std::cout << "Partition: Stage-4 passed." << std::endl;
    }
    MPI_Barrier(comm);
    /*
          std::cout<<rank<<": newSize: "<<nn<<" oldSize: "<<(nodeList.size())
          <<" send: "<<totSendToOthers<<" recv: "<<totRecvFromOthers<<std::endl;
        */
    MPI_Barrier(comm);
    #endif

    // allocate memory for the new arrays ...
    std::vector<T> newNodes(nn);

    #ifdef __DEBUG_PAR__
    MPI_Barrier(comm);
    if (!rank)
    {
      std::cout << "Partition: Final alloc successful." << std::endl;
    }
    MPI_Barrier(comm);
    #endif

    // perform All2All  ...
    T *nodeListPtr = NULL;
    T *newNodesPtr = NULL;
    if (!nodeList.empty())
    {
      nodeListPtr = &(*(nodeList.begin()));
    }
    if (!newNodes.empty())
    {
      newNodesPtr = &(*(newNodes.begin()));
    }
    par::Mpi_Alltoallv_sparse<T>(nodeListPtr, sendSz, sendOff,
                                  newNodesPtr, recvSz, recvOff, comm);

    #ifdef __DEBUG_PAR__
    MPI_Barrier(comm);
    if (!rank)
    {
      std::cout << "Partition: Stage-5 passed." << std::endl;
    }
    MPI_Barrier(comm);
    #endif

    // reset the pointer ...
    swap(nodeList, newNodes);
    newNodes.clear();

    // clean up...
    if (!nEmpty)
    {
      delete[] lscn;
      delete[] wts;
    }
    delete[] sendSz;
    sendSz = NULL;

    delete[] sendOff;
    sendOff = NULL;

    delete[] recvSz;
    recvSz = NULL;

    delete[] recvOff;
    recvOff = NULL;

    #ifdef __DEBUG_PAR__
    MPI_Barrier(comm);
    if (!rank)
    {
      std::cout << "Partition: Stage-6 passed." << std::endl;
    }
    MPI_Barrier(comm);
    #endif

    PROF_PARTW_END
} //end function


template <typename T>
int sampleSort(std::vector<T> &arr, std::vector<T> &SortedElem, std::vector<double> &stats, MPI_Comm comm)
{
  // std::cout<<"Sample Sort Execution Begin"<<std::endl;
  #ifdef __PROFILE_WITH_BARRIER__
    MPI_Barrier(comm);
  #endif
  PROF_SORT_BEGIN

  auto ss_start = std::chrono::system_clock::now();

  int npes;
  MPI_Comm_size(comm, &npes);
  //--
  int rank;
  MPI_Comm_rank(comm, &rank);

  //std::cout << rank << "Nodes Input : " << __func__ << arr.size() << std::endl;

  //--

  assert(arr.size());

  if (npes == 1)
  {
    std::cout << " have to use seq. sort"
              << " since npes = 1 . inpSize: " << (arr.size()) << std::endl;
    // std::sort(arr.begin(), arr.end());
    omp_par::merge_sort(arr.begin(), arr.end());
    SortedElem = arr;
    PROF_SORT_END
  }

  std::vector<T> splitters;
  std::vector<T> allsplitters;

  int myrank;
  MPI_Comm_rank(comm, &myrank);

  DendroIntL nelem = arr.size();
  DendroIntL nelemCopy = nelem;
  DendroIntL totSize;
  /*Calculating the total size of elements*/
  par::Mpi_Allreduce<DendroIntL>(&nelemCopy, &totSize, 1, MPI_SUM, comm);

  DendroIntL npesLong = npes;
  const DendroIntL FIVE = 5;

  if (totSize < (FIVE * npesLong * npesLong))
  {
    if (!myrank)
    {
      std::cout << " *Using bitonic sort since totSize < (5*(npes^2)). totSize: "
                << totSize << " npes: " << npes << std::endl;
    }

  #ifdef __DEBUG_PAR__
    MPI_Barrier(comm);
    if (!myrank)
    {
      std::cout << "SampleSort (small n): Stage-1 passed." << std::endl;
    }
    MPI_Barrier(comm);
  #endif

    SortedElem = arr;
    MPI_Comm new_comm;
    if (totSize < npesLong)
    {
      if (!myrank)
      {
        std::cout << " Input to sort is small. splittingComm: "
                  << npes << " -> " << totSize << std::endl;
      }
      par::splitCommUsingSplittingRank(static_cast<int>(totSize), &new_comm, comm);
    }
    else
    {
      new_comm = comm;
    }

  #ifdef __DEBUG_PAR__
    MPI_Barrier(comm);
    if (!myrank)
    {
      std::cout << "SampleSort (small n): Stage-2 passed." << std::endl;
    }
    MPI_Barrier(comm);
  #endif

    if (!SortedElem.empty())
    {
      par::bitonicSort<T>(SortedElem, new_comm);
      //std::cout<<"Rank:"<<rank<<" bitonic search complete"<<std::endl;
    }

  #ifdef __DEBUG_PAR__
    MPI_Barrier(comm);
    if (!myrank)
    {
      std::cout << "SampleSort (small n): Stage-3 passed." << std::endl;
    }
    MPI_Barrier(comm);
  #endif

    PROF_SORT_END
  } // end if

  // if(!rank)

  #ifdef __DEBUG_PAR__
  if (!myrank)
  {
    std::cout << "Using sample sort to sort nodes. n/p^2 is fine." << std::endl;
  }
  #endif

  //Re-part arr so that each proc. has atleast p elements.

  // if (!rank) std::cout << "[samplesort] repartitioning input" << std::endl;
  par::partitionW<T>(arr, NULL, comm);

  nelem = arr.size();

  // if (!rank) std::cout << RED "[samplesort] initial local sort" NRM << std::endl;
  // std::sort(arr.begin(),arr.end());
  auto splitterCalculation_start = std::chrono::system_clock::now();
  omp_par::merge_sort(arr.begin(), arr.end());

  // if (!rank) std::cout << GRN "[samplesort] initial local sort" NRM << std::endl;

  std::vector<T> sendSplits(npes - 1);
  splitters.resize(npes);
  std::vector<ot::TreeNode> treenode_split;
  ot::TreeNode temp;
  #pragma omp parallel for
  for (int i = 1; i < npes; i++)
  {

    sendSplits[i - 1] = arr[i * nelem / npes];

    //std::cout << BLU << "===============================================" << NRM << std::endl;
    // std::cout << RED " Splittersn in Sample Sort Rank:"<<rank<< NRM << std::endl;
    //std::cout << BLU << "===============================================" << NRM << std::endl;

    std::ostringstream convert;
    convert << sendSplits[i - 1];
    std::vector<std::string> results;
    std::stringstream s(convert.str());
    while (!s.eof())
    {
      std::string tmp;
      s >> tmp;
      results.push_back(tmp);
    }

    temp = ot::TreeNode(1, atoi(results[0].c_str()), atoi(results[1].c_str()), atoi(results[2].c_str()),
                        atoi(results[3].c_str()), 3, 8);
    treenode_split.push_back(temp);
    // std::cout <<"Rank: "<<rank<<"  "<<sendSplits[i-1]<<std::endl;

  } //end for i

  // treeNodesTovtk(treenode_split,rank,"SplittersbfBiTonicSort");
  //std::cout << myrank << ": snedSplit " << sendSplits[0] << "  ||||  " <<  sendSplits[1] << std::endl;
  //NOTE: @hari fails for npes=2
  // sort sendSplits using bitonic ...
  // if (!rank)
  // std::cout << rank << RED " [samplesort] bitonicsort " NRM << sendSplits.size() << std::endl;
  par::bitonicSort<T>(sendSplits, comm);
  // if (!rank)
  // std::cout << rank << GRN " [samplesort] done bitonicsort" NRM << std::endl;

  //treeNodesTovtk(treenode_spliters,rank,"SplitersafBitTonicSort");

  //std::cout << myrank << ": afterBitonic " << sendSplits[0] << "  ||||  " <<  sendSplits[1] << std::endl;

  // All gather with last element of splitters.
  T *sendSplitsPtr = NULL;
  T *splittersPtr = NULL;
  if (sendSplits.size() > static_cast<unsigned int>(npes - 2))
  {
    sendSplitsPtr = &(*(sendSplits.begin() + (npes - 2)));
  }
  if (!splitters.empty())
  {
    splittersPtr = &(*(splitters.begin()));
  }
  par::Mpi_Allgather<T>(sendSplitsPtr, splittersPtr, 1, comm);

  //       std::cout<<"Rank:"<<rank<<" :splittersPtr Size"<<sizeof(splittersPtr)<<std::endl;
  //       //for(int i=0;i<sizeof(splittersPtr);i++){
  // 	  std::cout<<"Rank:"<<rank<<":"<<(splittersPtr->m_uiX)<<std::endl;
  //splittersPtr[i].printTreeNode();
  //}
  sendSplits.clear();
  int *sendcnts = new int[npes];
  assert(sendcnts);

  int *recvcnts = new int[npes];
  assert(recvcnts);

  int *sdispls = new int[npes];
  assert(sdispls);

  int *rdispls = new int[npes];
  assert(rdispls);

  #pragma omp parallel for
  for (int k = 0; k < npes; k++)
  {
    sendcnts[k] = 0;
  }

  //To be parallelized
  /*      int k = 0;
      for (DendroIntL j = 0; j < nelem; j++) {
        if (arr[j] <= splitters[k]) {
          sendcnts[k]++;
        } else{
          k = seq::UpperBound<T>(npes-1, splittersPtr, k+1, arr[j]);
          if (k == (npes-1) ){
            //could not find any splitter >= arr[j]
            sendcnts[k] = (nelem - j);
            break;
          } else {
            assert(k < (npes-1));
            assert(splitters[k] >= arr[j]);
            sendcnts[k]++;
          }
        }//end if-else
      }//end for j
  */

  auto splitterCalculation_end = std::chrono::system_clock::now();
  auto all2all_start = std::chrono::system_clock::now();

  {
    int omp_p = omp_get_max_threads();
    int *proc_split = new int[omp_p + 1];
    DendroIntL *lst_split_indx = new DendroIntL[omp_p + 1];
    proc_split[0] = 0;
    lst_split_indx[0] = 0;
    lst_split_indx[omp_p] = nelem;
  #pragma omp parallel for
    for (int i = 1; i < omp_p; i++)
    {
      //proc_split[i] = seq::BinSearch(&splittersPtr[0],&splittersPtr[npes-1],arr[i*nelem/omp_p],std::less<T>());
      proc_split[i] =
          std::upper_bound(&splittersPtr[0], &splittersPtr[npes - 1], arr[i * nelem / omp_p], std::less<T>()) -
          &splittersPtr[0];
      if (proc_split[i] < npes - 1)
      {
        //lst_split_indx[i]=seq::BinSearch(&arr[0],&arr[nelem],splittersPtr[proc_split[i]],std::less<T>());
        lst_split_indx[i] =
            std::upper_bound(&arr[0], &arr[nelem], splittersPtr[proc_split[i]], std::less<T>()) - &arr[0];
      }
      else
      {
        proc_split[i] = npes - 1;
        lst_split_indx[i] = nelem;
      }
    }

    splitterCalculation_end = std::chrono::system_clock::now();
    all2all_start = std::chrono::system_clock::now();

  #pragma omp parallel for
    for (int i = 0; i < omp_p; i++)
    {
      int sendcnts_ = 0;
      int k = proc_split[i];
      for (DendroIntL j = lst_split_indx[i]; j < lst_split_indx[i + 1]; j++)
      {
        if (arr[j] <= splitters[k])
        {
          sendcnts_++;
        }
        else
        {
          if (sendcnts_ > 0)
            sendcnts[k] = sendcnts_;
          sendcnts_ = 0;
          k = seq::UpperBound<T>(npes - 1, splittersPtr, k + 1, arr[j]);
          if (k == (npes - 1))
          {
            //could not find any splitter >= arr[j]
            sendcnts_ = (nelem - j);
            break;
          }
          else
          {
            assert(k < (npes - 1));
            assert(splitters[k] >= arr[j]);
            sendcnts_++;
          }
        } //end if-else
      }   //end for j
      if (sendcnts_ > 0)
        sendcnts[k] = sendcnts_;
    }
    delete[] lst_split_indx;
    delete[] proc_split;
  }

  par::Mpi_Alltoall<int>(sendcnts, recvcnts, 1, comm);

  sdispls[0] = 0;
  rdispls[0] = 0;
  //      for (int j = 1; j < npes; j++){
  //        sdispls[j] = sdispls[j-1] + sendcnts[j-1];
  //        rdispls[j] = rdispls[j-1] + recvcnts[j-1];
  //      }
  omp_par::scan(sendcnts, sdispls, npes);
  omp_par::scan(recvcnts, rdispls, npes);

  DendroIntL totalRecv = rdispls[npes - 1] + recvcnts[npes - 1] - recvcnts[rank];
  DendroIntL totalSend = sdispls[npes - 1] + sendcnts[npes - 1] - sendcnts[rank];

  DendroIntL nsorted = rdispls[npes - 1] + recvcnts[npes - 1];
  SortedElem.resize(nsorted);

  T *arrPtr = NULL;
  T *SortedElemPtr = NULL;
  if (!arr.empty())
  {
    arrPtr = &(*(arr.begin()));
  }
  if (!SortedElem.empty())
  {
    SortedElemPtr = &(*(SortedElem.begin()));
  }

  par::Mpi_Alltoallv_Kway<T>(arrPtr, sendcnts, sdispls, SortedElemPtr, recvcnts, rdispls, comm);

  arr.clear();

  delete[] sendcnts;
  sendcnts = NULL;

  delete[] recvcnts;
  recvcnts = NULL;

  delete[] sdispls;
  sdispls = NULL;

  delete[] rdispls;
  rdispls = NULL;

  auto all2all_end = std::chrono::system_clock::now();

  DendroIntL localSz = SortedElem.size();

  // sort(SortedElem.begin(), SortedElem.end());
  // if (!rank) std::cout << "[samplesort] final local sort" << std::endl;
  auto localSort_start = std::chrono::system_clock::now();

  omp_par::merge_sort(&SortedElem[0], &SortedElem[nsorted]);
  // if (!rank) std::cout << "[samplesort] final local sort" << std::endl;

  auto localSort_end = std::chrono::system_clock::now();
  auto ss_end = std::chrono::system_clock::now();

  double splittercal = (std::chrono::duration_cast<std::chrono::milliseconds>((splitterCalculation_end - splitterCalculation_start)).count()) / (MILLISECOND_CONVERSION);
  double all2all = (std::chrono::duration_cast<std::chrono::milliseconds>((all2all_end - all2all_start)).count() / (MILLISECOND_CONVERSION));
  double localSort = (std::chrono::duration_cast<std::chrono::milliseconds>((localSort_end - localSort_start)).count() / (MILLISECOND_CONVERSION));

  double splitterCalculation_g[3];
  double all2all_g[3];
  double localSort_g[3];

  DendroIntL sendCounts_g[3];
  DendroIntL recvCounts_g[3];

  DendroIntL localSz_g[3];

  MPI_Reduce(&splittercal, splitterCalculation_g, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&splittercal, splitterCalculation_g + 1, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(&splittercal, splitterCalculation_g + 2, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  splitterCalculation_g[1] = splitterCalculation_g[1] / npes;

  MPI_Reduce(&all2all, all2all_g, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&all2all, all2all_g + 1, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(&all2all, all2all_g + 2, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  all2all_g[1] = all2all_g[1] / npes;

  MPI_Reduce(&localSort, localSort_g, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&localSort, localSort_g + 1, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(&localSort, localSort_g + 2, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  localSort_g[1] = localSort_g[1] / npes;

  MPI_Reduce(&totalSend, sendCounts_g, 1, MPI_LONG_LONG, MPI_MIN, 0, comm);
  MPI_Reduce(&totalSend, sendCounts_g + 1, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);
  MPI_Reduce(&totalSend, sendCounts_g + 2, 1, MPI_LONG_LONG, MPI_MAX, 0, comm);
  sendCounts_g[1] = sendCounts_g[1] / (double)npes;

  MPI_Reduce(&totalRecv, recvCounts_g, 1, MPI_LONG_LONG, MPI_MIN, 0, comm);
  MPI_Reduce(&totalRecv, recvCounts_g + 1, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);
  MPI_Reduce(&totalRecv, recvCounts_g + 2, 1, MPI_LONG_LONG, MPI_MAX, 0, comm);
  recvCounts_g[1] = recvCounts_g[1] / (double)npes;

  MPI_Reduce(&localSz, localSz_g, 1, MPI_LONG_LONG, MPI_MIN, 0, comm);
  //MPI_Reduce(&localSz,localSz_g+1,1,MPI_LONG_LONG,MPI_SUM,0,comm);
  MPI_Reduce(&localSz, localSz_g + 2, 1, MPI_LONG_LONG, MPI_MAX, 0, comm);
  localSz_g[1] = 0; //localSort_g[1]/(double)npes;

  if (!rank)
  {
    //                std::cout<<splitterCalculation_g[0]<<"\t"<<splitterCalculation_g[1]<<"\t"<<splitterCalculation_g[2];
    //                std::cout<<"\t"<<all2all_g[0]<<"\t"<<all2all_g[1]<<"\t"<<all2all_g[2];
    //                std::cout<<"\t"<<localSort_g[0]<<"\t"<<localSort_g[1]<<"\t"<<localSort_g[2]<<std::endl;

    stats.push_back(splitterCalculation_g[0]); // ss splitter min
    stats.push_back(splitterCalculation_g[1]); // ss splitter mean
    stats.push_back(splitterCalculation_g[2]); // ss splitter max

    stats.push_back(all2all_g[0]); // ss all2all min
    stats.push_back(all2all_g[1]); // ss all2all mean
    stats.push_back(all2all_g[2]); // ss all2all max

    stats.push_back(localSort_g[0]); // ss local sort min
    stats.push_back(localSort_g[1]); // ss local sort mean
    stats.push_back(localSort_g[2]); // ss local sort max

    stats.push_back(sendCounts_g[0]); // ss send counts min
    stats.push_back(sendCounts_g[1]); // ss send counts mean
    stats.push_back(sendCounts_g[2]); // ss send counts max

    stats.push_back(recvCounts_g[0]); // ss recv counts min
    stats.push_back(recvCounts_g[1]); // ss recv counts mean
    stats.push_back(recvCounts_g[2]); // ss recv counts max

    stats.push_back(localSz_g[0]); // ss local size min
    stats.push_back(localSz_g[1]); // ss local size mean
    stats.push_back(localSz_g[2]); // ss local size max
  }

  double ss_time = (std::chrono::duration_cast<std::chrono::milliseconds>((ss_end - ss_start)).count()) / (MILLISECOND_CONVERSION);
  double ss_time_g[3];

  MPI_Reduce(&ss_time, ss_time_g, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&ss_time, ss_time_g + 1, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(&ss_time, ss_time_g + 2, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  ss_time_g[1] = ss_time_g[1] / npes;

  if (!rank)
  {
    /*std::cout << "SS (min mean max) \t TS (min mean max)" << std::endl;
           std::cout << ss_time_g[0] << "\t" << ss_time_g[1] << "\t" << ss_time_g[2] <<std::endl;*/

    stats.push_back(ss_time_g[0]); // ss overall time min
    stats.push_back(ss_time_g[1]); // ss overall time mean
    stats.push_back(ss_time_g[2]); // ss overall time max

    //!!! Note : Quick fix for sample sort to match with the lengths for the ts_stats vector.
    stats.push_back(0);
    stats.push_back(0);
    stats.push_back(0);

    stats.push_back(0);
    stats.push_back(0);
    stats.push_back(0);
  }

  PROF_SORT_END

} //end function

 /********************************************************************/
 /*
  * which_keys is one of KEEP_HIGH or KEEP_LOW
  * partner    is the processor with which to Merge and Split.
  *
  */
template <typename T>
void MergeSplit(std::vector<T> &local_list, int which_keys, int partner, MPI_Comm comm)
{

  MPI_Status status;
  int send_size = local_list.size();
  int recv_size = 0;

  // first communicate how many you will send and how many you will receive ...

  par::Mpi_Sendrecv<int, int>(&send_size, 1, partner, 0,
                              &recv_size, 1, partner, 0, comm, &status);

  std::vector<T> temp_list(recv_size);

  T *local_listPtr = NULL;
  T *temp_listPtr = NULL;
  if (!local_list.empty())
  {
    local_listPtr = &(*(local_list.begin()));
  }
  if (!temp_list.empty())
  {
    temp_listPtr = &(*(temp_list.begin()));
  }

  par::Mpi_Sendrecv<T, T>(local_listPtr, send_size, partner,
                          1, temp_listPtr, recv_size, partner, 1, comm, &status);

  MergeLists<T>(local_list, temp_list, which_keys);

  temp_list.clear();
} // Merge_split

template <typename T>
void Par_bitonic_sort_incr(std::vector<T> &local_list, int proc_set_size, MPI_Comm comm)
{
  int eor_bit;
  int proc_set_dim;
  int stage;
  int partner;
  int my_rank;

  MPI_Comm_rank(comm, &my_rank);

  proc_set_dim = 0;
  int x = proc_set_size;
  while (x > 1)
  {
    x = x >> 1;
    proc_set_dim++;
  }

  eor_bit = (1 << (proc_set_dim - 1));
  for (stage = 0; stage < proc_set_dim; stage++)
  {
    partner = (my_rank ^ eor_bit);

    if (my_rank < partner)
    {
      MergeSplit<T>(local_list, KEEP_LOW, partner, comm);
    }
    else
    {
      MergeSplit<T>(local_list, KEEP_HIGH, partner, comm);
    }

    eor_bit = (eor_bit >> 1);
  }
} // Par_bitonic_sort_incr

template <typename T>
void Par_bitonic_sort_decr(std::vector<T> &local_list, int proc_set_size, MPI_Comm comm)
{
  int eor_bit;
  int proc_set_dim;
  int stage;
  int partner;
  int my_rank;

  MPI_Comm_rank(comm, &my_rank);

  proc_set_dim = 0;
  int x = proc_set_size;
  while (x > 1)
  {
    x = x >> 1;
    proc_set_dim++;
  }

  eor_bit = (1 << (proc_set_dim - 1));
  for (stage = 0; stage < proc_set_dim; stage++)
  {
    partner = my_rank ^ eor_bit;

    if (my_rank > partner)
    {
      MergeSplit<T>(local_list, KEEP_LOW, partner, comm);
    }
    else
    {
      MergeSplit<T>(local_list, KEEP_HIGH, partner, comm);
    }

    eor_bit = (eor_bit >> 1);
  }

} // Par_bitonic_sort_decr

template <typename T>
void Par_bitonic_merge_incr(std::vector<T> &local_list, int proc_set_size, MPI_Comm comm)
{
  int partner;
  int rank, npes;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &npes);

  unsigned int num_left = binOp::getPrevHighestPowerOfTwo(npes);
  unsigned int num_right = npes - num_left;

  // 1, Do merge between the k right procs and the highest k left procs.
  if ((static_cast<unsigned int>(rank) < num_left) &&
      (static_cast<unsigned int>(rank) >= (num_left - num_right)))
  {
    partner = static_cast<unsigned int>(rank) + num_right;
    MergeSplit<T>(local_list, KEEP_LOW, partner, comm);
  }
  else if (static_cast<unsigned int>(rank) >= num_left)
  {
    partner = static_cast<unsigned int>(rank) - num_right;
    MergeSplit<T>(local_list, KEEP_HIGH, partner, comm);
  }
}

template <typename T>
void bitonicSort_binary(std::vector<T> &in, MPI_Comm comm)
{
  int proc_set_size;
  unsigned int and_bit;
  int rank;
  int npes;

  MPI_Comm_size(comm, &npes);

  #ifdef __DEBUG_PAR__
  assert(npes > 1);
  assert(!(npes & (npes - 1)));
  assert(!(in.empty()));
  #endif

  MPI_Comm_rank(comm, &rank);

  for (proc_set_size = 2, and_bit = 2;
       proc_set_size <= npes;
       proc_set_size = proc_set_size * 2,
      and_bit = and_bit << 1)
  {

    if ((rank & and_bit) == 0)
    {
      Par_bitonic_sort_incr<T>(in, proc_set_size, comm);
    }
    else
    {
      Par_bitonic_sort_decr<T>(in, proc_set_size, comm);
    }
  } //end for
}

template <typename T>
void bitonicSort(std::vector<T> &in, MPI_Comm comm)
{
  int rank;
  int npes;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &rank);

  assert(!(in.empty()));

  //Local Sort first
  //std::sort(in.begin(),in.end());
  omp_par::merge_sort(in.begin(), in.end());

  if (npes > 1)
  {

    // check if npes is a power of two ...
    bool isPower = (!(npes & (npes - 1)));

    if (isPower)
    {
      bitonicSort_binary<T>(in, comm);
    }
    else
    {
      MPI_Comm new_comm;

      // Since npes is not a power of two, we shall split the problem in two ...
      //
      // 1. Create 2 comm groups ... one for the 2^d portion and one for the
      // remainder.
      unsigned int splitter = splitCommBinary(comm, &new_comm);

      if (static_cast<unsigned int>(rank) < splitter)
      {
        bitonicSort_binary<T>(in, new_comm);
      }
      else
      {
        bitonicSort<T>(in, new_comm);
      }

      // 3. Do a special merge of the two segments. (original comm).
      Par_bitonic_merge_incr(in, binOp::getNextHighestPowerOfTwo(npes), comm);
      MPI_Comm_free(&new_comm);
      splitter = splitCommBinaryNoFlip(comm, &new_comm);

      // 4. Now a final sort on the segments.
      if (static_cast<unsigned int>(rank) < splitter)
      {
        bitonicSort_binary<T>(in, new_comm);
      }
      else
      {
        bitonicSort<T>(in, new_comm);
      }
      MPI_Comm_free(&new_comm);
    } //end if isPower of 2
  }   //end if single processor
} //end function

template <typename T>
void MergeLists(std::vector<T> &listA, std::vector<T> &listB, int KEEP_WHAT)
{

  T _low, _high;

  assert(!(listA.empty()));
  assert(!(listB.empty()));

  // max ( min(A,B) )
  _low = ((listA[0] > listB[0]) ? listA[0] : listB[0]);
  // min ( max(A,B) )
  _high = ((listA[listA.size() - 1] < listB[listB.size() - 1]) ? listA[listA.size() - 1] : listB[listB.size() - 1]);

  // We will do a full merge first ...
  size_t list_size = listA.size() + listB.size();

  std::vector<T> scratch_list(list_size);

  unsigned int index1 = 0;
  unsigned int index2 = 0;

  for (size_t i = 0; i < list_size; i++)
  {
    //The order of (A || B) is important here,
    //so that index2 remains within bounds
    if ((index1 < listA.size()) &&
        ((index2 >= listB.size()) ||
         (listA[index1] <= listB[index2])))
    {
      scratch_list[i] = listA[index1];
      index1++;
    }
    else
    {
      scratch_list[i] = listB[index2];
      index2++;
    }
  }

  listA.clear();
  listB.clear();

  if (KEEP_WHAT == KEEP_LOW)
  {
    int ii = 0;

    while (((scratch_list[ii] < _low) ||
            (ii < (list_size / 2))) &&
           (scratch_list[ii] <= _high))
    {
      ii++;
    }
    if (ii)
    {
      listA.insert(listA.end(), scratch_list.begin(),
                   (scratch_list.begin() + ii));
    }
  }
  else
  {
    int ii = (list_size - 1);
    while (((ii >= (list_size / 2)) && (scratch_list[ii] >= _low)) || (scratch_list[ii] > _high))
    {
      ii--;
    }
    if (ii < (list_size - 1))
    {
      listA.insert(listA.begin(), (scratch_list.begin() + (ii + 1)),
                   (scratch_list.begin() + list_size));
    }
  }
  scratch_list.clear();
} //end function

template <typename T>
void parallel_rank(const T *in, unsigned int sz, DendroIntL *out, MPI_Comm comm)
{

  int rank, npes;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &npes);

  std::vector<_T<T>> key;
  std::vector<_T<T>> key_sorted;
  key.resize(sz);

  for (unsigned int i = 0; i < sz; i++)
  {

    key[i].p = rank;
    key[i].idx = i;
    key[i].val = in[i];
  }

  std::vector<double> stats;

  par::bitonicSort(key, comm);
  std::swap(key, key_sorted);
  key.clear();

  unsigned int localSz = sz;

  unsigned int *sorted_couts = new unsigned int[npes];
  unsigned int *sorted_offset = new unsigned int[npes];

  par::Mpi_Allgather(&localSz, sorted_couts, 1, comm);
  sorted_offset[0] = 0;
  omp_par::scan(sorted_couts, sorted_offset, npes);

  for (unsigned int i = 0; i < key_sorted.size(); i++)
    key_sorted[i].rank = sorted_offset[rank] + i;

  // for(unsigned int i=0;i<key_sorted.size();i++)
  // {
  //   std::cout<<"val: "<<key_sorted[i].val<<" p: "<<key_sorted[i].p<<" rank: "<<key_sorted[i].rank<<" idx : "<<key_sorted[i].idx<<std::endl;
  // }

  int *sCounts = new int[npes];
  int *rCounts = new int[npes];
  int *sOffset = new int[npes];
  int *rOffset = new int[npes];
  int *ccount = new int[npes];

  for (unsigned int i = 0; i < npes; i++)
  {
    sCounts[i] = 0;
    ccount[i] = 0;
  }

  for (unsigned int i = 0; i < key_sorted.size(); i++)
    sCounts[key_sorted[i].p]++;

  for (unsigned int i = 0; i < npes; i++)
    sCounts[i] = sizeof(_T<T>) * sCounts[i];

  par::Mpi_Alltoall(sCounts, rCounts, 1, comm);

  sOffset[0] = 0;
  rOffset[0] = 0;

  omp_par::scan(sCounts, sOffset, npes);
  omp_par::scan(rCounts, rOffset, npes);

  std::vector<_T<T>> sBuf;
  sBuf.resize((sOffset[npes - 1] + sCounts[npes - 1]) / sizeof(_T<T>));

  for (unsigned int i = 0; i < key_sorted.size(); i++)
  {
    sBuf[sOffset[key_sorted[i].p] + ccount[key_sorted[i].p]] = key_sorted[i];
    ccount[key_sorted[i].p]++;
  }

  key.resize((rOffset[npes - 1] + rCounts[npes - 1]) / sizeof(_T<T>));

  MPI_Alltoallv(&(*(sBuf.begin())), sCounts, sOffset, MPI_BYTE, &(*(key.begin())), rCounts, rOffset, MPI_BYTE, comm);

  for (unsigned int i = 0; i < key.size(); i++)
    out[key[i].idx] = key[i].rank;

  delete[] sCounts;
  delete[] rCounts;
  delete[] sOffset;
  delete[] rOffset;
  delete[] ccount;

  delete[] sorted_couts;
  delete[] sorted_offset;
}


template<typename T>
void computeOverallStats(T *stat, T *stat_g, MPI_Comm comm, const char* const prefix)
{
    int rank,npes;
    MPI_Comm_size(comm,&npes);
    MPI_Comm_rank(comm,&rank);

    par::Mpi_Reduce(stat,stat_g,1,MPI_MIN,0,comm);
    par::Mpi_Reduce(stat,stat_g+1,1,MPI_SUM,0,comm);
    par::Mpi_Reduce(stat,stat_g+2,1,MPI_MAX,0,comm);
    stat_g[1]/=(npes);

    if(!rank)
        printf("--- %s time(s) [min mean max]--- \n%.8E\t%.8E\t%.8E\n",prefix, stat_g[0],stat_g[1],stat_g[2]);

}


} // namespace par
