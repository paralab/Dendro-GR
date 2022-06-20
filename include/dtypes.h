#ifndef __DTYPES_H_
#define __DTYPES_H_

#include <mpi.h>
#include <complex>

/** 
  @file	dtypes.h
  @brief	Traits to determine MPI_DATATYPE from a C++ datatype
  @author	Hari Sundar, hsundar@gmail.com

  Traits to determine MPI_DATATYPE from a C++ datatype. For non standard
  C++ datatypes (like classes), we will need to define additional classes. An
  example is given for the case of the std. complex variable. Additional
  classes can be added as required.
 **/ 

namespace par {

/**
@class Mpi_datatype
@brief An abstract class used for communicating messages using user-defined datatypes.
 The user must implement the static member function "value()" that returns the MPI_Datatype corresponding 
to this user-defined datatype.
@author	Hari Sundar, hsundar@gmail.com
@see Mpi_datatype<bool>
*/
    template <typename T>
    class Mpi_datatype;

#define HS_MPIDATATYPE(CTYPE, MPITYPE)	\
	template <>  \
		class Mpi_datatype<CTYPE> \
		{  \
		    public: \
			static MPI_Datatype value() {\
                          return MPITYPE;\
                        } \
		};

    HS_MPIDATATYPE(short,          MPI_SHORT)
    HS_MPIDATATYPE(int,            MPI_INT)
    HS_MPIDATATYPE(long,           MPI_LONG)
    HS_MPIDATATYPE(unsigned short, MPI_UNSIGNED_SHORT)
    HS_MPIDATATYPE(unsigned int,   MPI_UNSIGNED)
    HS_MPIDATATYPE(unsigned long,  MPI_UNSIGNED_LONG)
    HS_MPIDATATYPE(float,          MPI_FLOAT)
    HS_MPIDATATYPE(double,         MPI_DOUBLE)
    HS_MPIDATATYPE(long double,    MPI_LONG_DOUBLE)
    HS_MPIDATATYPE(long long,    MPI_LONG_LONG_INT)
    HS_MPIDATATYPE(char,    MPI_CHAR)
    HS_MPIDATATYPE(unsigned char,    MPI_UNSIGNED_CHAR)

    //PetscScalar is simply a typedef for double. Hence no need to explicitly
    //define an mpi_datatype for it.

#undef HS_MPIDATATYPE

    template <typename T>
    class Mpi_datatype<std::complex<T> > {
    public:
        static MPI_Datatype value()
        {
            static bool         first = true;
            static MPI_Datatype datatype;

            if (first) {
                first = false;
                MPI_Type_contiguous(2, Mpi_datatype<T>::value(), &datatype);
                MPI_Type_commit(&datatype);
            }

            return datatype;
        }
    };

    /**
    @brief A template specialization of the abstract class Mpi_datatype. This can be used for communicating messages of type "bool"
    @author Rahul Sampath, rahul.sampath@gmail.com
    */
    template <>
    class Mpi_datatype<bool> {
        static void bool_LOR(void *in, void *inout, int* len, MPI_Datatype * dptr) {
            for (int i = 0; i < (*len); i++) {
                ((static_cast<bool*>(inout))[i]) =
                  ( ((static_cast<bool*>(in))[i]) || 
                    ((static_cast<bool*>(inout))[i]) );
            }//end for	
        }//end function


        static void bool_LAND(void *in, void *inout, int* len, MPI_Datatype * dptr) {
            for (int i = 0; i < (*len); i++) {
                ((static_cast<bool*>(inout))[i]) =
                  ( ((static_cast<bool*>(in))[i]) && ((static_cast<bool*>(inout))[i]));
            }//end for	
        }//end function

    public: 
        /** 
          @brief User defined MPI_Operation to perform: second[i] = (first[i] && second[i]), 
          @remark first and second are 2 arrays of type bool and '&&' represents the 'Logical AND' operation. 
         **/
        static MPI_Op LAND() {
            static bool         first = true;
            static MPI_Op land;
            if (first) {
                first = false;
                MPI_Op_create(Mpi_datatype<bool>::bool_LAND ,true ,&land);
            }

            return land;
        }

        /** 
          @brief User defined MPI_Operation to perform: second[i] = (first[i] || second[i]), 
          @remark first and second are 2 arrays of type bool and '||' represents the 'Logical OR' operation. 
         **/
        static MPI_Op LOR() {
            static bool         first = true;
            static MPI_Op lor;
            if (first) {
                first = false;
                MPI_Op_create(Mpi_datatype<bool>::bool_LOR ,true ,&lor);
            }

            return lor;
        }

	  /** 
          @return the MPI_Datatype for the C++ datatype "bool"
         **/
        static MPI_Datatype value() {
            static bool         first = true;
            static MPI_Datatype datatype;
            if (first) {
                first = false;
                MPI_Type_contiguous(sizeof(bool), MPI_BYTE, &datatype);
                MPI_Type_commit(&datatype);
            }
            return datatype;
        }
    };

} //end namespace

#endif

