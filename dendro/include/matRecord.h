//
// Created by milinda on 12/20/18.
//

/**
 * @brief Matrix record class to be used while setting values in the global matrix.
 * This is based from the old dendro MatRecord class, where now it supports  both petsc and non petsc matrices.
 *
 * @author Milinda Fernando, School of Computing, University of Utah.
 * */

#ifdef BUILD_WITH_PETSC
    #include "petsc.h"
#endif

#include <iostream>
#include <ostream>
#include "dendro.h"

#ifndef DENDRO_5_0_MATRECORD_H
#define DENDRO_5_0_MATRECORD_H

namespace ot
{

    class MatRecord {

        private:
            /**@brief local index of the node coresponding to the local row of the matrix */
            unsigned int m_uiRowID;
            /**@brief local index of the node corresponding to the local column of the matrix */
            unsigned int m_uiColID;
            /**@brief The degree of freedom (0-based) of the node corresponding to  the row of the matrix.*/
            unsigned int m_uiRowDim;
            /**@brief The degree of freedom (0-based) of the node corresponding to  the column of the matrix. */
            unsigned int m_uiColDim;

#ifdef BUILD_WITH_PETSC
            /**@brief  matrix value (petsc)*/
            PetscScalar m_uiVal;
#else
            /**@brief  matrix value (dendro)*/
            DendroScalar m_uiVal;
#endif

        public:

            /**
             *@brief The default constructor
             */
            MatRecord() {
                m_uiRowID = m_uiColID = m_uiRowDim = m_uiColDim = static_cast<unsigned int>(-1);
                m_uiVal = 0;
            }

            /**@brief returns the row ID*/
            inline unsigned int getRowID() const {return m_uiRowID;}
            /**@brief returns the col ID*/
            inline unsigned int getColID() const {return m_uiColID;}
            /**@brief return the row dof*/
            inline unsigned int getRowDim() const {return m_uiRowDim;}
            /**@brief return the col dof*/
            inline unsigned int getColDim() const {return m_uiColDim;}

            /**@brief sets the rowID value*/
            inline void setRowID(unsigned int rowID) {m_uiRowID=rowID;}

            /**@brief sets the rowID value*/
            inline void setColID(unsigned int colID) {m_uiColID=colID;}

            /**@brief sets the rowDim value*/
            inline void setRowDim(unsigned int rowDim) {m_uiRowDim=rowDim;}

            /**@brief sets the colDim value*/
            inline void setColDim(unsigned int colDim) {m_uiColDim=colDim;}


#ifdef BUILD_WITH_PETSC

            /**@breif MatRecord constructor
             * @param [in] rowID: row id
             * @param [in] colID: col id
             * @param [in] rowDim: row dim
             * @param [in] colDim: colDim
             * @param [in] value: value of for the entry.
             * */

            MatRecord(unsigned int rowID, unsigned int colID, unsigned int rowDim, unsigned int colDim, PetscScalar value)
            {
                m_uiRowID=rowID;
                m_uiColID=colID;

                m_uiRowDim=rowDim;
                m_uiColDim=colDim;

                m_uiVal=value;

            }



            /**@brief returns the entry value*/
            inline PetscScalar getMatVal() const {return m_uiVal;}

            /**@brief sets matrix value*/
            inline void setMatValue(PetscScalar value) {m_uiVal=value;}

#else

            /**@breif MatRecord constructor
             * @param [in] rowID: row id
             * @param [in] colID: col id
             * @param [in] rowDim: row dim
             * @param [in] colDim: colDim
             * @param [in] value: value of for the entry.
             * */

            MatRecord(unsigned int rowID, unsigned int colID, unsigned int rowDim, unsigned int colDim, DendroScalar value)
            {
                m_uiRowID=rowID;
                m_uiColID=colID;

                m_uiRowDim=rowDim;
                m_uiColDim=colDim;

                m_uiVal=value;

            }



            /**@brief returns the entry value*/
            inline DendroScalar getMatVal() const {return m_uiVal;}

            /**@brief sets matrix value*/
            inline void setMatValue(DendroScalar value) {m_uiVal=value;}

#endif

            /**@brief The copy constructor */
            MatRecord(const MatRecord & other) {
                m_uiRowID = other.getRowID();
                m_uiColID = other.getColID();
                m_uiRowDim = other.getRowDim();
                m_uiColDim = other.getColDim();
                m_uiVal = other.getMatVal();
            }


            /** @brief The assignment operator */
            MatRecord & operator = (MatRecord const  & other) {
                if(this == (&other)) {return *this;}
                m_uiRowID = other.getRowID();
                m_uiColID = other.getColID();
                m_uiRowDim = other.getRowDim();
                m_uiColDim = other.getColDim();
                m_uiVal = other.getMatVal();
                return *this;
            }

            /** @brief Overloaded == Operator */
            bool  operator == ( MatRecord const  &other) const {
                return ( (m_uiRowID == other.getRowID()) && (m_uiColID == other.getColID())
                         && (m_uiVal == other.getMatVal()) && (m_uiRowDim == other.getRowDim()) && (m_uiColDim == other.getColDim()) );
            }

            /** @brief Overloaded != Operator */
            bool  operator != (MatRecord const  &other) const {
                return (!((*this) == other));
            }

            /** @brief Overloaded < Operator */
            bool  operator < (MatRecord const  &other) const {
                if(m_uiRowID < other.getRowID()) {
                    return true;
                }else if(m_uiRowID == other.getRowID()) {
                    if(m_uiRowDim < other.getRowDim()) {
                        return true;
                    }else if(m_uiRowDim == other.getRowDim()) {
                        if(m_uiColID < other.getColID()) {
                            return true;
                        }else if (m_uiColID == other.getColID()) {
                            if(m_uiColDim < other.getColDim()) {
                                return true;
                            }else if(m_uiColDim == other.getColDim()) {
                                return (m_uiVal < other.getMatVal());
                            }else {
                                return false;
                            }
                        }else {
                            return false;
                        }
                    }else {
                        return false;
                    }
                }else {
                    return false;
                }
            }

            /** @brief Overloaded > Operator */
            bool  operator > (MatRecord const  &other) const {
                return ( (!((*this) < other)) && ((*this) != other) );
            }

            /** @brief Overloaded <= Operator */
            bool  operator <= (MatRecord const  &other) const {
                return ( ((*this) < other) || ((*this) == other) );
            }

            /** @brief Overloaded >= Operator */
            bool  operator >= (MatRecord const  &other) const {
                return (!((*this) < other)) ;
            }

            friend std::ostream & operator<< (std::ostream & os, MatRecord const & re) 
            {
                return (os << " row : "<<re.getRowID() << " col: "<<re.getColID()<<" rdim: "<<re.getRowDim()<<" cdim: "<<re.getColDim()<< "val: "<<re.getMatVal() );
            }

            

    };


}


namespace par
{
    template<typename T>
    class Mpi_datatype;

    template <>
    class Mpi_datatype< ot::MatRecord > 
    {
        public:
            static MPI_Datatype value()
            {
                static bool         first = true;
                static MPI_Datatype datatype;

                if (first)
                {
                    first = false;
                    MPI_Type_contiguous(sizeof(ot::MatRecord), MPI_BYTE, &datatype);
                    MPI_Type_commit(&datatype);
                }

                return datatype;
            }
    };

}


#endif //DENDRO_5_0_MATRECORD_H
