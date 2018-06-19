//
// Created by milinda on 12/25/16.
//

/**
 * @author Milinda Fernando
 * @author Hari Sundar
 * @breif Contains data structures to store the reference element information.
 *
 * @refference: Based of HOMG code written in matlab.
 * */
#ifndef SFCSORTBENCH_REFERENCEELEMENT_H
#define SFCSORTBENCH_REFERENCEELEMENT_H

#ifdef WITH_BLAS_LAPACK
    #include "basis.h"
#else
    #include "interpMatrices.h"
#endif

#include "tensor.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
template <typename T>
void printArray_1D(const T *a, int length)
{
    for (int i = 0; i < length; i++) { std::cout<<a[i]<<" "; }
    std::cout<<std::endl;
}


template <typename T>
void printArray_2D(const T *a, int length1,int length2)
{
    for (int i = 0; i < length1; i++) {
        for (int j = 0; j < length2; j++) {
            std::cout << a[i * length2 + j] << " ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}



class RefElement{


private :
    /** Dimension */
    int                 m_uiDimension;
    /** Polynomial Order */
    int                 m_uiOrder;
    /** Number of 3D interpolation points on the element */
    int                 m_uiNp;
    /** Number of 2D face interpolation points */
    int                 m_uiNfp;
    /** Number of 1D interpolation points */
    int                 m_uiNrp;

    /** reference element volume */
    unsigned int m_uiVol;

    /** 1D reference coordinates of the interpolation nodes (uniform nodal points) */
    std::vector<double> u;

    /** 1D reference coordinates of the interpolation nodes (gll points)*/
    std::vector<double>       r;

    /** 1D regular points corresponding to child 0 of u*/
    std::vector<double> u_0;

    /** 1D regular points corresponding to child 1 of u*/
    std::vector<double> u_1;

    /** 1D Gauss points (used for quadrature)*/
    std::vector<double>       g;

    /** 1D weights for gauss quadrature */
    std::vector<double> w;

    /** 1D weights for gll quadrature*/
    std::vector<double> wgll;


    /** 1D interpolation matrix for child 0 */
    std::vector<double> ip_1D_0;

    /** 1D interpolation matrix for child 1*/
    std::vector<double> ip_1D_1;


    /** 1D interpolation matrix for child 0 (transpose) */
    std::vector<double> ipT_1D_0;

    /** 1D interpolation matrix for child 1 (transpose)*/
    std::vector<double> ipT_1D_1;

    /**Vandermonde matrix for interpolation points r.   */
    std::vector<double> Vr;

    /**Vandermonde matrix for interpolation points u.   */
    std::vector<double> Vu;

    /**Vandermonde matrix for polynomial at gauss points */
    std::vector<double> Vg;

    /**gradient of the vandermonde for polynomial eval at points u*/
    std::vector<double> gradVu;

    /**gradient of the vandermonde for polynomial eval at points r*/
    std::vector<double> gradVr;

    /**gradient of the vandermonde for polynomial eval at points g*/
    std::vector<double> gradVg;

    /**derivative of the pol. eval at points r. */
    std::vector<double> Dr;

    /** derivative of the pol. eval at the gauss points. */
    std::vector<double> Dg;

    /** derivative of the pol. eval at the gauss points. (transpose) */
    std::vector<double> DgT;

    /** 1D quadrature matrix*/
    std::vector<double> quad_1D;

    /** 1D quadrature matrix transpose*/
    std::vector<double> quadT_1D;

    /**Vandermonde matrix for interpolation points of child 0   */
    std::vector<double> Vu_0;

    /**Vandermonde matrix for interpolation points of child 1   */
    std::vector<double> Vu_1;

    /**intermidiate vec 1 needed during interploation */
    std::vector<double> im_vec1;

    /**intermidiate vec 1 needed during interploation */
    std::vector<double> im_vec2;



public:
    RefElement();
    RefElement(unsigned int dim, unsigned int order);
    ~RefElement();
    // some getter methods to access required data.
    inline int getOrder() const {return m_uiOrder;}
    inline int getDim() const {return m_uiDimension;}
    inline int get1DNumInterpolationPoints(){return m_uiNrp;}

    inline const double * getIMChild0() const {return &(*(ip_1D_0.begin()));}
    inline const double * getIMChild1() const {return &(*(ip_1D_1.begin()));}

    inline const double * getQ1d() const {return &(*(quad_1D.begin()));}
    inline const double * getQT1d()const {return &(*(quadT_1D.begin()));}
    inline const double * getDg1d()const {return &(*(Dg.begin()));}
    inline const double * getDgT1d()const {return &(*(DgT.begin()));}
    inline const double * getDr1d()const {return &(*(Dr.begin()));}

    inline double * getImVec1() {return &(*(im_vec1.begin()));}
    inline double * getImVec2() {return &(*(im_vec2.begin()));}

    inline const double * getWgq()const {return &(*(w.begin()));}
    inline const double * getWgll()const {return &(*(wgll.begin()));}

    inline const double getElementSz()const {return (r.back()-r.front());}



    /**
     * @param[in] in: input function values.
     * @param[in] childNum: Morton ID of the child number where the interpolation is needed.
     * @param[out] out: interpolated values.
     *
     * @brief This is computed in way that 3d coordinates changes in the order of z, y, x
     * Which means first we need to fill all the z values in plane(x=0,y=0) then all the z values in plane (x=0,y=0+h) and so forth.
     *
     */
    inline void I3D_Parent2Child(const double * in, double* out, unsigned int childNum ) const
    {

         double * im1=(double *)&(*(im_vec1.begin()));
         double * im2=(double *)&(*(im_vec2.begin()));

        switch (childNum)
        {
            case 0:
                DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_0.begin())),in,im1); // along x
                DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_0.begin())),im1,im2); // along y
                DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_0.begin())),im2,out); // along z
                break;
            case 1:
                DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_1.begin())),in,im1); // along x
                DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_0.begin())),im1,im2); // along y
                DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_0.begin())),im2,out); // along z

                break;
            case 2:
                DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_0.begin())),in,im1);  // along x
                DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_1.begin())),im1,im2); // along y
                DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_0.begin())),im2,out); // along z
                break;
            case 3:
                DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_1.begin())),in,im1);  // along x
                DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_1.begin())),im1,im2); // along y
                DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_0.begin())),im2,out); // along z
                break;
            case 4:
                DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_0.begin())),in,im1);  // along x
                DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_0.begin())),im1,im2); // along y
                DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_1.begin())),im2,out); // along z
                break;
            case 5:
                DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_1.begin())),in,im1); // along x
                DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_0.begin())),im1,im2); // along y
                DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_1.begin())),im2,out); // along z
                break;
            case 6:
                DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_0.begin())),in,im1); // along x
                DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_1.begin())),im1,im2); // along y
                DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_1.begin())),im2,out); // along z
                break;
            case 7:
                DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_1.begin())),in,im1); // along x
                DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_1.begin())),im1,im2); // along y
                DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp,&(*(ip_1D_1.begin())),im2,out); // along z
                break;
            default:
                std::cout<<"[refel][error]: invalid child number specified for 3D interpolation."<<std::endl;
                break;

        }



    }

    /**
    * @param[in] in: input function values.
    * @param[in] childNum: Morton ID of the child number where the contribution needed to be computed.
    * @param[out] out: child to parent contribution values. (used in FEM integral ivaluation)
    *
    * @brief This is computed in way that 3d coordinates changes in the order of z, y, x
    * Which means first we need to fill all the z values in plane(x=0,y=0) then all the z values in plane (x=0,y=0+h) and so forth.
    *
    */

    inline void I3D_Child2Parent(const double * in, double* out, unsigned int childNum ) const
    {
        double * im1=(double *)&(*(im_vec1.begin()));
        double * im2=(double *)&(*(im_vec2.begin()));

        switch (childNum)
        {
            case 0:
                DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_0.begin())),in,im1); // along x
                DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_0.begin())),im1,im2); // along y
                DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_0.begin())),im2,out); // along z
                break;
            case 1:
                DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_1.begin())),in,im1); // along x
                DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_0.begin())),im1,im2); // along y
                DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_0.begin())),im2,out); // along z

                break;
            case 2:
                DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_0.begin())),in,im1);  // along x
                DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_1.begin())),im1,im2); // along y
                DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_0.begin())),im2,out); // along z
                break;
            case 3:
                DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_1.begin())),in,im1);  // along x
                DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_1.begin())),im1,im2); // along y
                DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_0.begin())),im2,out); // along z
                break;
            case 4:
                DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_0.begin())),in,im1);  // along x
                DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_0.begin())),im1,im2); // along y
                DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_1.begin())),im2,out); // along z
                break;
            case 5:
                DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_1.begin())),in,im1); // along x
                DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_0.begin())),im1,im2); // along y
                DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_1.begin())),im2,out); // along z
                break;
            case 6:
                DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_0.begin())),in,im1); // along x
                DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_1.begin())),im1,im2); // along y
                DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_1.begin())),im2,out); // along z
                break;
            case 7:
                DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_1.begin())),in,im1); // along x
                DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_1.begin())),im1,im2); // along y
                DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp,&(*(ipT_1D_1.begin())),im2,out); // along z
                break;
            default:
                std::cout<<"[refel][error]: invalid child number specified for 3D interpolation."<<std::endl;
                break;

        }
    }


    /**
     * @param[in] in: input function values.
     * @param[in] childNum: Morton ID of the child number where the interpolation is needed.
     * @param[out] out: interpolated values.
     * */

    inline void I2D_Parent2Child(const double * in, double* out, unsigned int childNum) const
    {


         double * im1=(double *)&(*(im_vec1.begin()));
         double * im2=(double *)&(*(im_vec2.begin()));



        switch (childNum)
        {

            case 0:
                DENDRO_TENSOR_IAX_APPLY_ELEM_2D(m_uiNrp,&(*(ip_1D_0.begin())),in,im1); // along x
                DENDRO_TENSOR_AIX_APPLY_ELEM_2D(m_uiNrp,&(*(ip_1D_0.begin())),im1,out); // along y (in 3d z)
                    break;
            case 1:
                DENDRO_TENSOR_IAX_APPLY_ELEM_2D(m_uiNrp,&(*(ip_1D_1.begin())),in,im1); // along x
                DENDRO_TENSOR_AIX_APPLY_ELEM_2D(m_uiNrp,&(*(ip_1D_0.begin())),im1,out); // along y (in 3d z)
                    break;

            case 2:
                DENDRO_TENSOR_IAX_APPLY_ELEM_2D(m_uiNrp,&(*(ip_1D_0.begin())),in,im1); // along x
                DENDRO_TENSOR_AIX_APPLY_ELEM_2D(m_uiNrp,&(*(ip_1D_1.begin())),im1,out); // along y (in 3d z)
                    break;

            case 3:
                DENDRO_TENSOR_IAX_APPLY_ELEM_2D(m_uiNrp,&(*(ip_1D_1.begin())),in,im1); // along x
                DENDRO_TENSOR_AIX_APPLY_ELEM_2D(m_uiNrp,&(*(ip_1D_1.begin())),im1,out); // along y (in 3d z)
                    break;
            default:
                std::cout<<"[refel][error]: invalid child number specified for 2D  interpolation."<<std::endl;
                break;

        }


    }

    /**
     * @param[in] in: input function values.
     * @param[in] childNum: Morton ID of the child number where the interpolation is needed.
     * @param[out] out: child to parent contribution values. (used in FEM integral ivaluation)
     * */

    inline void I2D_Child2Parent(const double * in, double* out, unsigned int childNum) const
    {


        double * im1=(double *)&(*(im_vec1.begin()));
        double * im2=(double *)&(*(im_vec2.begin()));

        switch (childNum)
        {

            case 0:
                DENDRO_TENSOR_IAX_APPLY_ELEM_2D(m_uiNrp,&(*(ipT_1D_0.begin())),in,im1); // along x
                DENDRO_TENSOR_AIX_APPLY_ELEM_2D(m_uiNrp,&(*(ipT_1D_0.begin())),im1,out); // along y (in 3d z)
                break;
            case 1:
                DENDRO_TENSOR_IAX_APPLY_ELEM_2D(m_uiNrp,&(*(ipT_1D_1.begin())),in,im1); // along x
                DENDRO_TENSOR_AIX_APPLY_ELEM_2D(m_uiNrp,&(*(ipT_1D_0.begin())),im1,out); // along y (in 3d z)
                break;

            case 2:
                DENDRO_TENSOR_IAX_APPLY_ELEM_2D(m_uiNrp,&(*(ipT_1D_0.begin())),in,im1); // along x
                DENDRO_TENSOR_AIX_APPLY_ELEM_2D(m_uiNrp,&(*(ipT_1D_1.begin())),im1,out); // along y (in 3d z)
                break;

            case 3:
                DENDRO_TENSOR_IAX_APPLY_ELEM_2D(m_uiNrp,&(*(ipT_1D_1.begin())),in,im1); // along x
                DENDRO_TENSOR_AIX_APPLY_ELEM_2D(m_uiNrp,&(*(ipT_1D_1.begin())),im1,out); // along y (in 3d z)
                break;
            default:
                std::cout<<"[refel][error]: invalid child number specified for 2D  interpolation."<<std::endl;
                break;

        }


    }



    /**
     * @param [in] in input function values
     * @param [in] childNum Morton ID of the child number
     * @param [out] interpolated values from parent to child.
     * */

    inline void I1D_Parent2Child(const double * in,double * out,unsigned int childNUm) const
    {

        switch (childNUm)
        {
            case 0:
                    for(unsigned int i=0;i<m_uiNrp;i++)
                    {
                        out[i]=0.0;
                        for(unsigned int j=0;j<m_uiNrp;j++)
                        {
                            out[i]+=ip_1D_0[j*m_uiNrp+i]*in[j];
                        }
                    }
                    break;
            case 1:
                for(unsigned int i=0;i<m_uiNrp;i++)
                {
                    out[i]=0.0;
                    for(unsigned int j=0;j<m_uiNrp;j++)
                    {
                        out[i]+=ip_1D_1[j*m_uiNrp+i]*in[j];
                    }
                }
                break;

            default:
                std::cout<<"[refel][error]: Invalid child number specified for 1D interpolation. "<<std::endl;
                break;

        }

    }


    /**
     * @param [in] in input function values
     * @param [in] childNum Morton ID of the child number
     * @param [out] child to parent contribution values. (used in FEM integral ivaluation)
     * */

    inline void I1D_Child2Parent(const double * in,double * out,unsigned int childNUm) const
    {

        switch (childNUm)
        {
            case 0:
                for(unsigned int i=0;i<m_uiNrp;i++)
                {
                    out[i]=0.0;
                    for(unsigned int j=0;j<m_uiNrp;j++)
                    {
                        out[i]+=ipT_1D_0[j*m_uiNrp+i]*in[j];
                    }
                }
                break;
            case 1:
                for(unsigned int i=0;i<m_uiNrp;i++)
                {
                    out[i]=0.0;
                    for(unsigned int j=0;j<m_uiNrp;j++)
                    {
                        out[i]+=ipT_1D_1[j*m_uiNrp+i]*in[j];
                    }
                }
                break;

            default:
                std::cout<<"[refel][error]: Invalid child number specified for 1D interpolation. "<<std::endl;
                break;

        }

    }


    void generateHeaderFile(char * fName);



};

#endif //SFCSORTBENCH_REFERENCEELEMENT_H
