//
// Created by milinda on 11/27/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief simplified templated  point class.
 *
*/
//

#ifndef SFCSORTBENCH_POINTT_H
#define SFCSORTBENCH_POINTT_H

template <class T>
class PointT
{

    public:
    /**@brief default constructor * */
    PointT() : _x((T)0), _y((T)0), _z((T)0){}

    /**@brief constructor with specific coords **/
    PointT(T x,T y, T z):_x(x), _y(y), _z(z){}

    /**@brief get x coord*/
    const T& x() const {return _x; };
    /**@brief get y coord*/
    const T& y() const {return _y; };
    /**@brief get z coord*/
    const T& z() const {return _z; };


    private:
        T _x;
        T _y;
        T _z;
};

#endif //SFCSORTBENCH_POINTT_H
