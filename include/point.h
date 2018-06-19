/**
  @file Point.h
  @brief A point class
  @author Hari Sundar
  */

/***************************************************************************
 *   Copyright (C) 2005 by Hari sundar   *
 *   hsundar@seas.upenn.edu   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
***************************************************************************/
#ifndef __POINT_H
#define __POINT_H

#include  <cmath>

/**
  @brief A point class
  @author Hari Sundar
  */
class Point{
  public:

    /** @name Constructors and Destructor */
    //@{
    Point();
    // virtual ~Point();

    Point(double newx, double newy, double newz);
    Point(int newx, int newy, int newz);
    Point(unsigned int newx, unsigned int newy, unsigned int newz);
    Point(const Point &newpoint);
    //@}

    /** @name Getters */
    //@{
    const double& x() const {return _x; };
    const double& y() const {return _y; };
    const double& z() const {return _z; };

    int xint() const {return static_cast<int>(_x); };
    int yint() const {return static_cast<int>(_y); };
    int zint() const {return static_cast<int>(_z); };
    //@}

    /** @name Overloaded Operators */
    //@{
    Point operator-() const;

    void operator += (const Point &other);
    void operator /= (const int divisor);
    void operator /= (const double divisor);
    void operator *= (const int factor);
    void operator *= (const double factor);

    Point& operator=(const Point &other);
    Point  operator+(const Point &other);
    Point  operator-(const Point &other);
    Point  operator-(const Point &other) const;

    Point  operator/(const double divisor);
    Point  operator*(const double factor);
    
    double magnitude();

    inline bool operator != (const Point &other) { return ( (xint() != other.xint() ) || (yint() != other.yint()) || (zint() != other.zint())); };

    inline bool operator == (const Point &other) { return ( ( xint() == other.xint() ) && ( yint() == other.yint()) && ( zint() == other.zint())); };
    //@}

    inline double dot(Point Other) { 
      return  (_x*Other._x+_y*Other._y+_z*Other._z);
    };				

    inline Point cross(Point  Other){
      return  Point(_y*Other._z-Other._y*_z, _z*Other._x-_x*Other._z, 
          _x*Other._y-_y*Other._x); 
    };

    inline double abs(){
      return sqrt(_x*_x+_y*_y+_z*_z); 
    };

    void normalize();

    static Point TransMatMultiply( double* transMat, Point inPoint);
  protected:
    double _x;
    double _y;
    double _z;
};

#endif // POINT_H
