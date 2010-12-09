/***************************************************************************
 * file:        CoordSpherical.h
 *
 * author:      Christian
 *
 *					 Inspired by NASA World Wind Java.
 ***************************************************************************
 * part of:     framework implementation developed by DLR/KN-FS
 **************************************************************************/


#ifndef _COORDSPHERICAL_H
#define _COORDSPHERICAL_H

#include <omnetpp.h>
#include "INETDefs.h"
#include "FWMath.h"
#include "CoordSpherical.h"

#ifndef M_PI
#define M_PI	3.141592653589793238462643383279502
#endif


#define	Angle	double
// earth radius in meters
#define	EARTH_RADIUS			6378137
// constants for calculations between degrees and radians
#define	DEGREES_TO_RADIANS	(Angle)(M_PI / 180)
#define	RADIANS_TO_DEGREES	(Angle)(180 / M_PI)
#define	PIOver2					(Angle)(M_PI / 2)
// degree constants
#define	ZERO						   0
#define	POS90						  90
#define	NEG90						 -90
#define	POS180					 180
#define	NEG180					-180
#define	POS360					 360


/**
 * @brief Class for storing host positions
 *
 * Class for a triple storing (latitude, longitude, elevation).
 * Some comparison and basic arithmetic operators on CoordSpherical
 * structures are implemented.
 *
 * @ingroup support
 * @author Christian
 */
class INET_API CoordSpherical : public cPolymorphic
{
public:
  /** @brief latitude, longitude and elevation of the position */
  Angle latitude;
  Angle longitude;
  double elevation;

public:
	CoordSpherical()
	{
		latitude;
		longitude;
		elevation = 0;
	}

	/** Initializes coordinates.*/
	CoordSpherical(Angle _latitude, Angle _longitude, float _elevation)
            : latitude(_latitude), longitude(_longitude), elevation(_elevation) {};

	/** Copy constructor: initializes coordinates.*/
	CoordSpherical(const CoordSpherical& pos)
	{
		latitude = pos.latitude;
		longitude = pos.longitude;
		elevation = pos.elevation;
	}

	std::string info() const;

	static CoordSpherical getZeroCoord()
	{
		return CoordSpherical(ZERO, ZERO, 0);
	}

	static CoordSpherical fromRadians(Angle latitude, Angle longitude, float elevation)
	{
		return CoordSpherical(RADIANS_TO_DEGREES * latitude, RADIANS_TO_DEGREES * longitude, elevation);
	}

	static CoordSpherical fromDegrees(Angle latitude, Angle longitude, float elevation)
	{
		return CoordSpherical(latitude, longitude, elevation);
	}

/*
	---------
	Operators
	---------
*/

	/** Adds two coordinate vectors.*/
	friend CoordSpherical operator+(CoordSpherical a, CoordSpherical b)
	{	
		  return CoordSpherical(a.latitude+b.latitude, a.longitude+b.longitude, a.elevation+b.elevation);
	}

	/** Subtracts two coordinate vectors.*/
	friend CoordSpherical operator-(CoordSpherical a, CoordSpherical b)
	{
		  return CoordSpherical(a.latitude-b.latitude, a.longitude-b.longitude, a.elevation-b.elevation);
	}

	/** Adds coordinate vector b to a.*/
	CoordSpherical operator+=(CoordSpherical a)
	{
		  latitude+=a.latitude;
		  longitude+=a.longitude;
		  elevation+=a.elevation;
		  return *this;
	}

	/** subtracts coordinate vector b from a.*/
	CoordSpherical operator-=(CoordSpherical a)
	{
		  latitude-=a.latitude;
		  longitude-=a.longitude;
		  elevation-=a.elevation;
		  return *this;
	}

	/** Subtracts two coordinate vectors.*/
//	friend Coord operator-(Coord a, Coord b) {
//        return Coord(a.x-b.x, a.y-b.y);
//  }

  /** Multiplies a coordinate vector by a real number.*/
	friend CoordSpherical operator*(CoordSpherical a, double f) {
		return CoordSpherical(a.latitude*f, a.longitude*f, a.elevation*f);
	}
  /** Divides a coordinate vector by a real number.*/
	friend CoordSpherical operator/(CoordSpherical a, double f) {
		return CoordSpherical(a.latitude/f, a.longitude/f, a.elevation/f);
	}

	/** Copy rhs to this. */
	CoordSpherical operator=(const CoordSpherical& rhs)
	{
		if (this != &rhs)
		{
			latitude = rhs.latitude;
			longitude = rhs.longitude;
			elevation = rhs.elevation;
		}
		return *this;
	}

	/**
   * Tests whether two coordinate vectors are equal. Because
   * coordinates are of type floating point, this is done through
	* the FWMath::close function.
   */
	friend bool operator==(CoordSpherical& lhs, CoordSpherical& rhs);

  /**
   * Tests whether two coordinate vectors are not equal. Negation of
   * the operator==.
   */
	friend bool operator!=(CoordSpherical& lhs, CoordSpherical& rhs)
	{
		return !(lhs==rhs);
	}

	/**
	 * Returns distance^2 to CoordSpherical a
	 */
	float sqrdist(const CoordSpherical& a) const;

    /**
     * Returns the distance to CoordSpherical a
     */
    float distance(const CoordSpherical& a) const
    {
        return (Angle) ::sqrt( sqrdist(a) );
    }

    /**
     * Returns the distance to CoordSpherical a in radians
     * Based on algorithm at http://williams.best.vwh.net/avform.htm#Crs
     * Note that differences in elevation are not considered!
     */
    float distanceRad(const CoordSpherical& a);

	/**
	 * Parses the content of the passed radian coordinate and writes the converted
	 * decimal lat/lon values to the respective passed arguments.
	 */
	static void convertCoordinate(std::string& coordinate, std::string& coord_lat, std::string& coord_lon);

	/**
	 * Parses the content of the passed radian coordinate and writes the converted
	 * decimal lat/lon values back to the passed arguments.
	 */
	static void convertCoordinate13And14Digit(std::string& coord_lat, std::string& coord_lon);

	friend inline std::ostream& operator<<(std::ostream& os, const CoordSpherical& coord)
	{
		return os << "(" << coord.latitude << "," << coord.longitude << "," << coord.elevation << ")";
	}

	/**
	* calculates a sequence of nrPoints equally spaced waypoints along the Great Circle route to destination
	* implementation of algorithm found at http://williams.best.vwh.net/avform.htm#Crs
	* differences in elevation between the endpoints are interpolated linearly, i.e. the calculated route is a Great Circle along the surface of the Earth!
	*/
	std::vector<CoordSpherical> calculateGreatCircleRoute(const CoordSpherical& destination, int nrPoints);
};

#endif
