/***************************************************************************
 * file:        CoordSpherical.h
 *
 * author:      Christian
 *
 *					 Inspired by NASA World Wind Java.
 ***************************************************************************
 * part of:     framework implementation developed by DLR/KN-FS
 **************************************************************************/


#include <omnetpp.h>
#include "INETDefs.h"
#include "FWMath.h"
#include "CoordSpherical.h"

std::string CoordSpherical::info() const
{
    std::stringstream os;
    os << "(" << latitude << "," << longitude << "," << elevation << " )";
    return os.str();
}

/**
* Tests whether two coordinate vectors are equal. Because
* coordinates are of type floating point, this is done through
* the FWMath::close function.
*/
bool operator==(CoordSpherical& lhs, CoordSpherical& rhs)
{
    // TODO elevation could probably be more coarse grained?
    return FWMath::close(lhs.latitude, rhs.latitude)
            && FWMath::close(lhs.longitude, rhs.longitude)
            && FWMath::close(lhs.elevation, rhs.elevation);
}

/*
float CoordSpherical::sqrdist(const CoordSpherical& a) const
{
    // distance based on great circle
    // Haversine formula
    // ref.: http://www.movable-type.co.uk/scripts/latlong.html
    double d_lat = DEGREES_TO_RADIANS * (a.latitude - this->latitude);
    double d_long = DEGREES_TO_RADIANS * (a.longitude - this->longitude);

    double b = ::sin(d_lat/2) * ::sin(d_lat/2) + ::cos(DEGREES_TO_RADIANS * this->latitude) * ::cos(DEGREES_TO_RADIANS * a.latitude) * sin(d_long/2) * sin(d_long/2);
    double c = EARTH_RADIUS * 2.0 * ::atan2(::sqrt(b),::sqrt(1.0-b));

    // now we use pythagoras to include the elevation into the distance
    return (Angle) (c * c + (this->elevation - a.elevation) * (this->elevation - a.elevation));
}
*/
float CoordSpherical::sqrdist(const CoordSpherical& a) const
{
    double r = EARTH_RADIUS + this->elevation;
    double theta = (90 - this->latitude) * DEGREES_TO_RADIANS;
    double phi = fmod( (this->longitude + 360) , 360 ) * DEGREES_TO_RADIANS;

    double r2 = EARTH_RADIUS + a.elevation;
    double theta2 = (90 - a.latitude) * DEGREES_TO_RADIANS;
    double phi2 = fmod( (a.longitude + 360) , 360 ) * DEGREES_TO_RADIANS;

    return (r*r + r2*r2 - 2*r*r2*( cos(theta)*cos(theta2) + sin(theta)*sin(theta2)*cos(phi-phi2) ));
}

float CoordSpherical::distanceRad(const CoordSpherical& a)
{
    double theta = (this->latitude) * DEGREES_TO_RADIANS;
    double phi = - this->longitude * DEGREES_TO_RADIANS;

    double theta2 = (a.latitude) * DEGREES_TO_RADIANS;
    double phi2 = - a.longitude * DEGREES_TO_RADIANS;

    double A = pow(sin((theta-theta2)/2),2);
    double B = pow(sin((phi-phi2)/2),2);

    return (2 * asin(sqrt(A+cos(theta)*cos(theta2)*B)));
}

void CoordSpherical::convertCoordinate(std::string& coordinate, std::string& coord_lat, std::string& coord_lon)
{
    if (coordinate.length() != 11)
        opp_error("Coordinate invalid - must have length 11!");

    float latitude, longitude;
    char c_lat[15], c_lon[15];

    // Decimal degree DD = dd + mm/60.0 + ss/3600.
    //
    // FIXME
    // make conversion faster by using c like strings: start from back and set '\0's appropriately.
    latitude = (float) ( atof( coordinate.substr(0, 2).c_str() ) + atof( coordinate.substr(2, 2).c_str() ) / 60 );
    longitude = (float) ( atof( coordinate.substr(5, 3).c_str() ) + atof( coordinate.substr(8, 2).c_str() ) / 60 );

    if (coordinate[4] != 'N')
        latitude *= -1.0;

    if (coordinate[10] != 'E')
        longitude *= -1.0;

    // convert numbers to strings
    sprintf(c_lat, "%f", latitude);
    sprintf(c_lon, "%f", longitude);

    coord_lat = c_lat;
    coord_lon = c_lon;
}

void CoordSpherical::convertCoordinate13And14Digit(std::string& coord_lat, std::string& coord_lon)
{
    if ( (coord_lat.length() != 13) && (coord_lon.length() != 14) )
        opp_error("Coordinate invalid - must have length 13 or 14!");

    float latitude, longitude;
    char c_lat[15], c_lon[15];

    // Decimal degree DD = dd + mm/60.0 + ss/3600.
    //
    // FIXME
    // make conversion faster by using c like strings: start from back and set '\0's appropriately.
    latitude = (float) ( atof( coord_lat.substr(2, 2).c_str() ) + atof( coord_lat.substr(5, 2).c_str() ) / 60 + atof( coord_lat.substr(8, 4).c_str() ) / 3600 );
    longitude = (float) ( atof( coord_lon.substr(2, 3).c_str() ) + atof( coord_lon.substr(6, 2).c_str() ) / 60 + atof( coord_lon.substr(9, 4).c_str() ) / 3600 );

    if (coord_lat[0] != 'N')
        latitude *= -1.0;

    if (coord_lon[0] != 'E')
        longitude *= -1.0;

    // convert numbers to strings
    sprintf(c_lat, "%f", latitude);
    sprintf(c_lon, "%f", longitude);

    coord_lat = c_lat;
    coord_lon = c_lon;
}

std::vector<CoordSpherical> CoordSpherical::calculateGreatCircleRoute(const CoordSpherical& destination, int nrPoints)
{
    float step = 1.0 / nrPoints;
    float f = 0; // fraction of distance between start and end point

    float currentLat, currentLon;
    CoordSpherical currentWaypoint;
    std::vector<CoordSpherical> waypoints;
    waypoints.push_back(*this);

    // calculate distance between endpoints in radians
    float d = distanceRad(destination);

    // pre-calculate some values to save time
    float sin_d = sin(d);
    float sin_lat_rad = sin(DEGREES_TO_RADIANS * latitude);
    float cos_lat_rad = cos(DEGREES_TO_RADIANS * latitude);
    float sin_lon_rad = sin(-DEGREES_TO_RADIANS * longitude);
    float cos_lon_rad = cos(-DEGREES_TO_RADIANS * longitude);
    float sin_dest_lon_rad = sin(-DEGREES_TO_RADIANS * destination.longitude); // algorithm assumes that Western longitudes are positive
    float cos_dest_lon_rad = cos(-DEGREES_TO_RADIANS * destination.longitude);
    float sin_dest_lat_rad = sin(DEGREES_TO_RADIANS * destination.latitude);
    float cos_dest_lat_rad = cos(DEGREES_TO_RADIANS * destination.latitude);
    float deltaElev = destination.elevation - elevation;

    for (int i=0; i<nrPoints-1; i++)
    {
        f = f + step;

        float A = sin((1-f)*d)/sin_d;
        float B = sin(f*d)/sin_d;

        float x = A*cos_lat_rad*cos_lon_rad +  B*cos_dest_lat_rad*cos_dest_lon_rad;
        float y = A*cos_lat_rad*sin_lon_rad +  B*cos_dest_lat_rad*sin_dest_lon_rad;
        float z = A*sin_lat_rad		       	 +  B*sin_dest_lat_rad;

        currentLat = atan2(z,sqrt(pow(x,2)+pow(y,2)));
        currentLon = atan2(y,x);

        // convert back to degrees
        currentLat = RADIANS_TO_DEGREES * currentLat;
        currentLon = -RADIANS_TO_DEGREES * currentLon;

        currentWaypoint.latitude = currentLat;
        currentWaypoint.longitude = currentLon;
        currentWaypoint.elevation = elevation + f*deltaElev;
        waypoints.push_back(currentWaypoint);
    }

    waypoints.push_back(destination);

    return waypoints;
}

