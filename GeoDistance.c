//
//  Distance.c
//  Fencing
//
//  Created by Cortex on 5/10/14.
//  Copyright (c) 2014 Zaheer Naby. All rights reserved.
//

#include <stdio.h>
#include "math.h"




inline double radians(const double degrees) {
    return degrees * (M_PI/180.0);
}


double haversineDistance(const double latitude1, const double longitude1, const double latitude2, const double longitude2)
{
    static double earthRadius = 6371;
    
    double deltaLatitude = radians(latitude2-latitude1);
    double deltaLongitude = radians(longitude2-longitude1);
 
    double angle = pow(sin(deltaLatitude/2.0f), 2.0f) + pow(sin(deltaLongitude/2),2.0f) * cos(latitude1) * cos(latitude2);
    double arc  = 2 * atan2(sqrt(angle), sqrt(1-angle));
    
    double distance =  earthRadius * arc;
    return distance*1000;
}

double vincentyDistance(const double latitude1, const double longitude1, const double latitude2, const double longitude2)
{

    const double majorSemiAxis = 6378137.0;
    const double minorSemiAxis = 6356752.314245;
    const double flattening = 1.0/298.257223563;

    const double L = radians(longitude2-longitude1);
    const double U1 = atan((1.0-flattening) * tan(radians(latitude1)));
    const double U2 = atan((1.0-flattening) * tan(radians(latitude2)));
    
    const double sinU1 = sin(U1);
    const double cosU1 = cos(U1);
    const double sinU2 = sin(U2);
    const double cosU2 = cos(U2);
    
    double cosSqAlpha = 0.0f;
    double sinSigma = 0.0f;
    double cosSigma = 0.0f;
    double cos2SigmaM = 0.0f;
    double sigma = 0.0f;
    double lambda = L;
    double lambdaP = 0;
    double iterLimit = 100.0;
    
    while(fabs(lambda-lambdaP) > 1e-12 && iterLimit > 0)
    {

        const double sinLambda = sin(lambda);
        const double cosLambda = cos(lambda);
    
        sinSigma = (sqrt((cosU2*sinLambda) * (cosU2*sinLambda) +
                         (cosU1*sinU2-sinU1*cosU2*cosLambda) *
                         (cosU1*sinU2-sinU1*cosU2*cosLambda)));
        
        if (sinSigma==0)
        {
            return 0; //co-incident points
        }
        
        cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda;
        sigma = atan2(sinSigma, cosSigma);
        
        const double sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
        cosSqAlpha = 1 - sinAlpha*sinAlpha;
        cos2SigmaM = cosSigma - 2*sinU1*sinU2/cosSqAlpha;
        
        if (isnan(cos2SigmaM))
        {
            cos2SigmaM = 0; //equatorial line: cosSqAlpha=0 (6)
        }
        
        
        const double C = flattening/16*cosSqAlpha*(4+flattening*(4-3*cosSqAlpha));
        lambdaP = lambda;
        lambda = (L + (1-C) * flattening * sinAlpha *
                (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM))));
        
        if (--iterLimit==0)
        {
            return -1; //formula failed to converge
        }
    }
    
    const double uSq = cosSqAlpha * (majorSemiAxis*majorSemiAxis - minorSemiAxis*minorSemiAxis) / (minorSemiAxis*minorSemiAxis);
    const double A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
    const double B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));
    const double deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)-
                                                          B/6*cos2SigmaM*(-3+4*sinSigma*sinSigma)*(-3+4*cos2SigmaM*cos2SigmaM)));
    const double s = minorSemiAxis*A*(sigma-deltaSigma);
    return s;
}
