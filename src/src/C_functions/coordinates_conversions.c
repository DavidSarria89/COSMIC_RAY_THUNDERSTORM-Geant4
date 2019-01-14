// *******************************************************************************

/* GEANT4 code for propagation of gamma-rays, electron and positrons in Earth's atmosphere */

//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************

/* ------------------------------------------------------------------- */
/* ------------------------------ INCLUDES --------------------------- */
/* ------------------------------------------------------------------- */

#include "coordinates_conversions.h"

/* ------------------------------------------------------------------- */
/* ------------------------------ FUCNTIONS -------------------------- */
/* ------------------------------------------------------------------- */

// All distances are in km
// All angles are in degrees

/* ------------------------------------------------------------------- */
void
geocentric_to_ecef(struct geocentric_coords *input, struct ecef_coords *output)
{
    static const double aa = 6378.137; // km,  WGS-84
    static const double pi = 3.14159265359;

    double sinlat, coslat, sinlong, coslong;

    double R = aa + input->alt;

    double lat_rad = (input->lat) * pi / 180.; // conversion to radians
    double long_rad = (input->lon) * pi / 180.; // conversion to radians

    sinlat = sin(lat_rad);
    coslat = cos(lat_rad);

    sinlong = sin(long_rad);
    coslong = cos(long_rad);

    //    sincos(lat_rad,  &sinlat,  &coslat);
    //    sincos(long_rad, &sinlong, &coslong);

    output->x = R * coslong * coslat;
    output->y = R * sinlong * coslat;
    output->z = R * sinlat;
}

/* ------------------------------------------------------------------- */
void
geodetic_to_ecef(struct geodetic_coords *input, struct ecef_coords *output)
{
    static const double aa = 6378.137;          // km,  WGS-84
    static const double e2 = 0.006694379990141; //  WGS-84
    static const double pi = 3.14159265359;

    double sinlat, coslat, sinlong, coslong;
    double localVertical[3];

    double lat_rad = (input->lat) * pi / 180.;
    double long_rad = (input->lon) * pi / 180.;

    sinlat = sin(lat_rad);
    coslat = cos(lat_rad);

    sinlong = sin(long_rad);
    coslong = cos(long_rad);

    //    sincos(lat_rad,  &sinlat,  &coslat);
    //    sincos(long_rad, &sinlong, &coslong);

    localVertical[0] = coslong * coslat;
    localVertical[1] = sinlong * coslat;
    localVertical[2] = sinlat;

    double Nphi = aa / sqrt(1. - e2 * pow(sinlat, 2));

    output->x = (Nphi + (input->alt)) * localVertical[0];
    output->y = (Nphi + (input->alt)) * localVertical[1];
    output->z = (Nphi * (1. - e2) + (input->alt)) * localVertical[2];
}

/* ------------------------------------------------------------------- */
void
ecef_to_geocentric(struct ecef_coords *input, struct geocentric_coords *output)
{
    static const double aa = 6378.137; // km,  WGS-84
    static const double pi = 3.14159265359;

    double R;
    double x, y, z;

    x = input->x;
    y = input->y;
    z = input->z;

    // longitude is easy
    R = sqrt(x * x + y * y + z * z); // sqrt(x*x+y*y)

    output->lon = atan2(y, x) * 180. / pi;

    // easy version for latitude (less accurate), of first guess for more accurate
    output->lat = asin(z / R) * 180. / pi;

    output->alt = R - aa;
}

/* ------------------------------------------------------------------- */
void
ecef_to_geodetic(struct ecef_coords *input, struct geodetic_coords *output)
{
    static const double aa = 6378.137;          // km,  WGS-84
    static const double e2 = 0.006694379990141; //  WGS-84
    static const double pi = 3.14159265359;

    double x, y, z;

    x = input->x;
    y = input->y;
    z = input->z;

    double Nphi = 0, R, lat;
    double coslon, sinlon, sinlat, coslat;

    // longitude is easy
    R = sqrt(x * x + y * y + z * z); // sqrt(x*x+y*y)
    output->lon = atan2(y, x) * 180. / pi;     // same as geocentric

    // first guess
    lat = asin(z / R); // radians

    double p = sqrt(x * x + y * y); //

    // Compute latitude recursively
    int i;

    for (i = 0; i < 6; i++) // 6 iterations is enough to get
    {
        // a precision of better than a centimeter
        sinlat = sin(lat);
        Nphi = aa / sqrt(1. - e2 * sinlat * sinlat);
        lat = atan((z + Nphi * e2 * sinlat) / p);
    }

    sinlat = sin(lat);
    coslat = cos(lat);

    //
    // // Get altitude from latitude
    output->alt = p * coslat + (z + e2 * Nphi * sinlat) * sinlat - Nphi;

    output->lat = lat * 180. / pi;
} /* ecef_to_geodetic */

/* ------------------------------------------------------------------- */
void
ecef_to_geodetic_olson(struct ecef_coords *input, struct geodetic_coords *output)
{
    static double a = 6378137.0;             // WGS-84 semi-major axis
    static double e2 = 6.6943799901377997e-3; // WGS-84 first eccentricity squared
    static double a1 = 4.2697672707157535e+4; // a1 = a*e2
    static double a2 = 1.8230912546075455e+9; // a2 = a1*a1
    static double a3 = 1.4291722289812413e+2; // a3 = a1*e2/2
    static double a4 = 4.5577281365188637e+9; // a4 = 2.5*a2
    static double a5 = 4.2840589930055659e+4; // a5 = a1+a3
    static double a6 = 9.9330562000986220e-1; // a6 = 1-e2
    static const double pi = 3.14159265359;
    static double zp, w2, w, r2, r, s2, c2, s, c, ss;
    static double g, rg, rf, u, v, m, f, p, x, y, z;

    double geo[3];

    x = (input->x) * 1000.; // km to m
    y = (input->y) * 1000.;
    z = (input->z) * 1000.;

    zp = fabs(z);
    w2 = x * x + y * y;
    w = sqrt(w2);
    r2 = w2 + z * z;
    r = sqrt(r2);
    geo[1] = atan2(y, x); // Lon (final)
    s2 = z * z / r2;
    c2 = w2 / r2;
    u = a2 / r;
    v = a3 - a4 / r;

    if (c2 > 0.3)
    {
        s = (zp / r) * (1.0 + c2 * (a1 + u + s2 * v) / r);
        geo[0] = asin(s); // Lat
        ss = s * s;
        c = sqrt(1.0 - ss);
    }
    else
    {
        c = (w / r) * (1.0 - s2 * (a5 - u - c2 * v) / r);
        geo[0] = acos(c); // Lat
        ss = 1.0 - c * c;
        s = sqrt(ss);
    }

    g = 1.0 - e2 * ss;
    rg = a / sqrt(g);
    rf = a6 * rg;
    u = w - rg * c;
    v = zp - rf * s;
    f = c * u + s * v;
    m = c * v - s * u;
    p = m / (rf / g + f);
    geo[0] = geo[0] + p;      // Lat
    geo[2] = f + m * p / 2.0; // Altitude

    if (z < 0.0)
    {
        geo[0] *= -1.0; // Lat
    }

    output->lat = geo[0] * 180. / pi;
    output->lon = geo[1] * 180. / pi;
    output->alt = geo[2] / 1000.;
} /* ecef_to_geodetic_olson */
