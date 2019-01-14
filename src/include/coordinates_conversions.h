// *******************************************************************************

// GEANT4 code for propagation of gamma-rays, electron and positrons in Earth's atmosphere

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

#include <math.h>   /* maths functions */
#include <stdio.h>  /* for error messages. */
#include <stdlib.h> /* for malloc/free */

// STRUCTURES

struct ecef_coords
{
  double x; // km
  double y; // km
  double z; // km
};

struct geocentric_coords
{
  double lat; // degree
  double lon; // degree
  double alt; // km
};

struct geodetic_coords
{             // = geographic
  double lat; // degree
  double lon; // degree
  double alt; // km
};

// PROTOTYPES

void geocentric_to_ecef(struct geocentric_coords *input,
                        struct ecef_coords       *output);

void geodetic_to_ecef(struct geodetic_coords *input,
                      struct ecef_coords     *output);

void ecef_to_geocentric(struct ecef_coords       *input,
                        struct geocentric_coords *output);

void ecef_to_geodetic(struct ecef_coords     *input,
                      struct geodetic_coords *output);

void ecef_to_geodetic_olson(struct ecef_coords     *input,
                            struct geodetic_coords *output);
