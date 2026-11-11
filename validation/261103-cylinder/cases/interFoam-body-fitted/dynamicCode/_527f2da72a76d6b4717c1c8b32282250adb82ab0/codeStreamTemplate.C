/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR AUTHOR, AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Template for use with codeStream.

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "Ostream.H"
#include "Pstream.H"
#include "pointField.H"
#include "tensor.H"
#include "unitConversion.H"

//{{{ begin codeInclude
#line 46 "/users/sandyz/slurm_of/waves2-3d/ichiBoundaryFoam/validation/261103-cylinder/cases/interFoam-body-fitted/system/blockMeshDict/#codeStream"
#include "pointField.H"
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C" void codeStream_527f2da72a76d6b4717c1c8b32282250adb82ab0(Foam::Ostream& os, const Foam::dictionary& dict)
{
//{{{ begin code
    #line 51 "/users/sandyz/slurm_of/waves2-3d/ichiBoundaryFoam/validation/261103-cylinder/cases/interFoam-body-fitted/system/blockMeshDict/#codeStream"
pointField points(19);
        points[0] = point(0.06, 0, -0.0075);
        points[1] = point(0.1432392, 0, -0.0075);
        points[2] = point(0.6, 0, -0.0075);
        points[3] = point(0.6, 0.1012854096517361, -0.0075);
        points[4] = point(0.1012854096517361, 0.1012854096517361, -0.0075);
        points[5] = point(0.04242640687119285, 0.04242640687119285, -0.0075);
        points[6] = point(0.6, 0.6, -0.0075);
        points[7] = point(0.1012854096517361, 0.6, -0.0075);

        // Mirror +x points to -x side
        points[11] = point(-points[0].x(), points[0].y(), points[0].z());
        points[12] = point(-points[1].x(), points[1].y(), points[1].z());
        points[13] = point(-points[2].x(), points[2].y(), points[2].z());
        points[14] = point(-points[3].x(), points[3].y(), points[3].z());
        points[15] = point(-points[4].x(), points[4].y(), points[4].z());
        points[16] = point(-points[5].x(), points[5].y(), points[5].z());
        points[17] = point(-points[6].x(), points[6].y(), points[6].z());
        points[18] = point(-points[7].x(), points[7].y(), points[7].z());

        // Points on the y-z-plane
        points[8] = point(0, 0.6, -0.0075);
        points[9] = point(0, 0.1432392, -0.0075);
        points[10] = point(0, 0.06, -0.0075);

        // Mirror -z points to +z side
        label sz = points.size();
        points.setSize(2*sz);
        for (label i = 0; i < sz; ++i)
        {
            const point& pt = points[i];
            points[i + sz] = point(pt.x(), pt.y(), 0.0075);
        }

        // Add an inner cylinder
        sz = points.size();
        label nAdd = 6;
        points.setSize(sz + nAdd);

        // Points within the inner cylinder
        points[sz] = point(0, 0, -0.0075);
        points[sz + 1] = point(0, 0.03, -0.0075);
        points[sz + 2] = point(0.03, 0.03, -0.0075);
        points[sz + 3] = point(0.03, 0, -0.0075);

        // Mirror points from +x side to -x side
        points[sz + 4] =
            point(-points[sz + 2].x(), points[sz + 2].y(), points[sz + 2].z());
        points[sz + 5] =
            point(-points[sz + 3].x(), points[sz + 3].y(), points[sz + 3].z());

        // Mirror -z points to +z side
        sz = points.size();
        points.setSize(sz + nAdd);
        for (label i = 0; i < nAdd; ++i)
        {
            const point& pt = points[i+sz-nAdd];
            points[i+sz] = point(pt.x(), pt.y(), 0.0075);
        }

        // Add downstream and upstream blocks
        sz = points.size();
        nAdd = 6;
        points.setSize(sz + nAdd);

        // Points on outlet
        points[sz] = point(1.120002, 0, -0.0075);
        points[sz + 1] = point(1.120002, 0.6, -0.0075);
        points[sz + 4] = point(1.120002, 0.1012854096517361, -0.0075);

        // Points on inlet
        points[sz + 2] = point(-0.6074999999999999, 0, -0.0075);
        points[sz + 3] = point(-0.6074999999999999, 0.6, -0.0075);
        points[sz + 5] = point(-0.6074999999999999, 0.1012854096517361, -0.0075);

        // Mirror -z points to +z side
        sz = points.size();
        points.setSize(sz + nAdd);
        for (label i = 0; i < nAdd; ++i)
        {
            const point& pt = points[i + sz - nAdd];
            points[i + sz] = point(pt.x(), pt.y(), 0.0075);
        }

        os  << points;
//}}} end code
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

