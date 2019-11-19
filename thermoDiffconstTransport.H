/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Class
    Foam::thermoDiffconstTransport

Group
    grpSpecieTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    thermoDiffconstTransportI.H
    thermoDiffconstTransport.C

\*---------------------------------------------------------------------------*/


#ifndef thermoDiffconstTransport_H
#define thermoDiffconstTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class thermoDiffconstTransport; 

template<class Thermo>
inline thermoDiffconstTransport<Thermo> operator+
(
    const thermoDiffconstTransport<Thermo>&,
    const thermoDiffconstTransport<Thermo>&
);

template<class Thermo>
inline thermoDiffconstTransport<Thermo> operator*
(
    const scalar,
    const thermoDiffconstTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const thermoDiffconstTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class thermoDiffconstTransport Declaration
\*---------------------------------------------------------------------------*/
/*****************************************/
#include "thermoDiffperfectGas.H" 
/******************************************/
template<class Thermo>
class thermoDiffconstTransport
:
    public Thermo          // public thermoDiffperfectGas<Specie> 
{
    // Private data

        //- Constant dynamic viscosity [Pa.s]
        scalar mu_;

        //- Reciprocal Prandtl Number []
        scalar rPr_;
/************************************************/
//NOTE: New Data members, this dimensionless quantiites are implemented as dictionary entries that are necessary to calculate mass diffusivity and thermodiffusion coefficient.
        //- Schmidt number []
        scalar Sc_;
        //- ThermoDiffusion ratio []
        scalar KT_;
/************************************************/

    // Private Member Functions

        //- Construct from components
        inline thermoDiffconstTransport
        (
            const Thermo& t,
            const scalar mu,
            const scalar Pr,
         /************************************************/
            const scalar Sc,
            const scalar KT
/************************************************/
        );


public:
    
        

    // Constructors

        //- Construct as named copy
        inline thermoDiffconstTransport(const word&, const thermoDiffconstTransport&);

        //- Construct from dictionary
        thermoDiffconstTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<thermoDiffconstTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<thermoDiffconstTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "thermoDiffconst<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;
/******************************************************************/
//NOTE: New member functions -Mass diffusivity needs to be calculated with density as a parameter, there is an inline function where it can be seen at the thermoDiffusionTransport.H file
        //- Mass diffusivity [m^2/s] 
        inline scalar Dab(const scalar p, const scalar T) const;
        //- thermodiffusion coefficient [m^2/s] 
        inline scalar DT(const scalar p, const scalar T) const;
/*******************************************************************/
        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline void operator=(const thermoDiffconstTransport&);

        inline void operator+=(const thermoDiffconstTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend thermoDiffconstTransport operator+ <Thermo>
        (
            const thermoDiffconstTransport&,
            const thermoDiffconstTransport&
        );

        friend thermoDiffconstTransport operator* <Thermo>
        (
            const scalar,
            const thermoDiffconstTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const thermoDiffconstTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoDiffconstTransportI.H"

#ifdef NoRepository
    #include "thermoDiffconstTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
