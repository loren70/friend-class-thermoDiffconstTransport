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
    Foam::thermoDiffperfectGas

Group
    grpSpecieEquationOfState

Description
    Perfect gas equation of state.

SourceFiles
    thermoDiffperfectGasI.H
    thermoDiffperfectGas.C

\*---------------------------------------------------------------------------*/

#ifndef thermoDiffperfectGas_H
#define thermoDiffperfectGas_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class thermoDiffperfectGas;

template<class Specie>
inline thermoDiffperfectGas<Specie> operator+
(
    const thermoDiffperfectGas<Specie>&,
    const thermoDiffperfectGas<Specie>&
);

template<class Specie>
inline thermoDiffperfectGas<Specie> operator*
(
    const scalar,
    const thermoDiffperfectGas<Specie>&
);

template<class Specie>
inline thermoDiffperfectGas<Specie> operator==
(
    const thermoDiffperfectGas<Specie>&,
    const thermoDiffperfectGas<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const thermoDiffperfectGas<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class thermoDiffperfectGas Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class thermoDiffperfectGas
:
    public Specie
{
/****************************************/
template<class Thermo> friend class thermoDiffconstTransport;//Friend class declaration
/****************************************/
public:

    // Constructors

        //- Construct from components
        inline thermoDiffperfectGas(const Specie& sp);

        //- Construct from dictionary
        thermoDiffperfectGas(const dictionary& dict);

        //- Construct as named copy
        inline thermoDiffperfectGas(const word& name, const thermoDiffperfectGas&);

        //- Construct and return a clone
        inline autoPtr<thermoDiffperfectGas> clone() const;

        // Selector from dictionary
        inline static autoPtr<thermoDiffperfectGas> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "thermoDiffperfectGas<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp departure [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return internal energy departure [J/kg]
            inline scalar E(const scalar p, const scalar T) const;

            //- Return Cv departure [J/(kg K]
            inline scalar Cv(scalar p, scalar T) const;

            //- Return entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const thermoDiffperfectGas&);
        inline void operator*=(const scalar);


    // Friend operators

        friend thermoDiffperfectGas operator+ <Specie>
        (
            const thermoDiffperfectGas&,
            const thermoDiffperfectGas&
        );

        friend thermoDiffperfectGas operator* <Specie>
        (
            const scalar s,
            const thermoDiffperfectGas&
        );

        friend thermoDiffperfectGas operator== <Specie>
        (
            const thermoDiffperfectGas&,
            const thermoDiffperfectGas&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const thermoDiffperfectGas&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoDiffperfectGasI.H"

#ifdef NoRepository
    #include "thermoDiffperfectGas.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
