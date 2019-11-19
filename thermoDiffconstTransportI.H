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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::thermoDiffconstTransport<Thermo>::thermoDiffconstTransport
(
    const Thermo& t,
    const scalar mu,
    const scalar Pr,
/*******************************************/
    const scalar Sc,
    const scalar KT
/******************************************/
)
:
    Thermo(t),
    mu_(mu),
    rPr_(1.0/Pr),
/*******************************************/
    Sc_(Sc),
    KT_(KT)
/*******************************************/
{}


template<class Thermo>
inline Foam::thermoDiffconstTransport<Thermo>::thermoDiffconstTransport
(
    const word& name,
    const thermoDiffconstTransport& ct
)
:
    Thermo(name, ct),
    mu_(ct.mu_),
    rPr_(ct.rPr_),
/*******************************************/
    Sc_(ct.Sc_),
    KT_(ct.KT_)
/*******************************************/
{}


template<class Thermo>
inline Foam::autoPtr<Foam::thermoDiffconstTransport<Thermo>>
Foam::thermoDiffconstTransport<Thermo>::clone() const
{
    return autoPtr<thermoDiffconstTransport<Thermo>>::New(*this);
}


template<class Thermo>
inline Foam::autoPtr<Foam::thermoDiffconstTransport<Thermo>>
Foam::thermoDiffconstTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<thermoDiffconstTransport<Thermo>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::thermoDiffconstTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return mu_;
}


template<class Thermo>
inline Foam::scalar Foam::thermoDiffconstTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return this->Cp(p, T)*mu(p, T)*rPr_;
}


template<class Thermo>
inline Foam::scalar Foam::thermoDiffconstTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{
    return mu(p, T)*rPr_;
}

/**************************************************************/
//NOTE: For the mass diffusivity calculation, I tried to emulate what is done for calculating kappa (at line 91) 
//The issue is that rho(p, T) is not a member of the Thermo class, but I wonder if rho can be accessed by other means.

//- Constant mass diffusivity for binary mixtures.
template<class Thermo> //
inline Foam::scalar Foam::thermoDiffconstTransport<Thermo>::Dab
(
    const scalar p,
    const scalar T
) const
{
    return mu(p, T)/(Foam::thermoDiffperfectGas<Specie>::rho(p, T)*Sc_);
}
//- Thermodiffusion coefficient.
template<class Thermo>
inline Foam::scalar Foam::thermoDiffconstTransport<Thermo>::DT
(
    const scalar p,
    const scalar T
) const
{
    return Dab(p, T)*KT_/T;
}
/****************************************************************/

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::thermoDiffconstTransport<Thermo>::operator=
(
    const thermoDiffconstTransport<Thermo>& ct
)
{
    Thermo::operator=(ct);

    mu_ = ct.mu_;
    rPr_ = ct.rPr_;
/****************************************************************/
    Sc_ = ct.Sc_;
    KT_ = ct.KT_;
/****************************************************************/
}


template<class Thermo>
inline void Foam::thermoDiffconstTransport<Thermo>::operator+=
(
    const thermoDiffconstTransport<Thermo>& st
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(st);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        scalar Y2 = st.Y()/this->Y();

        mu_ = Y1*mu_ + Y2*st.mu_;
        rPr_ = 1.0/(Y1/rPr_ + Y2/st.rPr_);
/****************************************************************/
        Sc_ = Y1*Sc_ + Y2*st.Sc_;
        KT_ = Y1*KT_ + Y2*st.KT_;
/****************************************************************/
    }
}


template<class Thermo>
inline void Foam::thermoDiffconstTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::thermoDiffconstTransport<Thermo> Foam::operator+
(
    const thermoDiffconstTransport<Thermo>& ct1,
    const thermoDiffconstTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) + static_cast<const Thermo&>(ct2)
    );

    if (mag(t.Y()) < SMALL)
    {
        return thermoDiffconstTransport<Thermo>
        (
            t,
            0,
            ct1.rPr_,
/****************************************************************/
            ct1.Sc_,
            ct1.KT_
/****************************************************************/
        );
    }
    else
    {
        scalar Y1 = ct1.Y()/t.Y();
        scalar Y2 = ct2.Y()/t.Y();

        return thermoDiffconstTransport<Thermo>
        (
            t,
            Y1*ct1.mu_ + Y2*ct2.mu_,
            1.0/(Y1/ct1.rPr_ + Y2/ct2.rPr_),
/****************************************************************/
            Y1*ct1.Sc_ + Y2*ct2.Sc_,
            Y1*ct1.KT_ + Y2*ct2.KT_
/****************************************************************/
        );
    }
}


template<class Thermo>
inline Foam::thermoDiffconstTransport<Thermo> Foam::operator*
(
    const scalar s,
    const thermoDiffconstTransport<Thermo>& ct
)
{
    return thermoDiffconstTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        ct.mu_,
        1.0/ct.rPr_,
/****************************************************************/
        ct.Sc_,
        ct.KT_
/****************************************************************/
    );
}


// ************************************************************************* //
