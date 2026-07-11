/*---------------------------------------------------------------------------*\
    libmyThermos  —  custom psiThermo instantiations for hypersonic ZPGFP work

    Uses makeThermos (PLURAL): registers each combo into
        basicThermo + fluidThermo + psiThermo
    so that solvers selecting via fluidThermo::New (e.g. HiSA) find them.
    (makeThermo, singular, would register into psiThermo ONLY.)

        polynomialTransport / hConstThermo / perfectGas / sensInternalEnergy
        polynomialTransport / janafThermo  / perfectGas / sensInternalEnergy

    Build:   wmake libso       (-> $FOAM_USER_LIBBIN/libmyThermos.so)
    Use:     libs ("libmyThermos.so");   in system/controlDict
\*---------------------------------------------------------------------------*/

#include "psiThermo.H"
#include "makeThermo.H"

#include "specie.H"
#include "perfectGas.H"
#include "hConstThermo.H"
#include "janafThermo.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"

#include "polynomialTransport.H"

#include "hePsiThermo.H"
#include "pureMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    polynomialTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    polynomialTransport,
    sensibleInternalEnergy,
    janafThermo,
    perfectGas,
    specie
);

} // End namespace Foam

// ************************************************************************* //