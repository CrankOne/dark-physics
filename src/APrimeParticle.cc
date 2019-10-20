/*
 * Copyright (c) 2016 Renat R. Dusaev <crank@qcrypt.org>
 * Author: Renat R. Dusaev <crank@qcrypt.org>
 * Author: Bogdan Vasilishin <togetherwithra@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

# include "dphmc-config.h"
# include "APrimeParticle.hh"

# include <G4ParticleTable.hh>
# include <G4SystemOfUnits.hh>
# include <G4PhaseSpaceDecayChannel.hh>
# include <G4DalitzDecayChannel.hh>
# include <G4DecayTable.hh>

namespace dphmc {

APrime::Parameters APrime::parameters = {
    /* mass ........... */ 10*MeV,
    /* mixingFactor ... */ 1e-4,
    /* enableDecays ... */ false,
    /* codePDG ........ */ 88
};

APrime * APrime::theInstance = nullptr;

APrime::APrime(){
    dphmc_dbg_msg("A' definition instantiated.");
}

APrime::~APrime(){}

APrime *
APrime::Definition() {
    if( theInstance ) {
        return theInstance;
    }
    const G4String name = "APrime";
    // search in particle table]
    G4ParticleTable * pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition * anInstance = pTable->FindParticle(name);
    if( !anInstance ) {
        // NOTE: this decay formula is given by eq. (11), PhysRevD.80.075018
        // (Bjorken et al. "New fixed-target experiments to search for dark
        // gauge forces"), and introduced here for m_{A'} < 2*m_{\mu}.
        // TODO: if( APrimeMass_GeV > 2*m_mu ) { emraise(...); }
        G4double lifeTime = 0.
               , decayWidth = 0.
               ;
        if( parameters.enableDecays ) {
            # if 0
            G4double lifeTime_s = 80.0/(3*pow(10., 8.)) *
                    pow(10., -6.) *
                    pow( pow(10.,-4.)/mixingFactor, 2. ) *
                    (100*pow(10.,6.)/( APrimeMass_GeV*pow(10.,9) ) );
            //std::cout << "lifeTime_s: " << lifeTime_s << std::endl;  // XXX
            G4double decayWidth_eV = 6.582*pow(10., -16.) / lifeTime_s;
            # else
            lifeTime =
                ((80.*CLHEP::micrometer)/1./*<- N_eff --- branching ratio */)
               *(1e-8/(parameters.mixingFactor * parameters.mixingFactor))
               *(0.1/*GeV*//(parameters.mass/GeV /*GeV*/))
               /CLHEP::c_light
               ;
            decayWidth = CLHEP::hbar_Planck / lifeTime;
            # endif
        }
        anInstance = new G4ParticleDefinition(
                /* Name ..................... */ name,
                /* Mass ..................... */ parameters.mass,
                /* Decay width .............. */ decayWidth,
                /* Charge ................... */ 0.*eplus,
		        /* 2*spin ................... */ 2,
                /* parity ................... */ 0,
                /* C-conjugation ............ */ 0,
		        /* 2*Isospin ................ */ 0,
                /* 2*Isospin3 ............... */ 0,
                /* G-parity ................. */ 0,
	            /* type ..................... */ "boson",
                /* lepton number ............ */ 0,
                /* baryon number ............ */ 0,
                /* PDG encoding ............. */ parameters.codePDG,
		        /* stable ................... */ parameters.enableDecays,
                /* lifetime.................. */ lifeTime,
                /* decay table .............. */ NULL,  // TODO
                /* shortlived ............... */ false,  // TODO
                /* subType .................. */ "geantino",
                /* anti particle encoding ... */ 90    // TODO: ???
            );

        if( parameters.enableDecays ) {
            // Create decay table
            G4DecayTable* decayTable = new G4DecayTable();
            // Create decay channel
            G4VDecayChannel* mode;
            // Create mode A'->e-e+
            mode = new G4PhaseSpaceDecayChannel("APrime", 1.0, 2, "e-", "e+");
            // Insert mode to the table
            decayTable->Insert(mode);
            // Set decay table to A'
            anInstance->SetDecayTable(decayTable);
        }
        // Bohr Magnetron
        G4double muB =  0 ;
        anInstance->SetPDGMagneticMoment( muB * 2.* 1.0011596521859 );
        dphmc_dbg_msg("A' particle type registered with parameters:");
        dphmc_dbg_msg("  mass ....... : %e", parameters.mass );
        dphmc_dbg_msg("  eps ........ : %e", parameters.mixingFactor );
        dphmc_dbg_msg("  decays ..... : %s", parameters.enableDecays ? "yes" : "no" );
        dphmc_dbg_msg("  PDG code ... : %d", parameters.codePDG );
    }
    theInstance = reinterpret_cast<APrime*>(anInstance);
    return theInstance;
}

APrime *
APrime::APrimeDefinition() {
    return Definition();
}

//APrime *
//APrime::APrime() {
//    return Definition();
//}

}  // namespace ::dphmc


