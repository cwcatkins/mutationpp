/**
 * @file OmegaVT.cpp
 *
 * @brief Implementation of OmegaVT.
 */

/*
 * Copyright 2014 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include "MillikanWhite.h"
#include "Mixture.h"
#include "TransferModel.h"
#include <cmath>

using namespace Mutation;

namespace Mutation {
    namespace Transfer {

/**
 * Represents a coupling between vibrational and translational energy modes.
 */
class OmegaVT : public TransferModel
{
public:

    OmegaVT(Mixture& mix)
        : TransferModel(mix), m_mw(mix)
    {
        m_const_Park_correction = std::sqrt(PI*KB/(8.E0*NA));
        m_ns              = m_mixture.nSpecies();
        m_transfer_offset = m_mixture.hasElectrons() ? 1 : 0;

        mp_Mw = new double [m_ns];
        for(int i = 0; i < m_ns; ++i)
            mp_Mw[i] = m_mixture.speciesMw(i);
        mp_hv = new double [m_ns];
        mp_hveq = new double [m_ns];
        mp_hel = new double [m_ns];
        mp_heleq = new double [m_ns];

    
	// Set the default source-term options:
	m_energy_diff = "vibrational";
	m_park_correction = "mutationpp";

	// Check to see if the options are used in the XML file:
	std::string filename =
	  Mutation::Utilities::getEnvironmentVariable("MPP_DATA_DIRECTORY") + "/transfer/VT.xml";
    
	// Open the VT file as an XML document
	Mutation::Utilities::IO::XmlDocument doc(filename);
	
	// Find the source-options tag if its present:
	Mutation::Utilities::IO::XmlElement::const_iterator iter = doc.root().findTag("source-options");
	if (iter != doc.root().end()){
	  const Mutation::Utilities::IO::XmlElement& options_element = *iter;

	  // Load the energy_difference and park_correction options:
	  options_element.getAttribute("energy_difference", m_energy_diff, m_energy_diff);
	  options_element.getAttribute("park_correction", m_park_correction, m_park_correction);

	  // Check that they match what we expect:
	  if ( m_energy_diff != "vibrational" && m_energy_diff != "vibrational-electronic" )
	    options_element.parseError("energy_difference not recognised!");
	  if ( m_park_correction != "mutationpp" && m_park_correction != "hy2foam")
	    options_element.parseError("park_correction not recognised!");
	}
    }

    /**
     * Computes the source terms of the Vibration-Translational energy transfer in \f$ [J/(m^3\cdot s)] \f$
     * using a Landau-Teller formula taking into account Park's correction (default; can be disabled by making zero the appropriate flag, see below):
     * \f[ \Omega^{VT}_m = \rho_m \frac{e^V_m\left(T\right)-e^V_m\left(T_{vm}\right)}
     *      {\tau^{VT}_m} \left|\frac{T_{sh}-T_{V}}{T_{sh}-T_{Vsh}}\right|^{s-1} \f]
     * with \f$ s = 3.5exp(-\frac{5000}{T_s}) \f$.
     *
     * The average relaxation time \f$ \tau^{VT}_m \f$ is given by the expression:
     *
     * \f[ \tau^{VT}_m = \frac{ \sum_{j \in \mathcal{H}} \rho_j / M_j}{ \sum_{j \in \mathcal{H}} \rho_j / (M_j \tau^{VT}_{mj}) } \f]
     *
     * More information about the above model can be found in @cite Park1993 and
     * @cite Schwarz1952 .
     *
     */

    double source()
    {
        const double * p_Y = m_mixture.Y();
        double rho = m_mixture.density();
        double T = m_mixture.T();
        double Tv = m_mixture.Tv();

        m_mixture.speciesHOverRT(T, T, T, T, T, NULL, NULL, NULL, mp_hveq, mp_heleq, NULL);
        m_mixture.speciesHOverRT(T, Tv, T, Tv, Tv, NULL, NULL, NULL, mp_hv, mp_hel, NULL);

        // int inv = 0;
	inv = 0;
        double src = 0.0;
        for (int iv = 0; iv-inv < m_mw.nVibrators(); ++iv){
            if(m_mixture.species(iv).type() != Mutation::Thermodynamics::MOLECULE){
                inv++;
            } else {
	      if (m_energy_diff == "vibrational-electronic" )
		src += p_Y[iv]*rho*RU*T/mp_Mw[iv]*(mp_hveq[iv]+mp_heleq[iv] - (mp_hv[iv]+mp_hel[iv]))/compute_tau_VT_m(iv-inv);
	      else
		src += p_Y[iv]*rho*RU*T/mp_Mw[iv]*(mp_hveq[iv] - mp_hv[iv])/compute_tau_VT_m(iv-inv);
            }
        }
        return src;
    }

private:
    //const Thermodynamics::Thermodynamics* mp_thermo;
    MillikanWhite m_mw;

    // Declaring a function to compute \tau_VT_mj for each species. It returns a pointer to tau_VT_m.
    inline double const compute_tau_VT_mj(int const, int const);

    /**
     * Computes the frequency avarage over heavy particles.
     */
    double compute_tau_VT_m(int const);

    /**
     * This function computes the Park correction for vibrational-translational energy trasfer.
     */
    inline double const compute_Park_correction_VT(int const, int const);

    /**
     * Necesseary variables
     */
    int m_ns;
    int m_transfer_offset;
    int inv;
    double* mp_Mw;
    double* mp_hv;
    double* mp_hveq;
    double* mp_hel;
    double* mp_heleq;
    double m_const_Park_correction;
  std::string m_energy_diff;
  std::string m_park_correction;
};
      
// Implementation of the Vibrational-Translational Energy Transfer.
      
inline double const OmegaVT::compute_Park_correction_VT(int const i_vibrator, int const i_partner)
{
    // Limiting cross section for Park's Correction
    double P = m_mixture.P();
    double T = m_mixture.T();
    double sigma;
    
    if ( m_park_correction == "mutationpp" ){
      if (T > 20000.0) {
	sigma = m_mw[i_vibrator].omega() * 6.25 ; // 6.25 = (50000/20000)^2
      } else {
	sigma = m_mw[i_vibrator].omega() *(2.5E9/(T*T));
      }
      
      return (m_const_Park_correction * sqrt(m_mw[i_vibrator][i_partner].mu()*T)/(sigma*P));
    }
    else {
      // hy2Foam relaxation time:
      sigma = m_mw[i_vibrator].omega() *(2.5E9/(T*T));
      double rho = m_mixture.density();
      double n_sr, c_s;
      const double * p_Y = m_mixture.Y();
      if ( i_vibrator+inv == i_partner+m_transfer_offset )
	n_sr = NA*rho*p_Y[i_vibrator+inv] / mp_Mw[i_vibrator+inv];
      else
	n_sr = NA*rho*( (p_Y[i_vibrator+inv]/mp_Mw[i_vibrator+inv]) +
    		      (p_Y[i_partner+m_transfer_offset]/mp_Mw[i_partner+m_transfer_offset]) );
      c_s = sqrt(8E0*T*RU/(mp_Mw[i_vibrator+inv]*PI));
      return (1.0/(c_s*n_sr*sigma));
    }
}
 
inline double const OmegaVT::compute_tau_VT_mj(int const i_vibrator, int const i_partner)
{
//    Enable in the future for multiple vibrational temperatures
      double P = m_mixture.P();
      double T = m_mixture.T();

      return( exp( m_mw[i_vibrator][i_partner].a() * (pow(T,-1.0/3.0) - m_mw[i_vibrator][i_partner].b()) -18.421) * ONEATM / P );
}
      
double OmegaVT::compute_tau_VT_m(int const i_vibrator)
{
    const double * p_Y = m_mixture.Y();
    
    double sum1 = 0.0;
    double sum2 = 0.0;
    
    // Partner offset
    for (int i_partner = m_transfer_offset; i_partner < m_ns; ++i_partner){
        double tau_j = compute_tau_VT_mj(i_vibrator, i_partner-m_transfer_offset) + compute_Park_correction_VT(i_vibrator, i_partner-m_transfer_offset);
        sum1 += p_Y[i_partner]/(mp_Mw[i_partner]);
        sum2 += p_Y[i_partner]/(mp_Mw[i_partner]*tau_j);
    }
  
    return(sum1/sum2);
}


// Register the transfer model
Utilities::Config::ObjectProvider<
    OmegaVT, TransferModel> omegaVT("OmegaVT");

    } // namespace Transfer
} // namespace Mutation 
