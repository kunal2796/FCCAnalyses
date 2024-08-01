
#ifndef JETCLUSTERINGUTILS_ANALYZERS_H
#define JETCLUSTERINGUTILS_ANALYZERS_H

#include "Math/Vector4D.h"
#include "ROOT/RVec.hxx"

#include "FastJet/JetClustering.h"
#include "fastjet/JetDefinition.hh"

#include "TRandom3.h"

namespace FCCAnalyses {
  /**
   * @brief Jet clustering interface utilities.
   *
   * This represents a set functions and utilities to perform jet clustering
   * from a list of particles.
   */
  namespace JetClusteringUtils {
    /**
     * Maximum number of d_{n, n+1} that are kept in FCCAnalysesJet
     */
    const int Nmax_dmerge = 10;

    /** Set fastjet pseudoJet for later reconstruction*/
    std::vector<fastjet::PseudoJet> set_pseudoJets(const ROOT::VecOps::RVec<float>& px,
                                                   const ROOT::VecOps::RVec<float>& py,
                                                   const ROOT::VecOps::RVec<float>& pz,
                                                   const ROOT::VecOps::RVec<float>& e);

    /** Set fastjet pseudoJet for later reconstruction using px, py, pz and m
     *
     * This version is to be preferred over the px,py,pz,E version when m is known
     * accurately, because it uses double precision to reconstruct the energy,
     * reducing the size of rounding errors on FastJet calculations (e.g. of
     * PseudoJet masses)
     *
    */
    std::vector<fastjet::PseudoJet> set_pseudoJets_xyzm(const ROOT::VecOps::RVec<float>& px,
                                                        const ROOT::VecOps::RVec<float>& py,
                                                        const ROOT::VecOps::RVec<float>& pz,
                                                        const ROOT::VecOps::RVec<float>& m);

    /** Get fastjet pseudoJet after reconstruction from FCCAnalyses jets*/
    ROOT::VecOps::RVec<fastjet::PseudoJet> get_pseudoJets(const JetClustering::FCCAnalysesJet& in);

    /** Get fastjet constituents after reconstruction from FCCAnalyses jets*/
    std::vector<std::vector<int>> get_constituents(const JetClustering::FCCAnalysesJet& in);

    /// return the dmin corresponding to the recombination that went from n+1 to n jets
    float get_exclusive_dmerge(const JetClustering::FCCAnalysesJet& in, int n);

    float get_exclusive_dmerge_max(const JetClustering::FCCAnalysesJet& in, int n);

    /** Get jet px. Details. */
    ROOT::VecOps::RVec<float> get_px(const ROOT::VecOps::RVec<fastjet::PseudoJet>& in);

    /** Get jet py. Details. */
    ROOT::VecOps::RVec<float> get_py(const ROOT::VecOps::RVec<fastjet::PseudoJet>& in);

    /** Get jet pz. Details. */
    ROOT::VecOps::RVec<float> get_pz(const ROOT::VecOps::RVec<fastjet::PseudoJet>& in);

    /** Get jet energy. Details. */
    ROOT::VecOps::RVec<float> get_e(const ROOT::VecOps::RVec<fastjet::PseudoJet>& in);

    /** Get jet pt. Details. */
    ROOT::VecOps::RVec<float> get_pt(const ROOT::VecOps::RVec<fastjet::PseudoJet>& in);

    /** Get jet p. Details. */
    ROOT::VecOps::RVec<float> get_p(const ROOT::VecOps::RVec<fastjet::PseudoJet>& in);

    /** Get jet mass. Details. */
    ROOT::VecOps::RVec<float> get_m(const ROOT::VecOps::RVec<fastjet::PseudoJet>& in);

    /** Get jet eta. Details. */
    ROOT::VecOps::RVec<float> get_eta(const ROOT::VecOps::RVec<fastjet::PseudoJet>& in);

    /** Get jet phi. Details (range [0,2*pi]). */
    ROOT::VecOps::RVec<float> get_phi(const ROOT::VecOps::RVec<fastjet::PseudoJet>& in);

    /** Get jet phi. Details (range [-pi,pi]). */
    ROOT::VecOps::RVec<float> get_phi_std(const ROOT::VecOps::RVec<fastjet::PseudoJet>& in);


    /** Get jet theta. Details. */
    ROOT::VecOps::RVec<float> get_theta(const ROOT::VecOps::RVec<fastjet::PseudoJet>& in);

    ///Select clustered jets with transverse momentum greader than a minimum value [GeV]
    struct sel_pt {
      sel_pt(float arg_min_pt);
      float m_min_pt = 1.; //> transverse momentum threshold [GeV]
      ROOT::VecOps::RVec<fastjet::PseudoJet>  operator() (ROOT::VecOps::RVec<fastjet::PseudoJet> in);
    };

    /// Internal methods
    JetClustering::FCCAnalysesJet initialise_FCCAnalysesJet();

    JetClustering::FCCAnalysesJet build_FCCAnalysesJet(const std::vector<fastjet::PseudoJet>& in,
                                                       const std::vector<float>& dmerge,
                                                       const std::vector<float>& dmerge_max);

    std::vector<fastjet::PseudoJet> build_jets(fastjet::ClusterSequence& cs, int exclusive, float cut, int sorted);

    bool check(unsigned int n, int exclusive, float cut);

    fastjet::RecombinationScheme recomb_scheme(int recombination);

    std::vector<float> exclusive_dmerge(fastjet::ClusterSequence& cs, int do_dmarge_max);

    // build the resonance from 2 <fastjet::PseudoJet> objects. Keep the closest to the mass given as input
    struct resonanceBuilder {
      float m_resonance_mass;
      resonanceBuilder(float arg_resonance_mass);
      ROOT::VecOps::RVec<fastjet::PseudoJet> operator()(ROOT::VecOps::RVec<fastjet::PseudoJet> legs);
    };

    struct recoilBuilder {
      recoilBuilder(float arg_sqrts);
      float m_sqrts = 240.0;
      double operator() (ROOT::VecOps::RVec<fastjet::PseudoJet> in);
    };

    /// --- From Edi --- ///

    /** Reshape the given variable (given as a flat vector for the event) into 2d vector per jet. */
    std::vector<std::vector<float>> reshape2jet(ROOT::VecOps::RVec<float> var, std::vector<std::vector<int>> constituents);
    std::vector<std::vector<int>> int_2d();
    std::vector<std::vector<float>> float_2d(std::vector<std::vector<int>> constituents);

    /** Get number of constituents per jet. */
    ROOT::VecOps::RVec<float> get_nConstituents(std::vector<std::vector<int>> in);

    //The below fncs take already the reshaped values as arguments... Different from above...  
    /** Get difference of jet theta and jet constituent theta. */
    std::vector<std::vector<float>> get_dTheta(ROOT::VecOps::RVec<float> jet_theta, std::vector<std::vector<float>> constituents_theta);

    /** Get difference of jet phi and jet constituent phi in [-pi, pi]. */
    std::vector<std::vector<float>> get_dPhi(ROOT::VecOps::RVec<float> jet_phi, std::vector<std::vector<float>> constituents_phi);
  
    /** Get difference dR in theta-phi space of jet and jet constituent. */
    std::vector<std::vector<float>> get_dR(std::vector<std::vector<float>> constituents_dTheta, std::vector<std::vector<float>> constituents_dPhi);
  
    /** Get ratio of jet constituent |p| and jet |p|. */
    std::vector<std::vector<float>> get_pRel(ROOT::VecOps::RVec<float> jet_p, std::vector<std::vector<float>> constituents_p);

    /** Get ratio of jet constituent e and jet e. */
    std::vector<std::vector<float>> get_eRel(ROOT::VecOps::RVec<float> jet_e, std::vector<std::vector<float>> constituents_e);
  
    /** Get angles of jet constituents and jet axis. */
    std::vector<std::vector<float>> get_dAngle(ROOT::VecOps::RVec<float> jet_px, 
					       ROOT::VecOps::RVec<float> jet_py, 
					       ROOT::VecOps::RVec<float> jet_pz, 
					       std::vector<std::vector<float>> constituents_px,
					       std::vector<std::vector<float>> constituents_py,
					       std::vector<std::vector<float>> constituents_pz
					       );
    std::vector<float> get_angularity(float kappa, float beta, ROOT::VecOps::RVec<float> jet_e, std::vector<std::vector<float>> constituents_e, std::vector<std::vector<float>> constituents_dAngle);

    /// --- From Edi --- ///

  }  // namespace JetClusteringUtils
}  // namespace FCCAnalyses

#endif
