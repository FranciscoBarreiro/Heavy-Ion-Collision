// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

#include "Rivet/Projections/HepMCHeavyIon.hh"
#include "Rivet/Projections/JetFinder.hh"
#include <fstream>

using namespace std;

namespace Rivet {


  /// Generic analysis looking at various distributions of final state particles
  /// @todo Rename as MC_HADRONS
  class AnalysisHadron : public Analysis {
  public:

    /// Constructor
    AnalysisHadron()
      : Analysis("AnalysisHadron")
    {    }


  public:

    /// @name Analysis methods
    //@{

    void init() {

      // Projections
      //????????ATENÇÃO: pode se alterar aqui o corte do soft?????
      //const FinalState soft(Cuts::abseta < 2.0);
	  const FinalState soft(Cuts::abseta < 2.0 && Cuts::pT > 30*GeV);
      declare(soft, "softplus");

      // get file name
        string name;
        std::cout << std::endl << std::endl << "Out file name: ";
        std::cin >> name;
        if(name == "e\n") name = "eventdata.dat";
        std::cout << std::endl;

      //Abrir um ficheiro para escrita
      Hadrondata.open(name);
      Hadrondata << "Ângulo	Jprodr	Weight	ParticleID	Energy	Px	Py	Pz	Pt	Ângulohadron\n";

    }



    /// Perform the per-event analysis
    void analyze(const Event& event) {

      //Obter o hadron mais energético
      const FinalState& soft = apply<FinalState>(event, "softplus");
      FourMomentum pi(0,0,0,0);
      HepMC::GenParticle gp(pi, PID::PROTON);
      for (const Particle& p1 : soft.particles()) {
        Particle gt(gp);
        if(p1.momentum().pT() > gt.momentum().pT() && p1.isHadron())
        {
          gp.set_momentum(p1.momentum());
          gp.set_pdg_id(p1.pid());
        }
      }


      /*cout << "Particle ID" << endl;
      cout << gp.pdg_id() << endl;
      cout << "Energy maximum" << endl;
      cout << gp.momentum().e() << endl;
      cout << "Px" << endl;
      cout << gp.momentum().px() << endl;
      cout << "Py" << endl;
      cout << gp.momentum().py() << endl;
      cout << "Pz" << endl;
      cout << gp.momentum().pz() << endl;*/


      //Imprimir o peso
      /*cout << "peso" << endl;
      cout << event.weights()[0] << endl;*/

      //Ler a linha H
      const HepMC::HeavyIon *hion = event.genEvent()->heavy_ion();
      if (hion) {
        /*cout << "phi" << "  " << "jprodr" << endl;
        cout << hion->event_plane_angle() << "  " << hion->eccentricity() << endl;
        cout << " " << endl;*/
        /*cout << "phi" << endl;
        cout << hi.Ncoll_hard() << endl;
        cout << hi.Npart_proj() << endl;
        cout << hi.Npart_targ() << endl;
        cout << hi.N_Nwounded_collisions() << endl;
        cout << hi.Nwounded_N_collisions() << endl;
        cout << hi.Nwounded_Nwounded_collisions() << endl;
        cout << hi.impact_parameter() << endl;
        cout << hi.event_plane_angle() << endl;
        cout << hi.sigma_inel_NN() << endl;
        cout << hi.Nspec_proj_n() << endl;
        cout << hi.Nspec_targ_n() << endl;
        cout << hi.Nspec_proj_p() << endl;
        cout << hi.Nspec_targ_p() << endl;
        std::map<int,double> hf=hi.participant_plane_angles();
        std::map<int,double>::iterator it1=hf.begin();
        cout << it1->first << endl;*/
      }
      //Escrever no ficheiro
      if(sqrt(gp.momentum().px()*gp.momentum().px()+gp.momentum().py()*gp.momentum().py()) > 0)
      {
        Hadrondata << hion->event_plane_angle() << "	" << hion->eccentricity() << "	" << event.weights()[0] << "	" << gp.pdg_id() << "	" << gp.momentum().e() << "	" << gp.momentum().px() << "	" << gp.momentum().py() << "	" << gp.momentum().pz() << "	" << sqrt(gp.momentum().px()*gp.momentum().px()+gp.momentum().py()*gp.momentum().py()) << "	" << gp.momentum().phi() <<"\n";
      //Hadrondata << hion->event_plane_angle() << " " << hion->eccentricity() << " " << event.weights()[0] << "\n";
      }
      
    }



    /// Finalize
    void finalize() {
      Hadrondata.close();
    }

    //@}


  private:

    /// @name Histograms
    //@{
    ofstream Hadrondata;

    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(AnalysisHadron);

}