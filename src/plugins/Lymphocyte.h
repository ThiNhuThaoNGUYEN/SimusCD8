// ****************************************************************************
//
//              SiMuScale - Multi-scale simulation framework
//
// ****************************************************************************
//
// Copyright: See the AUTHORS file provided with the package
// E-mail: simuscale-contact@lists.gforge.inria.fr
// Original Authors : Samuel Bernard, Carole Knibbe, David Parsons
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ****************************************************************************

// <TODO> Modify the include guard </TODO>
#ifndef SIMUSCALE_LYMPHOCYTE_H__
#define SIMUSCALE_LYMPHOCYTE_H__


// ============================================================================
//                                   Includes
// ============================================================================
#include <memory>
#include <string>
//using namespace std;

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "Cell.h"
#include "Alea.h"
#include "params.h"

/**
 * This is a cell formalism plugin boilerplate
 */
// <TODO> Modify the class name </TODO>

class Lymphocyte : public Cell {
public:

    typedef enum {
      APC_Encounters,TCC,T_Encounters // DifferentiationState, 
    } StateVariable;

    // =================================================================
    //                         Factory functions
    // =================================================================

    template <typename... Args>
    static std::shared_ptr<Cell> MakeCell(Args&&... args) {
        return std::shared_ptr<Lymphocyte>(new Lymphocyte(std::forward<Args>(args)...));
    }

    // ==========================================================================
    //                               Constructors
    // ==========================================================================
    Lymphocyte() = delete; //< no default Default ctor
    Lymphocyte(const Lymphocyte &model); //< Copy ctor
    Lymphocyte(Lymphocyte &&) = delete; //< Move ctor
    Lymphocyte(CellType cellType,
            const MoveBehaviour& move_behaviour,
            const Coordinates<double>& pos,
            double initial_volume,
            double volume_min,
            double doubling_time);
    Lymphocyte(gzFile& backup);
    
    

    // ==========================================================================
    //                                Destructor
    // ==========================================================================
    virtual ~Lymphocyte();

    // ==========================================================================
    //                                Operators
    // ==========================================================================
    /// Copy assignment
    Lymphocyte &operator=(const Lymphocyte &other) = delete;

    /// Move assignment
    Lymphocyte &operator=(Lymphocyte &&other) = delete;

    // ==========================================================================
    //                              Public Methods
    // ==========================================================================
    void InternalUpdate(const double& dt) override;
    Cell* Divide() override;
     
    void Save(gzFile backup) const override;
    void Load(gzFile& backup);

    double get_GeneParams(void);
//    double get_CyclingParams(void);
    double get_InitValues(void);
    double get_P(int i) const;
    double get_mRNA(int i) const;
    double getPhylogeny_t(int i) const;
    int getPhylogeny_id(int i) const;
   

    // ==========================================================================
    //                                Accessors
    // ==========================================================================
    double get_output(InterCellSignal signal) const override;
//    void UpdateCelltype() const override;
  //  bool isCycling() const override;
    bool isDividing() const override;
    bool isDying() const override;


protected:
    // ==========================================================================
    //                            Protected Method
    // ==========================================================================
    //double GetDifferentiationState();
//    double Get_initval(int i, double * P, const double APC_c, const double Tcc, const double x, double initval0);
//double KKon_bound_tau(int i, double * P, double * M, const double APC_c,const  double Tcc, const double S1);
//double GetSigma_i(double initval);
 void Get_Sigma( double * sigma, double * P, const double duration_APC, const double Tcc, const double x, bool AcceptNegative);
   // double GetSigma_i_new(int i, double initval);
    //int choice_i(std::vector<double> ProbaArray);
    void Intracellular_ExactEvol(double DeltaT, double * P, double * M, double * S1);
    //double random_mitotic_threshold(void){return Alea::exponential_random(mitotic_threshold_P0_);};

    void ODE_update(const double& dt);
  //  static int compute_dydt(double t, const double y[], double dydt[],
    //        void* cell);
//    int compute_dydt(double t, const double y[], double dydt[]) const;

    void Dying();
   // void UpdateMitoticStatus();
   // void UpdateCyclingStatus() override;
   // bool StopCycling() override;
   // bool StartCycling() override;
Coordinates<double> MotileDisplacement(const double& dt) override;
 
    void SetNewProteinLevel1(const double Ki, const double MotherProteins[],const  int sz)const ;
    void UpdateContactAPC(double Mother_APC_contact);
    void SetNewProteinLevel2(const double Ki,const double MotherProteins[],const int sz)const ;
     void SetNewRNALevel1(const double Km, const double MotherRNA[],const  int sz)const ;
    void SetNewRNALevel2(const double Km,const double MotherRNA[],const int sz)const ;
    void UpdatePhylogeny(std::vector<int> phy_id, std::vector<double> phy_t, int sz);
     void UpdateContactTCC(double Mother_TCC_contact);
 //    void Update_count_division(const double Mother_division[], const int s, const int z )const;
 void Update_count_division(double Mother_division, const int s )const;
 void Update_count_division_mother(double Mother_division);
    // ==========================================================================
    //                               Attributes
    // ==========================================================================
public:
    // <TODO> Modify the class id and keyword </TODO>
    // Class ID, uniquely identifies this class
    static constexpr CellFormalism classId_ = 0xdaa16ede; // CRC32 Hash
    // Keyword to be used in param files
    static constexpr char classKW_[] = "LYMPHOCYTE";

protected:
    // <TODO> Add your formalism-specific attributes here </TODO>
    static constexpr uint32_t odesystemsize_ = 3;

    static constexpr double sigma_  = 0.5;   // std of random movement

    // The internal state of the cell
    // Because we use the GSL, the internal state of the system must be an array
    // of double values.
  double* internal_state_;
  double* mRNA_array_;
  double* Protein_array_;
  int PhantomJumpCounts_= 0;
  int* TrueJumpCounts_array_;
  double* ThinningParam_array_;
  double* initialV_;
  double* initval1;

  double* GenesInteractionsMatrix_;
  double* KinParam_;
  double Duration_APC = 0.;
  double Duration_TCC = 0.;
  double Remember_encounter_APC = 0.;
  double Remember_encounter_T = 0.;
  double* Remember_division_ ;
  double count_division = 0. ;
  int Number_Tcells = 0 ;
  


 double count_death = 0.;
 
  
  int Number_Of_Genes_;
  int Number_Of_Parameters_ = 13;

  int ithGene_;
  double Time_NextJump_ = 0.;

  std::vector<double>  phylogeny_t_;
  std::vector<int>  phylogeny_id_;

//    std::vector<double>  initialV_;
  const char phylogeny_T_filename[16] = "phylogeny_T.txt";
  const char phylogeny_ID_filename[32] = "phylogeny_ID.txt";
  FILE *phylogeny_T_file_, *phylogeny_ID_file_;
    bool isMitotic_ = false; // distinguishes cells about to divide

    double mitotic_threshold_ = 0.;
    double cycling_threshold_ = 0.;
    double dividing_threshold_ = 0.;
    double death_threshold_ = 0.;
    
    double mitotic_threshold_P0_=0.;

    // Death/Survival parameters
    double k_survival = 10.0; // higher values = higher survival.

  

    // Cell cycle parameters
    // The values are taken from (Battogtokh & Tyson 2006 PRE)
    static constexpr double cy_j = 0.05;
    static constexpr double cy_p = 0.15;
    static constexpr double cy_k1 = 0.002;
    static constexpr double cy_k2 = 0.0795;
    static constexpr double cy_k3 = 0.01;
    static constexpr double cy_k4 = 2.0;
    static constexpr double cy_k5 = 0.05;
    static constexpr double cy_k6 = 0.04;
    static constexpr double cy_k7 = 1.5;
    static constexpr double cy_k8 = 0.19;
    static constexpr double cy_k9 = 0.64;
    static constexpr double cy_k10 = 0.0025;
    static constexpr double cy_k11 = 0.07;
    static constexpr double cy_k12 = 0.08;
    static constexpr double cy_A = 0.52;

private:
    /** dummy attribute - allows to register class in Simuscale statically */
    static bool registered_;
};

#endif // SIMUSCALE_Lymphocyte_H__
