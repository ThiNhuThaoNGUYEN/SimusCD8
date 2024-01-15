// ============================================================================
//                                   Includes
// ============================================================================
#include "Lymphocyte.h"
#include <vector>
using std::vector;
#include <iostream>
#include <random>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <memory>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "Simulation.h"
#include "Alea.h"
#include "Coordinates.h"
#include "Cell.h"
using std::unique_ptr;
using std::shared_ptr;
#include <sstream>
#include <string>
#include <fstream>
using std::ifstream;
using namespace std;

// ============================================================================
//                       Definition of static attributes
// ============================================================================
constexpr CellFormalism Lymphocyte::classId_;
constexpr char Lymphocyte::classKW_[];

constexpr uint32_t Lymphocyte::odesystemsize_;
// ============================================================================
//                                Constructors
// ============================================================================

Lymphocyte::Lymphocyte(const Lymphocyte &model) :
Cell(model),
isMitotic_(model.isMitotic_) {
    phylogeny_id_.push_back(this->id());
    phylogeny_t_.push_back(Simulation::sim_time());
    internal_state_ = new double[odesystemsize_];
    memcpy(internal_state_, model.internal_state_, odesystemsize_ * sizeof (*internal_state_));

    this->Number_Of_Genes_ = model.Number_Of_Genes_;
   
     
    mRNA_array_ = new double[Number_Of_Genes_];
    Protein_array_ = new double[Number_Of_Genes_+1];
    TrueJumpCounts_array_ = new int[Number_Of_Genes_];
    ThinningParam_array_ = new double[Number_Of_Genes_];
     memcpy(mRNA_array_, model.mRNA_array_, Number_Of_Genes_ * sizeof (*mRNA_array_));
     memcpy(Protein_array_, model.Protein_array_, (Number_Of_Genes_+1) * sizeof (*Protein_array_));
     memcpy(TrueJumpCounts_array_, model.TrueJumpCounts_array_, Number_Of_Genes_ * sizeof (*TrueJumpCounts_array_));
     memcpy(ThinningParam_array_, model.ThinningParam_array_, Number_Of_Genes_ * sizeof (*ThinningParam_array_));

    
    KinParam_ = new double[Number_Of_Parameters_];
    GenesInteractionsMatrix_ = new double[Number_Of_Genes_ * Number_Of_Genes_];

    memcpy(KinParam_, model.KinParam_, Number_Of_Parameters_ * sizeof (*KinParam_));
    memcpy(GenesInteractionsMatrix_, model.GenesInteractionsMatrix_, Number_Of_Genes_ * Number_Of_Genes_ * sizeof (*GenesInteractionsMatrix_));

}

/**
 * Create a new object with the provided constitutive elements and values
 */
Lymphocyte::Lymphocyte(CellType cellType,
        const MoveBehaviour& move_behaviour,
        const Coordinates<double>& pos,
        double initial_volume,
        double volume_min,
        double doubling_time) :
  Cell(cellType, std::move(move_behaviour),
pos, CellSize(initial_volume, volume_min),
doubling_time) {
    internal_state_ = new double[odesystemsize_];
     phylogeny_id_.push_back(this->id());
     phylogeny_t_.push_back(Simulation::sim_time());


    get_GeneParams();
    get_InitValues();
    
    
    mRNA_array_ = new double[Number_Of_Genes_];
    Protein_array_ = new double[Number_Of_Genes_+1];
    TrueJumpCounts_array_ = new int[Number_Of_Genes_];
    ThinningParam_array_ = new double[Number_Of_Genes_];
   
    
    for (u_int32_t j = 0; j < odesystemsize_; j++){
      internal_state_[j] = 0.;}

    double thenorm = 0.05;

    if (cellType == NICHE) {
          thenorm = 0.;
    }


//if there are n T-cells at day  0

if (this->id() <= Number_Tcells) { //id =1,..

 for (u_int32_t j = 0; j < Number_Of_Genes_; j++) {
if( j < Number_Of_Genes_ ){
        mRNA_array_[j] = initialV_[j+2+Number_Of_Genes_+(this->id()-1)*(Number_Of_Genes_*2+2)];
        Protein_array_[j] = initialV_[j+2+(this->id()-1)*(Number_Of_Genes_*2+2)];
} else {
  
        mRNA_array_[j] = 0.;
        Protein_array_[j] = 0.;
}

        TrueJumpCounts_array_[j] = 0;
        ThinningParam_array_[j] = thenorm * Alea::random();
     }
} 
 else{ 

    for (u_int32_t j = 0; j < Number_Of_Genes_; j++) {
        mRNA_array_[j] = 0.;
        Protein_array_[j] = 0.;
        TrueJumpCounts_array_[j] = 0;
        ThinningParam_array_[j] = thenorm * Alea::random();
     }

}

delete [] initialV_;



    Protein_array_[Number_Of_Genes_] = 0.;
    PhantomJumpCounts_ = 0;


    internal_state_[APC_Encounters] = 0.;
    internal_state_[TCC] = 0.;
    
  }

/**
 * Restore an object that was backed-up in backup_file
 */
Lymphocyte::Lymphocyte(gzFile& backup) :
Cell(backup) {

    internal_state_ = new double[odesystemsize_];
    Load(backup);
}

// ============================================================================
//                                 Destructor
// ============================================================================

Lymphocyte::~Lymphocyte() noexcept {
    delete [] internal_state_;

    delete [] mRNA_array_;
    delete [] Protein_array_;
    delete [] TrueJumpCounts_array_;
    delete [] ThinningParam_array_;

    delete [] GenesInteractionsMatrix_;
    delete [] KinParam_;
    
}

// ============================================================================
//                                   Methods
// ============================================================================

void Lymphocyte::InternalUpdate(const double& dt) {
    // Describe what happens inside the cell when taking a time step
    //std::cout << "Internal update" << std::endl;
    // Update the intracellular state of the cell

   if (cell_type_ != CellType::NICHE) {
    ODE_update(dt);
      
    }    // If the cell is dying, don't take any further actions
   if (isDying()) return;

    // Make the cell grow if it's cycling (and not dying)
    Grow(dt);
   }

void Lymphocyte::UpdatePhylogeny(vector<int> phy_id, vector<double> phy_t, int sz) {
    phylogeny_id_.insert(phylogeny_id_.begin(), phy_id.begin(), phy_id.end());
    phylogeny_t_.insert(phylogeny_t_.begin(), phy_t.begin(), phy_t.end());
}
void Lymphocyte::SetNewProteinLevel1(const double K,const double MotherProteins[],const int sz)const {
    for (u_int32_t j = 0; j < sz; j++) {
      Protein_array_[j] =  K * MotherProteins[j]/2.;
    }
}
void Lymphocyte::SetNewProteinLevel2(const  double K,const double MotherProteins[],const int sz)const {
    for (u_int32_t j = 0; j < sz; j++) {
      Protein_array_[j] = (2.0 - K) * MotherProteins[j]/2.;
    }
}

void Lymphocyte::SetNewRNALevel1(const double K,const double MotherRNA[],const int sz)const {
    for (u_int32_t j = 0; j < sz; j++) {
     mRNA_array_[j] =  K * MotherRNA[j]/2.;
    }
}
void Lymphocyte::SetNewRNALevel2(const  double K,const double MotherRNA[],const int sz)const {
    for (u_int32_t j = 0; j < sz; j++) {
      mRNA_array_[j] = (2.0 -  K) * MotherRNA[j]/2.;
    }
}
void Lymphocyte::UpdateContactAPC(double Mother_APC_contact) {
       Duration_APC =  Mother_APC_contact;    
}

void Lymphocyte::UpdateContactTCC(double Mother_TCC_contact) {
       Duration_TCC =  Mother_TCC_contact;    
}

void Lymphocyte::Update_count_division( double Mother_division,const int s )const {
   
  Protein_array_[s]  = Mother_division;
    
  }

void Lymphocyte::Update_count_division_mother( double Mother_division ) {
Mother_division += 1.;
} 


int num_division = 0;

Cell* Lymphocyte::Divide() {

    // Create a copy of the cell with no mechanical forces nor chemical stimuli
    // applied to it
    Lymphocyte* newCell = new Lymphocyte(*this);

    SeparateDividingCells(this, newCell);

    std::cout << "Cell " << this->id() << " is DIVIDING!!" <<  " at t=" << Simulation::sim_time() << std::endl;    ///remember dividin for the cell

   Protein_array_[Number_Of_Genes_] += 1.;
    newCell->Update_count_division(Protein_array_[Number_Of_Genes_], Number_Of_Genes_);
    newCell->UpdateContactAPC(Duration_APC);

    num_division += 1.; 
    
     double m = 20.;

     double a = 10.;
     double K_p, K_m;
     double randomCell;
   

    K_p = 1.0 - (m / 100.)* Alea::gaussian_random_0_2();
    K_m = 1.0 - (a / 100.)* Alea::gaussian_random_0_2();


 if ( K_p >= 2. || K_m >= 2. || K_p < 0. || K_m < 0.) { 
cout << "negative,**************************************** "  << K_p<< K_m << endl;
 K_p = 1.0 - (m / 100.)* Alea::gaussian_random_0_2();
 K_m = 1.0 - (a / 100.)* Alea::gaussian_random_0_2();

 }

 //partitioning of molecular content into daughter cells
 
    randomCell = Alea::random();
    if (randomCell < 0.5){
      
    newCell->SetNewProteinLevel1(K_p, Protein_array_, Number_Of_Genes_);
    newCell->SetNewRNALevel1(K_m, mRNA_array_, Number_Of_Genes_);

    this->SetNewProteinLevel2(K_p, Protein_array_, Number_Of_Genes_);
    this->SetNewRNALevel2(K_m, mRNA_array_, Number_Of_Genes_);
    } else
      
    {
    newCell->SetNewProteinLevel2(K_p, Protein_array_, Number_Of_Genes_);
    newCell->SetNewRNALevel2(K_m, mRNA_array_, Number_Of_Genes_);

    this->SetNewProteinLevel1(K_p, Protein_array_, Number_Of_Genes_);
    this->SetNewRNALevel1(K_m, mRNA_array_, Number_Of_Genes_);
    
      }


    
   newCell->UpdatePhylogeny(phylogeny_id_, phylogeny_t_, phylogeny_id_.size());

   //records the ID of the proliferating cell
   
    phylogeny_id_.push_back(this->id());
 
    phylogeny_t_.push_back(Simulation::sim_time());

  
  if ( ( phylogeny_ID_file_ = fopen(phylogeny_ID_filename,"a") ) != NULL ) {
    fprintf(phylogeny_ID_file_, "%d -> %d , t= %f \n", this->id(),  newCell-> id(),  Simulation::sim_time());
      fclose(phylogeny_ID_file_);
    }
    if ( ( phylogeny_T_file_ = fopen(phylogeny_T_filename,"a") ) != NULL ) {
      fprintf(phylogeny_T_file_, "%f \n ", Simulation::sim_time());
      fclose(phylogeny_T_file_);
    }
    
    return newCell;
}


  


void Lymphocyte::Save(gzFile backup) const {
    // Write my classId
    gzwrite(backup, &classId_, sizeof (classId_));

    // Write the generic cell stuff
    Cell::Save(backup);

    // Write the specifics of this class
    //  int8_t isMitotic = isMitotic_ ? 1 : 0;
    //  gzwrite(backup_file, &isMitotic, sizeof(isMitotic));
        // <TODO> write your specific attributes </TODO>

    for (u_int32_t i = 0; i < odesystemsize_; i++){
        gzwrite(backup, &internal_state_[i], sizeof (internal_state_[i]));
        
    }

}

void Lymphocyte::Load(gzFile& backup) {
    // <TODO> Load your specific attributes </TODO>
    
    for (u_int32_t i = 0; i < odesystemsize_; i++)
        gzread(backup, &internal_state_[i], sizeof (internal_state_[i]));
}


double Lymphocyte::get_GeneParams(void) {
    KinParam_ = new double[Number_Of_Parameters_];

    ////////////////////////////////////
    ifstream indatakin; // indata is like cin

string filename =Simulation::input_dir_+"kineticsparam.txt";

    indatakin.open(filename); // opens the file
    if (!indatakin) { // file couldn't be opened
        std::cerr << "Error: file kineticsparam.txt could not be opened" << std::endl;
        exit(1);
    }
    int iter = 0;
    std::string line;
    std::getline(indatakin, line);
    std::getline(indatakin, line);
    std::istringstream iss(line);
    double a;

    while (iter < Number_Of_Parameters_) {
        if (!(iss >> a)) {
            std::cerr << "Error: file kineticsparam.txt not readable - possible value missing" << std::endl;
            break;
        }// error
        else {
            KinParam_[iter] = a;
            iter = iter + 1;
        }
    }

    indatakin.close();

    ifstream indataf; // indata is like cin
    
    string filename2 =Simulation::input_dir_+"GeneInteractionsMatrix_3decoratings_9G.txt";

    indataf.open(filename2); // opens the file  

    if (!indataf) { // file couldn't be opened
        std::cerr << "Error: file GeneInteractionsMatrix.txt could not be opened" << std::endl;
        exit(1);
    }
    Number_Of_Genes_ = 0;
    std::getline(indataf, line);
    std::stringstream ss(line);
    std::string word;
    int i = 0;
    int j = 0;
   
    while (ss >> word) {
        Number_Of_Genes_ += 1;
    }

     Number_Of_Genes_ -= 1; // to include title column

    GenesInteractionsMatrix_ = new double[Number_Of_Genes_ * Number_Of_Genes_];


     while (std::getline(indataf, line)) {
        std::istringstream iss(line);
        iss >> word;
        j = 0;
        while (j < Number_Of_Genes_) {
            if ((!(iss >> a))) {
                std::cerr << "Error: possible value missing" << std::endl;
                break;
            }// error
            else {
                GenesInteractionsMatrix_[(j + Number_Of_Genes_ * i)] = a;
                j = j + 1;
            }
        }
    i = i + 1;

        // process pair (a,b)
    }

  
    indataf.close();

    return 0;
}



double Lymphocyte::get_InitValues() {
std::vector<double>  initialV0;

ifstream indataV; // indata

string filename3 =Simulation::input_dir_+"Day_1_3decoratings_22.txt";

 indataV.open(filename3); // opens the file
  if (!indataV) { // file couldn't be opened
      //  std::cerr << "Error: file Day_1.txt could not be opened" << std::endl;
  initialV_ =new double[Number_Of_Genes_*(Number_Of_Genes_)];

  for (u_int32_t i = 0; i < Number_Tcells*(Number_Of_Genes_*2+2); i++){
    initialV_[i]= 0.;
    indataV.close();
    }

  }
    std::string line2;
    double val;
    int c = 0;
    int l = 0;

while (std::getline(indataV, line2)) {
 Number_Tcells += 1;
        std::istringstream Vss(line2);
        l = 0;
        while ( l < (Number_Of_Genes_*2+2)) {
            if ((!(Vss >> val))) {
                std::cerr << "Error: Value missing in day_3/character" << std::endl;
                break;
            }// error
            else {
        initialV0.push_back(val);      
         l = l + 1;
            }
        }
        c = c + 1;
        }

    initialV_ =new double[Number_Tcells*(Number_Of_Genes_*2+2)];
    indataV.close();
    
    for (u_int32_t i = 0; i < Number_Tcells*(Number_Of_Genes_*2+2); i++){
    initialV_[i]=initialV0.at(i);
     }

    return 0.0;
}


double Lymphocyte::get_P(int i) const {
    if (i > Number_Of_Genes_ - 1) {
        i = Number_Of_Genes_ - 1;
    }
    return (Protein_array_[i]);
}


double Lymphocyte::getPhylogeny_t(int i) const {
    double ans;
    if (i <0) {
        ans = phylogeny_t_.size();
    }
    else {
        ans = phylogeny_t_[i];
    }
    return ans;
}
int Lymphocyte::getPhylogeny_id(int i) const {
    int ans;
    if (i <0) {
        ans = phylogeny_id_.size();
    }
    else {
        ans = phylogeny_id_[i];
    }
    return ans;
}




double Lymphocyte::get_mRNA(int i) const {
    if (i > Number_Of_Genes_ - 1) {
        i = Number_Of_Genes_ - 1;
    }
    return (mRNA_array_[i]);
}
double Lymphocyte::get_output(InterCellSignal signal) const {
    // <TODO> Determine the amount of <signal> the cell is secreting </TODO>
       switch (signal) {
	  case InterCellSignal::LYMPHOCYTE_TYPE:
            return 1.0 * (cell_type_ == CellType::LYMPHOCYTE)
                    + 2.0 * (cell_type_ == CellType::APC);
	  case InterCellSignal::APC_CONTACT:
            return internal_state_[APC_Encounters];
	  case InterCellSignal::APC_DURATION:
	      return  Duration_APC;
	  case InterCellSignal::T_DURATION:
            return  Duration_TCC;
      case InterCellSignal::VOLUME:
            return outer_volume();
	  case InterCellSignal::LYMPHOCYTE_P1:
            return Protein_array_[0];
      case InterCellSignal::LYMPHOCYTE_P2:
            return  Protein_array_[1];
      case InterCellSignal::LYMPHOCYTE_P3:
	        return Protein_array_[2];
      case InterCellSignal::LYMPHOCYTE_P4:
	      return Protein_array_[3];
	  case InterCellSignal::LYMPHOCYTE_P5:
            return Protein_array_[4];
      case InterCellSignal::LYMPHOCYTE_P6:
            return  Protein_array_[5];
      case InterCellSignal::LYMPHOCYTE_P7:
            return  Protein_array_[6];
      case InterCellSignal::LYMPHOCYTE_P8:
            return  Protein_array_[7];
    case InterCellSignal::LYMPHOCYTE_P9:
            return  Protein_array_[8];
case InterCellSignal::LYMPHOCYTE_P1_1:
       return  Protein_array_[9];
    case InterCellSignal::LYMPHOCYTE_P2_1:
       return  Protein_array_[10];
    case InterCellSignal::LYMPHOCYTE_P3_1:
       return  Protein_array_[11];
    case InterCellSignal::LYMPHOCYTE_P4_1:
       return  Protein_array_[12];
    case InterCellSignal::LYMPHOCYTE_P5_1:
       return  Protein_array_[13];
    case InterCellSignal::LYMPHOCYTE_P6_1:
       return  Protein_array_[14];
   case InterCellSignal::LYMPHOCYTE_P7_1:
        return  Protein_array_[15];
   case InterCellSignal::LYMPHOCYTE_P8_1:
       return  Protein_array_[16];
    case InterCellSignal::LYMPHOCYTE_P9_1:
       return  Protein_array_[17];
    case InterCellSignal::LYMPHOCYTE_P1_2:
       return  Protein_array_[18];
    case InterCellSignal::LYMPHOCYTE_P2_2:
       return  Protein_array_[19];
    case InterCellSignal::LYMPHOCYTE_P3_2:
       return  Protein_array_[20];
    case InterCellSignal::LYMPHOCYTE_P4_2:
       return  Protein_array_[21];
    case InterCellSignal::LYMPHOCYTE_P5_2:
       return  Protein_array_[22];
    case InterCellSignal::LYMPHOCYTE_P6_2:
       return  Protein_array_[23];
    case InterCellSignal::LYMPHOCYTE_P7_2:
       return  Protein_array_[24];
    case InterCellSignal::LYMPHOCYTE_P8_2:
       return  Protein_array_[25];
    case InterCellSignal::LYMPHOCYTE_P9_2:
       return  Protein_array_[26];
    case InterCellSignal::LYMPHOCYTE_P1_3:
       return  Protein_array_[27];
    case InterCellSignal::LYMPHOCYTE_P2_3:
       return  Protein_array_[28];
    case InterCellSignal::LYMPHOCYTE_P3_3:
       return  Protein_array_[29];
    case InterCellSignal::LYMPHOCYTE_P4_3:
       return  Protein_array_[30];
    case InterCellSignal::LYMPHOCYTE_P5_3:
       return  Protein_array_[31];
    case InterCellSignal::LYMPHOCYTE_P6_3:
       return  Protein_array_[32];
      case InterCellSignal::LYMPHOCYTE_P7_3:
       return  Protein_array_[33];
      case InterCellSignal::LYMPHOCYTE_P8_3:
       return  Protein_array_[34];
      case InterCellSignal::LYMPHOCYTE_P9_3:
       return  Protein_array_[35];
      case InterCellSignal::REMEMBER_DIVISION:
            return  Protein_array_[36]; 
      case InterCellSignal::LYMPHOCYTE_mRNA1:
            return mRNA_array_[0];
      case InterCellSignal::LYMPHOCYTE_mRNA2:
            return mRNA_array_[1];   
      case InterCellSignal::LYMPHOCYTE_mRNA3:
	        return  mRNA_array_[2];
      case InterCellSignal::LYMPHOCYTE_mRNA4:
         	return mRNA_array_[3];
      case InterCellSignal::LYMPHOCYTE_mRNA5:
            return mRNA_array_[4];
      case InterCellSignal::LYMPHOCYTE_mRNA6:
            return mRNA_array_[5];   
      case InterCellSignal::LYMPHOCYTE_mRNA7:
            return  mRNA_array_[6];
      case InterCellSignal::LYMPHOCYTE_mRNA8:
             return  mRNA_array_[7]; 
      case InterCellSignal::LYMPHOCYTE_mRNA9:
             return  mRNA_array_[8];
      case InterCellSignal::LYMPHOCYTE_mRNA1_1:
       return  mRNA_array_[9];
      case InterCellSignal::LYMPHOCYTE_mRNA2_1:
       return  mRNA_array_[10];
      case InterCellSignal::LYMPHOCYTE_mRNA3_1:
       return  mRNA_array_[11];
      case InterCellSignal::LYMPHOCYTE_mRNA4_1:
       return  mRNA_array_[12];
      case InterCellSignal::LYMPHOCYTE_mRNA5_1:
       return  mRNA_array_[13];
      case InterCellSignal::LYMPHOCYTE_mRNA6_1:
       return  mRNA_array_[14];
      case InterCellSignal::LYMPHOCYTE_mRNA7_1:
        return  mRNA_array_[15];
      case InterCellSignal::LYMPHOCYTE_mRNA8_1:
       return  mRNA_array_[16];
      case InterCellSignal::LYMPHOCYTE_mRNA9_1:
       return  mRNA_array_[17];
      case InterCellSignal::LYMPHOCYTE_mRNA1_2:
       return  mRNA_array_[18];
      case InterCellSignal::LYMPHOCYTE_mRNA2_2:
       return  mRNA_array_[19];
      case InterCellSignal::LYMPHOCYTE_mRNA3_2:
       return  mRNA_array_[20];
      case InterCellSignal::LYMPHOCYTE_mRNA4_2:
       return  mRNA_array_[21];
      case InterCellSignal::LYMPHOCYTE_mRNA5_2:
       return  mRNA_array_[22];
      case InterCellSignal::LYMPHOCYTE_mRNA6_2:
       return  mRNA_array_[23];
      case InterCellSignal::LYMPHOCYTE_mRNA7_2:
       return  mRNA_array_[24];
      case InterCellSignal::LYMPHOCYTE_mRNA8_2:
       return  mRNA_array_[25];
      case InterCellSignal::LYMPHOCYTE_mRNA9_2:
       return  mRNA_array_[26];
      case InterCellSignal::LYMPHOCYTE_mRNA1_3:
       return  mRNA_array_[27];
      case InterCellSignal::LYMPHOCYTE_mRNA2_3:
       return  mRNA_array_[28];
      case InterCellSignal::LYMPHOCYTE_mRNA3_3:
       return  mRNA_array_[29];
      case InterCellSignal::LYMPHOCYTE_mRNA4_3:
       return  mRNA_array_[30];
      case InterCellSignal::LYMPHOCYTE_mRNA5_3:
       return  mRNA_array_[31];
      case InterCellSignal::LYMPHOCYTE_mRNA6_3:
       return  mRNA_array_[32];
      case InterCellSignal::LYMPHOCYTE_mRNA7_3:
       return  mRNA_array_[33];
      case InterCellSignal::LYMPHOCYTE_mRNA8_3:
       return  mRNA_array_[34];
      case InterCellSignal::LYMPHOCYTE_mRNA9_3:
       return  mRNA_array_[35];

      case InterCellSignal::LYMPHOCYTE_CONTACT:
	        return internal_state_[T_Encounters];//internal_state_[TCC];
         
       default:
      return 0.0;  // if signal is unknown, return 0.0, not a failure
    }
}



bool Lymphocyte::isDividing() const { 

    if (Protein_array_[(Number_Of_Genes_)] > 0.){
        return (Protein_array_[1] >= 0.4 &&
	    size_.may_divide());}
    
       
  else { return (Protein_array_[0] >= 0.5 &&
           size_.may_divide());}
    
}

bool Lymphocyte::isDying() const {
    // <TODO> Should the cell die now ? </TODO>
return (cell_type_ == LYMPHOCYTE && (Protein_array_[2] >= 0.75)) ||(cell_type_ == APC && ((12.0*24.+ this->pos_x()*2. + this->pos_y()*2. + this->pos_z()/2.) < Simulation::sim_time() ))  ;}

void Lymphocyte::Get_Sigma( double *Sigma_i ,double * P, double Duration_APC, double Tcc, const double x,  bool AcceptNegative ) {

double Parr[Number_Of_Genes_];
 memcpy(Parr, P, Number_Of_Genes_ * sizeof (*P));

 double interIJ;
 double initval_i[Number_Of_Genes_];


  for (u_int32_t i = 0; i < Number_Of_Genes_; i++) {
// Basal activity for other genes 
   initval_i[i] = -x/2.;
//APC signal activates gene 1 
if ( i == 0){

   if (Duration_APC >= 10. && Protein_array_[(Number_Of_Genes_)] < 1.) {initval_i[i] += 2.*x; }
 
 }

// TCC activates gene3
 else if (i == 2){
   if (Tcc > 0. ){
     initval_i[i] += 2.*x;
}

}

 for (u_int32_t j = 0; j < Number_Of_Genes_; j++) {
   interIJ = GenesInteractionsMatrix_[i + Number_Of_Genes_ * j]; // J acts on I
    if (((interIJ > 0.)||AcceptNegative)) {      
      initval_i[i] += Parr[j] * interIJ;
      }
 }
}
 for (u_int32_t i = 0; i < Number_Of_Genes_; i++) {
Sigma_i[i] = (1. / (1. + exp(-initval_i[i])));
}
}


void Lymphocyte::Intracellular_ExactEvol(double thetime, double * P, double * M, double * S1) {
        // Default degradation rates
            double d0 = KinParam_[0];  // mRNA degradation rates
            double D2 = KinParam_[1];  //protein degradation rates


            double D1 = KinParam_[5];  //protein degradation rates * gene1

            double D3 = KinParam_[6];
            double D4 = KinParam_[7];  //protein degradation rates * gene4
            double D7 = KinParam_[8];  //protein degradation rates * gene7
            double D8 = KinParam_[9];
            double D5 = KinParam_[10];
            double D6 = KinParam_[11];
            double D9 = KinParam_[12]; //protein degradation rates * gene9
    double Marr[Number_Of_Genes_];
    double Parr[Number_Of_Genes_];

    memcpy(Marr, M, Number_Of_Genes_ * sizeof (*M));
    memcpy(Parr, P, Number_Of_Genes_ * sizeof (*P));
    // New values:
    for (u_int32_t i = 0; i < Number_Of_Genes_; i++) {
      
      if (i == 0 ){ 
        M[i] = Marr[i] * exp(-thetime * d0);
        P[i] = ((S1[i] / (d0 - D1)) * Marr[i]*(exp(-thetime * D1) - exp(-thetime * d0)) + Parr[i] * exp(-thetime * D1));
      }
     else if ( i == 3 ){
        M[i] = Marr[i] * exp(-thetime * d0);
        P[i] = ((S1[i] / (d0 - D4)) * Marr[i]*(exp(-thetime * D4) - exp(-thetime * d0))+ Parr[i] * exp(-thetime * D4));
      }
      else if (i == 6){
        M[i] = Marr[i] * exp(-thetime * d0);
        P[i] = ((S1[i] /(d0 - D7)) * Marr[i]*(exp(-thetime * D7) - exp(-thetime * d0))+ Parr[i] * exp(-thetime * D7));
      }
      
     else if (i == 7){
        M[i] = Marr[i] * exp(-thetime * d0);
        P[i] = ((S1[i] /(d0 - D8)) * Marr[i]*(exp(-thetime * D8) - exp(-thetime * d0)) + Parr[i] * exp(-thetime * D8));
      }
     else if (i == 8){
        M[i] = Marr[i] * exp(-thetime * d0);
        P[i] = ((S1[i] /(d0 - D9)) * Marr[i]*(exp(-thetime * D9) - exp(-thetime * d0)) + Parr[i] * exp(-thetime * D9));
      }
 else if (i == 4 ){
        M[i] = Marr[i] * exp(-thetime * d0);
        P[i] = ((S1[i] /(d0 - D5)) * Marr[i]*(exp(-thetime * D5) - exp(-thetime * d0)) + Parr[i] * exp(-thetime *  D5));
      }
else if  (i == 5) {
        M[i] = Marr[i] * exp(-thetime * d0);
        P[i] = ((S1[i] /(d0 - D6)) * Marr[i]*(exp(-thetime * D6) - exp(-thetime * d0)) + Parr[i] * exp(-thetime *  D6));
      }
 else if (i == 2 ){
        M[i] = Marr[i] * exp(-thetime * d0);
        P[i] = ((S1[i]/(d0 - D3)) * Marr[i]*(exp(-thetime * D3) - exp(-thetime * d0)) + Parr[i] * exp(-thetime *  D3));
      }

     else if (i == 1 ){
        M[i] = Marr[i] * exp(-thetime * d0);
        P[i] = ((S1[i] / (d0 - D2)) * Marr[i]*(exp(-thetime * D2) - exp(-thetime * d0))	+ Parr[i] * exp(-thetime * D2));}
    
else {
        M[i] = Marr[i] * exp(-thetime * d0);
        P[i] = ((S1[i] / (d0 - D1)) * Marr[i]*(exp(-thetime * D1) - exp(-thetime * d0)) + Parr[i] * exp(-thetime * D1));
}
}
}
/*
void Lymphocyte::UpdateMitoticStatus() {
    // Once a cell enters the mitotic state, it remains mitotic until it
    // dies or divides
  
    if (not isMitotic_)
    {
      isMitotic_ = (isCycling_ && ((Protein_array_[0] > 0. )));
	// && (cell_type_==APC || (get_P(0)>random_mitotic_threshold()));
        //if APC, condition for division is only on cycle variable CY_X
        //otherwise, type = Lymphocyte -> division depends on protein content P0
    }
}


void Lymphocyte::UpdateCyclingStatus() {
    if (isCycling_) {
        if (StopCycling()) isCycling_ = false;
    } else if (StartCycling()) {
        isCycling_ = true;
    }

}

bool Lymphocyte::StopCycling() {
    // <TODO> Should the cell stop cycling now ? </TODO>
    return (cell_type_ == CellType::NICHE or get_P(0) < 0.);
}

bool Lymphocyte::StartCycling() {
    // <TODO> Should the cell start cycling now ? </TODO>
    return cell_type_ != NICHE && // internal_state_[CY_X] > cycling_threshold_ &&
             (Protein_array_[0]>= 0.) &&
            ((cell_type_ == STEM && getInSignal(InterCellSignal::NICHE) > 0) ||
            cell_type_ == CANCER ||
            cell_type_ == LYMPHOCYTE ||
            cell_type_ == APC);

}
*/

// movement of T cell
Coordinates<double> Lymphocyte::MotileDisplacement(const double& dt) {
  
  std::vector<Cell*> neighb = neighbours();
  auto cell_it = neighb.begin();
  bool contact_APC = false;
  bool contact_T = false;
   while ( cell_it != neighb.end() ) {
     //if T cell contacts APC
    if ( (*cell_it)->cell_type() == APC ) {
      contact_APC = true;
      break;
    }
    //if T cell contacts other T cell
    if ( (*cell_it)->cell_type() == LYMPHOCYTE ) {
      contact_T = true;
      break;
    }
    cell_it++;
   }
   double sigma_t;
   if ( contact_APC ){
     if (Protein_array_[Number_Of_Genes_] < 1.) {
      sigma_t = 0.02;}
     else    sigma_t = 5.;

    }
   else if ( contact_T ) {
    sigma_t = 0.8;
    }
 else {
 sigma_t = 5.;
     }
  Coordinates<double> displ { sqrt(dt)*sigma_t*Alea::gaussian_random(),
                              sqrt(dt)*sigma_t*Alea::gaussian_random(),
                              sqrt(dt)*sigma_t*Alea::gaussian_random() };
  return displ;
}




void Lymphocyte::ODE_update(const double& dt){

    
double t = 0., t1 = dt;

   switch (cell_type_) {
        case LYMPHOCYTE:
        {

      internal_state_[APC_Encounters] = getInSignal(InterCellSignal::APC_CONTACT);   
        if ( internal_state_[APC_Encounters]> 0.){
         Duration_APC += dt;
         Remember_encounter_APC += internal_state_[APC_Encounters];
        }

    //Lymphocyte contacts effector

std::vector<Cell*> neighb = neighbours();
  auto cell_it = neighb.begin();
  bool contact_E = false;
  while ( cell_it != neighb.end() ) {
    if ( (*cell_it)->cell_type() == LYMPHOCYTE && (Duration_APC > 0. && Protein_array_[7] < 0.02 && Protein_array_[6] >= 0.01)) {

 contact_E = true;
    break;
    }
    cell_it++;
   }
  // have TCC signaling because it contacted the effector cell
   if ( contact_E ) {
     internal_state_[T_Encounters] = getInSignal(InterCellSignal::LYMPHOCYTE_CONTACT);
     internal_state_[TCC] = std::min(getInSignal(InterCellSignal::LYMPHOCYTE_CONTACT),1.);
    } 
	     
neighb.clear();

	 if ( internal_state_[TCC]> 0.){
         Duration_TCC +=  dt;
	 Remember_encounter_T +=  internal_state_[TCC];
	 }
	       
            // Here: PDMP (intracellular signalling)*
            // Default bursting parameters
            double a0 = KinParam_[2]; //0;
            double a1 = KinParam_[3];/// change to normalize protein
            double a2 = KinParam_[4]; 
            // Default degradation rates
            double d0 = KinParam_[0];  // mRNA degradation rates
            double D2 = KinParam_[1];  //protein degradation rates
	    
            //double d01 = KinParam_[5];  // mRNA degradation rates * gene1
            double D1 = KinParam_[5];  //protein degradation rates * gene1
	    
            double D3 = KinParam_[6]; 
	    //******test  model 4genes
	        double D4 = KinParam_[7];  //protein degradation rates * gene4
	        double D7 = KinParam_[8];  //protein degradation rates * gene7
            double D8 = KinParam_[9];
            double D5 = KinParam_[10];
            double D6 = KinParam_[11];
            double D9 = KinParam_[12]; //protein degradation rates * gene9
            const double K0 = a0*d0;
            const double K1 = a1*d0;

            bool jump = false;

double *S1 = new double[Number_Of_Genes_]; //init
 for (u_int32_t i = 0; i < Number_Of_Genes_; i++) {
           S1[i] = 0.;
         }	   
            const double r = a2;
            
	    
            double currentTime = 0.;
            double DeltaT = 0.;
            double *PPmax = new double[Number_Of_Genes_]; //init
          double *ProbaArray = new double[Number_Of_Genes_+1]; //init

    for (u_int32_t i = 0; i < Number_Of_Genes_; i++) {
	       PPmax[i] = 0.;
	     }
	     for (u_int32_t i = 0; i < Number_Of_Genes_+1; i++) {
	       ProbaArray[i] = 0.;
	     }

	  double *KKon = new double[Number_Of_Genes_]; //init
	  for (u_int32_t i = 0; i < Number_Of_Genes_; i++) {
	    KKon[i] = 0.;
	  }

for (u_int32_t i = 0; i < Number_Of_Genes_; i++) {
          if ( i == 0){
        S1[i] = d0 * D1 * a2 / K1;}
          else  if  ( i == 3){
            S1[i] = d0 * D4 * a2 / K1;}
          else if ( i == 6){
            S1[i] = d0 * D7 * a2 / K1;}
  else if ( i == 7){
            S1[i] = d0 * D8 * a2 / K1;}
 else if ( i == 8){
            S1[i] = d0 * D9 * a2 / K1;}
  else if ( i == 4){
            S1[i] =d0 * D5 * a2 / K1;}
  else if (i == 5){
            S1[i] =d0 * D6 * a2 / K1;}
  else if ( i==2 ){
            S1[i] =d0 * D3 * a2 / K1;}
else if ( i==1 ){
            S1[i] =d0 * D2 * a2 / K1;}     
else {
            S1[i] =d0 * D1 * a2 / K1;}
 }
            while (Time_NextJump_ < t1){

    double *Sigma = new double[Number_Of_Genes_]; //init
    
      for (u_int32_t i = 0; i < Number_Of_Genes_; i++) {
     Sigma[i] = 0.;
      }
 
 double Tau = 0.;
// -------------------------- ------------------------------- ------------------
//-- Advance until next jump, then solve ODE with new init values ---
// -------------------------- ------------------------------- ------------------        
	  
                Intracellular_ExactEvol(std::max(Time_NextJump_ - currentTime, 0.), Protein_array_, mRNA_array_, S1);
               
         if ((ithGene_ < Number_Of_Genes_)&&(ithGene_ >= 0)) {// i represents a real gene, it's a true jump
                jump = true;
		        TrueJumpCounts_array_[ithGene_] += 1;
		        mRNA_array_[ithGene_] += Alea::exponential_random(1./r);
		    } else {
                    PhantomJumpCounts_ += 1;
                }

                // ----------------- Calculate Kon & Tau ----------
double thetime = log(d0 / D2) / (d0 - D2);
 double thetime1 = log(d0 / D1) / (d0 - D1);
 double thetime4 = log(d0 / D4) / (d0 - D4);
 double thetime7 = log(d0 / D7) / (d0 - D7);
 double thetime8 = log(d0 / D8) / (d0 - D8);
 double thetime9 = log(d0 / D9) / (d0 - D9);
 double thetime5 = log(d0 / D5) / (d0 - D5);
 double thetime6 = log(d0 / D6) / (d0 - D6);
 double thetime3 = log(d0 / D3) / (d0 - D3);
                for (u_int32_t i = 0; i < Number_Of_Genes_; i++) {
		    if ( i == 0){
		   PPmax[i] = Protein_array_[i] + (S1[i] / (d0 - D1)) * mRNA_array_[i]*(exp(-thetime1 * D1) - exp(-thetime1 * d0));
		  	  }//////gene1
		   else if ( i == 3){
		   PPmax[i] = Protein_array_[i] + (S1[i] / (d0 - D4)) * mRNA_array_[i]*(exp(-thetime4 * D4) - exp(-thetime4 * d0));
		  	  }//////gene4
		     else if (i == 6){
		   PPmax[i] = Protein_array_[i] + (S1[i]/ (d0 - D7)) * mRNA_array_[i]*(exp(-thetime7 * D7) - exp(-thetime7 * d0));
		  	  }
             else if (i == 7){
                   PPmax[i] = Protein_array_[i] + (S1[i]/ (d0 - D8)) * mRNA_array_[i]*(exp(-thetime8 * D8) - exp(-thetime8 * d0));         
               }
             else if (i == 8){
                   PPmax[i] = Protein_array_[i] + (S1[i]/ (d0 - D9)) * mRNA_array_[i]*(exp(-thetime9 * D9) - exp(-thetime9 * d0));
               }
            else if (i == 4){
                   PPmax[i] = Protein_array_[i] + (S1[i]/ (d0 - D5)) * mRNA_array_[i]*(exp(-thetime5 * D5) - exp(-thetime5 * d0));
                          }
 else if  (i == 5){
                   PPmax[i] = Protein_array_[i] + (S1[i]/ (d0 - D6)) * mRNA_array_[i]*(exp(-thetime6 * D6) - exp(-thetime6 * d0));
                          }
  else if (i ==2 ){
                   PPmax[i] = Protein_array_[i] + (S1[i]/ (d0 - D3)) * mRNA_array_[i]*(exp(-thetime3 * D3) - exp(-thetime3 * d0));
       }
		   else if (i ==1 ){
		       PPmax[i] = Protein_array_[i] + (S1[i] / (d0 - D2)) * mRNA_array_[i]*(exp(-thetime * D2) - exp(-thetime * d0));
} 
else {
           PPmax[i] = Protein_array_[i] + (S1[i] / (d0 - D1)) * mRNA_array_[i]*(exp(-thetime1 * D1) - exp(-thetime1 * d0));

}
}


        Get_Sigma(Sigma, PPmax, Duration_APC, internal_state_[TCC], 10., false);
        
        for (u_int32_t i = 0; i < Number_Of_Genes_; i++) {

        KKon[i] = (1. - Sigma[i]) * K0 + Sigma[i] * K1 + exp(-10. * log(10.)); // Fix precision errors
		 
                    Tau += KKon[i];
		     
		}
        
        // ----------------- Calculate Delta = Next Jump Time ----------
	
                // Draw the waiting time before the next jump
	        DeltaT = Alea::exponential_random(1./Tau);//*Number_Of_Genes_));
		
        // ---------------------------------------- Select NEXT burst ------------------
                // ----------------- Construct probability array ----------
                double Proba_no_jump = 1.;

        Get_Sigma(Sigma, Protein_array_, Duration_APC, internal_state_[TCC], 10., true);
        
        for (u_int32_t i = 0; i < Number_Of_Genes_; i++) {
    
        KKon[i] = (1. - Sigma[i]) * K0 + Sigma[i]*K1; 
	    ProbaArray[i] = KKon[i] / Tau;
          
        Proba_no_jump = Proba_no_jump - ProbaArray[i];
         }

        
         ProbaArray[Number_Of_Genes_] = Proba_no_jump;
        
                // ----------------- Find gene # from probability array --------
         
		ithGene_ =  Alea::discrete_distribution(Number_Of_Genes_+1, ProbaArray);
	
        delete [] Sigma;
        
         currentTime = Time_NextJump_;

        Time_NextJump_ += DeltaT;
	
    }
	   
// -------------------------- ------------------------------- ------------------
// -------------------------- If next jump is too far away, solve ODE ----------
// -------------------------- ------------------------------- ------------------
	    Time_NextJump_ -= dt;
         
        Intracellular_ExactEvol(t1-currentTime, Protein_array_, mRNA_array_, S1);
        
        delete [] PPmax;
	    delete [] ProbaArray;
        delete [] KKon;    
        delete [] S1;
        break;

        
	}
        case APC:
        {
  
        internal_state_[TCC] = -5.;// to distinguish a APC and lymphocyte
	    internal_state_[APC_Encounters] = -5.; // to distinguish a APC and lymphocyte
        break;}

   }
    // Let y(t) = internal_state_(t)
    // The following loop computes f(y) = dy/dt in order to solve the system until time t + DT, where DT is the "big" time step.
    // DT is cut into several small time steps, with an adaptative small time step that may vary at each small time step:
    // t, t + h, t + h + h', t + h + h' + h'', ... t + DT. Each small timestep is computed by gsl_odeiv2_evolve_apply.
    // The function gsl_odeiv2_evolve_apply advances the system dy/dt = f(y) from time t and state y using
    // the function step with an adaptative method.

    // The last argument to gsl_odeiv2_evolve_apply , here internal_state_, is both an input (state at time t) and an output of the
    // routine (state at time t+h).
    // The h given as the 7th parameter is also both an input and an output: it is the initial guess for the
    // small time step.  The routine will make several calls to func before finding the best value of h.
    // Each call to func advances the system to time t+h. If the resulting error is too large, this candidate new state
    // is rejected, the system goes back to time t and a new, smaller, h is tested. Once the error is smaller
    // than the tolerance specified in argument c, the new state is accepted and the system advances to time t+h.

    // Note that the error is controlled at each timestep. Thus the "right" h value may change at each
    // time step: with the above notations, the successive "right" values would be h, then h', then h'', and so on.
    // With the RK-Cash-Karp method, there are roughly 6 calls to func to compute the candidate new state and its error.
    // If the candidate new state is rejected, in the worst case, 6 more calls to func will be necessary to find the
    // new candidate state corresponding to the smaller h. There are thus 6 * (number of candidate values for h) calls
    // to func at each small timestep.

    // When the right h is found, the candidate new state is accepted.
    // For the next timestep, a slightly larger h may be tried first, if the next error is likely to be smaller
    // than the tolerance.
    
    /*  while (t < t1) {
       int status = gsl_odeiv2_evolve_apply(e, c, s,
              &sys,
               &t, t1, // t will be updated several times until t1=t+DT is reached //
    &h, // the small time step can be updated several times too (adaptative time step) //
      internal_state_); // internal_state_ will be updated several times //
    

     if (status != GSL_SUCCESS) {

       fprintf(stderr, "Warning: error in ODE update\n");
       exit(EXIT_FAILURE);
     }
    */

     //   printf ("%f %f\n", t, h);
   //  }

      /* gsl_odeiv2_evolve_free(e);
 gsl_odeiv2_control_free(c);
 gsl_odeiv2_step_free(s);*/


}
/**
 * Piping method: wraps the non-static version of compute_dydt to be used in GSL
 */
/*int Lymphocyte::compute_dydt(double t,
        const double y[],
        double dydt[],
        void* cell) {
    return ((Lymphocyte*) cell)->compute_dydt(t, y, dydt);
    }*/

/**
 * Compute dy/dt given a time point and the state and parameters of the system
 *
 * \param t the initial time
 * \param y an estimate of the system state at time t
 * \param dydt computed dy/dt (output)
 * \return GSL_SUCCESS on success
 *
 * For a given small timestep (suppose we are at t0 and we want
 * to find a candidate new state for time t0+h), func will be called several times by gsl_odeiv2_evolve_apply,
 * with different time points, like, say,  t0, t0 + h/5, t0 + 3h/10 (actually 6 intermediate timepoints for the RK-Cash-Karp method,
 * for example). Then y(t+h) will be computed twice, with two different accuracies (4th and 5th order), with two different
 * linear combinations of the 6 values of f (a idea of Fehlberg's).
 * From the comparison between two orders, the local error can be estimated.
 * If this error is too large, then a smaller h must be used. If on the contrary, the error is smaller than the tolerance,
 * the candidate new state be will accepted, the system has really advanced to time t0+h.

 * In our case, y is internal_state_. func only uses it as an input to estimate the derivative (cf const keyword). func does
 * NOT update internal_state_. It is gsl_odeiv2_evolve_apply that updates internal_state_, with the right small timestep h,
 * by using the 5th order approximation.
 */
/*int Lymphocyte::compute_dydt(double t, const double y[], double dydt[]) const {
    // NOTA: t is unused (for now). It is left here because it might be used in
    // the future and it helps understand what's going on.
    assert(t == t); // Suppress warning

    // this is the intracellular network

    // INPUT FROM other cells
    // double F = getInSignal(InterCellSignal::CLOCK);
    // double death_signal = getInSignal(InterCellSignal::DEATH);
    double in_signal = getInSignal(InterCellSignal::CYCLE_Z);
    double death_signal = getInSignal(InterCellSignal::DEATH);

    // CELL CYCLE VARIABLES (ref: Introducing via reduced cell cycle model, Battogtokh & Tyson 2006 PRE)
    double cy_x = y[CY_X]; // y is a fictive intermediate internal state where the stepper wants to know the derivative
    double cy_z = y[CY_Z]; // we copy it into local variables to make the equations shorter and clearer

    // DEATH AND SURVIVAL VARIABLES
    double cd_death = y[CD_DEATH];

    // CELL CYCLE AUXILIARY FUNCTION
    double e_temp = cy_p - cy_x + cy_p * cy_j + cy_x*cy_j;
    double W0 = 2 * cy_x * cy_j / (e_temp + sqrt(e_temp * e_temp - 4 * cy_x * cy_j * (cy_p - cy_x)));
    double a_temp = cy_k6 + cy_k7*cy_z;
    double b_temp = cy_k8 *  (outer_volume() / size_.volume_min()) + cy_k9*cy_x;
    e_temp = b_temp - a_temp + b_temp * cy_j + a_temp*cy_j;
    double Y0 = 2 * a_temp * cy_j / (e_temp + sqrt(e_temp * e_temp - 4 * a_temp * cy_j * (b_temp - a_temp)));

    // CELL CYCLE EQUATIONS
    // cycling signal variables cy
    // dX/dt = m (k1+k2*W0)-(k3+k4*Y0+k5*Z)*X
    // dZ/dt = k10(1+A(1+sin(f*t)))+k11*X-k12*Z
    // ORIGINAL EQ  f[0] = cell_volume*(cy_k1+cy_k2*W0) - (cy_k3+cy_k4*Y0+cy_k5*cy_z)*cy_x;  (outer_volume() / size_.volume_min()
     dydt[CY_X] = ( (outer_volume() / size_.volume_min())* (cy_k1 + cy_k2 * W0) -
            (cy_k3 + cy_k4 * Y0 + cy_k5 * cy_z) * cy_x);
     dydt[CY_Z] = (cy_k10 * (1 + cy_A) + cy_k11 * cy_x - cy_k12 * cy_z);
    //dydt[CY_Z] = (cy_k10 * (1 + cy_A)  - cy_k12 * cy_z);
 
    // death signal variable cd
    dydt[CD_DEATH] = death_signal - k_survival * cd_death;
     dydt[TCC] = 0.0 * dydt[CY_Z];
     dydt[APC_Encounters] = 0.0 * dydt[CY_Z];
     dydt[DifferentiationState] = 0.0 * dydt[CY_Z];
    return GSL_SUCCESS;

    

}
*/


// Register this class in Cell
bool Lymphocyte::registered_ =
        Cell::RegisterClass(classId_, classKW_,
        [](const CellType type,
        const MoveBehaviour& move_behaviour,
        const Coordinates<double>& pos,
        double initial_volume,
        double volume_min,
        double doubling_time) {
            return static_cast<Cell*> (
                    new Lymphocyte(type,
                    move_behaviour,
                    pos,
                    initial_volume,
                    volume_min,
                    doubling_time));
        },
[](gzFile backup_file) {
    return static_cast<Cell*> (new Lymphocyte(backup_file));
}
);





