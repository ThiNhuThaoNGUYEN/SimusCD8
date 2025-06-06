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

#ifndef SIMULATION_H__
#define SIMULATION_H__

// =================================================================
//                              Includes
// =================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>

#include <string>

//using namespace std;

#include "params/SimulationParams.h"
#include "Grid.h"
#include "Population.h"
#include "params.h"
#include "Cell.h"
#include "OutputManager.h"

// =================================================================
//                          Class declarations
// =================================================================


using std::string;


/*!
  \brief TODO.
*/
class Simulation
{
 public :
  // =================================================================
  //                        Singleton management
  // =================================================================
 private:
  static Simulation instance_;
  Simulation();
  Simulation (const Simulation&) = delete;
  Simulation& operator= (const Simulation&) = delete;
  virtual ~Simulation();
 public:
  static Simulation& instance() {return instance_;};


  // =================================================================
  //                            Public Methods
  // =================================================================
  static void Setup(const SimulationParams& simParams,
                    const string& output_dir);

 static string output_dir_ ; // TODO add trailing / if missing
 static string input_dir_ ; // TODO add trailing / if missing

  static void Run();

  static void Save();
  static void Load(const string& input_dir,
                   int32_t timestep,
                   const string& output_dir);


  // =================================================================
  //                            Accessors
  // =================================================================
  static double sim_time() {return instance_.time_;};
  static double dt() {return instance_.dt_;};
  static int32_t timestep() {return instance_.timestep_;};
  static const OutputManager& output_manager() {
    return instance_.output_manager_;};
  // static double maxtime() {return instance_.maxtime_;};
  static const Population& pop() {return *(instance_.pop_);};
  static const std::list<InterCellSignal>& using_signals() { return instance_.using_signals_; };
  static const std::vector<std::string>& intercellular_signal_names() { return instance_.intercellular_signal_names_; };
  static const std::vector<std::string>& celltype_names() { return instance_.celltype_names_; };
  static const std::vector<double> max_signals() { return instance_.max_signals_; }
  static const std::vector<double> min_signals() { return instance_.min_signals_; }



 protected :
  // =================================================================
  //                           Protected Methods
  // =================================================================
  void DoSetup(const SimulationParams& simParams, const string& output_dir);
  void DoRun();
  void Finalize();
  void Update();
  void ComputeNeighbourhood();
  void CheckNeighbourhood();
  void ComputeInteractions();
  void ApplyUpdate();
  void UpdateMinMaxSignals(Cell* cell);

  void DoSave();
  void DoLoad(const string& input_dir,
              double time,
              const string& output_dir);


// void InitOutputDir(const string& output_dir);
  // =================================================================
  //                              Attributes
  // =================================================================
  /** Current time */
  double  time_;
  /** Current time-step */
  int32_t timestep_;
  /** Number of timesteps to be simulated */
  int32_t max_timestep_;
  /** Delta time between 2 time steps */
  double dt_;
  /** Number of timesteps between 2 backups */
  int32_t backup_dtimestep_;

  /** The cell population */
  Population* pop_;

  /** Intercellular signals to monitor */
  std::list<InterCellSignal> using_signals_;
  std::vector<std::string> intercellular_signal_names_;
  std::vector<std::string> celltype_names_;

  std::vector<double> max_signals_;
  std::vector<double> min_signals_;

  /** The spatial grid */
public:
  Grid* grid_;

  OutputManager output_manager_;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================

#endif // SIMULATION_H__
