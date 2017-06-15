/// \file ivoldual_main.cxx
/// generate interval volume using dual contouring algorithm
/// Version 0.2.1

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011-2017 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <iostream>

#include "isodual.h"
#include "ivoldualIO.h"
#include "ivoldual.h"

#include "ijkdual_triangulate.txx"

using namespace IJK;
using namespace IJKDUAL;
using namespace ISODUAL;
using namespace IVOLDUAL;

using namespace std;

// local subroutines
void memory_exhaustion();
void construct_interval_volume
(const IO_INFO & io_info, const DUALISO_DATA & dualiso_data,
 DUALISO_TIME & dualiso_time, IO_TIME & io_time);



// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
  time_t start_time;
  time(&start_time);

  DUALISO_TIME dualiso_time;
  IO_TIME io_time = {0.0, 0.0};
  IO_INFO io_info;
  IJK::ERROR error;

  try {

    std::set_new_handler(memory_exhaustion);

    // Set to interval volume.
    io_info.flag_interval_volume = true;

    parse_command_line(argc, argv, io_info);

    DUALISO_SCALAR_GRID full_scalar_grid, scalar_grid_4D;
    NRRD_HEADER nrrd_header;
    read_nrrd_file
      (io_info.input_filename, full_scalar_grid,  nrrd_header, io_time);

    if (!check_input(io_info, full_scalar_grid, error)) 
      { throw(error); };

    nrrd_header.GetSpacing(io_info.grid_spacing);

    // set DUAL datastructures and flags
    DUALISO_DATA dualiso_data;

    // Note: dualiso_data.SetScalarGrid must be called before set_mesh_data.
    // subsample and supersample parameters are hard-coded here.
    dualiso_data.SetScalarGrid
      (full_scalar_grid, io_info.flag_subsample, io_info.subsample_resolution, 
       io_info.flag_supersample, io_info.supersample_resolution);
    dualiso_data.Set(io_info);
    warn_non_manifold(io_info);
    report_num_cubes(full_scalar_grid, io_info, dualiso_data);

    construct_interval_volume
      (io_info, dualiso_data, dualiso_time, io_time);

    if (io_info.flag_report_time) {

      time_t end_time;
      time(&end_time);
      double total_elapsed_time = difftime(end_time, start_time);

      cout << endl;
      report_time(io_info, io_time, dualiso_time, total_elapsed_time);
    };

  } 
  catch (ERROR & error) {
    if (error.NumMessages() == 0) {
      cerr << "Unknown error." << endl;
    }
    else { error.Print(cerr); }
    cerr << "Exiting." << endl;
    exit(20);
  }
  catch (...) {
    cerr << "Unknown error." << endl;
    exit(50);
  };

}

// forward declaration
template <typename DUALISO_DATA_TYPE, typename DUAL_ISOSURFACE_TYPE>
void rescale_and_triangulate
(const IO_INFO & io_info, const DUALISO_DATA_TYPE & dualiso_data,
 const SCALAR_TYPE isovalue, DUAL_ISOSURFACE_TYPE & dual_isosurface);

void construct_interval_volume
(const IO_INFO & io_info, const DUALISO_DATA & dualiso_data,
 DUALISO_TIME & dualiso_time, IO_TIME & io_time)
{
  int dimension = dualiso_data.ScalarGrid().Dimension();
  const int num_cube_vertices = IJK::compute_num_cube_vertices(dimension);
  const int num_cubes = dualiso_data.ScalarGrid().ComputeNumCubes();

  io_time.write_time = 0;

  for (unsigned int i = 0; i+1 < io_info.isovalue.size(); i++) {

    const SCALAR_TYPE isovalue0 = io_info.isovalue[i];
    const SCALAR_TYPE isovalue1 = io_info.isovalue[i+1];

    DUALISO_INFO dualiso_info(dimension);
    dualiso_info.grid.num_cubes = num_cubes;

    // Dual contouring.  
    DUAL_INTERVAL_VOLUME dual_interval_volume(dimension, num_cube_vertices);

    dual_contouring_interval_volume
      (dualiso_data, isovalue0, isovalue1, dual_interval_volume,
       dualiso_info);

    // Time info
    dualiso_time.Add(dualiso_info.time);

    // Rescale 
    rescale_vertex_coord
      (dimension, dualiso_data.ScalarGrid().SpacingPtrConst(),
       dual_interval_volume.vertex_coord);

    OUTPUT_INFO output_info;
    set_output_info(io_info, i, output_info);
    output_info.SetDimension(dimension, num_cube_vertices);

    output_dual_isosurface
      (output_info, dualiso_data, dual_interval_volume, dualiso_info, io_time);
  }
}


void memory_exhaustion()
{
  cerr << "Error: Out of memory.  Terminating program." << endl;
  exit(10);
}

