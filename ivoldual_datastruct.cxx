/// \file ivoldual_datastruct.cxx
/// ivoldual data structures.

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2017 Rephael Wenger

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


#include "ivoldual_datastruct.h"

// **************************************************
// CLASS IVOLDUAL_DATA_FLAGS
// **************************************************

// Initialize
void IVOLDUAL::IVOLDUAL_DATA_FLAGS::Init()
{
  hex_tri_method = UNDEFINED_HEX_TRI;
  flag_split_ambig_pairs = false;
  flag_split_ambig_pairsB = false;
  flag_add_isov_dual_to_hexahedra = false;
  flag_orient_in = false;
  default_interior_code = 2;
}


// **************************************************
// CLASS IVOLDUAL_DATA MEMBER FUNCTIONS
// **************************************************

void IVOLDUAL::IVOLDUAL_DATA::SubdivideScalarGrid
(const DUALISO_SCALAR_GRID_BASE & scalar_grid2) 
{
  const int subdivide_resolution = 2;

  const int dimension = scalar_grid2.Dimension();
  IJK::ARRAY<COORD_TYPE> spacing(dimension);
  scalar_grid.Subdivide(scalar_grid2, subdivide_resolution);
  IJK::copy_coord(dimension, scalar_grid2.SpacingPtrConst(),
                  spacing.Ptr());
  IJK::divide_coord
    (dimension, subdivide_resolution, spacing.PtrConst(), spacing.Ptr());
  scalar_grid.SetSpacing(spacing.PtrConst()); 

  is_scalar_grid_set = true;
}


// Copy, subsample, supersample or subdivide scalar grid.
void IVOLDUAL::IVOLDUAL_DATA::SetScalarGrid
(const DUALISO_SCALAR_GRID_BASE & scalar_grid2, 
 const bool flag_subsample, const int subsample_resolution,
 const bool flag_supersample, const int supersample_resolution, 
 const bool flag_subdivide)
{
  IJK::PROCEDURE_ERROR error("IVOLDUAL_DATA::SetScalarGrid");

  if (flag_subsample + flag_supersample + flag_subdivide > 1) {
    error.AddMessage
      ("Scalar grid can only be one of subsampled or supersampled or subdivide.");
    throw error;
  }
  
  if (flag_subsample) {
    // subsample grid
    SubsampleScalarGrid(scalar_grid2, subsample_resolution);
  }
  else if (flag_supersample) {
    // supersample grid
    SupersampleScalarGrid(scalar_grid2, supersample_resolution);
  }
  else if (flag_subdivide) {
    // subdivide grid
    SubdivideScalarGrid(scalar_grid2);
  }
  else {
    CopyScalarGrid(scalar_grid2);
  };
}
