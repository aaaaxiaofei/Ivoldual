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
      ("At most one of -subsample, -supersample or -subdivide may be used.");
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

  void IVOLDUAL::IVOLDUAL_SCALAR_GRID::Subdivide
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid2, const int subdivide_period)
  {  
    const DTYPE dimension = scalar_grid2.Dimension();
    IJK::ARRAY<ATYPE> subdivide_axis_size(dimension);
    IJK::PROCEDURE_ERROR error("IVOLDUAL_SCALAR_GRID::Subdivide");

    for (DTYPE d = 0; d < dimension; d++) {
      subdivide_axis_size[d] =
        IJK::compute_supersample_size(scalar_grid2.AxisSize(d), subdivide_period);
    }

    this->SetSize(dimension, subdivide_axis_size.PtrConst());

    if (this->NumVertices() < 1) { return; };

    SupersampleCopy(scalar_grid2, subdivide_period);
    SubdivideInterpolate(subdivide_period);
  
  }

  void IVOLDUAL::IVOLDUAL_SCALAR_GRID::SubdivideInterpolate
  (const int subdivide_period)
  {
    const DTYPE dimension = this->Dimension();
    IJK::ARRAY<ATYPE> subgrid_axis_size(dimension);
    IJK::ARRAY<ATYPE> subsample_period(dimension);
    IJK::ARRAY<VTYPE> axis_increment(dimension);
    compute_increment(*this, axis_increment.Ptr());

    for (DTYPE d = 0; d < this->Dimension(); d++) {
      for (DTYPE j = 0; j < d; j++) { subsample_period[j] = 1; };
      for (DTYPE j = d; j < this->Dimension(); j++)
        { subsample_period[j] = subdivide_period; };
      for (DTYPE j = 0; j < this->Dimension(); j++)
        { subgrid_axis_size[j] = this->AxisSize(j); };
      subgrid_axis_size[d] = 1;

      NTYPE numv;
      IJK::compute_subsample_size
        (this->Dimension(), subgrid_axis_size.PtrConst(),
         subsample_period.PtrConst(), numv);

      IJK::ARRAY<VTYPE> vlist(numv);
      IJK::subsample_subgrid_vertices
        (*this, 0, subgrid_axis_size.PtrConst(),
         subsample_period.PtrConst(), vlist.Ptr());

      for (VTYPE x = 0; x+1 < this->AxisSize(d); x += subdivide_period) {

        VTYPE inc0 = x*axis_increment[d];
        VTYPE inc1 = inc0 + subdivide_period*axis_increment[d];
        for (VTYPE i = 0; i < numv; i++) {
          VTYPE v0 = vlist[i] + inc0;
          VTYPE v1 = vlist[i] + inc1;

          VTYPE v2 = v0;
          for (VTYPE j = 1; j < subdivide_period; j++) {
            v2 += axis_increment[d];
            STYPE s0 = this->scalar[v0];
            STYPE s1 = this->scalar[v1];
            this->scalar[v2] = std::max(s0, s1);
          }
        }
      }
    }
  }