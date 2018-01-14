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
#include "ijkgrid_macros.h"
#include "ivoldual_ivolpoly.txx"
#include "ivoldualtable.h"

// **************************************************
// CLASS IVOLDUAL_DATA_FLAGS
// **************************************************

// Initialize
void IVOLDUAL::IVOLDUAL_DATA_FLAGS::Init()
{
  hex_tri_method = UNDEFINED_HEX_TRI;
  flag_split_ambig_pairs = false;
  flag_split_ambig_pairsB = false;
  flag_rm_non_manifold = false;
  flag_add_isov_dual_to_hexahedra = false;
  flag_orient_in = false;
  flag_lsmooth_elength = false;
  flag_lsmooth_jacobian = false;
  flag_split_hex = false;
  default_interior_code = 2;
  lsmooth_elength_iter = 1;
  lsmooth_jacobian_iter = 1;
  lsmooth_elength_threshhold = 0.1;
  lsmooth_jacobian_threshhold = 0.0;
}

// **************************************************
// CLASS IVOLDUAL_POLY_INFO
// **************************************************

// Initialize
void IVOLDUAL::IVOLDUAL_POLY_INFO::Init()
{
  flag_subdivide_hex = false;
}

// **************************************************
// CLASS DUALISO INFO MEMBER FUNCTIONS
// **************************************************
IVOLDUAL::IVOLDUAL_INFO::IVOLDUAL_INFO() 
{
  Clear();
}

IVOLDUAL::IVOLDUAL_INFO::IVOLDUAL_INFO(const int dimension) 
{
  Clear();
}

void IVOLDUAL::IVOLDUAL_INFO::Clear() 
{
  non_manifold_changes = 0;
}

// **************************************************
// CLASS IVOLDUAL_DATA MEMBER FUNCTIONS
// **************************************************

void IVOLDUAL::IVOLDUAL_DATA::SubdivideScalarGrid
(const DUALISO_SCALAR_GRID_BASE & scalar_grid2, IVOLDUAL_SCALAR_GRID & scalar_grid3, 
 const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1) 
{
  const int subdivide_resolution = 2;

  const int dimension = scalar_grid2.Dimension();
  IJK::ARRAY<COORD_TYPE> spacing(dimension);
  scalar_grid3.Subdivide(scalar_grid2, subdivide_resolution, isovalue0, isovalue1);

  IJK::copy_coord(dimension, scalar_grid2.SpacingPtrConst(),
                  spacing.Ptr());
  IJK::divide_coord
    (dimension, subdivide_resolution, spacing.PtrConst(), spacing.Ptr());
  scalar_grid3.SetSpacing(spacing.PtrConst()); 

  is_scalar_grid_set = true;
}


// Copy, subsample, supersample or subdivide scalar grid.
void IVOLDUAL::IVOLDUAL_DATA::SetScalarGrid
(const DUALISO_SCALAR_GRID_BASE & scalar_grid2, 
 const bool flag_subsample, const int subsample_resolution,
 const bool flag_supersample, const int supersample_resolution, 
 const bool flag_subdivide, const SCALAR_TYPE isovalue0, 
 const SCALAR_TYPE isovalue1)
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
    IVOLDUAL_SCALAR_GRID scalar_grid3;
    SubdivideScalarGrid(scalar_grid2, scalar_grid3, isovalue0, isovalue1);
  }
  else {
    CopyScalarGrid(scalar_grid2);
  };
}

void IVOLDUAL::IVOLDUAL_DATA::EliminateAmbigFacets
(const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1, 
 VERTEX_INDEX & changes_of_ambiguity)
{
  int last_changes = 1;

  while (last_changes) {
    // Eliminate non-manifold cases caused by opposite diagonal plus vertices.
    last_changes = EliminateNonmanifold(isovalue0, isovalue1);

    // Eliminate non-manifold cases caused by ambiguous facets.
    last_changes += scalar_grid.EliminateAmbiguity(isovalue0, isovalue1);

    changes_of_ambiguity += last_changes;
  }

  scalar_grid.AddOuterLayer();
}

int IVOLDUAL::IVOLDUAL_DATA::EliminateNonmanifold
(const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1)
{
  int changes_of_cube = 0;
  const int DIM3(3);
  IVOLDUAL_CUBE_TABLE ivoldual_table(DIM3, true);
  const int interior_code(2);
  const int above_code(3);
  TABLE_INDEX table_index;

  IJK_FOR_EACH_GRID_CUBE(icube, scalar_grid, VERTEX_INDEX) {

    IVOLDUAL::compute_ivoltable_index_of_grid_cube<4>
      (scalar_grid, isovalue0, isovalue1, icube,
       interior_code, above_code, table_index);
    
    if (ivoldual_table.NumAmbiguousFacets(table_index) == 0 && 
        ivoldual_table.IsNonManifold(table_index)) {
      // Eliminate non-manifold caused by diagonal plus vertices
      scalar_grid.EliminateDiagonalNonmanifold(isovalue0, isovalue1, icube);
      changes_of_cube++;
    }
  }
  return changes_of_cube;
}

void IVOLDUAL::IVOLDUAL_SCALAR_GRID::Subdivide
(const DUALISO_SCALAR_GRID_BASE & scalar_grid2, const int subdivide_period, 
 const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1)
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
  SubdivideInterpolate(subdivide_period, isovalue0, isovalue1);

}

void IVOLDUAL::IVOLDUAL_SCALAR_GRID::SubdivideInterpolate
(const int subdivide_period, const SCALAR_TYPE isovalue0, 
 const SCALAR_TYPE isovalue1)
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

          // Linear interpolation.
          this->scalar[v2] = 0.5*(s0+s1);
        }
      }
    }
  }
}

/// Eliminate ambiguity facets in supersample grid.
int IVOLDUAL::IVOLDUAL_SCALAR_GRID::EliminateAmbiguity
(const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1) 
{
  const DTYPE dim = this->dimension;
  int changes_of_ambiguity = 0;

  ATYPE sizex = this->AxisSize(0);
  ATYPE sizey = this->AxisSize(1);
  ATYPE sizez = this->AxisSize(2);


  for (VTYPE i = 0; i < this->NumVertices(); i++) {
    

    for (DTYPE d = 0; d < dim; d++) {

      VTYPE v01_idx = this->NextVertex(i, d);
      VTYPE v10_idx = this->NextVertex(i, (d+1)%dim);
      VTYPE v11_idx = this->NextVertex(this->NextVertex(i, d), (d+1)%dim);

      // If current vertex is at the end of a direction, then continue.
      if ((d == 0 && v01_idx%sizex == 0) ||
          (d == 2 && v10_idx%sizex == 0) ||
          (d == 1 && (v01_idx/sizex)%sizey == 0) ||
          (d == 0 && (v10_idx/sizex)%sizey == 0) ||
          (d == 2 && (v01_idx/sizex/sizey)%sizez == 0) ||
          (d == 1 && (v10_idx/sizex/sizey)%sizez == 0) ) 
        { continue; } 

      // Get scalar values of a square.
      STYPE *s00 = &this->scalar[i];
      STYPE *s01 = &this->scalar[v01_idx];
      STYPE *s10 = &this->scalar[v10_idx];
      STYPE *s11 = &this->scalar[v11_idx];


      // In a unit cube before subdivide:
      // Original vertex is level 0
      // Vertex at unit cube edge center is level 1
      // Vertex at unit cube facet center is level 2
      // Vertex at unit cube body center is level 3
      if (*s00 < isovalue0 && *s01 > isovalue0 && *s10 > isovalue0 && *s11 < isovalue0) {
        
        VTYPE z_level = v11_idx/sizex/sizey;
        VTYPE y_level = (v11_idx%(sizex*sizey))/sizex;
        VTYPE x_level = v11_idx%sizex;
        VTYPE level = z_level%2 + y_level%2 + x_level%2;

        if (level == 0) *s00 = 0.5 * (isovalue0 + isovalue1); 
        else *s11 = 0.5 * (isovalue0 + isovalue1); 
        changes_of_ambiguity++; 
      }
      else if (*s00 > isovalue0 && *s01 < isovalue0 && *s10 < isovalue0 && *s11 > isovalue0) {
        
        VTYPE z_level = v01_idx/sizex/sizey;
        VTYPE y_level = (v01_idx%(sizex*sizey))/sizex;
        VTYPE x_level = v01_idx%sizex;
        VTYPE level = z_level%2 + y_level%2 + x_level%2;

        if (level == 0) *s10 = 0.5 * (isovalue0 + isovalue1); 
        else *s01 = 0.5 * (isovalue0 + isovalue1); 
        changes_of_ambiguity++; 
      }
      else if (*s00 > isovalue1 && *s01 < isovalue1 && *s10 < isovalue1 && *s11 > isovalue1) {
        
        VTYPE z_level = v11_idx/sizex/sizey;
        VTYPE y_level = (v11_idx%(sizex*sizey))/sizex;
        VTYPE x_level = v11_idx%sizex;
        VTYPE level = z_level%2 + y_level%2 + x_level%2;

        if (level == 0) *s00 = 0.5 * (isovalue0 + isovalue1); 
        else *s11 = 0.5 * (isovalue0 + isovalue1); 
        changes_of_ambiguity++; 
      }
      else if (*s00 < isovalue1 && *s01 > isovalue1 && *s10 > isovalue1 && *s11 < isovalue1) {
        
        VTYPE z_level = v01_idx/sizex/sizey;
        VTYPE y_level = (v01_idx%(sizex*sizey))/sizex;
        VTYPE x_level = v01_idx%sizex;
        VTYPE level = z_level%2 + y_level%2 + x_level%2;

        if (level == 0) *s10 = 0.5 * (isovalue0 + isovalue1); 
        else *s01 = 0.5 * (isovalue0 + isovalue1); 
        changes_of_ambiguity++; 
      }
    }  
  }

  return changes_of_ambiguity;
}

/// Add outer layer to subdivide grid
void IVOLDUAL::IVOLDUAL_SCALAR_GRID::AddOuterLayer() 
{

  const DTYPE dim = this->Dimension();
  
  int dim_x = this->AxisSize(0);
  int dim_y = this->AxisSize(1);
  int dim_z = this->AxisSize(2);

  AXIS_SIZE_TYPE axis_size[3] = {dim_x, dim_y, dim_z};

  IVOLDUAL_SCALAR_GRID scalar_grid3;
  scalar_grid3.SetSize(3, axis_size);
  scalar_grid3.CopyScalar(*this);

  ATYPE dim_x2 = dim_x + 2;
  ATYPE dim_y2 = dim_y + 2;
  ATYPE dim_z2 = dim_z + 2;

  AXIS_SIZE_TYPE axis_size2[3] = {dim_x2, dim_y2, dim_z2};
  this->SetSize(dim, axis_size2);

  VERTEX_INDEX iv = 0;
  VERTEX_INDEX iw = 0;

  for(int k = 0; k < dim_z; k++) {

    for(int j = 0; j < dim_y; j++) {

      for(int i = 0; i < dim_x; i++) {
        this->Set(iw++, scalar_grid3.Scalar(iv));
        if (i == 0 || i+1 == dim_x)
          this->Set(iw++, scalar_grid3.Scalar(iv));
        iv++;
      }
      if (j == 0 || j+1 == dim_y) {
        for (int i = 0; i < dim_x2; i++) {
          this->Set(iw, this->Scalar(iw-dim_x2));
          iw++;
        }
      }

    }

    if (k == 0 || k+1 == dim_z) {
      for (int j = 0; j < dim_y2; j++) {
        for (int i = 0; i < dim_x2; i++) {
          this->Set(iw, this->Scalar(iw-dim_x2*dim_y2));
          iw++;
        }
      }
    }
  }
}

void IVOLDUAL::IVOLDUAL_SCALAR_GRID::EliminateDiagonalNonmanifold
(const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
 VTYPE icube) 
{

  for (int k = 0; k < this->NumCubeVertices(); k++) {
    int idx = this->CubeVertex(icube, k);
    
    if (this->Scalar(idx) > isovalue1) {
      this->Set(idx, 0.5*(isovalue0 + isovalue1));
      break;
    }
  }
}


// **************************************************
// DUAL_IVOLVERT MEMBER FUNCTIONS
// **************************************************

void IVOLDUAL::DUAL_IVOLVERT::Init()
{
  num_incident_hex = 0;
  num_incident_iso_quad = 0;
  separation_vertex = 0;
  separation_edge_direction = 0;
  flag_missing_ivol_hexahedra = false;
  is_doubly_connected = false;
  in_loop = false;
  in_box = false;
  in_pseudobox = false;
}
