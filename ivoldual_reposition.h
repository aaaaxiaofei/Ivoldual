/// \file ivoldual_reposition.h
/// Reposition vertices for mesh quality optimization.

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

#ifndef _IVOLDUAL_REPOSITION_
#define _IVOLDUAL_REPOSITION_

#include "ijk.txx"

#include "ivoldual_types.h"
#include "ivoldual_datastruct.h"
#include "ivoldualtable.h"

/// ivoldual classes and routines.
namespace IVOLDUAL {

  // **************************************************
  // MESH OPTIMIZATION
  // **************************************************
  
  /// Eliminate non-manifold facets and cubes
  void eliminate_non_manifold_grid
  (IVOLDUAL_DATA & ivoldual_data, 
   const SCALAR_TYPE isovalue0, 
   const SCALAR_TYPE isovalue1, 
   IVOLDUAL_INFO & dualiso_info);

  /// Laplacian Smoothing for small edge length
  void laplacian_smooth_elength
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const IVOLDUAL_DATA_FLAGS & param,
   IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   const DUAL_IVOLVERT_ARRAY & ivolv_list,
   COORD_ARRAY & vertex_coord,
   float laplacian_smooth_limit, 
   int iteration);

  /// Laplacian Smoothing for negative Jacobian
  void laplacian_smooth_jacobian
  (const std::vector<VERTEX_INDEX> & ivolpoly_cube,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   const DUAL_IVOLVERT_ARRAY & ivolv_list,
   COORD_ARRAY & vertex_coord, 
   float jacobian_limit, 
   int iteration);

  void laplacian_smooth_jacobian
  (const std::vector<VERTEX_INDEX> & ivolpoly_cube,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   const DUAL_IVOLVERT_ARRAY & ivolv_list,
   COORD_ARRAY & vertex_coord, 
   std::vector<int> negative_jabocian_list);

  void gradient_move_vertex
	(const std::vector<VERTEX_INDEX> & ivolpoly_cube,
	 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
	 IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   IJK::VERTEX_POLY_INCIDENCE<int,int> vertex_poly_incidence,
	 const DUAL_IVOLVERT_ARRAY & ivolv_list,
	 COORD_ARRAY & vertex_coord, 
	 COORD_TYPE *ver_coord, int ver_index,	 
	 bool flag_onLower,  bool flag_onUpper);

}

#endif