/// \file ivoldual.h
/// Construct interval volume in arbitrary dimensions using dual contouring.

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

/*!
  \mainpage IVOLDUAL: EXTRACT INTERVAL VOLUME

  IVOLDUAL is a program for generating interval volumes
  using the dual contouring.
*/

#ifndef _IVOLDUAL_
#define _IVOLDUAL_

#include <string>

#include "ijk.txx"

#include "ivoldual_types.h"
#include "ivoldual_datastruct.h"

/// ivoldual classes and routines.
namespace IVOLDUAL {

// **************************************************
// DUAL CONTOURING INTERVAL VOLUME
// **************************************************

  /// Construct interval volume using dual contouring.
  void dual_contouring_interval_volume
    (const DUALISO_DATA & dualiso_data, 
     const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
     DUAL_INTERVAL_VOLUME & dual_interval_volume, DUALISO_INFO & dualiso_info);

  /// Construct interval volume using dual contouring.
  /// - Returns list of isosurface polytope vertices
  ///   and list of isosurface vertex coordinates
  void dual_contouring_interval_volume
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const DUALISO_DATA_FLAGS & param,
   std::vector<ISO_VERTEX_INDEX> & ivolpoly_vert,
   IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   COORD_ARRAY & vertex_coord,
   MERGE_DATA & merge_data, 
   DUALISO_INFO & dualiso_info);

  /// Construct interval volume using dual contouring.
  /// - Returns list of isosurface polytope vertices
  /// - Version which creates isodual_table.
  void dual_contouring_interval_volume
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
   const DUALISO_DATA_FLAGS & param,
   std::vector<ISO_VERTEX_INDEX> & ivolpoly_vert,
   IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   COORD_ARRAY & vertex_coord,
   MERGE_DATA & merge_data, 
   DUALISO_INFO & dualiso_info);

  /// Construct interval volume using dual contouring.
  /// Returns list of isosurface polytope vertices
  ///   and list of isosurface vertex coordinates
  void dual_contouring_interval_volume
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const DUALISO_DATA_FLAGS & param,
   std::vector<ISO_VERTEX_INDEX> & ivolpoly_vert,
   std::vector<GRID_CUBE_DATA> & cube_ivolv_list,
   std::vector<DUAL_ISOVERT> & ivolv_list,
   IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   COORD_ARRAY & vertex_coord,
   MERGE_DATA & merge_data, 
   DUALISO_INFO & dualiso_info);


  // **************************************************
  // ENCODE GRID VERTICES
  // **************************************************

  /// Encode grid vertices as 0, 1, 2 or 3.
  /// - 0: Below lower isovalue.
  /// - 1: Between lower and upper isovalue.
  /// - 2: Between lower and upper isovalue.
  /// - 3: Above upper isovalue.
  /// @param default_interior_code Default code (1 or 2) for vertices
  ///   between lower and upper isovalue.
  void encode_grid_vertices
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
   const GRID_VERTEX_ENCODING default_interior_code,
   IVOLDUAL_ENCODED_GRID & encoded_grid,
   DUALISO_INFO & dualiso_info);


  // **************************************************
  // COMPUTE IVOLTABLE INFO FOR EACH ACTIVE GRID CUBE
  // **************************************************

  void compute_cube_ivoltable_info
  (const IVOLDUAL_ENCODED_GRID & encoded_grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   std::vector<GRID_CUBE_DATA> & cube_ivolv_list);


  // **************************************************
  // EXTRACT DUAL INTERVAL VOLUME POLYTOPES
  // **************************************************

  void extract_dual_ivolpoly
  (const IVOLDUAL_ENCODED_GRID & encoded_grid,
   std::vector<ISO_VERTEX_INDEX> & ivolpoly,
   std::vector<POLY_VERTEX_INDEX> & poly_vertex,
   IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   DUALISO_INFO & dualiso_info);

  void extract_ivolpoly_dual_to_grid_edges
  (const IVOLDUAL_ENCODED_GRID & encoded_grid,
   std::vector<ISO_VERTEX_INDEX> & ivolpoly,
   std::vector<POLY_VERTEX_INDEX> & poly_vertex,
   IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   DUALISO_INFO & dualiso_info);

  void extract_ivolpoly_dual_to_grid_vertices
  (const IVOLDUAL_ENCODED_GRID & encoded_grid,
   std::vector<ISO_VERTEX_INDEX> & ivolpoly,
   std::vector<POLY_VERTEX_INDEX> & poly_vertex,
   IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   DUALISO_INFO & dualiso_info);


  // **************************************************
  // SPLIT DUAL INTERVAL VOLUME VERTICES
  // **************************************************

  void split_dual_ivolvert
  (const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const std::vector<VERTEX_INDEX> & ivolpoly_cube, 
   const std::vector<POLY_VERTEX_INDEX> & poly_vertex,
   const IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   std::vector<GRID_CUBE_DATA> & cube_list,
   std::vector<DUAL_ISOVERT> & ivolv_list, 
   std::vector<ISO_VERTEX_INDEX> & isopoly,
   int & num_split);


  // **************************************************
  // POSITION INTERVAL VOLUME VERTICES
  // **************************************************

  void position_all_dual_ivol_vertices
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const IJKDUAL::ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
   const SCALAR_TYPE isovalue0,
   const SCALAR_TYPE isovalue1,
   const std::vector<DUAL_ISOVERT> & ivolv_list, 
   COORD_TYPE * vertex_coord);

  void position_all_dual_ivol_vertices
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const IJKDUAL::ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
   const SCALAR_TYPE isovalue0,
   const SCALAR_TYPE isovalue1,
   const std::vector<DUAL_ISOVERT> & ivolv_list, 
   COORD_ARRAY & vertex_coord);

  void position_dual_ivolv_centroid_multi
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const SCALAR_TYPE isovalue0,
   const SCALAR_TYPE isovalue1,
   const DUAL_ISOVERT & ivolv_info,
   const CUBE_FACE_INFO & cube,
   COORD_TYPE * vcoord, 
   COORD_TYPE * temp_coord0, COORD_TYPE * temp_coord1, 
   COORD_TYPE * temp_coord2);

  void position_dual_ivolv_on_isosurface_centroid_multi
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const SCALAR_TYPE isovalue,
   const DUAL_ISOVERT & ivolv_info,
   const CUBE_FACE_INFO & cube,
   COORD_TYPE * vcoord, 
   COORD_TYPE * temp_coord0, COORD_TYPE * temp_coord1, 
   COORD_TYPE * temp_coord2);

  void position_dual_ivolv_in_interval_volume_centroid_multi
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const SCALAR_TYPE isovalue0,
   const SCALAR_TYPE isovalue1,
   const DUAL_ISOVERT & ivolv_info,
   const CUBE_FACE_INFO & cube,
   COORD_TYPE * vcoord, 
   COORD_TYPE * temp_coord0, COORD_TYPE * temp_coord1, 
   COORD_TYPE * temp_coord2);

}

#endif