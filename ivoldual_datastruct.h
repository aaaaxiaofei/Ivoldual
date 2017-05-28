/// \file ivoldual_datastruct.h
/// Data structure definitions for ivoldual.
/// Copies types from ijkdual.h.

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

#ifndef _IVOLDUAL_DATASTRUCT_
#define _IVOLDUAL_DATASTRUCT_

#include "ijkcube.txx"
#include "ijkdualtable.txx"

#include "ijkdual_datastruct.h"
#include "ivoldual_types.h"


namespace IVOLDUAL {

  // **************************************************
  // INTERVAL VOLUME GRID VERTEX ENCODING
  // **************************************************

  /// Type for encoding of grid vertices.
  /// - Encoded vertices have scalar values 0,1,2 or 3.
  typedef unsigned char GRID_VERTEX_ENCODING;


  // **************************************************
  // INTERVAL VOLUME LOOKUP TABLE ENTRY
  // **************************************************

  typedef typename IJKDUALTABLE::
  IVOLDUAL_TABLE_ENTRY<int,ISO_VERTEX_INDEX,TABLE_INDEX> IVOLDUAL_TABLE_ENTRY;


  // **************************************************
  // INTERVAL VOLUME LOOKUP TABLE
  // **************************************************

  typedef typename IJKDUALTABLE::
  IVOLDUAL_CUBE_TABLE<4,int,int,TABLE_INDEX,IVOLDUAL_TABLE_ENTRY>
  IVOLDUAL_CUBE_TABLE;


  // **************************************************
  // GRID DATA STRUCTURES
  // **************************************************

  typedef IJKDUAL::DUALISO_GRID DUALISO_GRID;
  typedef IJKDUAL::DUALISO_SCALAR_GRID_BASE DUALISO_SCALAR_GRID_BASE;
  typedef IJKDUAL::DUALISO_SCALAR_GRID DUALISO_SCALAR_GRID;

  /// Type of grid encoding grid vertices.
  /// - 0: Below lower isovalue.
  /// - 1: Between lower and upper isovalue.
  /// - 2: Between lower and upper isovalue.
  /// - 3: Above upper isovalue.
  typedef IJK::SCALAR_GRID<DUALISO_GRID, GRID_VERTEX_ENCODING> 
  IVOLDUAL_ENCODED_GRID;

  /// Type of cube with face info.
  typedef IJK::CUBE_FACE_INFO<int,int,int> CUBE_FACE_INFO;


  // **************************************************
  // DUAL ISOSURFACE VERTICES
  // **************************************************

  typedef IJKDUAL::DUAL_ISOVERT DUAL_ISOVERT;


  // **************************************************
  // GRID CUBE DATA
  // **************************************************

  typedef IJKDUAL::GRID_CUBE_DATA GRID_CUBE_DATA;


  // **************************************************
  // INTERVAL VOLUME POLY INFO
  // **************************************************

  class IVOLDUAL_POLY_INFO {

  public:

    /// If true, dual to grid edge.
    /// Else dual to grid vertex.
    bool flag_dual_to_edge;

    /// Grid vertex or lower/leftmost endoint of grid edge.
    ISO_VERTEX_INDEX v0;

    /// Direction of dual edge if dual to grid edge.
    int edge_direction;

    void SetDualToEdge
    (const ISO_VERTEX_INDEX iv0, const int edge_direction)
    {
      flag_dual_to_edge = true;
      v0 = iv0;
      this->edge_direction = edge_direction;
    }

    void SetDualToVertex(const ISO_VERTEX_INDEX iv0)
    {
      flag_dual_to_edge = false;
      v0 = iv0;
    }

  };


  /// Array of IVOLDUAL_POLY_INFO.
  typedef std::vector<IVOLDUAL_POLY_INFO> IVOLDUAL_POLY_INFO_ARRAY;


  // **************************************************
  // DUAL CONTOURING INTERVAL VOLUME
  // **************************************************

  typedef typename IJKDUAL::DUAL_ISOSURFACE_BASE<IVOLDUAL_POLY_INFO> 
  DUAL_INTERVAL_VOLUME;


  // **************************************************
  // DUAL CONTOURING INPUT DATA AND DATA STRUCTURES
  // **************************************************

  typedef IJKDUAL::DUALISO_DATA_FLAGS DUALISO_DATA_FLAGS;
  typedef IJKDUAL::DUALISO_DATA DUALISO_DATA;


  // **************************************************
  // DUALISO TIME
  // **************************************************

  typedef IJKDUAL::DUALISO_TIME DUALISO_TIME;


  // **************************************************
  // DUALISO INFO
  // **************************************************

  typedef IJKDUAL::DUALISO_INFO DUALISO_INFO;


  // **************************************************
  // MERGE DATA
  // **************************************************

  typedef IJKDUAL::MERGE_DATA MERGE_DATA;

}

#endif
