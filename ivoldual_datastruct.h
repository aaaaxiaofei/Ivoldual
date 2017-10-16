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
  IVOLDUAL_TABLE_ENTRY<int,ISO_VERTEX_INDEX,FACET_BITS_TYPE,TABLE_INDEX> 
  IVOLDUAL_TABLE_ENTRY;


  // **************************************************
  // INTERVAL VOLUME LOOKUP TABLE
  // **************************************************

  typedef typename IJKDUALTABLE::
  IVOLDUAL_CUBE_DOUBLE_TABLE<4,int,int,TABLE_INDEX,IVOLDUAL_TABLE_ENTRY>
  IVOLDUAL_CUBE_TABLE;


  // **************************************************
  // GRID DATA STRUCTURES
  // **************************************************

  typedef IJKDUAL::DUALISO_GRID DUALISO_GRID;
  typedef IJKDUAL::DUALISO_SCALAR_GRID_BASE DUALISO_SCALAR_GRID_BASE;
  typedef IJKDUAL::DUALISO_SCALAR_GRID DUALISO_SCALAR_GRID;

  class IVOLDUAL_SCALAR_GRID: public DUALISO_SCALAR_GRID {

  public:
    IVOLDUAL_SCALAR_GRID() {};

    // Subdivide \a scalar grid.
    void Subdivide
    (const DUALISO_SCALAR_GRID_BASE & scalar_grid2, const int subdivide_period, 
     const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1);

    // Interpolate between subdivide vertices.
    void SubdivideInterpolate
    (const int subdivide_period, const SCALAR_TYPE isovalue0, 
     const SCALAR_TYPE isovalue1); 

    void AddOuterLayer(IVOLDUAL_SCALAR_GRID & scalar_grid3);
  };

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

  /// Information about an interval volume vertex.
  class DUAL_IVOLVERT: public IJKDUAL::DUAL_ISOVERT {

  public:

    /// If true, the mesh is missing some interval volume hexahedra
    ///   which would normally be incident on the vertex.
    /// The hexahedra are missing because they are dual to grid edges
    ///   or vertices which are on the grid boundary.
    /// - NOT YET IMPLEMENTED.
    bool flag_missing_ivol_hexahedra;
  };
  typedef std::vector<DUAL_IVOLVERT> DUAL_IVOLVERT_ARRAY;
  

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

    /// If true, reverse the orientation of the polytope.
    bool flag_reverse_orient;

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
      flag_reverse_orient = false;
    }

    void SetDualToEdge
    (const ISO_VERTEX_INDEX iv0, const int edge_direction,
     const bool flag_reverse_orient)
    {
      flag_dual_to_edge = true;
      v0 = iv0;
      this->edge_direction = edge_direction;
      this->flag_reverse_orient = flag_reverse_orient;
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
  // DUAL CONTOURING INTERVAL VOLUME FLAGS
  // **************************************************

  /// Data structure for flags controlling interval volume construction.<br>
  /// To add a flag:
  /// - Add a field in this class;
  /// - Initialize the field in IVOLDUAL_DATA_FLAGS::Init()
  ///   in file ivoldual_datastruct.cxx;
  /// - Add a command line option type to OPTION_TYPE in file ivoldualIO.cxx;
  /// - Add a command line option using function option.AddOptionNoArg(),
  ///   or option.AddOption1Arg, in procedure set_command_line_options 
  ///   in file ivoldualIO.cxx.  Usage message and the help message are
  ///   passed as parameters in those commands;
  /// - Add a case to the switch statement in procedure process_option
  ///   that sets the appropriate IVOLDUAL_DATA_FLAGS field when the
  ///   option is invoked on the command line.
  class IVOLDUAL_DATA_FLAGS:public IJKDUAL::DUALISO_DATA_FLAGS {

  protected:
    void Init();

  public:

    /// If true, split vertices in ambiguous cubes sharing facets.
    bool flag_split_ambig_pairs;
    bool flag_split_ambig_pairsB;

    /// Hexahedral triangulation method.
    HEX_TRI_METHOD hex_tri_method;

    /// If true, add extra vertices dual to iso hex for triangulation.
    bool flag_add_isov_dual_to_hexahedra;

    /// Hexahedra orientation.
    /// - If true, orient hexahedra so that facet normals point in.
    bool flag_orient_in;

    /// Default encoded value (1 or 2) for vertices in the interior
    ///   of the interval volume.
    GRID_VERTEX_ENCODING default_interior_code;

  public:

    /// Constructor.
    IVOLDUAL_DATA_FLAGS() { Init(); };

    /// Set flags.
    void Set(const IVOLDUAL_DATA_FLAGS & data_flags)
    { *this = data_flags; };
  };


  // **************************************************
  // DUAL CONTOURING INTERVAL VOLUME
  // **************************************************

  /// Representation of dDual contouring interval volume created
  /// by dual contouring interval volume algorithm.
  class DUAL_INTERVAL_VOLUME:
    public IJKDUAL::DUAL_ISOSURFACE_BASE<IVOLDUAL_POLY_INFO> {

  public:
    /// List of interval volume vertices.
    DUAL_IVOLVERT_ARRAY ivolv_list;

  public:
    DUAL_INTERVAL_VOLUME
    (const int dimension, const VERTEX_INDEX numv_per_ivolpoly):
      DUAL_ISOSURFACE_BASE<IVOLDUAL_POLY_INFO>(dimension, numv_per_ivolpoly)
    {};
  };


  // **************************************************
  // IVOLDUAL INPUT DATA AND DATA STRUCTURES
  // **************************************************

  /// Input data to ivoldual.
  class IVOLDUAL_DATA:
    public IJKDUAL::DUALISO_DATA_BASE<IVOLDUAL_SCALAR_GRID,IVOLDUAL_DATA_FLAGS> 
  {
  public:
    IVOLDUAL_DATA() {};

    void SubdivideScalarGrid      /// Subdivide scalar_grid.
      (const DUALISO_SCALAR_GRID_BASE & scalar_grid2, IVOLDUAL_SCALAR_GRID & scalar_grid3, 
       const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1); 

    /// Copy, subsample, supersample or subdivide scalar grid.
    /// @pre At most one of flag_subsample, flag_supersample or
    ///   flag_subdivide may be true.
    void SetScalarGrid
      (const DUALISO_SCALAR_GRID_BASE & scalar_grid2, 
       const bool flag_subsample, const int subsample_resolution,
       const bool flag_supersample, const int supersample_resolution, 
       const bool flag_subdivide, const SCALAR_TYPE isovalue0, 
       const SCALAR_TYPE isovalue1);

  };

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
