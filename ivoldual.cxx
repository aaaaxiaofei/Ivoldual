/// \file ivoldual.cxx
/// Generate interval volume in arbitrary dimensions using dual contouring.
/// Version 0.1

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



#include "ijktime.txx"

#include "ijkdual.txx"
#include "ijkdual_extract.txx"
#include "ijkdual_position.txx"

#include "ivoldual.h"
#include "ivoldual_datastruct.h"

#include "ijkisopoly.txx"

using namespace IJK;
using namespace IVOLDUAL;


// **************************************************
// DUAL CONTOURING INTERVAL VOLUME
// **************************************************

// Construct interval volume using dual contouring.
void IVOLDUAL::dual_contouring_interval_volume
(const DUALISO_DATA & dualiso_data, 
 const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
 DUAL_INTERVAL_VOLUME & dual_interval_volume, DUALISO_INFO & dualiso_info)
{
  const int dimension = dualiso_data.ScalarGrid().Dimension();
  const AXIS_SIZE_TYPE * axis_size = dualiso_data.ScalarGrid().AxisSize();
  PROCEDURE_ERROR error("dual_contouring_interval_volume");

  clock_t t_start = clock();

  if (!dualiso_data.Check(error)) { throw error; };

  dual_interval_volume.Clear();
  dualiso_info.time.Clear();

  IJKDUAL::ISO_MERGE_DATA merge_data(dimension, axis_size);

  dual_contouring_interval_volume
    (dualiso_data.ScalarGrid(), isovalue0, isovalue1, dualiso_data,
     dual_interval_volume.isopoly_vert, dual_interval_volume.isopoly_info, 
     dual_interval_volume.vertex_coord, merge_data, dualiso_info);

  // store times
  clock_t t_end = clock();
  clock2seconds(t_end-t_start, dualiso_info.time.total);
}


// Construct interval volume using dual contouring.
// Returns list of isosurface polytope vertices
// Version which creates ivoldual_table.
void IVOLDUAL::dual_contouring_interval_volume
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
 const DUALISO_DATA_FLAGS & param,
 std::vector<ISO_VERTEX_INDEX> & ivolpoly_vert,
 IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
 COORD_ARRAY & vertex_coord,
 MERGE_DATA & merge_data, 
 DUALISO_INFO & dualiso_info)
{
  const int dimension = scalar_grid.Dimension();
  const bool flag_separate_neg = param.SeparateNegFlag();
  const bool flag_always_separate_opposite(true);

  IVOLDUAL_CUBE_TABLE
    ivoldual_table(dimension, flag_separate_neg, 
                   flag_always_separate_opposite);

  dual_contouring_interval_volume
    (scalar_grid, isovalue0, isovalue1, ivoldual_table, param, ivolpoly_vert, 
     ivolpoly_info, vertex_coord, merge_data, dualiso_info);
}


// Construct interval volume using dual contouring.
// Returns list of isosurface polytope vertices
//   and list of isosurface vertex coordinates
void IVOLDUAL::dual_contouring_interval_volume
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const DUALISO_DATA_FLAGS & param,
 std::vector<ISO_VERTEX_INDEX> & ivolpoly_vert,
 IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
 COORD_ARRAY & vertex_coord,
 MERGE_DATA & merge_data, 
 DUALISO_INFO & dualiso_info)
{
  std::vector<DUAL_ISOVERT> ivolv_list;
  std::vector<GRID_CUBE_DATA> cube_ivolv_list;

  dual_contouring_interval_volume
    (scalar_grid, isovalue0, isovalue1, ivoldual_table, param, ivolpoly_vert,
     cube_ivolv_list, ivolv_list, ivolpoly_info, vertex_coord, 
     merge_data, dualiso_info);
}


// Construct interval volume using dual contouring.
// Returns list of isosurface polytope vertices
//   and list of isosurface vertex coordinates
void IVOLDUAL::dual_contouring_interval_volume
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
 DUALISO_INFO & dualiso_info)
{
  const int dimension = scalar_grid.Dimension();
  const bool flag_separate_neg = param.SeparateNegFlag();
  const bool flag_always_separate_opposite(true);
  IJK::PROCEDURE_ERROR error("dual_contouring_interval_volume");
  clock_t t0, t1, t2;

  if (scalar_grid.Dimension() != ivoldual_table.Dimension()) {
    error.AddMessage
      ("Programming error.  Incorrect isodual table dimension.");
    error.AddMessage
      ("  Interval volume table dimension does not match scalar grid dimension.");
    error.AddMessage
      ("    Interval volume table dimension: ", ivoldual_table.Dimension(), "");
    error.AddMessage
      ("    Scalar grid dimension: ", scalar_grid.Dimension(), "");
    throw error;
  }

  t0 = clock();

  ivolpoly_vert.clear();
  dualiso_info.time.Clear();

  IVOLDUAL_ENCODED_GRID encoded_grid;
  const GRID_VERTEX_ENCODING default_interior_code = 2;
  encode_grid_vertices
    (scalar_grid, isovalue0, isovalue1, default_interior_code, 
     encoded_grid, dualiso_info);

  std::vector<ISO_VERTEX_INDEX> ivolpoly;
  std::vector<POLY_VERTEX_INDEX> poly_vertex;
  extract_dual_ivolpoly
    (encoded_grid, ivolpoly, poly_vertex, ivolpoly_info, dualiso_info);
  t1 = clock();

  std::vector<ISO_VERTEX_INDEX> cube_list;
  std::vector<ISO_VERTEX_INDEX> ivolpoly_cube;
  merge_identical(ivolpoly, cube_list, ivolpoly_cube, merge_data);
  t2 = clock();

  set_grid_cube_indices(scalar_grid, cube_list, cube_ivolv_list);

  compute_cube_ivoltable_info
    (encoded_grid, ivoldual_table, cube_ivolv_list);

  VERTEX_INDEX num_split;
  split_dual_ivolvert
    (ivoldual_table, ivolpoly_cube, poly_vertex, ivolpoly_info, 
     cube_ivolv_list, ivolv_list, ivolpoly_vert, num_split);

  IJKDUAL::ISODUAL_CUBE_TABLE_AMBIG isodual_table
    (dimension, flag_separate_neg, flag_always_separate_opposite);
  position_all_dual_ivol_vertices
    (scalar_grid, ivoldual_table, isodual_table, isovalue0, isovalue1, 
     ivolv_list, vertex_coord);

  dualiso_info.multi_isov.num_cubes_multi_isov = num_split;
  dualiso_info.multi_isov.num_cubes_single_isov =
    cube_list.size() - num_split;

  // store times
  IJK::clock2seconds(t1-t0, dualiso_info.time.extract);
  IJK::clock2seconds(t2-t1, dualiso_info.time.merge);

}


// **************************************************
// ENCODE VERTICES
// **************************************************

void IVOLDUAL::encode_grid_vertices
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
 const GRID_VERTEX_ENCODING default_interior_code,
 IVOLDUAL_ENCODED_GRID & encoded_grid,
 DUALISO_INFO & dualiso_info)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();

  encoded_grid.SetSize(dimension, axis_size);

  for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++) {
    if (scalar_grid.Scalar(iv) < isovalue0) 
      { encoded_grid.Set(iv, 0); }
    else if (scalar_grid.Scalar(iv) > isovalue1) 
      { encoded_grid.Set(iv, 3); }
    else {
      // Set all vertices with scalar value 
      //   in range [isovalue0,isovalue1] to default_interior_code
      encoded_grid.Set(iv,default_interior_code);
    }
  }
}


// **************************************************
// COMPUTE IVOLTABLE INFO FOR EACH ACTIVE GRID CUBE
// **************************************************

namespace {

  TABLE_INDEX compute_table_index_from_encoded_grid
  (const IVOLDUAL_ENCODED_GRID & encoded_grid,
   const int num_vertex_types,
   const VERTEX_INDEX icube)
  {
    const int num_cube_vertices = encoded_grid.NumCubeVertices();
    TABLE_INDEX table_index = 0;
    TABLE_INDEX factor = 1;
    for (int j = 0; j < num_cube_vertices; j++) {
      const int j2 = num_cube_vertices-1-j;
      const VERTEX_INDEX jv2 = encoded_grid.CubeVertex(icube, j2);
      table_index = 
        (table_index*num_vertex_types) + encoded_grid.Scalar(jv2);
    }

    return(table_index);
  }

}

void IVOLDUAL::compute_cube_ivoltable_info
(const IVOLDUAL_ENCODED_GRID & encoded_grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 std::vector<GRID_CUBE_DATA> & cube_ivolv_list)
{
  const int dimension = encoded_grid.Dimension();
  const int num_vertex_types = ivoldual_table.NumVertexTypes();

  for (int i = 0; i < cube_ivolv_list.size(); i++) {
    const VERTEX_INDEX cube_index = cube_ivolv_list[i].cube_index;
    const TABLE_INDEX table_index =
      compute_table_index_from_encoded_grid
      (encoded_grid, num_vertex_types, cube_index);

    cube_ivolv_list[i].table_index = table_index;
    cube_ivolv_list[i].num_isov = ivoldual_table.NumIsoVertices(table_index);
  }
}


// **************************************************
// EXTRACT DUAL INTERVAL VOLUME POLYTOPES
// **************************************************

void IVOLDUAL::extract_dual_ivolpoly
(const IVOLDUAL_ENCODED_GRID & encoded_grid,
 std::vector<ISO_VERTEX_INDEX> & ivolpoly,
 std::vector<POLY_VERTEX_INDEX> & poly_vertex,
 IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
 DUALISO_INFO & dualiso_info)
{
  dualiso_info.time.extract = 0;

  clock_t t0 = clock();

  // Initialize output
  ivolpoly.clear();

  if (encoded_grid.NumCubeVertices() < 1) { return; }

  extract_ivolpoly_dual_to_grid_edges
    (encoded_grid, ivolpoly, poly_vertex, ivolpoly_info, dualiso_info);
  extract_ivolpoly_dual_to_grid_vertices
    (encoded_grid, ivolpoly, poly_vertex, ivolpoly_info, dualiso_info);

  clock_t t1 = clock();
  IJK::clock2seconds(t1-t0, dualiso_info.time.extract);
}

namespace {

  void add_grid_facet_vertices
  (const DUALISO_GRID & grid, const VERTEX_INDEX iv0, const int edge_dir,
   std::vector<ISO_VERTEX_INDEX> & ivolpoly,
   std::vector<POLY_VERTEX_INDEX> & poly_vertex)
  {
    const int num_facet_vertices = grid.NumFacetVertices();

    for (int k = 0; k < num_facet_vertices; k++) {

      // jcube is the index of the cube containing the k'th vertex
      //   of the facet.
      const VERTEX_INDEX jcube = grid.FacetVertex(iv0, edge_dir, k);
      ivolpoly.push_back(jcube);
      poly_vertex.push_back(edge_dir*num_facet_vertices+k);
    }
  }

  bool does_grid_edge_have_dual_ivolpoly
  (const IVOLDUAL_ENCODED_GRID & encoded_grid,
   const VERTEX_INDEX iend0, const int edge_dir)
  {
    VERTEX_INDEX iend1 = encoded_grid.NextVertex(iend0, edge_dir);
    GRID_VERTEX_ENCODING s0 = encoded_grid.Scalar(iend0);
    GRID_VERTEX_ENCODING s1 = encoded_grid.Scalar(iend1);

    if (s0 > s1) { std::swap(s0,s1); }

    if (s0 == 0) {
      if (s1 >= 2) { return(true); }
    }
    else if (s0 == 1) {
      if (s1 == 3) { return(true); }
    }

    return(false);
  }

  bool does_grid_vertex_have_dual_ivolpoly
  (const IVOLDUAL_ENCODED_GRID & encoded_grid, const VERTEX_INDEX iv0)
  {
    GRID_VERTEX_ENCODING s0 = encoded_grid.Scalar(iv0);

    if (s0 == 1 || s0 == 2) { return(true); }
    return(false);
  }

}

void IVOLDUAL::extract_ivolpoly_dual_to_grid_edges
(const IVOLDUAL_ENCODED_GRID & encoded_grid,
 std::vector<ISO_VERTEX_INDEX> & ivolpoly,
 std::vector<POLY_VERTEX_INDEX> & poly_vertex,
 IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
 DUALISO_INFO & dualiso_info)
{
  const int num_facet_vertices = encoded_grid.NumFacetVertices();

  if (num_facet_vertices == 0) { return; }

  IJK_FOR_EACH_INTERIOR_GRID_EDGE(iend0, edge_dir, encoded_grid, VERTEX_INDEX) {

    if (does_grid_edge_have_dual_ivolpoly(encoded_grid, iend0, edge_dir)) {

      const VERTEX_INDEX increment = 
        encoded_grid.FacetVertexIncrement(edge_dir,num_facet_vertices-1);
      const VERTEX_INDEX iv0 = iend0 - increment;

      for (int j = 0; j < 2; j++) {
        add_grid_facet_vertices
          (encoded_grid, iv0, edge_dir, ivolpoly, poly_vertex);
      }

      IVOLDUAL_POLY_INFO info;
      info.SetDualToEdge(iend0, edge_dir);
      ivolpoly_info.push_back(info);
    }
  }

}


void IVOLDUAL::extract_ivolpoly_dual_to_grid_vertices
(const IVOLDUAL_ENCODED_GRID & encoded_grid,
 std::vector<ISO_VERTEX_INDEX> & ivolpoly,
 std::vector<POLY_VERTEX_INDEX> & poly_vertex,
 IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
 DUALISO_INFO & dualiso_info)
{
  const int num_facet_vertices = encoded_grid.NumFacetVertices();
  const int num_cube_vertices = encoded_grid.NumCubeVertices();

  if (num_facet_vertices == 0) { return; }

  IJK_FOR_EACH_INTERIOR_GRID_VERTEX(iv0, encoded_grid, VERTEX_INDEX) {

    if (does_grid_vertex_have_dual_ivolpoly(encoded_grid, iv0)) {
      const VERTEX_INDEX increment = 
        encoded_grid.CubeVertexIncrement(num_cube_vertices-1);
      const VERTEX_INDEX iv1 = iv0 - increment;


      for (int k = 0; k < num_cube_vertices; k++) {
        ivolpoly.push_back(encoded_grid.CubeVertex(iv1, k));
        poly_vertex.push_back(k);
      }

      IVOLDUAL_POLY_INFO info;
      info.SetDualToVertex(iv0);
      ivolpoly_info.push_back(info);
    }
  }
}


// **************************************************
// SPLIT DUAL INTERVAL VOLUME VERTICES
// **************************************************

namespace {

  void set_dual_ivolv_vertices
  (const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const std::vector<GRID_CUBE_DATA> & cube_list,
   const std::vector<VERTEX_INDEX> & ivolpoly_cube, 
   const std::vector<POLY_VERTEX_INDEX> & poly_vertex,
   const IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   std::vector<ISO_VERTEX_INDEX> & ivolpoly_vert)
  {
    const int dimension = ivoldual_table.Dimension();
    IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);
    const int num_cube_vertices = cube.NumVertices();
    const int num_facet_vertices = cube.NumFacetVertices();

    ivolpoly_vert.resize(ivolpoly_cube.size());

    for (ISO_VERTEX_INDEX i = 0; i < ivolpoly_cube.size(); i++) {
      const VERTEX_INDEX k = ivolpoly_cube[i];
      const TABLE_INDEX it = cube_list[k].table_index;
      const int ipoly = i/num_cube_vertices;
      const int poly_vertex_i = poly_vertex[i];
      const int icorner = i%num_cube_vertices;

      // Compute index of vertex opposite to poly_vertex[i]
      if (ivolpoly_info[ipoly].flag_dual_to_edge) {
        const int ifacet = cube.FacetIndex(poly_vertex_i);
        const int j = poly_vertex_i - ifacet*num_facet_vertices;
        int opposite_vertex = (num_facet_vertices-1) - j;
        opposite_vertex += (ifacet*num_facet_vertices);

        if (icorner < num_cube_vertices/2) {
          // Vertex on lower facet.
          ivolpoly_vert[i] = cube_list[k].first_isov +
            ivoldual_table.LowerIncident(it, opposite_vertex);
        }
        else {
          // Vertex on upper facet.
          ivolpoly_vert[i] = cube_list[k].first_isov +
            ivoldual_table.UpperIncident(it, opposite_vertex);
        }
      }
      else {
        const int opposite_vertex = cube.OppositeVertex(poly_vertex_i);

        ivolpoly_vert[i] = cube_list[k].first_isov +
          ivoldual_table.IncidentIVolVertex(it, opposite_vertex);
      }
    }
  }

}

void IVOLDUAL::split_dual_ivolvert
(const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const std::vector<VERTEX_INDEX> & ivolpoly_cube, 
 const std::vector<POLY_VERTEX_INDEX> & poly_vertex,
 const IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
 std::vector<GRID_CUBE_DATA> & cube_list,
 std::vector<DUAL_ISOVERT> & ivolv_list, 
 std::vector<ISO_VERTEX_INDEX> & ivolpoly_vert,
 int & num_split)
{
  IJK::construct_dual_isovert_list(ivoldual_table, cube_list, ivolv_list);

  set_dual_ivolv_vertices
    (ivoldual_table, cube_list, ivolpoly_cube, poly_vertex, ivolpoly_info, 
     ivolpoly_vert);

  compute_num_splitB(ivoldual_table, cube_list, num_split);
}


// **************************************************
// POSITION INTERVAL VOLUME VERTICES
// **************************************************

void IVOLDUAL::position_all_dual_ivol_vertices
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const IJKDUAL::ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
 const SCALAR_TYPE isovalue0,
 const SCALAR_TYPE isovalue1,
 const std::vector<DUAL_ISOVERT> & ivolv_list, 
 COORD_TYPE * vertex_coord)
{
  const int dimension = scalar_grid.Dimension();
  DUAL_ISOVERT dual_isovert;

  if (dimension < 1) { return; }

  IJK::ARRAY<COORD_TYPE> coord0(dimension);
  IJK::ARRAY<COORD_TYPE> coord1(dimension);
  IJK::ARRAY<COORD_TYPE> coord2(dimension);
  CUBE_FACE_INFO cube(dimension);

  for (ISO_VERTEX_INDEX ivolv = 0; ivolv < ivolv_list.size(); ivolv++) {
    position_dual_ivolv_centroid_multi
      (scalar_grid, ivoldual_table, isovalue0, isovalue1, 
       ivolv_list[ivolv], cube, 
       vertex_coord+ivolv*dimension, 
       coord0.Ptr(), coord1.Ptr(), coord2.Ptr());
  }

}

void IVOLDUAL::position_all_dual_ivol_vertices
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const IJKDUAL::ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
 const SCALAR_TYPE isovalue0,
 const SCALAR_TYPE isovalue1,
 const std::vector<DUAL_ISOVERT> & ivolv_list, 
 COORD_ARRAY & vertex_coord)
{
  const int dimension = scalar_grid.Dimension();

  vertex_coord.resize(ivolv_list.size()*dimension);
  position_all_dual_ivol_vertices
    (scalar_grid, ivoldual_table, isodual_table, isovalue0, isovalue1, 
     ivolv_list, &(vertex_coord.front()));
}


void IVOLDUAL::position_dual_ivolv_centroid_multi
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const SCALAR_TYPE isovalue0,
 const SCALAR_TYPE isovalue1,
 const DUAL_ISOVERT & ivolv_info,
 const CUBE_FACE_INFO & cube,
 COORD_TYPE * vcoord, 
 COORD_TYPE * temp_coord0, COORD_TYPE * temp_coord1, COORD_TYPE * temp_coord2)
{
  const int dimension = scalar_grid.Dimension();
  const int ivolv = ivolv_info.patch_index;
  const TABLE_INDEX it = ivolv_info.table_index;

  if (ivoldual_table.OnLowerIsosurface(it, ivolv)) {
    position_dual_ivolv_on_isosurface_centroid_multi
      (scalar_grid, ivoldual_table, isovalue0, ivolv_info, cube,
       vcoord, temp_coord0, temp_coord1, temp_coord2);
  }
  else if (ivoldual_table.OnUpperIsosurface(it, ivolv)) {
    position_dual_ivolv_on_isosurface_centroid_multi
      (scalar_grid, ivoldual_table, isovalue1, ivolv_info, cube,
       vcoord, temp_coord0, temp_coord1, temp_coord2);
  }
  else {
    position_dual_ivolv_in_interval_volume_centroid_multi
      (scalar_grid, ivoldual_table, isovalue0, isovalue1, ivolv_info, cube,
       vcoord, temp_coord0, temp_coord1, temp_coord2);
  }

}


void IVOLDUAL::position_dual_ivolv_on_isosurface_centroid_multi
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const SCALAR_TYPE isovalue,
 const DUAL_ISOVERT & ivolv_info,
 const CUBE_FACE_INFO & cube,
 COORD_TYPE * vcoord, 
 COORD_TYPE * temp_coord0, COORD_TYPE * temp_coord1, COORD_TYPE * temp_coord2)
{
  const int dimension = scalar_grid.Dimension();
  const VERTEX_INDEX icube = ivolv_info.cube_index;
  const int ivolv = ivolv_info.patch_index;
  const TABLE_INDEX table_index = ivolv_info.table_index;
  int num_intersected_edges = 0;
  IJK::set_coord(dimension, 0.0, vcoord);

  for (int ie = 0; ie < cube.NumEdges(); ie++) {
    int k0 = cube.EdgeEndpoint(ie, 0);
    int k1 = cube.EdgeEndpoint(ie, 1);
    VERTEX_INDEX iend0 = scalar_grid.CubeVertex(icube, k0);
    VERTEX_INDEX iend1 = scalar_grid.CubeVertex(icube, k1);

    bool flag_intersect = false;
    if (ivoldual_table.EdgeHasDualIVolPoly(table_index, ie)) {

      if (ivoldual_table.UpperIncident(table_index, ie) == ivolv ||
          ivoldual_table.LowerIncident(table_index, ie) == ivolv) 
        { flag_intersect = true; }
    }

    if (ivoldual_table.IsInIntervalVolume(table_index, k0)) {
      if (ivoldual_table.IncidentIVolVertex(table_index, k0) == ivolv) {
        if (ivoldual_table.IsBelowIntervalVolume(table_index, k1)) {
          if (ivoldual_table.OnLowerIsosurface(table_index, ivolv)) 
            { flag_intersect = true;  }
        }
        else if (ivoldual_table.IsAboveIntervalVolume(table_index, k1)) {
          if (ivoldual_table.OnUpperIsosurface(table_index, ivolv))
            { flag_intersect = true; }
        }
      }
    }
    else if (ivoldual_table.IsInIntervalVolume(table_index, k1)) {
      if (ivoldual_table.IncidentIVolVertex(table_index, k1) == ivolv) {
        if (ivoldual_table.IsBelowIntervalVolume(table_index, k0)) {
          if (ivoldual_table.OnLowerIsosurface(table_index, ivolv)) 
            { flag_intersect = true; }
        }
        else if (ivoldual_table.IsAboveIntervalVolume(table_index, k0)) {
          if (ivoldual_table.OnUpperIsosurface(table_index, ivolv))
            { flag_intersect = true; }
        }
      }
    }

    if (flag_intersect) {

      SCALAR_TYPE s0 = scalar_grid.Scalar(iend0);
      SCALAR_TYPE s1 = scalar_grid.Scalar(iend1);

      scalar_grid.ComputeCoord(iend0, temp_coord0);
      scalar_grid.ComputeCoord(iend1, temp_coord1);


      if ((s0 < isovalue && s1 < isovalue) || 
          (s0 > isovalue && s1 > isovalue)) {
        // Use edge midpoint.
        IJK::linear_interpolate_coord
          (dimension, 0.5, temp_coord0, temp_coord1, temp_coord2);
      }
      else {
        IJK::linear_interpolate_coord
          (dimension, s0, temp_coord0, s1, temp_coord1, isovalue, 
           temp_coord2);
      }

      IJK::add_coord(dimension, vcoord, temp_coord2, vcoord);
      num_intersected_edges++;
    }
  }

  if (num_intersected_edges > 0) {
    IJK::multiply_coord
      (dimension, 1.0/num_intersected_edges, vcoord, vcoord);
  }
  else {
    scalar_grid.ComputeCubeCenterCoord(icube, vcoord);
  }

}


namespace {

  /// Return true if interval volume polytope dual to iv0 and
  ///   incident on ivolv intersects grid edge (iv0,iv1) 
  ///   and some vertex of (iv0,iv1) is not in the interval volume.
  /// @param isovalueX Isovalue of intersection point.
  bool does_poly_dual_to_grid_vertex_intersect_grid_edge
  (const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const SCALAR_TYPE isovalue0,
   const SCALAR_TYPE isovalue1,
   const TABLE_INDEX table_index,
   const int ivolv,
   const int ie,
   const int iv0,
   const int iv1,
   const SCALAR_TYPE s0,
   const SCALAR_TYPE s1,
   SCALAR_TYPE & isovalueX)
  {
    if (ivoldual_table.IsInIntervalVolume(table_index, iv0)) {
      if (ivoldual_table.IncidentIVolVertex(table_index, iv0) == ivolv) {
        if (ivoldual_table.IsBelowIntervalVolume(table_index, iv1)) {
          if (ivoldual_table.EdgeHasDualIVolPoly(table_index, ie)) 
            { isovalueX = (s0+isovalue0)/2.0; }
          else 
            { isovalueX = isovalue0; }
          return(true);
        }
        else if (ivoldual_table.IsAboveIntervalVolume(table_index, iv1)) {
          if (ivoldual_table.EdgeHasDualIVolPoly(table_index, ie)) 
            { isovalueX = (s0+isovalue1)/2.0; }
          else 
            { isovalueX = isovalue1; }
          return(true);
        }
        else {
          return(false);
        }
      }
    }

    return(false);
  }
   

  /// Return true if grid edge ie intersects interval volume polytope incident 
  /// on ivolv and some vertex of ie is not in the interval volume.
  /// @param isovalueX Isovalue of intersection point.
  bool does_grid_edge_intersect_incident_poly
  (const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const SCALAR_TYPE isovalue0,
   const SCALAR_TYPE isovalue1,
   const TABLE_INDEX table_index,
   const int ivolv,
   const int ie,
   const int iv0,
   const int iv1,
   const SCALAR_TYPE s0,
   const SCALAR_TYPE s1,
   SCALAR_TYPE & isovalueX)
  {

    if (does_poly_dual_to_grid_vertex_intersect_grid_edge
        (ivoldual_table, isovalue0, isovalue1, table_index, ivolv,
         ie, iv0, iv1, s0, s1, isovalueX))
      { return(true); }
    else if (does_poly_dual_to_grid_vertex_intersect_grid_edge
        (ivoldual_table, isovalue0, isovalue1, table_index, ivolv,
         ie, iv1, iv0, s1, s0, isovalueX))
      { return(true); }
    else if (ivoldual_table.EdgeHasDualIVolPoly(table_index, ie)) {
      if (ivoldual_table.UpperIncident(table_index, ie) == ivolv) {
        isovalueX = isovalue1;
        return(true);
      }
      else if (ivoldual_table.LowerIncident(table_index, ie) == ivolv) {
        isovalueX = isovalue0;
        return(true);
      }
    }

    return(false);
  }

}

void IVOLDUAL::position_dual_ivolv_in_interval_volume_centroid_multi
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const SCALAR_TYPE isovalue0,
 const SCALAR_TYPE isovalue1,
 const DUAL_ISOVERT & ivolv_info,
 const CUBE_FACE_INFO & cube,
 COORD_TYPE * vcoord, 
 COORD_TYPE * temp_coord0, COORD_TYPE * temp_coord1, COORD_TYPE * temp_coord2)
{
  const int dimension = scalar_grid.Dimension();
  const VERTEX_INDEX icube = ivolv_info.cube_index;
  const int ivolv = ivolv_info.patch_index;
  const TABLE_INDEX table_index = ivolv_info.table_index;
  SCALAR_TYPE isovalueX;

  int num_intersected_edges = 0;
  IJK::set_coord(dimension, 0.0, vcoord);

  for (int ie = 0; ie < cube.NumEdges(); ie++) {
    const int k0 = cube.EdgeEndpoint(ie, 0);
    const int k1 = cube.EdgeEndpoint(ie, 1);
    const VERTEX_INDEX iend0 = scalar_grid.CubeVertex(icube, k0);
    const VERTEX_INDEX iend1 = scalar_grid.CubeVertex(icube, k1);
    const SCALAR_TYPE s0 = scalar_grid.Scalar(iend0);
    const SCALAR_TYPE s1 = scalar_grid.Scalar(iend1);

    if (does_grid_edge_intersect_incident_poly
        (ivoldual_table, isovalue0, isovalue1, table_index, ivolv,
         ie, k0, k1, s0, s1, isovalueX)) {

      scalar_grid.ComputeCoord(iend0, temp_coord0);
      scalar_grid.ComputeCoord(iend1, temp_coord1);

      if ((s0 < isovalueX && s1 < isovalueX) || 
          (s0 > isovalueX && s1 > isovalueX)) {
        // Use edge midpoint.
        IJK::linear_interpolate_coord
          (dimension, 0.5, temp_coord0, temp_coord1, temp_coord2);
      }
      else {
        IJK::linear_interpolate_coord
          (dimension, s0, temp_coord0, s1, temp_coord1, isovalueX, temp_coord2);
      }

      IJK::add_coord(dimension, vcoord, temp_coord2, vcoord);
      num_intersected_edges++;
    }
  }

  if (num_intersected_edges > 0) {
    IJK::multiply_coord
      (dimension, 1.0/num_intersected_edges, vcoord, vcoord);
  }
  else {
    scalar_grid.ComputeCubeCenterCoord(icube, vcoord);
  }

}