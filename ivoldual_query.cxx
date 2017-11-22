/// \file ivoldual_query.cxx
/// Query routines for ivoldual.


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


#include "ivoldual_query.h"

// Return true if vertex is lower left vertex of an isosurface box.
bool IVOLDUAL::is_ivolv_lower_left_vertex_of_isosurface_box
(const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 const IVOL_VERTEX_INDEX ivolv0,
 IVOL_VERTEX_INDEX box_corner[8])
{
  typedef IVOLDUAL_TABLE_VERTEX_INFO::CUBE_VERTEX_TYPE CUBE_VERTEX_TYPE;

  const int DIM3(3);
  const int NUM_CUBE_VERTICES(8);
  const int NUM_CUBE_FACET_VERTICES(NUM_CUBE_VERTICES/2);
  const CUBE_VERTEX_TYPE UNDEFINED_CUBE_VERTEX =
    IVOLDUAL_TABLE_VERTEX_INFO::UNDEFINED_CUBE_VERTEX;
  const CUBE_VERTEX_TYPE separation_vertex = 
    ivolv_list[ivolv0].SeparationVertex();

  if (separation_vertex == UNDEFINED_CUBE_VERTEX) 
    { return(false); }

  if (ivolv_list[ivolv0].NumIncidentIsoQuad() != 3) { return(false); }

  box_corner[0] = ivolv0;
  if (!vertex_adjacency_list.GetAdjacentVertexInOrientedDirection
      (ivolv0, 0, 1, box_corner[1])) { return(false); }
  if (!vertex_adjacency_list.GetAdjacentVertexInOrientedDirection
      (ivolv0, 1, 1, box_corner[2])) { return(false); }
  if (!vertex_adjacency_list.GetAdjacentVertexInOrientedDirection
      (box_corner[1], 1, 1, box_corner[3])) { return(false); }

  for (int j = 0; j < NUM_CUBE_FACET_VERTICES; j++) {
    if (!vertex_adjacency_list.GetAdjacentVertexInOrientedDirection
        (box_corner[j], 2, 1, box_corner[4+j])) { return(false); }
  }

  for (int i = 1; i < NUM_CUBE_VERTICES; i++) {
    const VERTEX_INDEX ivolv = box_corner[i];
    if (ivolv_list[ivolv].separation_vertex != separation_vertex)
      { return(false); }
    if (ivolv_list[ivolv].NumIncidentIsoQuad() != 3)
      { return(false); }
  }

  return(true);
}


// Return true if vertex is lower left vertex of an isosurface pseudobox.
bool IVOLDUAL::is_ivolv_lower_left_vertex_of_isosurface_pseudobox
(const DUALISO_GRID & grid,
 const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 const IVOL_VERTEX_INDEX ivolv0,
 IVOL_VERTEX_INDEX box_corner[8])
{
  typedef IVOLDUAL_TABLE_VERTEX_INFO::CUBE_VERTEX_TYPE CUBE_VERTEX_TYPE;

  const int DIM3(3);
  const int NUM_CUBE_VERTICES(8);
  const int NUM_CUBE_FACET_VERTICES(NUM_CUBE_VERTICES/2);
  const CUBE_VERTEX_TYPE UNDEFINED_CUBE_VERTEX =
    IVOLDUAL_TABLE_VERTEX_INFO::UNDEFINED_CUBE_VERTEX;
  const VERTEX_INDEX cube_index = ivolv_list[ivolv0].cube_index;

  // Index of vertex surrounded by pseudobox.
  const CUBE_VERTEX_TYPE grid_vertex = grid.CubeVertex(cube_index, 7);

  /* OBSOLETE
  box_corner[0] = ivolv0;
  if (!vertex_adjacency_list.GetAdjacentVertexInOrientedDirection
      (ivolv0, 0, 1, box_corner[1])) { return(false); }
  if (!vertex_adjacency_list.GetAdjacentVertexInOrientedDirection
      (ivolv0, 1, 1, box_corner[2])) { return(false); }
  if (!vertex_adjacency_list.GetAdjacentVertexInOrientedDirection
      (box_corner[1], 1, 1, box_corner[3])) { return(false); }

  for (int j = 0; j < NUM_CUBE_FACET_VERTICES; j++) {
    if (!vertex_adjacency_list.GetAdjacentVertexInOrientedDirection
        (box_corner[j], 2, 1, box_corner[4+j])) { return(false); }
  }

  for (int i = 1; i < NUM_CUBE_VERTICES; i++) {
    const VERTEX_INDEX ivolv = box_corner[i];
    if (ivolv_list[ivolv].separation_vertex != separation_vertex)
      { return(false); }
    if (ivolv_list[ivolv].NumIncidentIsoQuad() != 3)
      { return(false); }
  }
  */

  return(true);
}

