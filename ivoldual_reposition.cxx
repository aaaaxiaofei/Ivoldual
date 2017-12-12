/// \file ivoldual_reposition.cxx
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

#include "ivoldual_compute.h"
#include "ivoldual_reposition.h"

#include "ijktriangulate.txx"

using namespace IJK;
using namespace IVOLDUAL;

void IVOLDUAL::eliminate_non_manifold_grid
(IVOLDUAL_DATA & ivoldual_data, 
 const SCALAR_TYPE isovalue0, 
 const SCALAR_TYPE isovalue1, 
 IVOLDUAL_INFO & dualiso_info) 
{
  VERTEX_INDEX changes_of_non_manifold = 0;
  ivoldual_data.EliminateAmbigFacets
    (isovalue0, isovalue1, changes_of_non_manifold);
  dualiso_info.non_manifold_changes = changes_of_non_manifold;
}

void IVOLDUAL::laplacian_smooth_elength
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const IVOLDUAL_DATA_FLAGS & param,
 IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 COORD_ARRAY & vertex_coord, 
 int iteration)
{
  const int d = 3;
  float laplacian_smooth_limit = 0.1;
  float dist;
  COORD_TYPE * vcoord = &(vertex_coord.front());

  for (int it = 0; it < 2*iteration+1; it++) {

    bool skipSurfaceVert = (it % 2 == 0);

    // Loop over all vertices
    for (int cur = 0; cur < vertex_adjacency_list.NumVertices(); cur++) {
      
      // Current node coordinates.
      COORD_TYPE *cur_coord = vcoord + cur*d;

      // Check if current node is on isosurface.
      const int ivolv_cur = ivolv_list[cur].patch_index;
      const TABLE_INDEX table_cur = ivolv_list[cur].table_index;
      bool curOnLower = ivoldual_table.OnLowerIsosurface(table_cur, ivolv_cur);
      bool curOnUpper = ivoldual_table.OnUpperIsosurface(table_cur, ivolv_cur);
      bool isOnSurface = curOnLower || curOnUpper;

      if (isOnSurface == skipSurfaceVert) continue;

      // Store sum of neighbor coordinates.
      COORD_TYPE neigh_sum[d]; 
      IJK::set_coord(d, 0.0, neigh_sum);
      bool needMoving = false;
      int adj_count = 0;

      // Loop over adjacent vertices of the current vertex
      for (int  k = 0; k < vertex_adjacency_list.NumAdjacent(cur); k++) {

        // Neighbor node coordinates
        int adj = vertex_adjacency_list.AdjacentVertex(cur, k);
        COORD_TYPE *neigh_coord = vcoord + adj*d;

        // Check if neighbor node is on isosurface.
        const int ivolv_adj = ivolv_list[adj].patch_index;
        const TABLE_INDEX table_adj = ivolv_list[adj].table_index;
        bool adjOnLower = ivoldual_table.OnLowerIsosurface(table_adj, ivolv_adj);
        bool adjOnUpper = ivoldual_table.OnUpperIsosurface(table_adj, ivolv_adj);

        // Skip if a vertex and its adjacent vertex are not on the same surface.
        if ((curOnLower && !adjOnLower) || (curOnUpper && !adjOnUpper))
            continue;

        IJK::add_coord(d, neigh_sum, neigh_coord, neigh_sum);
        IJK::compute_distance(d, cur_coord, neigh_coord, dist);

        // Check if minimum distance is valid.
        if (dist < laplacian_smooth_limit) {
          needMoving = true;
        }

        adj_count++;
      }

      // Update current node coordinate.
      if (needMoving) {
        IJK::divide_coord(d, adj_count, neigh_sum, neigh_sum);
        IJK::copy_coord(d, neigh_sum, cur_coord);
      }
    }
  }
}

void IVOLDUAL::laplacian_smooth_jacobian
(const std::vector<VERTEX_INDEX> & ivolpoly_cube,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 COORD_ARRAY & vertex_coord, 
 int iteration)
 {

  using namespace std;

  const int DIMENSION(3);
  const float step_base(0.02);
  const int NUM_VERT_PER_HEXAHEDRON(8); 
  COORD_TYPE * vcoord = &(vertex_coord.front());
  float jacobian_limit = 0.0;

  // Polytopes dual to vertex.
  IJK::POLYMESH_DATA<VERTEX_INDEX,int, 
    IJK::HEX_TRIANGULATION_INFO<char,char>> hex_data;
  hex_data.AddPolytopes(ivolpoly_cube, NUM_VERT_PER_HEXAHEDRON);
  IJK::VERTEX_POLY_INCIDENCE<int,int> vertex_poly_incidence(hex_data);

  for (int it = 0; it < iteration; it++) {

    std::vector<int> negative_jabocian_list;

    // Find all vertices with negative Jacobian.
    for (int ihex = 0; ihex < ivolpoly_cube.size()/8; ihex++) {
      for (int i = 0; i < 8; i++) {
        // Compute Jacobian at current vertex
        COORD_TYPE jacob;        
        compute_hexahedron_Jacobian_determinant
          (ivolpoly_cube, ihex, vertex_coord, i, jacob);
        if (jacob < jacobian_limit) {
          negative_jabocian_list.push_back(ivolpoly_cube[ihex * 8 + i]);
        }
      }
    }

    for (int cur : negative_jabocian_list) {
      // Current node coordinates.
      COORD_TYPE *cur_coord = vcoord + cur*DIMENSION;

      // Check if current node is on isosurface.
      const int ivolv_cur = ivolv_list[cur].patch_index;
      const TABLE_INDEX table_cur = ivolv_list[cur].table_index;
      const VERTEX_INDEX cube_cur = ivolv_list[cur].cube_index;
      bool curOnLower = ivoldual_table.OnLowerIsosurface(table_cur, ivolv_cur);
      bool curOnUpper = ivoldual_table.OnUpperIsosurface(table_cur, ivolv_cur);

      // Loop over adjacent vertices of the current vertex
      for (int j = 0; j < vertex_adjacency_list.NumAdjacent(cur); j++) {

        // Neighbor node coordinates
        int adj = vertex_adjacency_list.AdjacentVertex(cur, j);
        COORD_TYPE *neigh_coord = vcoord + adj*DIMENSION;
        
        // Check if neighbor node is on isosurface.
        const int ivolv_adj = ivolv_list[adj].patch_index;
        const TABLE_INDEX table_adj = ivolv_list[adj].table_index;
        const VERTEX_INDEX cube_adj = ivolv_list[adj].cube_index;
        bool adjOnLower = ivoldual_table.OnLowerIsosurface(table_adj, ivolv_adj);
        bool adjOnUpper = ivoldual_table.OnUpperIsosurface(table_adj, ivolv_adj);

        // Neighbor vertex is in the same cube
        if (cube_cur == cube_adj) {
          COORD_TYPE target[DIMENSION];
          float pre_jacobian = 0.0;
          for (int d = 0; d < DIMENSION; d++) {
            target[d] = neigh_coord[d];
          }
          for (int k = 0; k < vertex_adjacency_list.NumAdjacent(adj); k++) {

            int adj2 = vertex_adjacency_list.AdjacentVertex(adj, k);
            COORD_TYPE *neigh_coord2 = vcoord + adj2*DIMENSION;

            const int ivolv_adj2 = ivolv_list[adj2].patch_index;
            const TABLE_INDEX table_adj2 = ivolv_list[adj2].table_index;
            const VERTEX_INDEX cube_adj2 = ivolv_list[adj2].cube_index;
            bool adjOnLower2 = ivoldual_table.OnLowerIsosurface(table_adj2, ivolv_adj2);
            bool adjOnUpper2 = ivoldual_table.OnUpperIsosurface(table_adj2, ivolv_adj2);

            // Skip if two vertices are not on same surface.
            if ((adjOnLower && !adjOnLower2) || (adjOnUpper && !adjOnUpper2))
              continue;

            // copy adj coord to a temp_adj_coord
            COORD_TYPE neigh_temp[DIMENSION];

            
            for (int step = 1; step <= 20; step++) {
              // Backup Coord
              for (int d = 0; d < DIMENSION; d++) {
                neigh_temp[d] = neigh_coord[d];
              }

              float step_size = step_base * step;
              for (int d = 0; d < DIMENSION; d++) {
                neigh_coord[d] = (1.0-step_size)*neigh_coord[d] + step_size*neigh_coord2[d];
              }

              COORD_TYPE min_jacobian = 1.0;

              for (int ipoly = 0; ipoly < vertex_poly_incidence.NumIncidentPoly(adj); ipoly++) {
                const int ihex = vertex_poly_incidence.IncidentPoly(adj, ipoly);
                COORD_TYPE min_jacob, max_jacob;
                compute_min_max_hexahedron_Jacobian_determinant
                  (ivolpoly_cube, ihex, vertex_coord, min_jacob, max_jacob);
                min_jacobian = std::min(min_jacob, min_jacobian);
              }

              if (min_jacobian > pre_jacobian) {
                pre_jacobian = min_jacobian;
                for (int d = 0; d < DIMENSION; d++) {
                  target[d] = neigh_coord[d];
                }
              }

              // Restore neighbor coordinate to backup coord
              for (int d = 0; d < DIMENSION; d++) {
                neigh_coord[d] = neigh_temp[d];
              }
            }
          }
          for (int d = 0; d < DIMENSION; d++) {
            neigh_coord[d] = target[d];
          }
        }
      }
    }
  }
}