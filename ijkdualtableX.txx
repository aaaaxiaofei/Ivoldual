/// \file ijkdualtableX.txx
/// Extended dual lookup tables of isosurface and interval volume vertices.
/// Version 0.2.0

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

#ifndef _IJKDUALTABLE_X_
#define _IJKDUALTABLE_X_

#include "ijkdualtable.txx"
#include "ijklist.txx"

#include <vector>

namespace IJKDUALTABLE {

  // **************************************************
  // DUAL TABLE VERTEX INFORMATION
  // **************************************************

  /// Information about each of the isosurface/interval volume vertices
  ///   for each table entry.
  template <typename NTYPE, typename VINFO_TYPE>
  class DUAL_TABLE_VERTEX_INFO:
    public IJK::LIST_OF_LISTS<VINFO_TYPE,NTYPE> {

  public:
    typedef VINFO_TYPE VERTEX_INFO_TYPE;
    typedef NTYPE NUMBER_TYPE;

  public:
    DUAL_TABLE_VERTEX_INFO() {};
    template <typename DUAL_TABLE_TYPE>
    DUAL_TABLE_VERTEX_INFO(const DUAL_TABLE_TYPE & table);

    template <typename TI_TYPE>
    NTYPE NumVertices(const TI_TYPE ientry) const
    { return(this->ListLength(ientry)); };

    template <typename DUAL_TABLE_TYPE>
    void Set(const DUAL_TABLE_TYPE & table);

    /// Return reference to vertex info.
    template <typename TI_TYPE>
    const VINFO_TYPE & VertexInfo
    (const TI_TYPE ientry, const NTYPE isov) const
    { return(this->ElementRefConst(ientry, isov)); }

    /// Return non-constant reference to vertex info.
    template <typename TI_TYPE>
    VINFO_TYPE & VertexInfoNC
    (const TI_TYPE ientry, const NTYPE isov)
    { return(this->ElementRef(ientry, isov)); }
  };


  // **************************************************
  // CLASS DUAL_TABLE_VERTEX_INFO MEMBER FUNCTIONS
  // **************************************************

  template <typename NTYPE, typename VINFO_TYPE>
  template <typename DUAL_TABLE_TYPE>
  DUAL_TABLE_VERTEX_INFO<NTYPE,VINFO_TYPE>::
  DUAL_TABLE_VERTEX_INFO(const DUAL_TABLE_TYPE & table)
  {
    Set(table);
  }

  template <typename NTYPE, typename VINFO_TYPE>
  template <typename DUAL_TABLE_TYPE>
  void DUAL_TABLE_VERTEX_INFO<NTYPE,VINFO_TYPE>::
  Set(const DUAL_TABLE_TYPE & table)
  {
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    std::vector<NTYPE> num_isov(table.NumTableEntries());

    for (TABLE_INDEX i = 0; i < table.NumTableEntries(); i++) 
      { num_isov[i] = table.Entry(i).NumVertices(); }
      
    this->CreateLists(num_isov);
  }


  // **************************************************
  // COMPUTE VERTEX DEGREE
  // **************************************************

  /// Compute number of isosurface polytopes incident on each vertex.
  template <typename DUAL_TABLE_TYPE, 
            typename DUAL_TABLE_VINFO_TYPE>
  void compute_dual_isotable_vertex_degrees
  (const DUAL_TABLE_TYPE & table,
   DUAL_TABLE_VINFO_TYPE & dual_table_vertex_info)
  {
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_TABLE_VINFO_TYPE::NUMBER_TYPE NTYPE;

    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {

      // Initialize
      for (NTYPE j = 0; j < table.Entry(it).NumVertices(); j++) 
        { dual_table_vertex_info.VertexInfoNC(it, j).degree = 0; }
        

      for (NTYPE ie = 0; ie < table.NumPolyEdges(); ie++) {
        if (table.IsBipolar(it, ie)) {
          const NTYPE isov = table.IncidentIsoVertex(it, ie);
          dual_table_vertex_info.VertexInfoNC(it, isov).degree++;
        }
      }
    }

  }


  /// Compute number of interval volume polytopes incident on each vertex.
  template <typename DUAL_TABLE_TYPE, 
            typename DUAL_TABLE_VINFO_TYPE>
  void compute_ivoldual_table_num_incident_poly
  (const DUAL_TABLE_TYPE & table,
   DUAL_TABLE_VINFO_TYPE & dual_table_vertex_info)
  {
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_TABLE_VINFO_TYPE::NUMBER_TYPE NTYPE;

    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {

      // Initialize
      for (NTYPE j = 0; j < table.Entry(it).NumVertices(); j++) 
        { dual_table_vertex_info.VertexInfoNC(it, j).num_incident_poly = 0; }
        

      for (NTYPE ie = 0; ie < table.NumPolyEdges(); ie++) {
        if (table.EdgeHasDualIVolPoly(it, ie)) {
          const NTYPE isov0 = table.LowerIncident(it, ie);
          const NTYPE isov1 = table.UpperIncident(it, ie);
          dual_table_vertex_info.VertexInfoNC(it, isov0).num_incident_poly++;
          dual_table_vertex_info.VertexInfoNC(it, isov1).num_incident_poly++;
        }
      }

      for (NTYPE iv = 0; iv < table.NumPolyVertices(); iv++) {
        if (table.IsInIntervalVolume(it, iv)) {
          const NTYPE isov = table.IncidentIVolVertex(it, iv);
          dual_table_vertex_info.VertexInfoNC(it, isov).num_incident_poly++;
        }
      }
    }

  }


  /// Compute number of isosurface polytopes incident on each vertex.
  /// - Note: An isosurface polytope is a facet of an interval volume polytope.
  template <typename DUAL_TABLE_TYPE, 
            typename DUAL_TABLE_VINFO_TYPE>
  void compute_ivoldual_table_num_incident_isosurface_poly
  (const DUAL_TABLE_TYPE & table,
   DUAL_TABLE_VINFO_TYPE & vinfo)
  {
    typedef typename DUAL_TABLE_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_TABLE_VINFO_TYPE::NUMBER_TYPE NTYPE;

    NTYPE iend[2];

    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {

      // Initialize
      for (NTYPE j = 0; j < table.Entry(it).NumVertices(); j++) 
        { vinfo.VertexInfoNC(it, j).num_incident_isopoly = 0; }

      for (NTYPE ie = 0; ie < table.NumPolyEdges(); ie++) {

        iend[0] = table.Cube().EdgeEndpoint(ie, 0);
        iend[1] = table.Cube().EdgeEndpoint(ie, 1);

        if (table.EdgeHasDualIVolPoly(it, ie)) {
          const NTYPE ivolv0 = table.LowerIncident(it, ie);
          const NTYPE ivolv1 = table.UpperIncident(it, ie);

          if (table.IsBelowIntervalVolume(it, iend[0]) ||
              table.IsBelowIntervalVolume(it, iend[1])) {

            if (table.OnLowerIsosurface(it, ivolv0))
              { vinfo.VertexInfoNC(it, ivolv0).num_incident_isopoly++; }
            if (table.OnLowerIsosurface(it, ivolv1))
              { vinfo.VertexInfoNC(it, ivolv1).num_incident_isopoly++; }
          }

          if (table.IsAboveIntervalVolume(it, iend[0]) ||
              table.IsAboveIntervalVolume(it, iend[1])) {

            if (table.OnUpperIsosurface(it, ivolv0))
              { vinfo.VertexInfoNC(it, ivolv0).num_incident_isopoly++; }
            if (table.OnUpperIsosurface(it, ivolv1))
              { vinfo.VertexInfoNC(it, ivolv1).num_incident_isopoly++; }
          }
        }
        else if (table.IsInIntervalVolume(it, iend[0]) ||
                 table.IsInIntervalVolume(it, iend[1])) {

          if (table.IsInIntervalVolume(it, iend[0]))
            { std::swap(iend[0], iend[1]); }

          if (table.IsInIntervalVolume(it, iend[0])) { continue; }

          const NTYPE ivolv = table.IncidentIVolVertex(it, iend[1]);

          if (table.IsBelowIntervalVolume(it, iend[0])) {
            if (table.OnLowerIsosurface(it, ivolv))
              { vinfo.VertexInfoNC(it, ivolv).num_incident_isopoly++; }
          }

          if (table.IsAboveIntervalVolume(it, iend[0])) {
            if (table.OnUpperIsosurface(it, ivolv))
              { vinfo.VertexInfoNC(it, ivolv).num_incident_isopoly++; }
          }
        }
      }

    }
  }


  /// Determine separation cube vertices which are below/above 
  ///   the interval volume.
  /// @pre Call compute_ivoldual_table_num_incident_poly
  ///   and compute_ivoldual_table_num_incident_isosurface_poly
  ///   before calling this routine.
  template <typename DUAL_TABLE_TYPE, 
            typename DUAL_TABLE_VINFO_TYPE>
  void determine_separation_vertices_below_above_ivol
  (const DUAL_TABLE_TYPE & table,
   DUAL_TABLE_VINFO_TYPE & vinfo)
  {
    typedef typename DUAL_TABLE_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_TABLE_VINFO_TYPE::NUMBER_TYPE NTYPE;
    typedef typename DUAL_TABLE_VINFO_TYPE::VERTEX_INFO_TYPE::CUBE_VERTEX_TYPE
      CUBE_VERTEX_TYPE;

    const CUBE_VERTEX_TYPE UNDEFINED_CUBE_VERTEX =
      DUAL_TABLE_VINFO_TYPE::VERTEX_INFO_TYPE::UNDEFINED_CUBE_VERTEX;
    const NTYPE CUBE_VERTEX_DEGREE(table.Dimension());
    NTYPE ivolv[2], iend[2];
    std::vector<CUBE_VERTEX_TYPE> candidate_separation_vertex;
    std::vector<NTYPE> num_incident_edges_dual_to_isosurface;
    IJK::PROCEDURE_ERROR error ("determine_separation_vertices");

    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {

      const NTYPE num_ivolv = table.Entry(it).NumVertices();

      candidate_separation_vertex.resize(num_ivolv);
      num_incident_edges_dual_to_isosurface.resize(num_ivolv);

      // Initialize
      for (NTYPE j = 0; j < num_ivolv; j++) {
        candidate_separation_vertex[j] = UNDEFINED_CUBE_VERTEX;
        num_incident_edges_dual_to_isosurface[j] = 0;
      }

      for (NTYPE ie = 0; ie < table.NumPolyEdges(); ie++) {
        if (table.EdgeHasDualIVolPoly(it, ie)) {
          ivolv[0] = table.LowerIncident(it, ie);
          ivolv[1] = table.UpperIncident(it, ie);

          iend[0] = table.Cube().EdgeEndpoint(ie, 0);
          iend[1] = table.Cube().EdgeEndpoint(ie, 1);

          for (NTYPE i = 0; i < 2; i++) {
            for (NTYPE j = 0; j < 2; j++) {

              const NTYPE iend_i = iend[i];
              const NTYPE ivolv_j = ivolv[j];
              if (vinfo.VertexInfo(it, ivolv_j).num_incident_isopoly != 
                  CUBE_VERTEX_DEGREE)
                { continue; }

              if (table.OnLowerIsosurface(it, ivolv_j) &&
                  table.IsBelowIntervalVolume(it,iend_i)) {
                if (candidate_separation_vertex[ivolv_j] == 
                    UNDEFINED_CUBE_VERTEX) {
                  candidate_separation_vertex[ivolv_j] = iend_i;
                  num_incident_edges_dual_to_isosurface[ivolv_j] = 1;
                }
                else if (candidate_separation_vertex[ivolv_j] == iend_i) {
                  num_incident_edges_dual_to_isosurface[ivolv_j]++;
                }
              }

              if (table.OnUpperIsosurface(it, ivolv_j) &&
                  table.IsAboveIntervalVolume(it,iend_i)) {
                if (candidate_separation_vertex[ivolv_j] == 
                    UNDEFINED_CUBE_VERTEX) {
                  candidate_separation_vertex[ivolv_j] = iend_i;
                  num_incident_edges_dual_to_isosurface[ivolv_j] = 1;
                }
                else if (candidate_separation_vertex[ivolv_j] == iend_i) {
                  num_incident_edges_dual_to_isosurface[ivolv_j]++;
                }
              }
            }
          }
        }
      }

      for (NTYPE iv0 = 0; iv0 < table.NumPolyVertices(); iv0++) {
        if (table.IsInIntervalVolume(it, iv0)) {
          const NTYPE ivolv = table.IncidentIVolVertex(it, iv0);

          for (DTYPE d = 0; d < table.Dimension(); d++) {

            const NTYPE iv1 = table.Cube().VertexNeighbor(iv0, d);
            if (table.OnLowerIsosurface(it, ivolv) &&
                table.IsBelowIntervalVolume(it,iv1)) {
              if (candidate_separation_vertex[ivolv] == 
                  UNDEFINED_CUBE_VERTEX) {
                candidate_separation_vertex[ivolv] = iv1;
                num_incident_edges_dual_to_isosurface[ivolv] = 1;
              }
              else if (candidate_separation_vertex[ivolv] == iv1) {
                num_incident_edges_dual_to_isosurface[ivolv]++;
              }
            }

            if (table.OnUpperIsosurface(it, ivolv) &&
                table.IsAboveIntervalVolume(it,iv1)) {
              if (candidate_separation_vertex[ivolv] == 
                  UNDEFINED_CUBE_VERTEX) {
                candidate_separation_vertex[ivolv] = iv1;
                num_incident_edges_dual_to_isosurface[ivolv] = 1;
              }
              else if (candidate_separation_vertex[ivolv] == iv1) {
                num_incident_edges_dual_to_isosurface[ivolv]++;
              }
            }
          }
        }
      }

      for (NTYPE j = 0; j < num_ivolv; j++) {
        if (candidate_separation_vertex[j] != UNDEFINED_CUBE_VERTEX &&
            num_incident_edges_dual_to_isosurface[j] == CUBE_VERTEX_DEGREE) {
          vinfo.VertexInfoNC(it,j).separation_vertex = 
            candidate_separation_vertex[j];
        }
      }

      // Set separation vertices which are separated from other vertices
      //   by an isosurface below the interval volume and an isosurface
      //   above the interval volume.
      for (NTYPE j = 0; j < num_ivolv; j++) {
        candidate_separation_vertex[j] = UNDEFINED_CUBE_VERTEX;
        num_incident_edges_dual_to_isosurface[j] = 0;
      }

      for (NTYPE ie = 0; ie < table.NumPolyEdges(); ie++) {
        if (table.EdgeHasDualIVolPoly(it, ie)) {
          ivolv[0] = table.LowerIncident(it, ie);
          ivolv[1] = table.UpperIncident(it, ie);

          iend[0] = table.Cube().EdgeEndpoint(ie, 0);
          iend[1] = table.Cube().EdgeEndpoint(ie, 1);

          if (table.IsInIntervalVolume(it, iend[0]) ||
              table.IsInIntervalVolume(it, iend[1])) { continue; }

          for (NTYPE i = 0; i < 2; i++) {
            for (NTYPE j = 0; j < 2; j++) {

              const NTYPE iendA = iend[i];
              const NTYPE ivolvC = ivolv[j];
              const NTYPE ivolvD = ivolv[1-j];
              if (vinfo.VertexInfo(it, ivolvC).num_incident_isopoly != 
                  CUBE_VERTEX_DEGREE)
                { continue; }

              if (vinfo.VertexInfo(it, ivolvD).separation_vertex != iendA)
                { continue; }

              if (table.OnLowerIsosurface(it, ivolvC) &&
                  table.IsAboveIntervalVolume(it,iendA)) {
                if (candidate_separation_vertex[ivolvC] == 
                    UNDEFINED_CUBE_VERTEX) {
                  candidate_separation_vertex[ivolvC] = iendA;
                  num_incident_edges_dual_to_isosurface[ivolvC] = 1;
                }
                else if (candidate_separation_vertex[ivolvC] == iendA) {
                  num_incident_edges_dual_to_isosurface[ivolvC]++;
                }
              }

              if (table.OnUpperIsosurface(it, ivolvC) &&
                  table.IsBelowIntervalVolume(it,iendA)) {
                if (candidate_separation_vertex[ivolvC] == 
                    UNDEFINED_CUBE_VERTEX) {
                  candidate_separation_vertex[ivolvC] = iendA;
                  num_incident_edges_dual_to_isosurface[ivolvC] = 1;
                }
                else if (candidate_separation_vertex[ivolvC] == iendA) {
                  num_incident_edges_dual_to_isosurface[ivolvC]++;
                }
              }
            }
          }
        }
      }

      for (NTYPE j = 0; j < num_ivolv; j++) {
        if (candidate_separation_vertex[j] != UNDEFINED_CUBE_VERTEX &&
            num_incident_edges_dual_to_isosurface[j] == CUBE_VERTEX_DEGREE) {
          vinfo.VertexInfoNC(it,j).separation_vertex = 
            candidate_separation_vertex[j];
        }
      }

    }
  }


  /// Determine separation cube vertices which are in the interval volume.
  /// @pre Call compute_ivoldual_table_num_incident_poly
  ///   and compute_ivoldual_table_num_incident_isosurface_poly
  ///   before calling this routine.
  template <typename DUAL_TABLE_TYPE, 
            typename DUAL_TABLE_VINFO_TYPE>
  void determine_separation_vertices_in_ivol
  (const DUAL_TABLE_TYPE & table,
   DUAL_TABLE_VINFO_TYPE & vinfo)
  {
    typedef typename DUAL_TABLE_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_TABLE_VINFO_TYPE::NUMBER_TYPE NTYPE;
    typedef typename DUAL_TABLE_VINFO_TYPE::VERTEX_INFO_TYPE::CUBE_VERTEX_TYPE
      CUBE_VERTEX_TYPE;

    const CUBE_VERTEX_TYPE UNDEFINED_CUBE_VERTEX =
      DUAL_TABLE_VINFO_TYPE::VERTEX_INFO_TYPE::UNDEFINED_CUBE_VERTEX;
    const NTYPE CUBE_VERTEX_DEGREE(table.Dimension());
    NTYPE ivolv[2], iend[2];
    std::vector<CUBE_VERTEX_TYPE> candidate_separation_vertex;
    std::vector<NTYPE> num_incident_edges_dual_to_isosurface;
    IJK::PROCEDURE_ERROR error ("determine_separation_vertices");

    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {

      const NTYPE num_ivolv = table.Entry(it).NumVertices();

      candidate_separation_vertex.resize(num_ivolv);
      num_incident_edges_dual_to_isosurface.resize(num_ivolv);

      // Initialize
      for (NTYPE j = 0; j < num_ivolv; j++) {
        candidate_separation_vertex[j] = UNDEFINED_CUBE_VERTEX;
        num_incident_edges_dual_to_isosurface[j] = 0;
      }

      // Set separation vertices contained in a single isolated
      //   interval volume polyhedron.
      for (NTYPE iv0 = 0; iv0 < table.NumPolyVertices(); iv0++) {
        if (table.IsInIntervalVolume(it, iv0)) {
          const NTYPE ivolv = table.IncidentIVolVertex(it, iv0);

          if (vinfo.VertexInfo(it,ivolv).num_incident_poly == 1) {
            vinfo.VertexInfoNC(it,ivolv).separation_vertex = iv0;
          }
        }
      }

      for (NTYPE j = 0; j < num_ivolv; j++) {
        candidate_separation_vertex[j] = UNDEFINED_CUBE_VERTEX;
        num_incident_edges_dual_to_isosurface[j] = 0;
      }

      // Set separation vertices contained in an interval volume polyhedron
      //   surrounded by a "layer" of polyhedra dual to cube edges.
      for (NTYPE ie = 0; ie < table.NumPolyEdges(); ie++) {
        if (table.EdgeHasDualIVolPoly(it, ie)) {
          ivolv[0] = table.LowerIncident(it, ie);
          ivolv[1] = table.UpperIncident(it, ie);

          iend[0] = table.Cube().EdgeEndpoint(ie, 0);
          iend[1] = table.Cube().EdgeEndpoint(ie, 1);

          for (NTYPE i = 0; i < 2; i++) {
            for (NTYPE j = 0; j < 2; j++) {

              const NTYPE iendA = iend[i];
              const NTYPE iendB = iend[1-i];
              const NTYPE ivolv_j = ivolv[j];
              if (vinfo.VertexInfo(it, ivolv_j).num_incident_isopoly != 
                  CUBE_VERTEX_DEGREE)
                { continue; }

              if (table.IsInIntervalVolume(it,iendB)) {

                if (table.OnLowerIsosurface(it, ivolv_j) &&
                    table.IsBelowIntervalVolume(it,iendA)) {
                  if (candidate_separation_vertex[ivolv_j] == 
                      UNDEFINED_CUBE_VERTEX) {
                    candidate_separation_vertex[ivolv_j] = iendB;
                    num_incident_edges_dual_to_isosurface[ivolv_j] = 1;
                  }
                  else if (candidate_separation_vertex[ivolv_j] == iendB) {
                    num_incident_edges_dual_to_isosurface[ivolv_j]++;
                  }
                }

                if (table.OnUpperIsosurface(it, ivolv_j) &&
                    table.IsAboveIntervalVolume(it,iendA)) {
                  if (candidate_separation_vertex[ivolv_j] == 
                      UNDEFINED_CUBE_VERTEX) {
                    candidate_separation_vertex[ivolv_j] = iendB;
                    num_incident_edges_dual_to_isosurface[ivolv_j] = 1;
                  }
                  else if (candidate_separation_vertex[ivolv_j] == iendB) {
                    num_incident_edges_dual_to_isosurface[ivolv_j]++;
                  }
                }
              }
            }
          }
        }
      }

      for (NTYPE j = 0; j < num_ivolv; j++) {
        if (candidate_separation_vertex[j] != UNDEFINED_CUBE_VERTEX &&
            num_incident_edges_dual_to_isosurface[j] == CUBE_VERTEX_DEGREE) {

          const NTYPE sep_vert = candidate_separation_vertex[j];
          vinfo.VertexInfoNC(it,j).separation_vertex = sep_vert;

          // Also set candidate_separation_vertex[j] to be separation vertex
          //   of the incident interval volume vertex associated with j.
          const NTYPE ivolv = table.IncidentIVolVertex(it, sep_vert);
          vinfo.VertexInfoNC(it,ivolv).separation_vertex = sep_vert;
        }
      }

    }

  }


  /// Determine separation cube vertex for ivol vertex incident
  ///   on 3 isosurface polytopes.
  /// @pre Call compute_ivoldual_table_num_incident_poly
  ///   and compute_ivoldual_table_num_incident_isosurface_poly
  ///   before calling this routine.
  template <typename DUAL_TABLE_TYPE, 
            typename DUAL_TABLE_VINFO_TYPE>
  void determine_separation_vertices
  (const DUAL_TABLE_TYPE & table,
   DUAL_TABLE_VINFO_TYPE & vinfo)
  {
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_TABLE_VINFO_TYPE::NUMBER_TYPE NTYPE;
    typedef typename DUAL_TABLE_VINFO_TYPE::VERTEX_INFO_TYPE::CUBE_VERTEX_TYPE
      CUBE_VERTEX_TYPE;

    const CUBE_VERTEX_TYPE UNDEFINED_CUBE_VERTEX =
      DUAL_TABLE_VINFO_TYPE::VERTEX_INFO_TYPE::UNDEFINED_CUBE_VERTEX;

    // Initialize
    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {
      const NTYPE num_ivolv = table.Entry(it).NumVertices();

      for (NTYPE j = 0; j < num_ivolv; j++) {
        vinfo.VertexInfoNC(it, j).separation_vertex = UNDEFINED_CUBE_VERTEX;
      }
    }

    determine_separation_vertices_below_above_ivol(table, vinfo);
    determine_separation_vertices_in_ivol(table, vinfo);
  }


  /// Determine interval volume vertices which are doubly connected
  ///   to some adjacent cube.
  template <typename DUAL_TABLE_TYPE, 
            typename DUAL_TABLE_VINFO_TYPE>
  void determine_doubly_connected_ivol3D_vertices
  (const DUAL_TABLE_TYPE & table,
   DUAL_TABLE_VINFO_TYPE & vinfo)
  {
    typedef typename DUAL_TABLE_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_TABLE_VINFO_TYPE::NUMBER_TYPE NTYPE;
    typedef typename DUAL_TABLE_VINFO_TYPE::VERTEX_INFO_TYPE::CUBE_VERTEX_TYPE
      CUBE_VERTEX_TYPE;

    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {
      const NTYPE num_ivolv = table.Entry(it).NumVertices();

      for (NTYPE j = 0; j < num_ivolv; j++) {
        NTYPE dc_facet;
        if (is_ivol3D_vertex_doubly_connected(table, it, j, dc_facet)) {
          vinfo.VertexInfoNC(it,j).is_doubly_connected = true; 
          vinfo.VertexInfoNC(it,j).doubly_connected_facet = dc_facet;
        }
        else {
          vinfo.VertexInfoNC(it,j).is_doubly_connected = false; 
          vinfo.VertexInfoNC(it,j).doubly_connected_facet = 0;
        }
      }
    }
  }

}

#endif

