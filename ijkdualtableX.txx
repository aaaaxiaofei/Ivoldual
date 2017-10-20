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

  template <typename DUAL_TABLE_TYPE, 
            typename DUAL_TABLE_VINFO_TYPE>
  void compute_dual_ivoltable_vertex_degrees
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
        if (table.EdgeHasDualIVolPoly(it, ie)) {
          const NTYPE isov0 = table.LowerIncident(it, ie);
          const NTYPE isov1 = table.UpperIncident(it, ie);
          dual_table_vertex_info.VertexInfoNC(it, isov0).degree++;
          dual_table_vertex_info.VertexInfoNC(it, isov1).degree++;
        }
      }

      for (NTYPE iv = 0; iv < table.NumPolyVertices(); iv++) {
        if (table.IsInIntervalVolume(it, iv)) {
          const NTYPE isov = table.IncidentIVolVertex(it, iv);
          dual_table_vertex_info.VertexInfoNC(it, isov).degree++;
        }
      }
    }

  }

}

#endif

