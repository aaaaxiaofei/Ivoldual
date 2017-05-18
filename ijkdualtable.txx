/// \file ijkdualtable.txx
/// Class containing dual lookup table of isosurface vertices.
/// Version 0.1.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2012-2017 Rephael Wenger

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

#ifndef _IJKDUALTABLE_
#define _IJKDUALTABLE_

#include <climits>
#include <limits>

#include "ijk.txx"
#include "ijkcube.txx"

/// Classes and routines for storing and manipulating 
///   dual isosurface lookup table.
namespace IJKDUALTABLE {

  // **************************************************
  // COMPUTE FUNCTIONS
  // **************************************************

  /// Compute complement index.
  template <typename ITYPE, typename NTYPE>
  inline ITYPE compute_complement
  (const ITYPE ival, const NTYPE num_table_entries)
  { return(num_table_entries-1-ival); }
    

  // **************************************************
  // ISODUAL CUBE TABLE ROUTINES
  // **************************************************

  /// Create isodual cube table entry ientry.
  /// @tparam TI_TYPE Table index type.
  template <typename TI_TYPE, typename NTYPE, typename FIND_TYPE,
            typename CUBE_TYPE, typename ENTRY_TYPE>
  void create_isodual_cube_table_entry
  (const TI_TYPE ientry, const NTYPE num_poly_vertices,
   const bool flag_separate_neg, const bool flag_separate_opposite,
   const CUBE_TYPE & cube,
   FIND_TYPE & find_component, 
   ENTRY_TYPE & table_entry)
  {
    const bool flag_separate_pos = (!flag_separate_neg);

    find_component.ClearAll();
    find_component.SetVertexFlags(ientry);

    if (flag_separate_opposite) {
      if (flag_separate_neg) {
        if (IJK::is_two_opposite_ones(ientry, num_poly_vertices))
          { find_component.NegateVertexFlags(); };
      }
      else {
        if (IJK::is_two_opposite_zeros(ientry, num_poly_vertices))
          { find_component.NegateVertexFlags(); };
      }
    }


    NTYPE num_components(0);
    for (NTYPE i = 0; i < num_poly_vertices; i++) {
      if (find_component.Component(i) == 0) {

        if (find_component.VertexFlag(i) == flag_separate_pos) {
          num_components++;
          find_component.Search(i, num_components);
        }
      }
    }

    NTYPE num_zeros, num_ones;
    IJK::count_bits(ientry, num_poly_vertices, num_zeros, num_ones);
    if (num_zeros == 0 || num_ones == 0) {
      // All vertices are negative or all vertices are positive.
      table_entry.num_vertices = 0;
    }
    else {
      table_entry.num_vertices = num_components;
    }

    for (NTYPE ie = 0; ie < cube.NumEdges(); ie++) {
      NTYPE iv0 = cube.EdgeEndpoint(ie, 0);
      NTYPE iv1 = cube.EdgeEndpoint(ie, 1);

      if (find_component.VertexFlag(iv0) == find_component.VertexFlag(iv1)) {
        table_entry.is_bipolar[ie] = false;
        table_entry.incident_isovertex[ie] = 0;
      }
      else {
        table_entry.is_bipolar[ie] = true;
        int icomp = find_component.Component(iv0);
        if (find_component.VertexFlag(iv1) == flag_separate_pos) {
          // Vertex iv1 is negative.
          icomp = find_component.Component(iv1);
        }
        table_entry.incident_isovertex[ie] = icomp-1;
      }
    }

  }


  // **************************************************
  // AMBIGUITY ROUTINES
  // **************************************************

  template <typename TI_TYPE, typename FIND_TYPE>
  bool is_cube_ambiguous
  (const TI_TYPE ientry, FIND_TYPE & find_component)
  {
    typedef typename FIND_TYPE::NUMBER_TYPE NTYPE;

    NTYPE num_pos_components = 
      find_component.ComputeNumComponents(ientry, true);
    NTYPE num_neg_components = 
      find_component.ComputeNumComponents(ientry, false);

    if (num_pos_components > 1 || num_neg_components > 1)
      { return(true); }
    else
      { return(false); }
  }

  template <typename TI_TYPE, typename FTYPE, typename FIND_TYPE>
  bool is_cube_facet_ambiguous
  (const TI_TYPE ientry, const FTYPE kf, 
   FIND_TYPE & find_component)
  {
    typedef typename FIND_TYPE::NUMBER_TYPE NTYPE;

    NTYPE num_pos_components = 
      find_component.ComputeNumComponentsInFacet(ientry, kf, true);
    NTYPE num_neg_components = 
      find_component.ComputeNumComponentsInFacet(ientry, kf, false);

    if (num_pos_components > 1 || num_neg_components > 1)
      { return(true); }
    else
      { return(false); }
  }

  template <typename TI_TYPE, typename NTYPE, typename NTYPE2,
            typename FSET_TYPE, typename FIND_TYPE>
  void compute_ambiguous_cube_facets
  (const TI_TYPE ientry, 
   const NTYPE num_facets,
   FSET_TYPE & facet_set,  
   NTYPE2 & num_ambiguous_facets,
   FIND_TYPE & find_component)
  {
    facet_set = 0;
    num_ambiguous_facets = 0;
    FSET_TYPE mask = FSET_TYPE(1);

    for (NTYPE kf = 0; kf < num_facets; kf++) {
      if (is_cube_facet_ambiguous(ientry, kf, find_component)) {
        num_ambiguous_facets++;
        facet_set = (facet_set | mask);
      }
      mask = (mask << FSET_TYPE(1));
    }
  }


  // **************************************************
  // IVOLDUAL CUBE TABLE ROUTINES
  // **************************************************

  /// Compute lifted ivol indices.
  template <typename TI_TYPE, typename NTYPE, 
            typename TI_TYPE_L, typename TI_TYPE_U>
  void compute_lifted_ivol_indices4
  (const TI_TYPE ientry, const NTYPE num_cube_vertices, 
   TI_TYPE_L & ilower, TI_TYPE_U & iupper)
  {
    const int NUM_CUBE_VERTEX_TYPES = 4;
    TI_TYPE factorA, factorB;

    TI_TYPE x = ientry;
    factorA = 1;
    factorB = (TI_TYPE(1) << num_cube_vertices);
    ilower = 0;
    iupper = 0;

    while (x != 0) {
      TI_TYPE vertex_type = x%NUM_CUBE_VERTEX_TYPES;
      if (vertex_type == 1) {
        iupper += factorB;
      }
      else if (vertex_type == 2) {
        iupper += factorA;
        iupper += factorB;
        ilower += factorB;
      }
      else if (vertex_type == 3) {
        iupper += factorA;
        iupper += factorB;
        ilower += factorA;
        ilower += factorB;
      }

      factorA = (factorA << 1);
      factorB = (factorB << 1);

      x = x/NUM_CUBE_VERTEX_TYPES;
    }
  }


  /// Create ivoldual cube table entry ientry.
  /// @tparam TI_TYPE Table index type.
  template <typename TI_TYPE, typename NTYPE, typename FIND_TYPE,
            typename CUBE_TYPE, typename ISODUAL_ENTRY_TYPE,
            typename ENTRY_TYPE>
  void create_ivoldual_cube_table_entry
  (const TI_TYPE ientry, const NTYPE num_poly_vertices,
   const bool flag_separate_neg, const bool flag_separate_opposite,
   const CUBE_TYPE & cube,
   const CUBE_TYPE & lifted_cube,
   FIND_TYPE & find_component,
   ISODUAL_ENTRY_TYPE & lower_isodual,
   ISODUAL_ENTRY_TYPE & upper_isodual,
   ENTRY_TYPE & table_entry)
  {
    typedef typename CUBE_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE lifted_dimension = lifted_cube.Dimension();
    const NTYPE num_lifted_poly_vertices = lifted_cube.NumVertices();
    TI_TYPE ilower, iupper;

    compute_lifted_ivol_indices4(ientry, num_poly_vertices, ilower, iupper);
    create_isodual_cube_table_entry
      (ilower, num_lifted_poly_vertices, 
       flag_separate_neg, flag_separate_opposite,
       lifted_cube, find_component, lower_isodual);
    create_isodual_cube_table_entry
      (iupper, num_lifted_poly_vertices, 
       flag_separate_neg, flag_separate_opposite,
       lifted_cube, find_component, upper_isodual);

    table_entry.num_vertices = 
      lower_isodual.NumVertices() + upper_isodual.NumVertices();

    for (NTYPE ie = 0; ie < lifted_cube.NumEdges(); ie++) {
      const DTYPE edge_dir = lifted_cube.EdgeDir(ie);
      NTYPE iend0 = lifted_cube.EdgeEndpoint(ie,0);
      NTYPE iend1 = lifted_cube.EdgeEndpoint(ie,1);
      if (iend0 > iend1) { std::swap(iend0, iend1); }

      if (edge_dir+1 == lifted_dimension) {
        if (lower_isodual.IsBipolar(ie)) {
          const NTYPE iv = lower_isodual.IncidentIsoVertex(ie);
          table_entry.vertex_info[iend0].SetIncident(iv);
        }
        else if (upper_isodual.IsBipolar(ie)) {
          const NTYPE iv = upper_isodual.IncidentIsoVertex(ie);
          table_entry.vertex_info[iend0].SetIncident(iv);
        }
      }
      else {
        if (iend1 < cube.NumVertices()) {
          if (upper_isodual.IsBipolar(ie)) {
            const NTYPE je = cube.IncidentEdge(iend0, edge_dir);
            const NTYPE iv = lower_isodual.IncidentIsoVertex(ie);
            table_entry.edge_info[je].SetUpperIncident(iv);
          }
        }
        else {
          if (lower_isodual.IsBipolar(ie)) {
            const NTYPE jend0 = iend0 - cube.NumVertices();
            const NTYPE je = cube.IncidentEdge(jend0, edge_dir);
            const NTYPE iv = lower_isodual.IncidentIsoVertex(ie);
            table_entry.edge_info[je].SetLowerIncident(iv);
          }
        }
      }
    }
  }


  // **************************************************
  // UTILITY FUNCTIONS
  // **************************************************

  /// Calculate number of entries required in ISOSURFACE_TABLE.
  template <typename NUM_ENTRIES_TYPE,
            typename NTYPEV, typename NTYPEC>
  NUM_ENTRIES_TYPE calculate_num_entries
  (const NTYPEV num_vert, const NTYPEC num_colors)
  {
    const char * procname = "calculate_num_entries";

    NUM_ENTRIES_TYPE num_table_entries = 0;

    if (num_colors < 1)
      throw IJK::PROCEDURE_ERROR
        (procname, "Number of colors must be positive.");

    const NUM_ENTRIES_TYPE max2 = 
      std::numeric_limits<NUM_ENTRIES_TYPE>::max()/num_colors;

    num_table_entries = 1;
    for (NTYPEV iv = 0; iv < num_vert; iv++) {
      if (num_table_entries > max2)
        throw IJK::PROCEDURE_ERROR
          (procname, "Number of entries is too large.");

      num_table_entries = num_table_entries * num_colors;
    };

    return(num_table_entries);
  }


  /// Convert integer to boolean flags.
  /// @param[out] flag[] Array of boolean flags.
  /// @pre flag[] is pre-allocated to size at least num_flags.
  template <typename TI_TYPE, typename NTYPE>
  void convert2bool
  (const TI_TYPE ival, bool flag[], const NTYPE num_flags)
  {
    const TI_TYPE base = 2;

    TI_TYPE jval = ival;
    for (NTYPE i = 0; i < num_flags; i++) {
      flag[i] = bool(jval % base);
      jval = jval/base;
    };

    if (jval != 0) {
      IJK::PROCEDURE_ERROR error("convert2bool");
      error.AddMessage("Error converting ", ival, " to base 2.");
      error.AddMessage("Output has more than ", num_flags, " digits.");

      throw error;
    };
  }


  // **************************************************
  // ISODUAL TABLE ENTRY
  // **************************************************

  /// Entry in the dual isosurface lookup table.
  template <typename NTYPE, typename ISOV_TYPE>
  class ISODUAL_TABLE_ENTRY {

  public:

    typedef ISOV_TYPE ISO_VERTEX_INDEX_TYPE;

    NTYPE num_vertices;         ///< Number of dualiso vertices in cube.
    ISODUAL_TABLE_ENTRY();      ///< constructor
    ~ISODUAL_TABLE_ENTRY();     ///< destructor

    /// incident_isovertex[kf] = Isosurface vertex incident on face kf.
    ///       Face kf is dual to polytope edge kf.
    ISOV_TYPE * incident_isovertex;

    /// is_bipolar[ke] = True if polytope edge ke is bipolar.
    ///       Cube edge ke is dual to isosurface face kf.
    bool * is_bipolar;


    /// Return number of dual isosurface vertices in cube.
    NTYPE NumVertices() const
    { return(num_vertices); }

    /// Return true if edge ke is bipolar.
    template <typename NTYPE2>
    bool IsBipolar(const NTYPE2 ke) const
    { return(is_bipolar[ke]); }

    /// Return incident isovertex.
    template <typename NTYPE2>
    ISOV_TYPE IncidentIsoVertex(const NTYPE2 ke) const
    { return(incident_isovertex[ke]); }

    /// Allocate incident_isovert[] and is_bipolar[].
    /// - num_poly_vertices is not used in this Allocate()
    ///   but other table entries may use it.
    template <typename NUMV_TYPE, typename NUME_TYPE>
    void Allocate(const NUMV_TYPE num_poly_vertices,
                  const NUME_TYPE num_poly_edges);

    bool Check(IJK::ERROR & error_msg) const;
    void FreeAll();                            // free all memory
  };


  // **************************************************
  // IVOLDUAL TABLE ENTRY
  // **************************************************

  template <typename ISOV_TYPE>
  class IVOLDUAL_EDGE_INFO {

  protected:
    void Init();

  public:
    /// Isosurface vertex incident on isosurface polytope dual to edge
    ///   on lower isosurface.
    ISOV_TYPE lower_incident_isovertex;

    /// Isosurface vertex incident on isosurface polytope dual to edge
    ///   on upper isosurface.
    ISOV_TYPE upper_incident_isovertex;

    /// If flag_dual_isopoly is true, then some isosurface polytope
    ///   is dual to the edge.
    /// - If true, then both lower_incident_isovertex
    ///   and upper_incident_isovertex are set.
    bool flag_dual_isopoly;

    IVOLDUAL_EDGE_INFO() 
    { Init(); }

    /// Set lower_incident_isovertex.
    template <typename ISOV_TYPE2>
    void SetLowerIncident(const ISOV_TYPE2 ilower);

    /// Set upper_incident_isovertex.
    template <typename ISOV_TYPE2>
    void SetUpperIncident(const ISOV_TYPE2 iupper);
  };

  template <typename ISOV_TYPE>
  class IVOLDUAL_VERTEX_INFO {

  protected:
    void Init();

  public:
    /// Isosurface vertex incident on isosurface polytope dual to edge.
    ISOV_TYPE incident_isovertex;

    /// If flag_dual_isopoly is true, then some isosurface polytope
    ///   is dual to the edge.
    /// - If true, then both lower_incident_isovertex
    ///   and upper_incident_isovertex are set.
    bool flag_dual_isopoly;

    /// Set incident_isovertex.
    template <typename ISOV_TYPE2>
    void SetIncident(const ISOV_TYPE2 isov);

    IVOLDUAL_VERTEX_INFO() 
    { Init(); }
  };

  /// Entry in the dual interval volume lookup table.
  template <typename NTYPE, typename ISOV_TYPE>
  class IVOLDUAL_TABLE_ENTRY {

  public:

    typedef ISOV_TYPE IVOL_VERTEX_INDEX_TYPE;

    NTYPE num_vertices;         ///< Number of dualiso vertices in cube.
    IVOLDUAL_TABLE_ENTRY();      ///< constructor
    ~IVOLDUAL_TABLE_ENTRY();     ///< destructor

    IVOLDUAL_VERTEX_INFO<ISOV_TYPE> * vertex_info;
    IVOLDUAL_EDGE_INFO<ISOV_TYPE> * edge_info;

    /// Allocate vertex_info[] and edge_info[].
    template <typename NUMV_TYPE, typename NUME_TYPE>
    void Allocate(const NUMV_TYPE num_poly_vertices,
                  const NUME_TYPE num_poly_edges);

    bool Check(IJK::ERROR & error_msg) const;
    void FreeAll();                            // free all memory
  };


  // **************************************************
  // ISODUAL TABLE BASE
  // **************************************************

  /// Dual isosurface lookup table.
  /// Stores isosurface vertices and incident faces for each configuration 
  ///   of +/- labels at polytope vertex.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  class ISODUAL_TABLE_BASE {

  public:

  /// Index of entry in isosurface lookup table.
  /// Define within ISODUAL_TABLE_BASE for use in templates.
  typedef TI_TYPE TABLE_INDEX;    


  protected:

    DTYPE dimension;               ///< Dimension.
    NTYPE num_poly_vertices;       ///< Number of polytope vertices.
    NTYPE num_poly_edges;          ///< Number of polytope edges;
    ENTRY_TYPE * entry;            ///< Array of dual isosurface table entries.
    TI_TYPE num_table_entries;     ///< Number of entries in table.

    /// Maximum number of vertices allowed for cube.
    NTYPE max_num_vertices; 

    /// True, if array entry[] is allocated.
    bool is_table_allocated;  

    /// Initialization routine.
    template <typename DTYPE2>
    void Init(const DTYPE2 dimension);


  public:
    ISODUAL_TABLE_BASE();
    template <typename DTYPE2>
    ISODUAL_TABLE_BASE(const DTYPE2 d);
    ~ISODUAL_TABLE_BASE();                ///< Destructor

    // Get functions.
    DTYPE Dimension() const            ///< Return dimension.
    { return(dimension); };
    NTYPE NumPolyVertices() const      ///< Return number of polytope vertices.
    { return(num_poly_vertices); };
    NTYPE NumPolyEdges() const         ///< Return number of polytope edges.
    { return(num_poly_edges); };

    /// Return number of lookup table entries.
    TI_TYPE NumTableEntries() const { return(num_table_entries); };

    /// Return complement of table index it
    template <typename TI_TYPE2>
    TI_TYPE2 Complement(const TI_TYPE2 it) const
    { return(compute_complement(it, num_table_entries)); }

    /// Return number of vertices in isosurface patch for table entry \a it.
    template <typename TI_TYPE2>
    NTYPE NumIsoVertices(const TI_TYPE2 it) const
    { return(entry[it].num_vertices); }; 

    /// Return index of isosurface vertex incident on face kf.
    /// Undefined if polytope edge k is not bipolar.
    /// @param it Index of table entry.
    /// @param kf Isosurface face kf, dual to polytope edge kf.
    template <typename TI_TYPE2>
    NTYPE IncidentIsoVertex
    (const TI_TYPE2 it, const NTYPE kf) const
    { return(entry[it].incident_isovertex[kf]); };

    /// Return true if vertex iv is positive.
    /// @param iv Vertex index.
    template <typename TI_TYPE2>
    bool IsPositive(const TI_TYPE2 it, const NTYPE iv) const;

    /// Return maximum number of polytope vertices permitted in any table.
    /// Note: Even tables for polytopes of this size are probably impossible 
    ///   to compute/store.
    NTYPE MaxNumVertices() const { return(max_num_vertices); };

    /// Return true if table memory is allocated.
    bool IsTableAllocated() const
    { return(is_table_allocated); };

    // Set functions.
    template <typename DTYPE2>
    void SetDimension(const DTYPE2 d);
    void SetNumPolyVertices(const NTYPE num_vertices);
    void SetNumPolyEdges(const NTYPE num_edges);

    /// Set polytope to be a cube.
    template <typename DTYPE2>
    void SetToCube(const DTYPE2 d);

    /// Allocate table
    template <typename NTYPE2>
    void SetNumTableEntries(const NTYPE2 num_table_entries);

    // Check functions
    template <typename DTYPE2>
    bool CheckDimension(const DTYPE2 d) const;
    bool CheckDimension() const
    { return(CheckDimension(Dimension())); };
    bool CheckTable(IJK::ERROR & error_msg) const;
    bool Check(IJK::ERROR & error_msg) const;

    virtual void FreeAll();                     /// Free all memory.
  };


  // **************************************************
  // ISODUAL CUBE TABLE
  // **************************************************

  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  class ISODUAL_CUBE_TABLE:
    public ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE> {

  protected:
    //! If true, separate negative vertices.
    bool flag_separate_neg;

    //! If true, always separate two diagonally opposite
    //!   positive or negative vertices.
    bool flag_always_separate_opposite;     

    /// Create table entries.
    /// @param flag_separate_opposite If true, always separate two
    ///        diagonally opposite positive or negative vertices.
    void CreateTableEntries
      (const bool flag_separate_neg, const bool flag_separate_opposite);

  public:
    ISODUAL_CUBE_TABLE() {};

    template <typename DTYPE2>
    ISODUAL_CUBE_TABLE(const DTYPE2 dimension) :
      ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE> (dimension)
    { Create(dimension); }
      
    template <typename DTYPE2>
    ISODUAL_CUBE_TABLE
    (const DTYPE2 dimension, const bool flag_separate_opposite) :
      ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>(dimension)
    { Create(dimension, flag_separate_opposite); }

    template <typename DTYPE2>
    ISODUAL_CUBE_TABLE
      (const DTYPE2 dimension, const bool flag_separate_neg,
       const bool flag_separate_opposite) :
        ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>(dimension)
    { Create(dimension, flag_separate_neg, flag_separate_opposite); }


    // Get functions.

    /// Return true if edge ke is bipolar.
    /// @param it Index of table entry.
    /// @param ke Polytope edge ke, dual to isosurface face ke.
    template <typename TI_TYPE2>
    bool IsBipolar(const TI_TYPE2 it, const NTYPE ke) const
    { return(this->entry[it].is_bipolar[ke]); };

    template <typename DTYPE2>
    void Create(const DTYPE2 dimension);
    template <typename DTYPE2>
    void Create(const DTYPE2 dimension, const bool flag_separate_opposite);
    template <typename DTYPE2>
    void Create(const DTYPE2 dimension, const bool flag_separate_neg, 
                const bool flag_separate_opposite);
  };


  // **************************************************
  // ISODUAL TABLE AMBIG BASE
  // **************************************************

  //! Isodual table plus ambiguity information.
  //! @tparam TI_TYPE Table index type.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  class ISODUAL_TABLE_AMBIG_BASE:
    public ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE> {

  public:
    typedef unsigned char FACET_INDEX;          ///< Index of facet.
    typedef int FACET;          ///< Bits representing vertices in facet.
    typedef int FACET_SET;      ///< Bits representing set of facets.

  protected:
    bool * is_ambiguous;            ///< True for ambiguous configurations.
    FACET_SET * ambiguous_facet;    ///< k'th bit is 1 if facet k is ambiguous

    /// Number of ambiguous facets.
    FACET_INDEX * num_ambiguous_facets;           

    /// Number of active facets (with both positive and negative vertices.)
    FACET_INDEX * num_active_facets;

    void Init();                    ///< Initialization routine.
    void Alloc();                   ///< Allocate memory.
    void FreeAll();                 ///< Free all memory.

  public:
    
    // constructors
    ISODUAL_TABLE_AMBIG_BASE()
    { Init(); };

    template <typename DTYPE2>
    ISODUAL_TABLE_AMBIG_BASE(const DTYPE2 dimension) :
      ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>(dimension)
    { Init(); };

    ~ISODUAL_TABLE_AMBIG_BASE();                // destructor

    // get functions
    bool IsAmbiguous(const TI_TYPE it) const
    { return(is_ambiguous[it]); };
    bool IsFacetAmbiguous
      (const TI_TYPE it, const FACET_INDEX jf) const
    { return((ambiguous_facet[it] & ((1L) << jf)) != 0); };
    FACET_SET AmbiguousFacetBits(const TI_TYPE it) const
    { return(ambiguous_facet[it]); };
    FACET_INDEX NumAmbiguousFacets(const TI_TYPE it) const
    { return(num_ambiguous_facets[it]); };
    FACET_INDEX NumActiveFacets(const TI_TYPE it) const
    { return(num_active_facets[it]); };

    /// Allocate table with ambiguity information.
    template <typename NTYPE2>
    void SetNumTableEntries(const NTYPE2 num_table_entries);
  };


  // **************************************************
  // ISODUAL CUBE TABLE AMBIG
  // **************************************************

  //! Isodual cube table plus ambiguity information.
  //! @tparam TI_TYPE Table index type.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  class ISODUAL_CUBE_TABLE_AMBIG:
    public ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE> {

  protected:
    //! If true, separate negative vertices.
    bool flag_separate_neg;

    //! If true, always separate two diagonally opposite
    //!   positive or negative vertices.
    bool flag_always_separate_opposite;     

    /// Create table entries.
    /// @param flag_separate_opposite If true, always separate two
    ///        diagonally opposite positive or negative vertices.
    void CreateTableEntries
    (const bool flag_separate_neg, const bool flag_separate_opposite);

    /// Compute ambiguity information.
    void ComputeAmbiguityInformation();

    /// Compute number of active facets.
    void ComputeNumActiveFacets();

  public:

    ISODUAL_CUBE_TABLE_AMBIG() {};

    template <typename DTYPE2>
    ISODUAL_CUBE_TABLE_AMBIG(const DTYPE2 dimension) :
      ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>
      (dimension)
    { Create(dimension); }
      
    template <typename DTYPE2>
    ISODUAL_CUBE_TABLE_AMBIG
    (const DTYPE2 dimension, const bool flag_separate_opposite) :
      ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>(dimension)
    { Create(dimension, flag_separate_opposite); }

    template <typename DTYPE2>
    ISODUAL_CUBE_TABLE_AMBIG
    (const DTYPE2 dimension, const bool flag_separate_neg,
     const bool flag_separate_opposite) :
      ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>(dimension)
    { Create(dimension, flag_separate_neg, flag_separate_opposite); }

    // Get functions.

    /// Return true if edge ke is bipolar.
    /// @param it Index of table entry.
    /// @param ke Polytope edge ke, dual to isosurface face ke.
    template <typename TI_TYPE2>
    bool IsBipolar(const TI_TYPE2 it, const NTYPE ke) const
    { return(this->entry[it].is_bipolar[ke]); };

    /// Return true if separate negative vertices.
    bool FlagSeparateNeg() const
    { return(flag_separate_neg); }

    template <typename DTYPE2>
    void Create(const DTYPE2 dimension);
    template <typename DTYPE2>
    void Create(const DTYPE2 dimension, const bool flag_separate_opposite);
    template <typename DTYPE2>
    void Create(const DTYPE2 dimension, const bool flag_separate_neg, 
                const bool flag_separate_opposite);
  };



  // **************************************************
  // IVOLDUAL CUBE TABLE
  // **************************************************

  //! Isodual cube table plus ambiguity information.
  //! @tparam TI_TYPE Table index type.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  class IVOLDUAL_CUBE_TABLE:
    public ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE> {


  protected:
    //! If true, separate negative vertices.
    bool flag_separate_neg;

    //! If true, always separate two diagonally opposite
    //!   positive or negative vertices.
    bool flag_always_separate_opposite;     

    /// Create table entries.
    /// @param flag_separate_opposite If true, always separate two
    ///        diagonally opposite positive or negative vertices.
    void CreateTableEntries
    (const bool flag_separate_neg, const bool flag_separate_opposite);

  public:

    IVOLDUAL_CUBE_TABLE() {};

    template <typename DTYPE2>
    IVOLDUAL_CUBE_TABLE(const DTYPE2 dimension) :
      ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>
      (dimension) 
    { Create(dimension); };
      
    template <typename DTYPE2>
    IVOLDUAL_CUBE_TABLE
    (const DTYPE2 dimension, const bool flag_separate_opposite) :
      ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>
      (dimension) 
    { Create(dimension, flag_separate_opposite); }

    template <typename DTYPE2>
    IVOLDUAL_CUBE_TABLE
    (const DTYPE2 dimension, const bool flag_separate_neg,
     const bool flag_separate_opposite) :
      ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>
      (dimension)
    { Create(dimension, flag_separate_neg, flag_separate_opposite); }


    // Get functions.

    /// Return true if edge ie has dual isosurface polytope.
    template <typename TI_TYPE2, typename NTYPE2>
    bool EdgeHasDualIsoPoly(const TI_TYPE2 ientry, const NTYPE2 ie) const
    { return(this->entry[ientry].edge_info[ie].flag_dual_isopoly); }

    /// Return lower incident isovertex on polytope dual to edge ie.
    template <typename TI_TYPE2, typename NTYPE2>
    bool LowerIncident(const TI_TYPE2 ientry, const NTYPE2 ie) const
    { return(this->entry[ientry].edge_info[ie].lower_incident_isovertex); }

    /// Return upper incident isovertex on polytope dual to edge ie.
    template <typename TI_TYPE2, typename NTYPE2>
    bool UpperIncident(const TI_TYPE2 ientry, const NTYPE2 ie) const
    { return(this->entry[ientry].edge_info[ie].upper_incident_isovertex); }

    /// Return true if vertex iv has dual isosurface polytope.
    template <typename TI_TYPE2, typename NTYPE2>
    bool VertexHasDualIsoPoly(const TI_TYPE2 ientry, const NTYPE2 iv) const
    { return(this->entry[ientry].vertex_info[iv].flag_dual_isopoly); }

    /// Return incident isovertex on polytope dual to cube vertex iv.
    template <typename TI_TYPE2, typename NTYPE2>
    bool IncidentIsoVertex(const TI_TYPE2 ientry, const NTYPE2 iv) const
    { return(this->entry[ientry].vertex_info[iv].incident_isovertex); }



    template <typename DTYPE2>
    void Create(const DTYPE2 dimension);
    template <typename DTYPE2>
    void Create(const DTYPE2 dimension, const bool flag_separate_opposite);
    template <typename DTYPE2>
    void Create(const DTYPE2 dimension, const bool flag_separate_neg, 
                const bool flag_separate_opposite);
  };


  // **************************************************
  // CLASS FIND_COMPONENT
  // **************************************************

  /// Find connected component among cube vertices.
  template <typename DTYPE, typename NTYPE>
  class FIND_COMPONENT {

  protected:
    DTYPE dimension;
    NTYPE num_cube_vertices;
    bool * vertex_flag;
    NTYPE * component;

  public:
    typedef DTYPE DIMENSION_TYPE;         ///< Dimension type.
    typedef NTYPE NUMBER_TYPE;            ///< Number type.

  public:
    template <typename DTYPE2>
    FIND_COMPONENT(const DTYPE2 dimension);
    ~FIND_COMPONENT();

    // set functions
    template <typename TI_TYPE2>
    void SetVertexFlags(const TI_TYPE2 ival);
    void NegateVertexFlags();
    void ClearAll();

    // get functions
    DTYPE Dimension() const
    { return(dimension); }
    bool VertexFlag(const int i) const
    { return(vertex_flag[i]); }
    template <typename ITYPE>
    NTYPE Component(const ITYPE i) const
    { return(component[i]); }
    NTYPE NumCubeVertices() const
    { return(num_cube_vertices); }

    /// Search starting at vertex i.
    /// @pre icomp is not zero.
    template <typename ITYPE0, typename ITYPE1>
    void Search(const ITYPE0 i, const ITYPE1 icomp);

    /// Search facet starting at vertex i.
    /// @pre Facet kf contains vertex i.
    /// @pre icomp is not zero.
    template <typename ITYPE0, typename ITYPE1, typename ITYPE2>
    void SearchFacet(const ITYPE0 kf, const ITYPE1 i, const ITYPE2 icomp);

    /// Compute number of components.
    /// @param flag_positive If true, compute components of positive vertices.
    ///                      If false, compute components of negative vertices.
    template <typename TI_TYPE2>
    NTYPE ComputeNumComponents
    (const TI_TYPE2 ientry, const bool flag_positive);

    /// Compute number of components in facet.
    /// @param flag_positive If true, compute components of positive vertices.
    ///                      If false, compute components of negative vertices.
    template <typename TI_TYPE2, typename ITYPE>
    NTYPE ComputeNumComponentsInFacet
    (const TI_TYPE2 ientry, const ITYPE kf, const bool flag_positive);
  };


  // **************************************************
  // CLASS ISODUAL_CUBE_FACE_INFO
  // *************************************************

  /// Class with routines to query number of positive or negative vertices
  ///   in a cube facet and number of active facets.
  template <typename DTYPE, typename NTYPE, typename VTYPE>
  class ISODUAL_CUBE_FACE_INFO:public IJK::CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE> {

  public:
    ISODUAL_CUBE_FACE_INFO() {};
    ISODUAL_CUBE_FACE_INFO(const DTYPE dimension):
      IJK::CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>(dimension) {};

    template <typename TI_TYPE>
    void ComputeNumCubeFacetBits
    (const TI_TYPE ientry, const NTYPE ifacet,
     NTYPE & num_zeros, NTYPE & num_ones) const;

    template <typename TI_TYPE>
    bool IsCubeFacetActive(const TI_TYPE ientry, const NTYPE ifacet) const;

    template <typename TI_TYPE>
    NTYPE ComputeNumActiveCubeFacets(const TI_TYPE ientry) const;
  };


  // **************************************************
  // ISODUAL TABLE ENTRY MEMBER FUNCTIONS
  // **************************************************

  // constructor
  template <typename NTYPE, typename ISOV_TYPE>
  ISODUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE>::ISODUAL_TABLE_ENTRY()
  {
    num_vertices = 0;
    incident_isovertex = NULL;
    is_bipolar = NULL;
  }

  // destructor
  template <typename NTYPE, typename ISOV_TYPE>
  ISODUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE>::~ISODUAL_TABLE_ENTRY()
  {
    FreeAll();
  }

  template <typename NTYPE, typename ISOV_TYPE>
  template <typename NUMV_TYPE, typename NUME_TYPE>
  void ISODUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE>::
  Allocate(const NUMV_TYPE num_poly_vertices,
           const NUME_TYPE num_poly_edges)
  {
    FreeAll();

    incident_isovertex = new ISOV_TYPE[num_poly_edges];
    is_bipolar = new bool[num_poly_edges];
  }

  template <typename NTYPE, typename ISOV_TYPE>
  bool ISODUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE>::
  Check(IJK::ERROR & error_msg) const
  {
    if (num_vertices < 0) {
      error_msg.AddMessage
        ("Error.  Dual isosurface table entry contains negative number of isosurface vertices.");
      return(false);
    }

    if (incident_isovertex == NULL) {
      error_msg.AddMessage
        ("Memory for incident isosurface vertices not allocated.");
      return(false);
    }

    if (is_bipolar == NULL) {
      error_msg.AddMessage("Memory for bipolar edge flags not allocated.");
      return(false);
    }

    return(true);
  }

  // Free all memory.
  template <typename NTYPE, typename ISOV_TYPE>
  void ISODUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE>::FreeAll()
  {
    if (incident_isovertex != NULL) {
      delete [] incident_isovertex;
      incident_isovertex = NULL;
    }

    if (is_bipolar != NULL) {
      delete [] is_bipolar;
      is_bipolar = NULL;
    }

    num_vertices = 0;
  }


  // **************************************************
  // IVOLDUAL TABLE ENTRY MEMBER FUNCTIONS
  // **************************************************
 
  // IVOLDUAL_EDGE_INFO initalization routine
  template <typename ISOV_TYPE>
  void IVOLDUAL_EDGE_INFO<ISOV_TYPE>::Init()
  {
    flag_dual_isopoly = false;
  }

  template <typename ISOV_TYPE>
  template <typename ISOV_TYPE2>
  void IVOLDUAL_EDGE_INFO<ISOV_TYPE>::
  SetLowerIncident(const ISOV_TYPE2 ilower)
  {
    lower_incident_isovertex = ilower;
    flag_dual_isopoly = true;
  }

  template <typename ISOV_TYPE>
  template <typename ISOV_TYPE2>
  void IVOLDUAL_EDGE_INFO<ISOV_TYPE>::
  SetUpperIncident(const ISOV_TYPE2 iupper)
  {
    upper_incident_isovertex = iupper;
    flag_dual_isopoly = true;
  }

  // IVOLDUAL_VERTEX_INFO initalization routine
  template <typename ISOV_TYPE>
  void IVOLDUAL_VERTEX_INFO<ISOV_TYPE>::Init()
  {
    flag_dual_isopoly = false;
  }

  // Set incident_isovertex.
  template <typename ISOV_TYPE>
  template <typename ISOV_TYPE2>
  void IVOLDUAL_VERTEX_INFO<ISOV_TYPE>::
  SetIncident(const ISOV_TYPE2 isov)
  {
    incident_isovertex = isov;
    flag_dual_isopoly = true;
  }

  // constructor
  template <typename NTYPE, typename ISOV_TYPE>
  IVOLDUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE>::IVOLDUAL_TABLE_ENTRY()
  {
    num_vertices = 0;
    vertex_info = NULL;
    edge_info = NULL;
  }

  // destructor
  template <typename NTYPE, typename ISOV_TYPE>
  IVOLDUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE>::~IVOLDUAL_TABLE_ENTRY()
  {
    FreeAll();
  }

  template <typename NTYPE, typename ISOV_TYPE>
  template <typename NUMV_TYPE, typename NUME_TYPE>
  void IVOLDUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE>::
  Allocate(const NUMV_TYPE num_poly_vertices,
           const NUME_TYPE num_poly_edges)
  {
    FreeAll();

    vertex_info = new IVOLDUAL_VERTEX_INFO<ISOV_TYPE>[num_poly_vertices];
    edge_info = new IVOLDUAL_EDGE_INFO<ISOV_TYPE>[num_poly_edges];
  }

  template <typename NTYPE, typename ISOV_TYPE>
  bool IVOLDUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE>::
  Check(IJK::ERROR & error_msg) const
  {
    if (num_vertices < 0) {
      error_msg.AddMessage
        ("Programming error.  Dual interval volume table entry contains negative number");
      error_msg.AddMessage("  of interval volume vertices.");
      return(false);
    }

    if (vertex_info == NULL) {
      error_msg.AddMessage
        ("Programming error.  Memory for isosurface elements dual to polytope vertices not allocated.");
      return(false);
    }

    if (edge_info == NULL) {
      error_msg.AddMessage
        ("Programming errror.  Memory for isosurface elements dual to polytope edges not allocated.");
      return(false);
    }

    return(true);
  }

  // Free all memory.
  template <typename NTYPE, typename ISOV_TYPE>
  void IVOLDUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE>::FreeAll()
  {
    if (vertex_info != NULL) {
      delete [] vertex_info;
      vertex_info = NULL;
    }

    if (edge_info != NULL) {
      delete [] edge_info;
      edge_info = NULL;
    }

    num_vertices = 0;
  }


  // **************************************************
  // ISODUAL TABLE BASE MEMBER FUNCTIONS
  // **************************************************

  // default constructor. dimension = 3
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::ISODUAL_TABLE_BASE()
  {
    const DTYPE DIM3(3);

    Init(DIM3);
  }

  // Constructor.
  // d = dimension of space containing isosurface.  Should be 2, 3 or 4.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  ISODUAL_TABLE_BASE(const DTYPE2 d)
  {
    Init(d);
  }

  // Initialize
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Init(const DTYPE2 dimension)
  {
    const char * procname = "ISODUAL_TABLE_BASE::Init";

    max_num_vertices = 20;
    // Note: Even tables for polytopes of this size are probably impossible 
    //   to compute/store

    num_table_entries = 0;
    entry = NULL;
    is_table_allocated = false;

    SetDimension(dimension);
  }

  // destructor
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::~ISODUAL_TABLE_BASE()
  {
    FreeAll();
  }

  // Set dimension.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  SetDimension(const DTYPE2 dimension)
  {
    const char * procname = "ISODUAL_TABLE_BASE::SetDimension";

    if (IsTableAllocated()) {
      throw IJK::PROCEDURE_ERROR
        (procname, "Dual isosurface table already allocated.");
    }

    this->dimension = dimension;

    if (!CheckDimension())
      throw IJK::PROCEDURE_ERROR(procname, "Illegal polytope dimension.");
  }

  // Set number of polytope vertices.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  void ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  SetNumPolyVertices(const NTYPE num_vertices)
  {
    const char * procname = "ISODUAL_TABLE_BASE::SetNumVertices";

    if (IsTableAllocated()) {
      throw IJK::PROCEDURE_ERROR
        (procname, "Dual isosurface table already allocated.");
    }

    num_poly_vertices = num_vertices;
  }

  // Set number of polytope edges
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  void ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  SetNumPolyEdges(const NTYPE num_edges)
  {
    const char * procname = "ISODUAL_TABLE_BASE::SetNumEdges";

    if (IsTableAllocated()) {
      throw IJK::PROCEDURE_ERROR
        (procname, "Dual isosurface table already allocated.");
    }

    num_poly_edges = num_edges;
  }

  // Set polytope to be a cube.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  SetToCube(const DTYPE2 dimension)
  {
    const char * procname = "ISODUAL_TABLE_BASE::SetToCube";

    if (IsTableAllocated()) {
      throw IJK::PROCEDURE_ERROR
        (procname, "Dual isosurface table already allocated.");
    }

    this->SetDimension(dimension);

    NTYPE num_vertices = IJK::compute_num_cube_vertices(dimension);
    this->SetNumPolyVertices(num_vertices);

    NTYPE num_edges = IJK::compute_num_cube_edges(dimension);
    this->SetNumPolyEdges(num_edges);
  }


  // Allocate table
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename NTYPE2>
  void ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  SetNumTableEntries(const NTYPE2 num_table_entries)
  {
    const char * procname = "ISODUAL_TABLE_BASE::SetNumTableEntries";

    if (entry != NULL) delete [] entry;
    entry = NULL;
    this->num_table_entries = 0;
    is_table_allocated = false;

    entry = new ENTRY_TYPE[num_table_entries];
    if (entry == NULL && num_table_entries > 0)
      throw IJK::PROCEDURE_ERROR
        (procname, "Unable to allocate memory for dual isosurface table.");

    for (int i = 0; i < num_table_entries; i++) {
      entry[i].Allocate(NumPolyVertices(), NumPolyEdges());
    }

    this->num_table_entries = num_table_entries;
    is_table_allocated = true;
  }

  // Return true if vertex iv is positive.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename TI_TYPE2>
  bool ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  IsPositive(const TI_TYPE2 it, const NTYPE iv) const
  {
    const TI_TYPE2 mask = (TI_TYPE2(1) << iv);

    if ((it & mask) == 0) { return(false); }
    else { return(true); }
  }

  // check dimension
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  bool ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  CheckDimension(const DTYPE2 d) const
  {
    if (d < 1)
      return(false);
    else
      return(true);
  }


  // check table
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  bool ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  CheckTable(IJK::ERROR & error_msg) const
  {
    if (NumPolyVertices() >= CHAR_BIT*sizeof(TI_TYPE)) {
      error_msg.AddMessage("Too many polytope vertices");
      return(false);
    }

    if (NumPolyVertices() > MaxNumVertices()) {
      error_msg.AddMessage("Too many polytope vertices");
      return(false);
    }

    if (NumPolyVertices() < 1) {
      error_msg.AddMessage("Polytope must have at least one vertex.");
      return(false);
    }

    if (entry == NULL) {
      error_msg.AddMessage("Memory for dual isosurface table not allocated.");
      return(false);
    }

    for (TI_TYPE it = 0; it < NumTableEntries(); it++)
      if (!entry[it].Check(error_msg)) {
        error_msg.AddMessage
          ("Error detected at isosurface table entry ", it, ".");
        return(false);
      }

    return(true);
  }


  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  bool ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Check(IJK::ERROR & error_msg) const
  {
    if (!CheckTable(error_msg)) return(false);
    return(true);
  }


  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  void ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::FreeAll()
  // free all memory
  {
    if (entry != NULL) {
      for (TI_TYPE i = 0; i < num_table_entries; i++)
        { entry[i].FreeAll(); }
      delete [] entry;
      entry = NULL;
    };
    num_table_entries = 0;
    is_table_allocated = false;
  }


  // **************************************************
  // ISODUAL CUBE TABLE MEMBER FUNCTIONS
  // **************************************************

  // Create table.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void ISODUAL_CUBE_TABLE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Create(const DTYPE2 dimension, const bool flag_separate_neg,
         const bool flag_separate_opposite)
  {
    this->SetToCube(dimension);

    TI_TYPE n = calculate_num_entries<TI_TYPE>(this->NumPolyVertices(), 2);
    this->SetNumTableEntries(n);
    CreateTableEntries(flag_separate_neg, flag_separate_opposite);
  }


  // Create table.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void ISODUAL_CUBE_TABLE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Create(const DTYPE2 dimension, const bool flag_separate_opposite)
  {
    Create(dimension, true, flag_separate_opposite);
  }


  // Create table.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void ISODUAL_CUBE_TABLE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Create(const DTYPE2 dimension)
  {
    Create(dimension, true);
  }


  // Create table entries.
  // @param flag_separate_neg  If true, separate negative vertices.
  // @param flag_separate_opposite If true, always separate two diagonally 
  //        opposite negative or positive vertices 
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  void ISODUAL_CUBE_TABLE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  CreateTableEntries
  (const bool flag_separate_neg, const bool flag_separate_opposite)
  {
    this->flag_separate_neg = flag_separate_neg;
    this->flag_always_separate_opposite = flag_separate_opposite;

    const DTYPE dimension = this->Dimension();

    FIND_COMPONENT<DTYPE,NTYPE> find_component(dimension);
    IJK::CUBE_FACE_INFO<NTYPE,NTYPE,NTYPE> cube(dimension);

    for (TI_TYPE ientry = 0; ientry < this->NumTableEntries(); ientry++) {
      create_isodual_cube_table_entry
        (ientry, this->NumPolyVertices(), 
         flag_separate_neg, flag_separate_opposite, cube,
         find_component, this->entry[ientry]);
    }
  }


  // **************************************************
  // ISODUAL TABLE AMBIG BASE MEMBER FUNCTIONS
  // **************************************************

  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  void ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Init()
  {
    is_ambiguous = NULL;
    ambiguous_facet = NULL;
    num_ambiguous_facets = NULL;
    num_active_facets = NULL;
  }


  // Destructor.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  ~ISODUAL_TABLE_AMBIG_BASE()
  {
    FreeAll();
  }

  // Free all memory.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  void ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  FreeAll()
  {
    delete [] is_ambiguous;
    is_ambiguous = NULL;
    delete [] ambiguous_facet;
    ambiguous_facet = NULL;
    delete [] num_ambiguous_facets;
    num_ambiguous_facets = NULL;
    delete [] num_active_facets;
    num_active_facets = NULL;
    this->num_table_entries = 0;
  }

  // Allocate memory.
  // Free any previously allocated memory.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  void ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Alloc()
  {
    if (is_ambiguous != NULL) { FreeAll(); }

    is_ambiguous = new bool[this->NumTableEntries()];
    ambiguous_facet = new FACET_SET[this->NumTableEntries()];
    num_ambiguous_facets = new FACET_INDEX[this->NumTableEntries()];
    num_active_facets = new FACET_INDEX[this->NumTableEntries()];
  }


  // Set number of table entries
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename NTYPE2>
  void ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  SetNumTableEntries(const NTYPE2 num_table_entries)
  {
    ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
      SetNumTableEntries(num_table_entries);
    Alloc();
  }


  // **************************************************
  // ISODUAL CUBE TABLE AMBIG MEMBER FUNCTIONS
  // **************************************************

  // Create table.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void ISODUAL_CUBE_TABLE_AMBIG<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Create(const DTYPE2 dimension, const bool flag_separate_neg,
         const bool flag_separate_opposite)
  {
    this->SetToCube(dimension);

    TI_TYPE n = calculate_num_entries<TI_TYPE>(this->NumPolyVertices(), 2);
    this->SetNumTableEntries(n);
    CreateTableEntries(flag_separate_neg, flag_separate_opposite);
  }


  // Create table.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void ISODUAL_CUBE_TABLE_AMBIG<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Create(const DTYPE2 dimension, const bool flag_separate_opposite)
  {
    Create(dimension, true, flag_separate_opposite);
  }


  // Create table.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void ISODUAL_CUBE_TABLE_AMBIG<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Create(const DTYPE2 dimension)
  {
    Create(dimension, true);
  }


  // Create table entries.
  // @param flag_separate_neg  If true, separate negative vertices.
  // @param flag_separate_opposite If true, always separate two diagonally 
  //        opposite negative or positive vertices 
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  void ISODUAL_CUBE_TABLE_AMBIG<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  CreateTableEntries
  (const bool flag_separate_neg, const bool flag_separate_opposite)
  {
    this->flag_separate_neg = flag_separate_neg;
    this->flag_always_separate_opposite = flag_separate_opposite;

    const DTYPE dimension = this->Dimension();

    FIND_COMPONENT<DTYPE,NTYPE> find_component(dimension);
    IJK::CUBE_FACE_INFO<NTYPE,NTYPE,NTYPE> cube(dimension);

    for (TI_TYPE ientry = 0; ientry < this->NumTableEntries(); ientry++) {
      create_isodual_cube_table_entry
        (ientry, this->NumPolyVertices(), 
         flag_separate_neg, flag_separate_opposite, cube,
         find_component, this->entry[ientry]);
    }

    ComputeAmbiguityInformation();
    ComputeNumActiveFacets();
  }

  // Compute ambiguity information.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  void ISODUAL_CUBE_TABLE_AMBIG<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  ComputeAmbiguityInformation()
  {
    const DTYPE dimension = this->Dimension();
    const NTYPE num_cube_facets = IJK::compute_num_cube_facets(dimension);
    FIND_COMPONENT<DTYPE,NTYPE> find_component(dimension);

    for (TI_TYPE ientry = 0; ientry < this->NumTableEntries(); ientry++) {
      this->is_ambiguous[ientry] = is_cube_ambiguous(ientry, find_component);

      compute_ambiguous_cube_facets
        (ientry, num_cube_facets, this->ambiguous_facet[ientry], 
         this->num_ambiguous_facets[ientry], find_component);
    }
  }

  // Compute number of active facets.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  void ISODUAL_CUBE_TABLE_AMBIG<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  ComputeNumActiveFacets()
  {
    const DTYPE dimension = this->Dimension();
    ISODUAL_CUBE_FACE_INFO<NTYPE,NTYPE,NTYPE> cube_face_info(dimension);

    for (TI_TYPE ientry = 0; ientry < this->NumTableEntries(); ientry++) {
      this->num_active_facets[ientry] = 
        cube_face_info.ComputeNumActiveCubeFacets(ientry);
    }
  }


  // **************************************************
  // IVOLDUAL CUBE TABLE MEMBER FUNCTIONS
  // **************************************************

  // Create table.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void IVOLDUAL_CUBE_TABLE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Create(const DTYPE2 dimension, const bool flag_separate_neg,
         const bool flag_separate_opposite)
  {
    const int NUM_CUBE_VERTEX_TYPES = 4;

    this->SetToCube(dimension);

    TI_TYPE n = calculate_num_entries<TI_TYPE>
      (this->NumPolyVertices(), NUM_CUBE_VERTEX_TYPES);
    this->SetNumTableEntries(n);
    CreateTableEntries(flag_separate_neg, flag_separate_opposite);
  }


  // Create table.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void IVOLDUAL_CUBE_TABLE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Create(const DTYPE2 dimension, const bool flag_separate_opposite)
  {
    Create(dimension, true, flag_separate_opposite);
  }


  // Create table.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void IVOLDUAL_CUBE_TABLE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Create(const DTYPE2 dimension)
  {
    Create(dimension, true);
  }


  // Create table entries.
  // @param flag_separate_neg  If true, separate negative vertices.
  // @param flag_separate_opposite If true, always separate two diagonally 
  //        opposite negative or positive vertices 
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  void IVOLDUAL_CUBE_TABLE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  CreateTableEntries
  (const bool flag_separate_neg, const bool flag_separate_opposite)
  {
    typedef typename ENTRY_TYPE::IVOL_VERTEX_INDEX_TYPE ISOV_TYPE;

    this->flag_separate_neg = flag_separate_neg;
    this->flag_always_separate_opposite = flag_separate_opposite;

    const DTYPE dimension = this->Dimension();

    FIND_COMPONENT<DTYPE,NTYPE> find_component(dimension+1);
    IJK::CUBE_FACE_INFO<NTYPE,NTYPE,NTYPE> cube(dimension);
    IJK::CUBE_FACE_INFO<NTYPE,NTYPE,NTYPE> lifted_cube(dimension+1);

    ISODUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE> lower_isodual, upper_isodual;

    lower_isodual.Allocate(lifted_cube.NumVertices(), lifted_cube.NumEdges());
    upper_isodual.Allocate(lifted_cube.NumVertices(), lifted_cube.NumEdges());

    for (TI_TYPE ientry = 0; ientry < this->NumTableEntries(); ientry++) {
      // Create ivoldual cube table entry
      create_ivoldual_cube_table_entry
        (ientry, this->NumPolyVertices(), 
         flag_separate_neg, flag_separate_opposite, cube, lifted_cube,
         find_component, lower_isodual, upper_isodual, this->entry[ientry]);
    }

  }


  // **************************************************
  // CLASS FIND_COMPONENT MEMBER FUNCTIONS
  // **************************************************

  template <typename DTYPE, typename NTYPE>
  template <typename DTYPE2>
  FIND_COMPONENT<DTYPE,NTYPE>::FIND_COMPONENT(const DTYPE2 dimension)
  {
    this->dimension = dimension;
    num_cube_vertices = NTYPE(1) << Dimension();

    vertex_flag = new bool[num_cube_vertices];
    component = new NTYPE[num_cube_vertices];
  }


  template <typename DTYPE, typename NTYPE>
  FIND_COMPONENT<DTYPE,NTYPE>::~FIND_COMPONENT()
  {
    if (component != NULL) { delete [] component; }
    component = NULL;

    if (vertex_flag != NULL) { delete [] vertex_flag; }
    vertex_flag = NULL;
  }


  template <typename DTYPE, typename NTYPE>
  void FIND_COMPONENT<DTYPE,NTYPE>::ClearAll()
  {
    for (NTYPE i = 0; i < NumCubeVertices(); i++) {
      vertex_flag[i] = false;
      component[i] = 0;
    }
  }


  template <typename DTYPE, typename NTYPE>
  template <typename TI_TYPE>
  void FIND_COMPONENT<DTYPE,NTYPE>::SetVertexFlags(const TI_TYPE ival)
  {
    convert2bool(ival, vertex_flag, num_cube_vertices);
  }


  template <typename DTYPE, typename NTYPE>
  void FIND_COMPONENT<DTYPE,NTYPE>::NegateVertexFlags()
  {
    for (NTYPE i = 0; i < num_cube_vertices; i++) 
      { vertex_flag[i] = (!vertex_flag[i]); }
  }


  template <typename DTYPE, typename NTYPE>
  template <typename ITYPE0, typename ITYPE1>
  void FIND_COMPONENT<DTYPE,NTYPE>::
  Search(const ITYPE0 i, const ITYPE1 icomp)
  {
    std::vector<NTYPE> vlist;

    vlist.reserve(num_cube_vertices);

    if (icomp == 0) {
      IJK::PROCEDURE_ERROR error("FIND_COMPONENT::Search");
      error.AddMessage("Programming error. Component number cannot be zero.");
      throw error;
    }

    bool flag = vertex_flag[i];
    vlist.push_back(i);
    component[i] = icomp;
    while(vlist.size() > 0) {
      NTYPE j = vlist.back();
      vlist.pop_back();
      for (int d = 0; d < dimension; d++) {
        int j2;
        IJK::compute_cube_vertex_neighbor(j, d, j2);
        if (vertex_flag[j2] == flag && component[j2] == 0) {
          vlist.push_back(j2);
          component[j2] = icomp;
        }
      }
    }
  }


  template <typename DTYPE, typename NTYPE>
  template <typename ITYPE0, typename ITYPE1, typename ITYPE2>
  void FIND_COMPONENT<DTYPE,NTYPE>::
  SearchFacet(const ITYPE0 kf, const ITYPE1 i, const ITYPE2 icomp)
  {
    std::vector<NTYPE> vlist;

    vlist.reserve(num_cube_vertices);

    if (icomp == 0) {
      IJK::PROCEDURE_ERROR error("FIND_COMPONENT::Search");
      error.AddMessage("Programming error. Component number cannot be zero.");
      throw error;
    }

    if (!IJK::cube_facet_contains(dimension, kf, i)) {
      IJK::PROCEDURE_ERROR error("FIND_COMPONENT::Search");
      error.AddMessage("Programming error. Facet ", kf,
                       " does not contain vertex ", i, ".");
      throw error;
    }

    bool flag = vertex_flag[i];
    vlist.push_back(i);
    component[i] = icomp;
    while(vlist.size() > 0) {
      NTYPE j = vlist.back();
      vlist.pop_back();
      for (int d = 0; d < dimension; d++) {
        int j2;
        IJK::compute_cube_vertex_neighbor(j, d, j2);
        if (vertex_flag[j2] == flag && component[j2] == 0 &&
            IJK::cube_facet_contains(dimension, kf, i)) {
          vlist.push_back(j2);
          component[j2] = icomp;
        }
      }
    }
  }


  // Compute number of components of vertices.
  template <typename DTYPE, typename NTYPE>
  template <typename TI_TYPE>
  NTYPE FIND_COMPONENT<DTYPE,NTYPE>::ComputeNumComponents
  (const TI_TYPE ientry, const bool flag_positive)
  {
    ClearAll();
    SetVertexFlags(ientry);
    if (!flag_positive) { NegateVertexFlags(); }

    NTYPE num_components(0);
    for (NTYPE i = 0; i < NumCubeVertices(); i++) {
      if (Component(i) == 0) {

        if (VertexFlag(i)) {
          num_components++;
          Search(i, num_components);
        }
      }
    }

    return(num_components);
  }


  // Compute number of components of vertices in facet.
  template <typename DTYPE, typename NTYPE>
  template <typename TI_TYPE, typename ITYPE>
  NTYPE FIND_COMPONENT<DTYPE,NTYPE>::ComputeNumComponentsInFacet
  (const TI_TYPE ientry, const ITYPE kf, const bool flag_positive)
  {
    ClearAll();
    SetVertexFlags(ientry);
    if (!flag_positive) { NegateVertexFlags(); }

    NTYPE num_components(0);
    for (NTYPE i = 0; i < NumCubeVertices(); i++) {
      if (Component(i) == 0 &&
          IJK::cube_facet_contains(dimension, kf, i)) {

        if (VertexFlag(i)) {
          num_components++;
          SearchFacet(kf, i, num_components);
        }
      }
    }

    return(num_components);
  }


  // **************************************************
  // CLASS ISODUAL_CUBE_FACE_INFO MEMBER FUNCTIONS
  // **************************************************


  template <typename DTYPE, typename NTYPE, typename VTYPE>
  template <typename TI_TYPE>
  void ISODUAL_CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::
  ComputeNumCubeFacetBits(const TI_TYPE ientry, const NTYPE ifacet,
                          NTYPE & num_zeros, NTYPE & num_ones) const
  {
    num_ones = 0;
    num_zeros = 0;

    for (NTYPE j = 0; j < this->NumFacetVertices(); j++) {
      const NTYPE jv = this->FacetVertex(ifacet, j);
      const TI_TYPE mask = (TI_TYPE(1) << jv);
      if ((mask & ientry) == 0) {
        num_zeros++;
      }
      else {
        num_ones++;
      }
    }
  }


  template <typename DTYPE, typename NTYPE, typename VTYPE>
  template <typename TI_TYPE>
  bool ISODUAL_CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::
  IsCubeFacetActive(const TI_TYPE ientry, const NTYPE ifacet) const
  {
    NTYPE num_zeros, num_ones;

    ComputeNumCubeFacetBits(ientry, ifacet, num_zeros, num_ones);
    if (num_zeros > 0 && num_ones > 0) 
      { return(true); }
    else 
      { return(false); }
  }

  template <typename DTYPE, typename NTYPE, typename VTYPE>
  template <typename TI_TYPE>
  NTYPE ISODUAL_CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::
  ComputeNumActiveCubeFacets(const TI_TYPE ientry) const
  {
    NTYPE num_active_facets = 0;
    for (NTYPE ifacet = 0; ifacet < this->NumFacets(); ifacet++) {
      if (IsCubeFacetActive(ientry, ifacet)) 
        { num_active_facets++; }
    }

    return(num_active_facets);
  }

}

#endif
