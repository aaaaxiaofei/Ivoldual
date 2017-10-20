/// \file ivoldualtable.h
/// ivoldual lookup table.

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

#ifndef _IVOLDUALTABLE_
#define _IVOLDUALTABLE_

#include "ijkcube.txx"
#include "ijkdualtableX.txx"

#include "ivoldual_types.h"


namespace IVOLDUAL {

  // **************************************************
  // INTERVAL VOLUME LOOKUP TABLE ENTRY
  // **************************************************

  class IVOLDUAL_TABLE_ENTRY:public 
  IJKDUALTABLE::IVOLDUAL_TABLE_ENTRY
  <int,ISO_VERTEX_INDEX,FACET_BITS_TYPE,TABLE_INDEX> 
  {
  public:
    // ADDITIONAL DATA ABOUT THE IVOLDUAL CONFIGURATION.
  };

  // **************************************************
  // INTERVAL VOLUME VERTEX INFORMATION
  // **************************************************

  class IVOLDUAL_TABLE_VERTEX_INFO {

  public:
    int degree;

    int Degree() const
    { return(degree); }
  };


  // **************************************************
  // INTERVAL VOLUME LOOKUP TABLE
  // **************************************************

  class IVOLDUAL_CUBE_TABLE:
    public IJKDUALTABLE::IVOLDUAL_CUBE_DOUBLE_TABLE
  <4,int,int,TABLE_INDEX,IVOLDUAL_TABLE_ENTRY>
  {
  
  protected:
    void Init();

  public:
    IJKDUALTABLE::DUAL_TABLE_VERTEX_INFO<int,IVOLDUAL_TABLE_VERTEX_INFO>
    vertex_info;

  public:
    IVOLDUAL_CUBE_TABLE(const int dimension, const bool flag_separate_neg):
      IVOLDUAL_CUBE_DOUBLE_TABLE(dimension, flag_separate_neg)
    { Init(); }; 

    const IVOLDUAL_TABLE_VERTEX_INFO & VertexInfo
    (const TABLE_INDEX ientry, const int ivolv) const
    { return(vertex_info.VertexInfo(ientry, ivolv)); }
  };

}

#endif
