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
#include "ijkdualtable.txx"

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
  // INTERVAL VOLUME LOOKUP TABLE
  // **************************************************

  typedef typename IJKDUALTABLE::
  IVOLDUAL_CUBE_DOUBLE_TABLE<4,int,int,TABLE_INDEX,IVOLDUAL_TABLE_ENTRY>
  IVOLDUAL_CUBE_TABLE;
}

#endif
