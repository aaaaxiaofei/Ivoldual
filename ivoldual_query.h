/// \file ivoldual_query.h
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

#ifndef _IVOLDUAL_QUERY_
#define _IVOLDUAL_QUERY_

#include "ivoldual_types.h"
#include "ivoldual_datastruct.h"
#include "ivoldualtable.h"

/// ivoldual classes and routines.
namespace IVOLDUAL {

  /// Return true if vertex is lower left vertex of an isosurface box.
  /// @param[out] box_corners[] Interval volume vertices which form the
  ///    8 corners of the isosurface box.
  bool is_ivolv_lower_left_vertex_of_isosurface_box
  (const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   const DUAL_IVOLVERT_ARRAY & ivolv_list,
   const IVOL_VERTEX_INDEX ivolv,
   IVOL_VERTEX_INDEX box_corner[8]);

}

#endif
