/**
 * Copyright (C) 2019 Dean De Leo, email: dleo[at]cwi.nl
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

/**
 * The library requires to be set up with a few macros to enable/disable its features. Which is
 * what we do here.
 * This header acts a wrapper before including the actual llama.h header, ensuring all translation
 * units are compiled in a consistent way.
 */

// Benchmark counters?
//#define LL_COUNTERS

// Profile the cost of ll_writable_graph#add_edge_if_not_exists and #build
//#define LL_PROFILE_UPDATES


// We're finally ready to include the actual library
#include "EdgeMapVertexMap/internal/EdgeMap.hpp"
#include "EdgeMapVertexMap/internal/VertexMap.hpp"
#include "EdgeMapVertexMap/internal/VertexSubset.hpp"
#include "SSTGraph/SparseMatrix.hpp"
