/*
 *
 * lhsdefines.h: A C++ header used in the LHS package
 * Copyright (C) 2012  Robert Carnell
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 */

#ifndef LHSDEFINES
#define LHSDEFINES

void improvedLHS_C(int* N, int* K, int* DUP, int* result);
void maximinLHS_C(int* N, int* K, int* DUP, int* result);
void optimumLHS_C(int* N, int* K, int* MAXSWEEPS, double* EPS, int* pOld, int* JLen, int* bVerbose);
void optSeededLHS_C(int* N, int* K, int* MAXSWEEPS, double* EPS, double* pOld, int* JLen, int* bVerbose);

#endif
