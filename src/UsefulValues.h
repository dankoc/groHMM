/***************************************************************************
**
**   Copyright 2009, 2010, 2011 Charles Danko.
**
**   This program is part of the GRO-seq R package
**
**   groHMM is free software: you can redistribute it and/or modify it 
**   under the terms of the GNU General Public License as published by 
**   the Free Software Foundation, either version 3 of the License, or  
**   (at your option) any later version.
**
**   This program is distributed in the hope that it will be useful, but 
**   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
**   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
**   for more details.
**
**   You should have received a copy of the GNU General Public License along 
**   with this program.  If not, see <http://www.gnu.org/licenses/>.
**
***************************************************************************/

/*************************************************************
 *
 *  This file is comprised of prototypes and structure 
 *  definitions for the components of the HMM class.
 *
 *  Started: 11-25-2009.
 *
 *************************************************************/

// Assume that: exp(a) + exp(b) == exp(a); if(a-b > APPROX_EXP_VALUE_THRESHOLD)
#ifndef APPROXEXPVALUETHRESHOLD
#define APPROXEXPVALUETHRESHOLD
static double APPROX_EXP_VALUE_THRESHOLD = 700;
#endif

#ifndef VERYLARGEDOUBLEVALUE
#define VERYLARGEDOUBLEVALUE
static double VERY_LARGE_DOUBLE_VALUE = 1e20;
#endif

/*************************************
 *
 * MAX and MIN macros ... 
 *
 *************************************/
#ifndef max
    #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
    #define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif
