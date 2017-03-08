# auxiliary.R
# 
# Copyright (c) 2017, Aarhus University (Denmark)
# 
# Software written by Carlos P. Roca, as research funded by the European Union.
# 
# This software may be modified and distributed under the terms of the MIT 
# license. See the LICENSE file for details.


# Auxiliary functions


# Checks if x is a single atomic value

is.atomic.singleton <- function( x )
{
    is.atomic( x ) && length( x ) == 1
}


# Checks if x is a vector of atomic elements
# Note that is.vector( list( 1:10 ) ) == TRUE

is.atomic.vector <- function( x )
{
    is.atomic( x ) && ( is.null( dim( x ) ) || length( dim( x ) ) == 1 )
}

