# watson_u2.R
# 
# Copyright (c) 2017, Aarhus University (Denmark)
# 
# Software written by Carlos P. Roca, as research funded by the European Union.
# 
# This software may be modified and distributed under the terms of the MIT 
# license. See the LICENSE file for details.


# Watson U2 statistic

# See section 6.5 of Durbin, Distribution Theory for Tests Based on the 
# Sample Distribution Function, SIAM, Philadelphia (1973).


watson.u2 <- function( x, y )
{
    n <- length( x )
    m <- length( y )
    
    r <- c( sort( x ), sort( y ) )
    r.rank <- rank( r, ties.method = "average" )
    
    z <- ( r.rank[ 1:n ] - 1:n ) / m - ( 1:n - 1/2 ) / n 
    
    ( m / (n+m) ) * sum( ( z - mean( z ) )^2 ) + 
        ( m*(m+2*n) ) / ( 12*n*m*(n+m) ) 
}

