# normalize_svcd_aux.R
# 
# Copyright (c) 2017, Aarhus University (Denmark)
# 
# Software written by Carlos P. Roca, as research funded by the European Union.
# 
# This software may be modified and distributed under the terms of the MIT 
# license. See the LICENSE file for details.


# Auxiliary functions for standard-vector condition-decompositon ormalization


normalize.svcd.check.argument <- function( expression.data, 
    expression.condition, restrict.feature, search.h0.feature, 
    convergence.threshold, stdvec.graph, p.value.graph, verbose )
{
    # check expression.data
    if ( ! ( is.matrix( expression.data ) && is.numeric( expression.data ) ) )
        stop( "normalize.svcd: expression.data must be a numeric matrix" )
    
    # check expression.condition
    if ( ! ( ( ( is.atomic.vector( expression.condition ) && 
                ( is.character( expression.condition ) || 
                    is.numeric( expression.condition ) ) ) || 
            is.factor( expression.condition ) ) && 
        length( expression.condition ) == ncol( expression.data ) ) )
        stop( "normalize.svcd: bad expression.condition, see help" )
    
    # check restrict.feature
    if ( ! ( is.null( restrict.feature ) || 
        ( is.atomic.vector( restrict.feature ) && 
            ( is.character( expression.condition ) || 
                is.numeric( expression.condition ) ) ) ) )
        stop( "normalize.svcd: bad restrict.feature, see help" )
    
    # check search.h0.feature
    if ( ! ( is.atomic.singleton( search.h0.feature ) && 
        is.logical( search.h0.feature ) ) )
        stop( "normalize.svcd: search.h0.feature must be a logical value" )
    
    # check convergence.threshold
    if ( ! ( is.atomic.vector( convergence.threshold ) && 
        is.numeric( convergence.threshold ) && 
        length( convergence.threshold ) == 4 ) )
        stop( "normalize.svcd: bad convergence.threshold, see help" )
    
    # check stdvec.graph
    if ( ! ( is.null( stdvec.graph ) || 
        ( is.atomic.singleton( stdvec.graph ) && 
            is.character( stdvec.graph ) ) ) )
        stop( "normalize.svcd: stdvec.graph must be NULL or a character string" )

    # check p.value.graph
    if ( ! ( is.null( p.value.graph ) || 
        ( is.atomic.singleton( p.value.graph ) && 
            is.character( p.value.graph ) ) ) )
        stop( "normalize.svcd: p.value.graph must be NULL or a character string" )
    
    # check verbose
    if ( ! ( is.atomic.singleton( verbose ) && is.logical( verbose ) ) )
        stop( "normalize.svcd: verbose must be a logical value" )
}


normalize.stdvec.within.condition <- function( edata, condition, 
    norm.feature.idx, convergence.threshold, stdvec.graph, verbose )
{
    # calculate normalization
    norm.stdvec.result <- normalize.standard.vector( 
        edata[ norm.feature.idx, ], condition, FALSE, NULL, NULL, NULL, 
        convergence.threshold, stdvec.graph, NULL, verbose )
    
    # normalize with obtained offset
    edata <- sweep( edata, 2, norm.stdvec.result$offset )
    
    list( data = edata, 
        offset = norm.stdvec.result$offset, 
        convergence = norm.stdvec.result$convergence )
}


normalize.stdvec.between.condition <- function( edata, norm.cond, 
    norm.sample.cond, norm.feature.idx, search.h0.feature, 
    convergence.threshold, stdvec.graph, p.value.graph, verbose )
{
    # identify samples available per condition
    within.cond.n <- as.vector( table( norm.sample.cond )[ norm.cond ] )
    names( within.cond.n ) <- norm.cond
    
    # calculate balanced within-condition means
    bal.mean.n <- min( within.cond.n )
    
    expr.bal.mean.data <- sapply( norm.cond, function( cond ) 
        if ( within.cond.n[ cond ] == bal.mean.n )
            rowMeans( edata[ norm.feature.idx, cond == norm.sample.cond ] )
        else  # within.cond.n[ cond ] > bal.mean.n
            rowMeans( t( apply( 
                t( edata[ norm.feature.idx, cond == norm.sample.cond ] ), 
                2, function( ed ) sample( ed, bal.mean.n ) ) ) )
    )
    
    if ( search.h0.feature )
    {
        # calculate normalization while looking for h0 features
        
        # calculate within-condition means for F-statistic
        expr.mean.data <- sapply( norm.cond, function( cond ) 
            rowMeans( edata[ norm.feature.idx, cond == norm.sample.cond ] ) )
        
        # calculate within-condition variances for F-statistic
        within.cond.var <- sweep( sapply( norm.cond, function( cond ) 
            apply( edata[ norm.feature.idx, cond == norm.sample.cond ], 1, 
                var ) ), 
            2, within.cond.n - 1, "*" )
        within.cond.var <- rowSums( within.cond.var ) / 
            ( sum( within.cond.n ) - length( norm.cond ) )
        
        norm.stdvec.search.h0.feature.result <- normalize.standard.vector( 
            expr.bal.mean.data, "between.condition.search.h0.feature", 
            TRUE, expr.mean.data, within.cond.var, within.cond.n, 
            convergence.threshold, stdvec.graph, p.value.graph, verbose )
        
        h0.feature <- norm.stdvec.search.h0.feature.result$h0.feature
        h0.feature.convergence <- 
            norm.stdvec.search.h0.feature.result$convergence
        
        norm.feature <- h0.feature
    }
    else
    {
        # use all features given by norm.feature.idx for normalization
        norm.feature <- rownames( expr.bal.mean.data )
        
        h0.feature <- NULL
        h0.feature.convergence <- NULL
    }
    
    # calculate normalization
    norm.stdvec.result <- normalize.standard.vector( 
        expr.bal.mean.data[ norm.feature, ], "between.condition", FALSE, 
        NULL, NULL, NULL, convergence.threshold, stdvec.graph, NULL, verbose )
    
    # normalize each condition with obtained offset
    for( cond in norm.cond )
        edata[ , cond == norm.sample.cond ] <- 
            edata[ , cond == norm.sample.cond ] - 
            norm.stdvec.result$offset[ cond ]
    
    list( data = edata, 
        offset = norm.stdvec.result$offset, 
        convergence = norm.stdvec.result$convergence, 
        h0.feature = h0.feature, 
        h0.feature.convergence = h0.feature.convergence )
}


normalize.standard.vector <- function( edata, condition, search.h0.feature, 
    edata.fstat, within.cond.var, within.cond.n, convergence.threshold, 
    stdvec.graph, p.value.graph, verbose )
{
    if ( ! search.h0.feature )
    {
        iter.max <- 100
        single.threshold <- convergence.threshold[ 1 ]
        accum.threshold <- convergence.threshold[ 2 ]
        offset.accum.step.threshold <- 10
    }
    else
    {
        iter.max <- 200
        single.threshold <- convergence.threshold[ 3 ]
        accum.threshold <- convergence.threshold[ 4 ]
        offset.accum.step.threshold <- 10
        common.h0.feature.accum.step.max <- 10
    }
    
    stdvec.graph.feature.num <- 10000
    
    if ( verbose )
        cat( paste0( condition, "\n" ) )
    
    stdvec.offset <- rep( 0, ncol( edata ) )
    
    norm.stdvec.offset <- matrix( nrow=0, ncol = ncol( edata ) )
    norm.stdvec.offset.sd <- vector( "numeric" )
    norm.stdvec.offset.stderr <- vector( "numeric" )
    norm.stdvec.offset.delta.sd <- vector( "numeric" )
    norm.stdvec.offset.accum.step <- vector( "numeric" )
    norm.stdvec.numerical.demand <- vector( "numeric" )
    
    if ( ncol( edata ) >= 3 )
        norm.stdvec.watson.u2 <- matrix( nrow=0, 
            ncol = ( ncol( edata ) - 1 ) %/% 3 + 1 )
    else
        norm.stdvec.watson.u2 <- NULL
    
    if ( search.h0.feature )
    {
        norm.stdvec.common.h0.feature.accum.step <- vector( "numeric" )
        norm.stdvec.h0.feature.num <- vector( "numeric" )
        norm.stdvec.common.h0.feature.num <- vector( "numeric" )
    }
    
    if ( ! is.null( stdvec.graph ) )
    {
        # select a sample of features for plotting standardized sample vectors
        edata.feature <- rownames( edata )
        edata.feature.num <- length( edata.feature )
        
        if ( edata.feature.num <= stdvec.graph.feature.num )
            stdvec.graph.feature <- edata.feature
        else
            stdvec.graph.feature <- 
                sample( edata.feature, stdvec.graph.feature.num )
    }
    else
        stdvec.graph.feature <- NULL
    
    last.norm.data <- edata
    stdvec.offset.accum.step <- 0
    
    if ( search.h0.feature )
    {
        last.norm.data.fstat <- edata.fstat
        stdvec.common.h0.feature.accum.step <- 0
    }
    
    iter <- 0
    overall.convergence <- FALSE
    
    while ( iter < iter.max && ! overall.convergence )
    {
        iter <- iter + 1
        
        # obtain next step of standard vector offset
        if ( ! search.h0.feature )
            stdvec.offset.step <- calculate.stdvec.offset( last.norm.data, 
                condition, FALSE, NULL, NULL, NULL, stdvec.graph, 
                stdvec.graph.feature, NULL, iter )
        else
            stdvec.offset.step <- calculate.stdvec.offset( last.norm.data, 
                condition, TRUE, last.norm.data.fstat, within.cond.var, 
                within.cond.n, stdvec.graph, stdvec.graph.feature, 
                p.value.graph, iter )
        
        stdvec.offset.delta <- stdvec.offset.step$value
        stdvec.offset.stderr <- stdvec.offset.step$stderr
        stdvec.numerical.demand <- stdvec.offset.step$numerical.demand
        stdvec.watson.u2 <- stdvec.offset.step$watson.u2
        
        if ( search.h0.feature )
        {
            stdvec.h0.feature <- stdvec.offset.step$h0.feature
            stdvec.h0.feature.num <- length( stdvec.h0.feature )
        }
        
        # check errors
        if ( any( is.nan( stdvec.offset.delta ) ) || 
                is.nan( stdvec.offset.stderr ) )
            stop( "normalize.svcd: NaN error" )
        
        if ( stdvec.numerical.demand < .Machine$double.eps * 10^3 )
            stop( "normalize.svcd: Numerical error" )
        
        # update total standard vector offset
        stdvec.offset <- stdvec.offset + stdvec.offset.delta
        stdvec.offset <- stdvec.offset - mean( stdvec.offset )
        
        # update data at once with total offset
        last.norm.data <- sweep( edata, 2, stdvec.offset )
        if ( search.h0.feature )
            last.norm.data.fstat <- sweep( edata.fstat, 2, stdvec.offset )
        
        # check convergence
        stdvec.offset.sd <- sd( stdvec.offset )
        stdvec.offset.delta.sd <- sd( stdvec.offset.delta )
        
        stdvec.offset.stderr.ratio <- ifelse( stdvec.offset.sd > 0, 
            stdvec.offset.stderr / stdvec.offset.sd, 1 )
        stdvec.offset.delta.sd.ratio <- ifelse( stdvec.offset.sd > 0, 
            stdvec.offset.delta.sd / stdvec.offset.sd, 1 )
        
        if ( stdvec.offset.delta.sd == 0 )
            stdvec.offset.ratio <- 0
        else if ( stdvec.offset.stderr == 0 )
            stdvec.offset.ratio <- 1
        else
            stdvec.offset.ratio <- stdvec.offset.delta.sd / stdvec.offset.stderr
        
        stdvec.offset.accum.step <- ifelse( 
            stdvec.offset.ratio < accum.threshold, 
            stdvec.offset.accum.step + 1, 0 )
        
        stdvec.convergence <- stdvec.offset.ratio < single.threshold || 
            stdvec.offset.accum.step > offset.accum.step.threshold
        
        if ( ! search.h0.feature )
            overall.convergence <- stdvec.convergence
        else
        {
            stdvec.common.h0.feature.accum.step <- ifelse( stdvec.convergence, 
                stdvec.common.h0.feature.accum.step + 1, 0 )
            
            overall.convergence <- stdvec.common.h0.feature.accum.step >= 
                common.h0.feature.accum.step.max
        }
        
        # store last results
        last.stdvec.offset <- stdvec.offset
        
        if ( search.h0.feature )
        {
            if ( stdvec.common.h0.feature.accum.step == 0 )
                stdvec.common.h0.feature <- NULL
            else if ( stdvec.common.h0.feature.accum.step == 1 )
                stdvec.common.h0.feature <- stdvec.h0.feature
            else
                stdvec.common.h0.feature <- intersect( stdvec.h0.feature, 
                    stdvec.common.h0.feature )
            
            stdvec.common.h0.feature.num <- ifelse( 
                is.null( stdvec.common.h0.feature ), NA, 
                length( stdvec.common.h0.feature ) )
        }
        
        # store step results
        norm.stdvec.offset <- rbind( norm.stdvec.offset, stdvec.offset )
        norm.stdvec.offset.sd <- c( norm.stdvec.offset.sd, stdvec.offset.sd )
        norm.stdvec.offset.stderr <- c( norm.stdvec.offset.stderr, 
            stdvec.offset.stderr )
        norm.stdvec.offset.delta.sd <- c( norm.stdvec.offset.delta.sd, 
            stdvec.offset.delta.sd )
        norm.stdvec.offset.accum.step <- c( norm.stdvec.offset.accum.step, 
            stdvec.offset.accum.step )
        norm.stdvec.numerical.demand <- c( norm.stdvec.numerical.demand, 
            stdvec.numerical.demand )
        norm.stdvec.watson.u2 <- rbind( norm.stdvec.watson.u2, 
            stdvec.watson.u2 )
        
        if ( search.h0.feature )
        {
            norm.stdvec.common.h0.feature.accum.step <- 
                c( norm.stdvec.common.h0.feature.accum.step, 
                    stdvec.common.h0.feature.accum.step )
            norm.stdvec.h0.feature.num <- c( norm.stdvec.h0.feature.num, 
                stdvec.h0.feature.num )
            norm.stdvec.common.h0.feature.num <- 
                c( norm.stdvec.common.h0.feature.num, 
                    stdvec.common.h0.feature.num )
        }
        
        if ( verbose )
        {
            if ( ! is.null( stdvec.watson.u2 ) )
                stdvec.watson.u2.char <- paste0( signif( stdvec.watson.u2, 6 ), 
                    collapse=" " )
            else
                stdvec.watson.u2.char <- ""
            
            cat( sprintf( "  %2d %g %g %g %02d %g [%s]", iter, 
                stdvec.offset.sd, stdvec.offset.stderr.ratio, 
                stdvec.offset.delta.sd.ratio, stdvec.offset.accum.step, 
                stdvec.numerical.demand, stdvec.watson.u2.char ) )
            
            if ( search.h0.feature )
            {
                cat( sprintf( " %02d %d", stdvec.common.h0.feature.accum.step, 
                    stdvec.h0.feature.num ) )
                
                if ( ! is.na( stdvec.common.h0.feature.num ) )
                    cat( paste0( " ", stdvec.common.h0.feature.num ) )
            }
            
            cat( "\n" )
        }
    }
    
    if ( verbose )
        cat( "\n" )
    
    if ( ! overall.convergence )
        stop( "normalize.svcd: No convergence" )
    
    # remove sample or condition names from step results
    dimnames( norm.stdvec.offset ) <- NULL
    if ( ! is.null( norm.stdvec.watson.u2 ) )
        dimnames( norm.stdvec.watson.u2 ) <- NULL
    
    norm.stdvec.convergence <- list( offset = norm.stdvec.offset, 
        offset.sd = norm.stdvec.offset.sd, 
        offset.stderr = norm.stdvec.offset.stderr, 
        offset.delta.sd = norm.stdvec.offset.delta.sd, 
        offset.accum.step = norm.stdvec.offset.accum.step, 
        numerical.demand = norm.stdvec.numerical.demand, 
        watson.u2 = norm.stdvec.watson.u2 )
    
    if ( search.h0.feature )
    {
        norm.stdvec.convergence$common.h0.feature.accum.step <- 
            norm.stdvec.common.h0.feature.accum.step
        norm.stdvec.convergence$h0.feature.num <- 
            norm.stdvec.h0.feature.num
        norm.stdvec.convergence$common.h0.feature.num <- 
            norm.stdvec.common.h0.feature.num
        
        norm.stdvec.h0.feature <- stdvec.common.h0.feature
    }
    else
        norm.stdvec.h0.feature <- NULL

    list( offset = last.stdvec.offset, 
        convergence = norm.stdvec.convergence, 
        h0.feature = norm.stdvec.h0.feature )
}


calculate.stdvec.offset <- function( edata, condition, search.h0.feature, 
    edata.fstat, within.cond.var, within.cond.n, stdvec.graph, 
    stdvec.graph.feature, p.value.graph, iter )
{
    stdvec.trim <- 0.01
    ks.test.alpha <- 1e-3
    
    if ( ! search.h0.feature )
    {
        # identify features for normalization
        expr.var <- apply( edata, 1, var )
        expr.var[ expr.var < max( expr.var ) * .Machine$double.eps  ] <- NA
        
        stdvec.feature <- names( which( 
            expr.var > quantile( expr.var, stdvec.trim/2, na.rm=TRUE ) & 
            expr.var < quantile( expr.var, 1 - stdvec.trim/2, na.rm=TRUE ) ) )
        
        h0.feature <- NULL
    }
    else
    {
        # calculate f statistics for each feature
        expr.k <- length( within.cond.n )
        expr.n <- sum( within.cond.n )
        
        expr.grand.mean <- apply( edata.fstat, 1, function( ef ) 
            sum( ef * within.cond.n ) ) / expr.n
        
        between.cond.var <- apply( ( edata.fstat - expr.grand.mean )^2, 1, 
            function( ef2 ) sum( ef2 * within.cond.n ) ) / ( expr.k - 1 )
        
        expr.f <- between.cond.var / within.cond.var
        
        expr.f <- na.omit( expr.f )  # in case of 0/0
        attr( expr.f, "na.action" ) <- NULL
        
        expr.p.value <- pf( expr.f, df1 = expr.k-1, df2 = expr.n-expr.k, 
            lower.tail = FALSE )
        
        # identify h0 features with one-sided up Kolmogorov-Smirnov test
        ks.test.d <- sqrt( - log( ks.test.alpha ) / 2 )
        
        epv <- sort( expr.p.value )
        epv.n <- length( expr.p.value )
        
        epv.i <- 1
        ks.test.D.up <- 1:epv.n / epv.n - epv
        ks.test.reject <- any( ks.test.D.up > ks.test.d / sqrt( epv.n ) )
        
        while ( ks.test.reject )
        {
            epv.i <- epv.i + 1
            
            ks.test.D.up <- 
                ( epv.i : epv.n - epv.i + 1 ) / ( epv.n - epv.i + 1 ) - 
                ( epv[ epv.i : epv.n ] - epv[ epv.i - 1 ] ) / 
                    ( 1 - epv[ epv.i - 1 ] )
            
            ks.test.reject <- 
                any( ks.test.D.up > ks.test.d / sqrt( epv.n - epv.i + 1 ) )
        }
        
        epv.h0.i <- epv.i
        epv.h0.n <- epv.n - epv.i + 1
        epv.h0.p <-  ifelse( epv.i == 1, 0, epv[ epv.i - 1 ] )
        epv.h0.q <- ( epv.i - 1 ) / epv.n
        
        h0.feature.idx <- order( expr.p.value, decreasing=TRUE )[ 1 : epv.h0.n ]
        h0.feature <- names( expr.p.value )[ h0.feature.idx ]
        
        pi0.est <- ( 1 - epv.h0.q ) / ( 1 - epv.h0.p )
        
        # plot graph of p-values
        if ( ! is.null( p.value.graph ) )
        {
            if ( p.value.graph == "" )
                p.value.graph <- "."
            else if ( ! file.exists( p.value.graph ) )
                dir.create( p.value.graph, recursive = TRUE )
            
            png.filename <- sprintf( "%s/p_value_iter%02d.png", 
                p.value.graph, iter )
            
            png( png.filename, width = 1280, height = 720 )
            
            par( mfrow = c(1,2), pty = "s", mar = c( 3.4, 4.0, 0, 1.6 ), 
                oma = c( 0, 4.2, 3.4, 0 ), mgp = c( 4.0, 1.2, 0 ) )
            
            epv.quant <- 1:epv.n / epv.n
            epv.h0.idx <- epv.h0.i : epv.n
            
            epv.x0 <- epv.h0.p
            epv.y0 <- epv.h0.q
            epv.yd <- ks.test.d * ( sqrt( epv.h0.n ) / epv.n )
            
            xylim <- list( c(0,0), c( epv.h0.p, epv.h0.q ) )
            
            for ( i in 1:2  )
            {
                plot( 0, type = "n", 
                    xlim = c( xylim[[i]][1], 1 ), ylim = c( xylim[[i]][2], 1 ), 
                    xlab = "p-value", ylab = "", cex.axis = 2.5, cex.lab = 2.5 )
                
                segments( x0 = c( 0, epv[ - epv.h0.idx ] ), 
                    y0 = c( 0, epv.quant[ - epv.h0.idx ] ), 
                    x1 = c( epv[ - epv.h0.idx ], epv[ epv.h0.idx ][ 1 ] ), 
                    lwd = 2.5, col = "black" )
                segments( x0 = epv[ epv.h0.idx ], y0 = epv.quant[ epv.h0.idx ], 
                    x1 = c( epv[ epv.h0.idx ][ -1 ], 1 ), lwd = 2.5, 
                    col = "red" )
                points( epv.h0.p, epv.h0.q, pch = 20, cex = 3 )
                
                segments( 0, 1 - pi0.est, 1, 1, lwd = 2, lty = 2, col = "blue" )
                segments( 0, 1 - pi0.est + epv.yd, 1, 1 + epv.yd, lwd = 2, 
                    lty = 3, col = "blue" )
            }
            
            mtext( "F( p-value )", side = 2, line = 1.0, outer = TRUE, 
                cex = 2.5 )
            
            graph.title <- substitute( 
                paste( "iter=", iter, "    -    #", H[0], "=", h0.n ), 
                list( iter = sprintf( "%02d", iter ), 
                    h0.n = sprintf( "%5d", epv.h0.n ) ) )
            mtext( graph.title, 3, line = -1.2, outer = TRUE, cex = 2.7 )
            
            dev.off()
        }
        
        # identify features for normalization
        expr.var <- apply( edata.fstat[ h0.feature, ], 1, var )
        expr.var[ expr.var < max( expr.var ) * .Machine$double.eps  ] <- NA
        
        stdvec.feature <- h0.feature[ 
            expr.var > quantile( expr.var, stdvec.trim/2, na.rm=TRUE ) & 
            expr.var < quantile( expr.var, 1 - stdvec.trim/2, na.rm=TRUE ) ]
    }
    
    # center and scale expression data
    expr.mean <- rowMeans( edata[ stdvec.feature, ] )
    expr.centered <- sweep( edata[ stdvec.feature, ], 1, expr.mean, "-" )
    
    expr.sd.inv <- 1 / apply( expr.centered, 1, sd )
    expr.scaled <- sweep( expr.centered, 1, expr.sd.inv, "*" )
    
    # calculate offset
    expr.sd.inv.sum <- sum( expr.sd.inv )
    stdvec.offset <- apply( expr.scaled, 2, sum ) / expr.sd.inv.sum
    stdvec.offset <- stdvec.offset - mean( stdvec.offset )
    
    # estimate error and numerical demand
    stdvec.offset.stderr <- sqrt( length( stdvec.feature ) ) / expr.sd.inv.sum
    expr.sd.inv.min <- min( expr.sd.inv )
    stdvec.numerical.demand <- expr.sd.inv.min / expr.sd.inv.sum
    
    # calculate and plot density distribution of standard vector angles
    
    dimension.num <- ncol( expr.scaled )
    
    if ( dimension.num < 3 )
        theta.watson.u2 <- NULL
    else
    {
        # identify condition groups
        dimension.group.num <- ( dimension.num - 1 ) %/% 3 + 1
        
        dimension.group <- lapply( 1 : dimension.group.num, function ( g ) 
            if ( g < dimension.group.num )
                ( 3*g - 2 ) : ( 3*g )
            else
                ( dimension.num - 2 ) : dimension.num )
        
        theta.watson.u2 <- vector( "numeric", dimension.group.num )
        
        uv <- matrix( c( 0, -1/sqrt(2), 1/sqrt(2), 
            2/sqrt(6), -1/sqrt(6), -1/sqrt(6) ), nrow = 3 )
        
        for ( dim.group.idx in 1:dimension.group.num )
        {
            # select expression values for each condition group
            expr.dim.group <- expr.scaled[ , 
                dimension.group[[ dim.group.idx ]] ]
            
            if ( dimension.group.num > 1 )
            {
                # re-standardize again for this group
                expr.dim.group.mean <- rowMeans( expr.dim.group )
                expr.dim.group.centered <- sweep( expr.dim.group, 1, 
                    expr.dim.group.mean, "-" )
                
                expr.dim.group.sd <- apply( expr.dim.group.centered, 1, sd )
                expr.dim.group.sel <- expr.dim.group.sd != 0
                expr.dim.group <- sweep( 
                    expr.dim.group.centered[ expr.dim.group.sel, ], 1, 
                    expr.dim.group.sd[ expr.dim.group.sel ], "/" )
            }
            
            expr.uv <- expr.dim.group %*% uv
            expr.u <- expr.uv[ , 1 ]
            expr.v <- expr.uv[ , 2 ]
            
            # calculate density distribution of angles
            theta.density.n <- 2^11
            theta.density.adjust <- 0.5
            
            expr.theta <- atan2( expr.v, expr.u )
            expr.theta.ex <- c( expr.theta, 
                expr.theta + ifelse( expr.theta > 0, -2*pi, 2*pi ) )
            
            expr.theta.density <- density( expr.theta.ex, 
                adjust = theta.density.adjust, n = theta.density.n )
            expr.theta.density.sel <- expr.theta.density$x > -pi & 
                expr.theta.density$x <= pi
            
            # calculate density distribution of angles after permutations
            expr.theta.permu <- cbind( expr.theta.ex, 
                expr.theta.ex + (2*pi)/3, 
                expr.theta.ex - (2*pi)/3, 
                - expr.theta.ex + pi, 
                - expr.theta.ex + pi/3, 
                - expr.theta.ex - pi/3 )
            
            invar.theta <- expr.theta.permu[ expr.theta.permu > -pi & 
                expr.theta.permu <= pi ]
            invar.theta.ex <- c( invar.theta, 
                invar.theta + ifelse( invar.theta > 0, -2*pi, 2*pi ) )
            
            invar.theta.density <- density( invar.theta.ex, 
                adjust = theta.density.adjust, n = theta.density.n )
            invar.theta.density.sel <- invar.theta.density$x > -pi & 
                invar.theta.density$x <= pi
            
            # calculate Watson U2 statistic
            theta.watson.u2[ dim.group.idx ] <- 
                watson.u2( expr.theta, sample( invar.theta, 
                    length( expr.theta ), replace = TRUE ) )
            
            if ( ! is.null( stdvec.graph ) && 
                requireNamespace( "plotrix", quietly = TRUE ) )
            {
                if ( stdvec.graph == "" )
                    stdvec.graph <- "."
                else if ( ! file.exists( stdvec.graph ) )
                    dir.create( stdvec.graph, recursive = TRUE )
                
                png.filename <- sprintf( 
                    "%s/stdvec_%s%s_iter%02d.png", stdvec.graph, condition, 
                    ifelse( dimension.group.num > 1, 
                        sprintf( "_dg%02d", dim.group.idx ), "" ), 
                    iter )
                
                png( png.filename, width = 1280, height = 720 )
                
                par( mfrow = c(1,2), pty = "s", xpd = FALSE, cex.lab = 2, 
                    mar = c( 1, 0.5, 2.2, 0.5 ), oma = c( 0, 0, 0, 0 ) )
                
                # select offset values for condition group
                stdvec.offset.uv <- stdvec.offset[ 
                    dimension.group[[ dim.group.idx ]] ] %*% uv
                stdvec.offset.u <- stdvec.offset.uv[ , 1 ]
                stdvec.offset.v <- stdvec.offset.uv[ , 2 ]
                
                # plot a sample of standardized sample vectors
                uv.lim <- c( -1.5, 1.5 )
                plot( 0, type = "n", xlim = uv.lim, ylim = uv.lim, 
                    axes = FALSE, ann = FALSE, frame.plot = FALSE, asp = 1 )
                
                expr.uv.feature <- names( expr.u )
                if ( length( expr.uv.feature ) > 
                    length( stdvec.graph.feature ) )
                    expr.uv.feature <- intersect( expr.uv.feature, 
                        stdvec.graph.feature )
                expr.uv.feature.num <- length( expr.uv.feature )
                
                expr.uv.factor <- 1.06
                expr.uv.color <- gray( 0.3 )
                expr.uv.width <- ifelse( expr.uv.feature.num > 1000, 0.1, 0.2 )
                segments( 0, 0, expr.uv.factor * expr.u[ expr.uv.feature ], 
                    expr.uv.factor * expr.v[ expr.uv.feature ], 
                    lwd = expr.uv.width, col = expr.uv.color )
                
                grid.pos.x <- c( 0, -sqrt(3/4), sqrt(3/4) )
                grid.pos.y <- c( 1, -1/2, -1/2 )
                
                grid.length <- expr.uv.factor * sqrt( 2 )
                segments( 0, 0, grid.length * grid.pos.x, 
                    grid.length * grid.pos.y, lwd = 2.5, col = "black" )
                
                stdvec.offset.uv.factor <- 10
                stdvec.offset.uv.color <- "red"
                segments( 0, 0, stdvec.offset.uv.factor * stdvec.offset.u, 
                    stdvec.offset.uv.factor * stdvec.offset.v, lwd = 3, 
                    col = "red" )
                
                par( xpd = TRUE )
                
                grid.label.length <- 1.63
                grid.labels <- c( "s1", "s2", "s3" )
                text( grid.label.length * grid.pos.x, 
                    grid.label.length * grid.pos.y, grid.labels, cex = 2.5 )
                
                par( xpd = FALSE )
                
                stdvec.offset.mag <- sqrt( sum( stdvec.offset[ 
                    dimension.group[[ dim.group.idx ]] ]^2 ) )
                mtext( paste0( "||offset|| = ", sprintf( "%.3e", 
                    stdvec.offset.mag ) ), side = 1, line = 0.7, cex = 2.5 )
                
                # plot polar distributions of standard vector angles
                polar.expr.theta <- expr.theta.density$x[ 
                    expr.theta.density.sel ]
                polar.expr.rho <- 2 * expr.theta.density$y[ 
                    expr.theta.density.sel ]
                
                polar.invar.theta <- invar.theta.density$x[ 
                    invar.theta.density.sel ]
                polar.invar.rho <- 2 * invar.theta.density$y[ 
                    invar.theta.density.sel ]
                
                theta.labels <- c( "", "", "" )
                rho.grid <- seq( 0, 3/(2*pi), length.out=4 )
                rho.labels <- c( "", "", expression(1/pi), "" )
                
                plotrix::radial.plot( polar.expr.rho, polar.expr.theta - pi/2, 
                    rp.type = "p", start = pi/2, radial.lim = rho.grid, 
                    show.grid.labels = length( theta.labels ), 
                    labels = theta.labels, radial.labels = rho.labels, 
                    lwd = 2.5, mar = par( "mar" ) )
                
                plotrix::radial.plot( polar.invar.rho, 
                    polar.invar.theta - pi/2, rp.type = "p", start=pi/2, 
                    radial.lim = rho.grid, lwd = 2, lty = 2, line.col = "blue", 
                    add = TRUE )
                
                grid.label.length <- 0.52
                text( grid.label.length * grid.pos.x, 
                    grid.label.length * grid.pos.y, grid.labels, cex = 2.5 )
                
                mtext( substitute( paste( "Watson U"^"2", " = ", wu2 ), 
                    list( wu2 = sprintf( "%.3e", 
                        theta.watson.u2[ dim.group.idx ] ) ) ), 
                    side = 1, line = 0.8, cex = 2.5 )
                
                graph.title <- sprintf( "%s%s    -    iter=%02d", condition, 
                    ifelse( dimension.group.num > 1, 
                        sprintf( ":%02d", dim.group.idx ), "" ), 
                    iter )
                title( main = graph.title, outer = TRUE, line = -3, 
                    font.main = 1, cex.main = 2.7 )
                
                dev.off()
            }
        }
    }
    
    list( value = stdvec.offset, 
        stderr = stdvec.offset.stderr, 
        numerical.demand = stdvec.numerical.demand, 
        watson.u2 = theta.watson.u2, 
        h0.feature = h0.feature )
}

