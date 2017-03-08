# normalize_mediancd_aux.R
# 
# Copyright (c) 2017, Aarhus University (Denmark)
# 
# Software written by Carlos P. Roca, as research funded by the European Union.
# 
# This software may be modified and distributed under the terms of the MIT 
# license. See the LICENSE file for details.


# Auxiliary fuctions for median condition-decomposition normalization


normalize.mediancd.check.argument <- function( expression.data, 
    expression.condition, normalization.probability, restrict.feature, 
    search.h0.feature, convergence.threshold, p.value.graph, verbose )
{
    # check expression.data
    if ( ! ( is.matrix( expression.data ) && is.numeric( expression.data ) ) )
        stop( "normalize.mediancd: expression.data must be a numeric matrix" )
    
    # check expression.condition
    if ( ! ( ( ( is.atomic.vector( expression.condition ) && 
                ( is.character( expression.condition ) || 
                    is.numeric( expression.condition ) ) ) || 
            is.factor( expression.condition ) ) && 
        length( expression.condition ) == ncol( expression.data ) ) )
        stop( "normalize.mediancd: bad expression.condition, see help" )
    
    # check normalization.probability
    if ( ! ( is.atomic.singleton( normalization.probability ) && 
        is.numeric( normalization.probability ) && 
        0 <= normalization.probability && normalization.probability <= 1 ) )
        stop( "normalize.mediancd: bad normalization.probability, see help" )
    
    # check restrict.feature
    if ( ! ( is.null( restrict.feature ) || 
        ( is.atomic.vector( restrict.feature ) && 
            ( is.character( expression.condition ) || 
                is.numeric( expression.condition ) ) ) ) )
        stop( "normalize.mediancd: bad restrict.feature, see help" )
    
    # check search.h0.feature
    if ( ! ( is.atomic.singleton( search.h0.feature ) && 
        is.logical( search.h0.feature ) ) )
        stop( "normalize.mediancd: search.h0.feature must be a logical value" )
    
    # check convergence.threshold
    if ( ! ( is.atomic.vector( convergence.threshold ) && 
        is.numeric( convergence.threshold ) && 
        length( convergence.threshold ) == 2 ) )
        stop( "normalize.mediancd: bad convergence.threshold, see help" )
    
    # check p.value.graph
    if ( ! ( is.null( p.value.graph ) || 
        ( is.atomic.singleton( p.value.graph ) && 
            is.character( p.value.graph ) ) ) )
        stop( "normalize.mediancd: p.value.graph must be NULL or a character string" )
    
    # check verbose
    if ( ! ( is.atomic.singleton( verbose ) && is.logical( verbose ) ) )
        stop( "normalize.mediancd: verbose must be a logical value" )
}


normalize.median.within.condition <- function( edata, condition, norm.prob, 
    norm.feature.idx, verbose )
{
    if ( verbose )
        cat( paste0( condition, "\n" ) )
    
    # calculate normalization
    norm.median.offset <- apply( edata[ norm.feature.idx, ], 2, quantile, 
        probs = norm.prob )
    norm.median.offset <- norm.median.offset - mean( norm.median.offset )
    
    # normalize with obtained offset
    edata <- sweep( edata, 2, norm.median.offset )
    
    norm.within.cond.result <- list( data = edata, offset = norm.median.offset )
    
    norm.within.cond.result
}


normalize.median.between.condition <- function( edata, norm.cond, 
    norm.sample.cond, norm.prob, norm.feature.idx, search.h0.feature, 
    convergence.threshold, p.value.graph, verbose )
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
        
        if ( verbose )
            cat( "between.condition.search.h0.feature\n" )
        
        norm.median.search.h0.feature.result <- normalize.median.selection( 
            expr.bal.mean.data, expr.mean.data, within.cond.var, within.cond.n, 
            norm.prob, convergence.threshold, p.value.graph, verbose )
        
        h0.feature <- norm.median.search.h0.feature.result$h0.feature
        h0.feature.convergence <- 
            norm.median.search.h0.feature.result$convergence
        
        norm.feature <- h0.feature
    }
    else
    {
        # use all features given by norm.feature.idx for normalization
        norm.feature <- rownames( expr.bal.mean.data )
        
        h0.feature <- NULL
        h0.feature.convergence <- NULL
    }
    
    if ( verbose )
        cat( "between.condition\n" )
    
    # calculate normalization
    norm.median.offset <- apply( expr.bal.mean.data[ norm.feature, ], 2, 
        quantile, probs=norm.prob )
    norm.median.offset <- norm.median.offset - mean( norm.median.offset )
    
    # normalize each condition with obtained offset
    for( cond in norm.cond )
        edata[ , cond == norm.sample.cond ] <- 
            edata[ , cond == norm.sample.cond ] - norm.median.offset[ cond ]
    
    list( data = edata, 
        offset = norm.median.offset, 
        h0.feature = h0.feature, 
        h0.feature.convergence = h0.feature.convergence )
}


normalize.median.selection <- function( edata, edata.fstat, within.cond.var, 
    within.cond.n, norm.prob, convergence.threshold, p.value.graph, verbose )
{
    iter.max <- 200
    single.threshold <- convergence.threshold[ 1 ]
    accum.threshold <- convergence.threshold[ 2 ]
    offset.accum.step.threshold <- 10
    common.h0.feature.accum.step.max <- 10
    
    median.offset <- rep( 0, ncol( edata ) )
    
    norm.median.offset <- matrix( nrow=0, ncol = ncol( edata ) )
    norm.median.offset.sd <- vector( "numeric" )
    norm.median.offset.delta.sd <- vector( "numeric" )
    norm.median.offset.accum.step <- vector( "numeric" )
    norm.median.common.h0.feature.accum.step <- vector( "numeric" )
    norm.median.h0.feature.num <- vector( "numeric" )
    norm.median.common.h0.feature.num <- vector( "numeric" )
    
    last.norm.data <- edata
    last.norm.data.fstat <- edata.fstat
    
    median.offset.accum.step <- 0
    median.common.h0.feature.accum.step <- 0
    
    iter <- 0
    overall.convergence <- FALSE
    
    while ( iter < iter.max && ! overall.convergence )
    {
        iter <- iter + 1
        
        # obtain next step of median offset
        median.offset.step <- calculate.median.offset( last.norm.data, 
            last.norm.data.fstat, within.cond.var, within.cond.n, norm.prob, 
            p.value.graph, iter )
        
        median.offset.delta <- median.offset.step$value
        
        median.h0.feature <- median.offset.step$h0.feature
        median.h0.feature.num <- length( median.h0.feature )
        
        # check errors
        if ( any( is.nan( median.offset.delta ) ) )
            stop( "normalize.mediancd: NaN error" )
        
        # update total median offset
        median.offset <- median.offset + median.offset.delta
        median.offset <- median.offset - mean( median.offset )
        
        # update data at once with total offset
        last.norm.data <- sweep( edata, 2, median.offset )
        last.norm.data.fstat <- sweep( edata.fstat, 2, median.offset )
        
        # check convergence
        median.offset.sd <- sd( median.offset )
        median.offset.delta.sd <- sd( median.offset.delta )
        
        median.offset.delta.sd.ratio <- ifelse( median.offset.sd > 0, 
            median.offset.delta.sd / median.offset.sd, 1 )
        
        if ( median.offset.delta.sd == 0 )
            median.offset.ratio <- 0
        else if ( median.offset.sd == 0 )
            median.offset.ratio <- 1
        else
            median.offset.ratio <- median.offset.delta.sd / median.offset.sd
        
        median.offset.accum.step <- ifelse( 
            median.offset.ratio < accum.threshold, 
            median.offset.accum.step + 1, 0 )
        
        median.convergence <- median.offset.ratio < single.threshold || 
            median.offset.accum.step > offset.accum.step.threshold
        
        median.common.h0.feature.accum.step <- ifelse( median.convergence, 
            median.common.h0.feature.accum.step + 1, 0 )
        
        overall.convergence <- median.common.h0.feature.accum.step >= 
            common.h0.feature.accum.step.max
        
        # store last results
        last.median.offset <- median.offset
        
        if ( median.common.h0.feature.accum.step == 0 )
            median.common.h0.feature <- NULL
        else if ( median.common.h0.feature.accum.step == 1 )
            median.common.h0.feature <- median.h0.feature
        else
            median.common.h0.feature <- intersect( median.h0.feature, 
                median.common.h0.feature )
        
        median.common.h0.feature.num <- ifelse( 
            is.null( median.common.h0.feature ), NA, 
            length( median.common.h0.feature ) )
        
        # store step results
        norm.median.offset <- rbind( norm.median.offset, median.offset )
        norm.median.offset.sd <- c( norm.median.offset.sd, median.offset.sd )
        norm.median.offset.delta.sd <- c( norm.median.offset.delta.sd, 
            median.offset.delta.sd )
        norm.median.offset.accum.step <- c( norm.median.offset.accum.step, 
            median.offset.accum.step )
        norm.median.common.h0.feature.accum.step <- 
            c( norm.median.common.h0.feature.accum.step, 
                median.common.h0.feature.accum.step )
        norm.median.h0.feature.num <- c( norm.median.h0.feature.num, 
            median.h0.feature.num )
        norm.median.common.h0.feature.num <- 
            c( norm.median.common.h0.feature.num, median.common.h0.feature.num )
        
        if ( verbose )
        {
            cat( sprintf( "  %2d %g %g %02d %02d %d", iter, 
                median.offset.sd, median.offset.delta.sd.ratio, 
                median.offset.accum.step, median.common.h0.feature.accum.step, 
                median.h0.feature.num ) )
            
            if ( ! is.na( median.common.h0.feature.num ) )
                cat( paste0( " ", median.common.h0.feature.num ) )
            
            cat( "\n" )
        }
    }
    
    if ( ! overall.convergence )
        stop( "normalize.mediancd: No convergence" )
    
    # remove sample or condition names from step results
    dimnames( norm.median.offset ) <- NULL
    
    norm.median.convergence <- list( offset = norm.median.offset, 
        offset.sd = norm.median.offset.sd, 
        offset.delta.sd = norm.median.offset.delta.sd, 
        offset.accum.step = norm.median.offset.accum.step, 
        common.h0.feature.accum.step = 
            norm.median.common.h0.feature.accum.step, 
        h0.feature.num = norm.median.h0.feature.num, 
        common.h0.feature.num = norm.median.common.h0.feature.num )
    
    list( offset = last.median.offset, 
        convergence = norm.median.convergence, 
        h0.feature = median.common.h0.feature )
}


calculate.median.offset <- function( edata, edata.fstat, within.cond.var, 
    within.cond.n, norm.prob, p.value.graph, iter )
{
    ks.test.alpha <- 1e-3
    
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
                x1 = c( epv[ epv.h0.idx ][ -1 ], 1 ), lwd = 2.5, col = "red" )
            points( epv.h0.p, epv.h0.q, pch = 20, cex = 3 )
            
            segments( 0, 1 - pi0.est, 1, 1, lwd = 2, lty = 2, col = "blue" )
            segments( 0, 1 - pi0.est + epv.yd, 1, 1 + epv.yd, lwd = 2, lty = 3, 
                col = "blue" )
        }
        
        mtext( "F( p-value )", side = 2, line = 1.0, outer = TRUE, cex = 2.5 )
        
        graph.title <- substitute( 
            paste( "iter=", iter, "    -    #", H[0], "=", h0.n ), 
            list( iter = sprintf( "%02d", iter ), 
                h0.n = sprintf( "%5d", epv.h0.n ) ) )
        mtext( graph.title, 3, line = -1.2, outer = TRUE, cex = 2.7 )
        
        dev.off()
    }
    
    # calculate offset
    median.offset <- apply( edata[ h0.feature, ], 2, quantile, 
        probs = norm.prob )
    median.offset <- median.offset - mean( median.offset )
    
    list( value = median.offset, h0.feature = h0.feature )
}

