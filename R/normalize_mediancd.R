# normalize_mediancd.R
# 
# Copyright (c) 2017, Aarhus University (Denmark)
# 
# Software written by Carlos P. Roca, as research funded by the European Union.
# 
# This software may be modified and distributed under the terms of the MIT 
# license. See the LICENSE file for details.


#' MedianCD (Median Condition-Decomposition) Normalization
#' 
#' Normalizes gene expression data following the MedianCD algorithm. 
#' It also provides the estimated normalization factors and the no-variation 
#' features (e.g. genes) detected. 
#' 
#' Only features (e.g. genes) with no missing values in any sample are used 
#' for normalizing. 
#' First, the samples of each experimental condition are normalized separately. 
#' Then, conditions means are normalized while iteratively searching for 
#' no-variation features, that is, for a subset of features that do not show 
#' evidence of being differentially expressed accross conditions. 
#' After this, conditions means are actually normalized using only the detected 
#' no-variation features. 
#' For \eqn{C} experimental conditions, \eqn{C+2} normalizations are carried 
#' out, all of them using Median normalization. 
#' Finally, the normalization factors obtained in the normalizations 
#' \eqn{1, 2, \dots, C, C+2} are combined to obtain the overall normalization 
#' factors and normalize the expression data. 
#' See reference below for more details. 
#' 
#' If \code{expression.condition} indicates that all samples correspond to the 
#' same experimental condition, then only one Median normalization is 
#' performed, without searching for no-variation features. 
#' 
#' @keywords MedianCD Median Condition-Decomposition normalization gene 
#'     expression high-throughput assay
#' 
#' @param expression.data Numeric matrix with expression data. 
#'     Rows correspond to features, for example genes. 
#'     Columns correspond to samples. 
#'     If rows and/or columns do not have names, they are assigned names with 
#'     format \code{"feature.[number]"} and/or \code{"sample.[number]"}, 
#'     respectively. 
#' @param expression.condition Character or numeric vector defining the 
#'     experimental conditions. 
#'     It can also be a factor. 
#'     The length of \code{expression.condition} must be equal to the number of 
#'     columns of \code{expression.data}, that is, to the number of samples. 
#'     It can contain \code{NA} values, meaning samples to be ignored in the 
#'     normalization. 
#' @param normalization.probability Probabilty whose quantile is used to 
#'     normalize at. 
#'     It must be a numerical value in the interval \code{[0,1]}. 
#' @param restrict.feature When \code{restrict.feature} is \code{NULL}, all 
#'     features can be used for normalizing. 
#'     Otherwise, \code{restrict.feature} is expected to be a character or 
#'     numeric vector identifying the features, by name or index respectively, 
#'     that can be used for normalizing. 
#' @param search.h0.feature Logical value indicating whether no-variation 
#'     features should be searched for the final between-condition 
#'     normalization, thus restricting the set of features used in this 
#'     normalization. 
#'     When \code{search.h0.feature} is \code{FALSE}, the complete set of 
#'     features used in the within-condition normalizations is also used for 
#'     normalizing between conditions. 
#' @param convergence.threshold Numeric vector with two elements, defining 
#'     convergence parameters for the algorithm. 
#'     The format is \emph{(single.step.search.h0.feature, 
#'     multiple.step.search.h0.feature)}. 
#'     Using values different from the default ones is only advisable when the 
#'     implementation of convergence is understood in detail. 
#' @param p.value.graph When not \code{NULL}, it provides a character string 
#'     with the name of a directory to save graphs displaying the distributions 
#'     of p-values used for the identification of no-variation features. 
#'     The directory can be entered as an absolute or relative path, and it is 
#'     created if it does not exist. 
#'     To save in the current working directory, use \code{p.value.graph="."} 
#'     or \code{p.value.graph=""}. 
#' @param verbose Logical value indicating whether convergence information 
#'     should be printed to the console. 
#' 
#' @return List with the elements
#' \item{data}{Matrix with normalized data.}
#' \item{offset}{Vector of detected normalization factors, with one factor per 
#'     sample.}
#' \item{h0.feature}{Vector of detected no-variation features.}
#' \item{within.condition.offset}{Normalization factors detected in the 
#'     within-condition normalizations.}
#' \item{between.condition.offset}{Normalization factors detected in the 
#'     between-condition normalization.}
#' \item{h0.feature.convergence}{List with convergence information for the 
#'     detection of no-variation features.}
#' 
#' @author Carlos P. Roca, \email{carlosproca@@gmail.com}
#' 
#' @references Roca, Gomes, Amorim & Scott-Fordsmand: 
#'     Variation-preserving normalization unveils blind spots in gene 
#'     expression profiling. 
#'     \emph{Sci. Rep.} \bold{7}, 42460; 
#'     \href{http://dx.doi.org/10.1038/srep42460}{doi:10.1038/srep42460} 
#'     (2017). 
#' 
#' @seealso \code{\link{normalize.svcd}} Implements SVCD normalization.
#' 
#' @examples
#' # no offset
#' gene.n <- 1000
#' sample.n <- 9
#' expr.data <- matrix( rnorm( gene.n * sample.n ), nrow = gene.n )
#' expr.condition <- rep( c( 1, 2, 3 ), each = 3 )
#' normalize.result <- normalize.mediancd( expr.data, expr.condition )
#' sd( normalize.result$offset )
#' length( normalize.result$h0.feature )
#' 
#' \dontrun{
#' # with offset
#' gene.n <- 10000
#' sample.n <- 9
#' expr.data <- matrix( rnorm( gene.n * sample.n ), nrow = gene.n )
#' expr.condition <- rep( c( "treatment.a", "treatment.b", "control" ), 
#'     each = 3 )
#' offset.added <- rnorm( sample.n )
#' expr.data <- sweep( expr.data, 2, offset.added, "+" )
#' normalize.result <- normalize.mediancd( expr.data, expr.condition, 
#'     p.value.graph = "mediancd_p_value", verbose = TRUE )
#' sd( normalize.result$offset - offset.added )
#' length( normalize.result$h0.feature )
#' }
#' 
#' @export

normalize.mediancd <- function( expression.data, expression.condition, 
    normalization.probability = 0.5, restrict.feature = NULL, 
    search.h0.feature = TRUE, convergence.threshold = c( 0.001, 0.1 ), 
    p.value.graph = NULL, verbose = FALSE )
{
    # check arguments
    normalize.mediancd.check.argument( expression.data, expression.condition, 
        normalization.probability, restrict.feature, search.h0.feature, 
        convergence.threshold, p.value.graph, verbose )
    
    # add row names to expression.data if missing
    if ( is.null( rownames( expression.data ) ) )
    {
        feature.num <- nrow( expression.data )
        feature.name.width <- floor( log10( feature.num ) ) + 1
        feature.name <- sprintf( "feature.%0*d", feature.name.width, 
            1 : feature.num )
        rownames( expression.data ) <- feature.name
    }
    
    # add column names to expression.data if missing
    if ( is.null( colnames( expression.data ) ) )
    {
        sample.num <- ncol( expression.data )
        sample.name.width <- floor( log10( sample.num ) ) + 1
        sample.name <- sprintf( "sample.%0*d", sample.name.width, 
            1 : sample.num )
        colnames( expression.data ) <- sample.name
    }
    
    # identify samples and conditions to normalize
    normalize.sample <- 
        colnames( expression.data )[ ! is.na( expression.condition ) ]
    normalize.sample.condition <- 
        as.character( na.omit( expression.condition ) )
    normalize.condition <- unique( normalize.sample.condition )
    
    if ( length( normalize.condition ) == 0 )
        stop( "normalize.mediancd: No condition to normalize" )
    else if ( min( table( normalize.sample.condition ) ) < 2 )
        stop( "normalize.mediancd: There must be 2 or more samples for each condition" )
    
    # select features for normalization
    expression.feature <- rownames( expression.data )
    
    if ( is.null( restrict.feature ) )
        restrict.feature.idx <- 1 : length( expression.feature )
    else
    {
        restrict.feature.idx <- match( restrict.feature, expression.feature )
        if ( any( is.na( restrict.feature.idx ) ) )
            stop( "normalize.mediancd: Bad restrict.feature argument" )
    }
    
    # enforce no missing values in any normalization sample
    all.sample.feature.idx <- which( rowSums( 
        is.na( expression.data[ , normalize.sample ] ) ) == 0 )
    
    restrict.feature.idx <- 
        intersect( restrict.feature.idx, all.sample.feature.idx )
    
    # normalize within conditions

    normalize.expr.data <- matrix( nrow = length( expression.feature ), 
        ncol = length( normalize.sample ) )
    rownames( normalize.expr.data ) <- expression.feature
    colnames( normalize.expr.data ) <- normalize.sample
    
    normalize.within.cond.offset <- rep( 0, length( normalize.sample ) )
    names( normalize.within.cond.offset ) <- normalize.sample
    
    condition.sample.idx <- split( 1 : length( expression.condition ), 
        expression.condition )
    
    for ( sample.idx in condition.sample.idx )
    {
        condition <- as.character( expression.condition[ sample.idx[ 1 ] ] )
        norm.sample.idx <- which( condition == normalize.sample.condition )
        
        within.cond.norm.result <- normalize.median.within.condition( 
            expression.data[ , sample.idx ], condition, 
            normalization.probability, restrict.feature.idx, verbose )
        
        normalize.expr.data[ , norm.sample.idx ] <- within.cond.norm.result$data
        
        normalize.within.cond.offset[ norm.sample.idx ] <- 
            within.cond.norm.result$offset
    }
    
    if ( length( normalize.condition ) == 1 )
    {
        # no additional normalization needed
        normalize.result <- list( data = normalize.expr.data, 
            offset = normalize.within.cond.offset, 
            h0.feature = NULL, 
            within.condition.offset = normalize.within.cond.offset, 
            between.condition.offset = NULL, 
            h0.feature.convergence = NULL )
        
        return( normalize.result )
    }
    
    # normalize between conditions
    
    between.cond.norm.result <- normalize.median.between.condition( 
        normalize.expr.data, normalize.condition, normalize.sample.condition, 
        normalization.probability, restrict.feature.idx, search.h0.feature, 
        convergence.threshold, p.value.graph, verbose )
    
    normalize.expr.data <- between.cond.norm.result$data
    
    normalize.between.cond.offset <- between.cond.norm.result$offset
    normalize.between.cond.h0.feature <- between.cond.norm.result$h0.feature
    normalize.between.cond.h0.feature.convergence <- 
        between.cond.norm.result$h0.feature.convergence
    
    normalize.offset <- normalize.within.cond.offset + 
        normalize.between.cond.offset[ normalize.sample.condition ]
    
    list( data = normalize.expr.data, 
        offset = normalize.offset, 
        h0.feature = normalize.between.cond.h0.feature, 
        within.condition.offset = normalize.within.cond.offset, 
        between.condition.offset = normalize.between.cond.offset, 
        h0.feature.convergence = normalize.between.cond.h0.feature.convergence )
}

