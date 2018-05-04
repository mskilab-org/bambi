
#' bambi 
#'
#' motivation explanation
#' queries WGS and SS 10x BAMs
#'
#'
#' @section Arguments:
#' \describe{
#'   \item{p}{A \code{process} object.}
#'   \item{command}{Character scalar, the command to run. It will be
#'     escaped via \code{\link[base]{shQuote}}.}
#'   \item{args}{Character vector, arguments to the command. The will be
#'     escaped via \code{\link[base]{shQuote}}.}
#'   \item{commandline}{A character scalar, a full command line.
#'     No escaping will be performed on it.}
#'   \item{stdout}{What to do with the standard output. Possible values:
#'     \code{FALSE}: discard it; a string, redirect it to this file,
#'     \code{TRUE}: redirect it to a temporary file.}
#'   \item{stdout}{What to do with the standard error. Possible values:
#'     \code{FALSE}: discard it; a string, redirect it to this file,
#'     \code{TRUE}: redirect it to a temporary file.}
#'   \item{grace}{Grace pediod between the TERM and KILL signals, in
#'     seconds.}
#'   \item{...}{Extra arguments are passed to the
#'     \code{\link[base]{readLines}} function.}
#' }
#'
#' @section Details:
#' \code{$grab_bx()} grab BX tags (if exist)
#' 
#' grab_bx = function(query, data.table = FALSE, verbose = FALSE, mc.cores = 1)
#'
#' \code{$grab_cb()} grab CB tags (if exist)
#'
#' grab_cb = function(query, data.table = FALSE, verbose = FALSE, mc.cores = 1)
#'
#' \code{$grab_ub()} grab CB tags (if exist)
#'
#' grab_ub = function(query, data.table = FALSE, verbose = FALSE, mc.cores = 1)
#'
#' \code{$fetch_by_tag()} grab CB tags (if exist)
#'
#' fetch_by_tag = function(tag, query, data.table = FALSE, verbose = FALSE, mc.cores = 1)
#'
#' @import R6
#' @importFrom R6 R6Class
#' @import gUtils
#' @name bambi
#'
NULL

#' @export

bambi = R6::R6Class('bambi',

    public = list(

        bam_file = NULL,
        bamdb_path = NULL,

        initialize = function(bam_file, bamdb_path){  

            if (!file.exists(bam_file) | is.null(bam_file)){
                stop("BAM file not found. A valid BAM for 'bam_file' must be provided.")
            }

            check_valid_bam = readChar(gzfile(bam_file, 'r'), 4)
            if (!identical(check_valid_bam, 'BAM\1')){
                stop("Cannot open BAM. A valid BAM for 'bam_file' must be provided.")
            }

            if (is.null(bamdb_path)){
                if (file.exists(bamdb_path <- gsub('.bam$', '_lmdb', bam_file))){
                    bamdb_path = bamdb_path 
                } else{
                    stop('Error: Pathname to bamdb LMDB subdirectory must be provided')
                }
            }

            ## checks for LMDB path? Check it is a subdirectory at least

            self$bam_file = bam_file
            self$bamdb_path = bamdb_path  
        },


        grab_bx = function(barcodes, query=NULL, data.table = FALSE, verbose = FALSE, mc.cores = 1){     

            if (!inherits(barcodes, "character")){
                stop("Invalid barcodes. Input 'barcodes' must be a character vector. Must provide barcode with input BAM.")
            }

            if (!is.null(query)){
                if (inherits(query, 'GRanges')){

                    if (is.null(query$BX)){

                        if (verbose){
                            message('Using read.bam to pull reads under GRanges query from BAM file and find their BX tags.')
                            now = Sys.time()
                        }

                        query = read.bam(bam,  gr = query, tag = c('MD'), pairs.grl = FALSE)   ### what if there are no MD tags? Use RSamtools to check for MD tags?

                        if (verbose){
                            message(sprintf('Retrieved %s reads with %s unique barcodes:', length(query), length(unique(query$BX))))
                            print(Sys.time() - now)
                        }
                    } 
    
                    barcodes = query$BX  ## vector of BX's 
    
                } else if(inherits(query, 'data.frame') | inherits(query, 'data.table')){
             
                    if (is.null(query$BX)){

                        if (verbose){
                            message('Using read.bam to pull reads under data.table/data.frame query from BAM file and find their BX tags.')
                            now = Sys.time()
                        }

                        query = read.bam(bam,  gr = dt2gr(query), tag = c('MD'), pairs.grl = FALSE)   ### what if there are no MD tags? Use RSamtools to check for MD tags?

                        if (verbose){
                            message(sprintf('Retrieved %s reads with %s unique barcodes:', length(query), length(unique(query$CB))))
                            print(Sys.time() - now)
                        }
                    } 

                    barcodes = query$BX  ## vector of BX's 

                } else{
                    stop("Invalid query. Input 'query' must be a data.table, data.frame or GRanges.")
                }
            
            }
            
            barcodes = as.character(barcodes)

            out = query_bam_index(self$bam_file, self$bamdb_path, "BX", barcodes)

            out = as.data.table(out)
            out[, BX := barcodes]  ## for one
        
            if (any(nnix <<- out$cigar=='*')){
                out$cigar[nnix] = NA
            }

            if (data.table == TRUE){
                return(as.data.table(out)) ### check format
            } else{
                parse_outputs(out)
            }
        }, 

        grab_cb = function(barcodes, query=NULL, data.table = FALSE, verbose = FALSE, mc.cores = 1){     

            if (!inherits(barcodes, "character")){
                stop("Invalid barcodes. Input 'barcodes' must be a character vector. Must provide barcode with input BAM.")
            }

            if (!is.null(query)){
                if (inherits(query, 'GRanges')){

                    if (is.null(query$CB)){

                        if (verbose){
                            message('Using read.bam to pull reads under GRanges query from BAM file and find their BX tags.')
                            now = Sys.time()
                        }

                        query = read.bam(bam,  gr = query, tag = c('MD'), pairs.grl = FALSE)   ### what if there are no MD tags? Use RSamtools to check for MD tags?

                        if (verbose){
                            message(sprintf('Retrieved %s reads with %s unique barcodes:', length(query), length(unique(query$BX))))
                            print(Sys.time() - now)
                        }
                    } 
    
                    barcodes = query$CB  ## vector of CB's 
    
                } else if(inherits(query, 'data.frame') | inherits(query, 'data.table')){
             
                    if (is.null(query$UB)){

                        if (verbose){
                            message('Using read.bam to pull reads under data.table/data.frame query from BAM file and find their BX tags.')
                            now = Sys.time()
                        }

                        query = read.bam(bam,  gr = dt2gr(query), tag = c('MD'), pairs.grl = FALSE)   ### what if there are no MD tags? Use RSamtools to check for MD tags?

                        if (verbose){
                            message(sprintf('Retrieved %s reads with %s unique barcodes:', length(query), length(unique(query$CB))))
                            print(Sys.time() - now)
                        }
                    } 

                    barcodes = query$CB ## vector of CB's 

                } else{
                    stop("Invalid query. Input 'query' must be a data.table, data.frame or GRanges.")
                }
            
            }
            
            barcodes = as.character(barcodes)

            out = query_bam_index(self$bam_file, self$bamdb_path, "CB", barcodes)

            out = as.data.table(out)
            out[, CB := barcodes]  ## for one
        
            if (any(nnix <<- out$cigar=='*')){
                out$cigar[nnix] = NA
            }

            if (data.table == TRUE){
                return(as.data.table(out)) ### check format
            } else{
                parse_outputs(out)
            }
        }, 

        grab_ub = function(barcodes, query=NULL, data.table = FALSE, verbose = FALSE, mc.cores = 1){     

            if (!inherits(barcodes, "character")){
                stop("Invalid barcodes. Input 'barcodes' must be a character vector. Must provide barcode with input BAM.")
            }

            if (!is.null(query)){
                if (inherits(query, 'GRanges')){

                    if (is.null(query$UB)){

                        if (verbose){
                            message('Using read.bam to pull reads under GRanges query from BAM file and find their BX tags.')
                            now = Sys.time()
                        }

                        query = read.bam(bam,  gr = query, tag = c('MD'), pairs.grl = FALSE)   ### what if there are no MD tags? Use RSamtools to check for MD tags?

                        if (verbose){
                            message(sprintf('Retrieved %s reads with %s unique barcodes:', length(query), length(unique(query$BX))))
                            print(Sys.time() - now)
                        }
                    } 
    
                    barcodes = query$UB  ## vector of UB's 
    
                } else if(inherits(query, 'data.frame') | inherits(query, 'data.table')){
             
                    if (is.null(query$UB)){

                        if (verbose){
                            message('Using read.bam to pull reads under data.table/data.frame query from BAM file and find their BX tags.')
                            now = Sys.time()
                        }

                        query = read.bam(bam,  gr = dt2gr(query), tag = c('MD'), pairs.grl = FALSE)   ### what if there are no MD tags? Use RSamtools to check for MD tags?

                        if (verbose){
                            message(sprintf('Retrieved %s reads with %s unique barcodes:', length(query), length(unique(query$CB))))
                            print(Sys.time() - now)
                        }
                    } 

                    barcodes = query$UB ## vector of UB's 

                } else{
                    stop("Invalid query. Input 'query' must be a data.table, data.frame or GRanges.")
                }
            
            }
            
            barcodes = as.character(barcodes)

            out = query_bam_index(self$bam_file, self$bamdb_path, "UB", barcodes)

            out = as.data.table(out)
            out[, UB := barcodes]  ## for one
        
            if (any(nnix <<- out$cigar=='*')){
                out$cigar[nnix] = NA
            }

            if (data.table == TRUE){
                return(as.data.table(out)) ### check format
            } else{
                parse_outputs(out)
            }
        }, 


        fetch_by_tag = function(tag, tag_queries, query=NULL, data.table = FALSE, verbose = FALSE, mc.cores = 1){     

            if (!inherits(tag, "character")){
                stop("Invalid tag. Input 'tag' must be a character vector. Must provide valid BAM field.")
            }

            if (!inherits(tag_queries, "character")){
                stop("Invalid tag_queries. Input 'tag_queries' must be a character vector. Must provide valid BAM queries by BAM field 'tag'.")
            }

            ## currently there's no infrastructure for multiple queries of a different type
            if (length(unique(tag)) != 1){
                stop("Invalid tag. Currently, only one unique type of tag at a time supported. ")
            }

            if (!is.null(query)){
                if (inherits(query, 'GRanges')){

                    if (is.null(query$BX)){

                        if (verbose){
                            message('Using read.bam to pull reads under GRanges query from BAM file and find their BX tags.')
                            now = Sys.time()
                        }

                        query = read.bam(bam,  gr = query, tag = c('MD'), pairs.grl = FALSE)   ### what if there are no MD tags? Use RSamtools to check for MD tags?

                        if (verbose){
                            message(sprintf('Retrieved %s reads with %s unique barcodes:', length(query), length(unique(query$BX))))
                            print(Sys.time() - now)
                        }
                    } 
    
                    tag_queries = query$tag  ## vector of tag_queries by BAM field 'tag'
    
                } else if(inherits(query, 'data.frame') | inherits(query, 'data.table')){
             
                    if (is.null(query$BX)){

                        if (verbose){
                            message('Using read.bam to pull reads under data.table/data.frame query from BAM file and find their BX tags.')
                            now = Sys.time()
                        }

                        query = read.bam(bam,  gr = dt2gr(query), tag = c('MD'), pairs.grl = FALSE)   ### what if there are no MD tags? Use RSamtools to check for MD tags?

                        if (verbose){
                            message(sprintf('Retrieved %s reads with %s unique barcodes:', length(query), length(unique(query$CB))))
                            print(Sys.time() - now)
                        }
                    } 

                    tag_queries = query$tag  ## vector of tag_queries by BAM field 'tag'

                } else{
                    stop("Invalid query. Input 'query' must be a data.table, data.frame or GRanges.")
                }
            
            }
            
            tag_queries = as.character(tag_queries)

            out = query_bam_index(self$bam_file, self$bamdb_path, tag, tag_queries)

            out = as.data.table(out)
            out[, tag := tag_queries]  ## for one
        
            if (any(nnix <<- out$cigar=='*')){
                out$cigar[nnix] = NA
            }

            if (data.table == TRUE){
                return(as.data.table(out)) ### check format
            } else{
                parse_outputs(out)
            }
        }
    )
)


