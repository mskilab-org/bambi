
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
#' @import data.table
#' @name bambi
#'
NULL

#' @export

bambi = R6::R6Class('bambi',

    public = list(

        bam_file = NULL,
        bamdb_path = NULL,

        initialize = function(bam_file, bamdb_path=NULL){  

            if (!file.exists(bam_file) | is.null(bam_file)){
                stop("BAM file not found. A valid BAM for 'bam_file' must be provided.")
            }

            check_gz = gzfile(bam_file, 'r')
            check_valid_bam = readChar(check_gz, 4)
            if (!identical(check_valid_bam, 'BAM\1')){
                stop("Cannot open BAM. A valid BAM for 'bam_file' must be provided.")
            }
            on.exit(close(check_gz))

            if (is.null(bamdb_path)){
                if (dir.exists(bamdb_path <- gsub('.bam$', '_lmdb', bam_file))){
                    bamdb_path = bamdb_path 
                } else{
                    stop('Error: Pathname to bamdb LMDB subdirectory must be provided')
                }
            }

            self$bam_file = bam_file
            self$bamdb_path = bamdb_path  
        },


        grab_bx = function(barcodes=NULL, query=NULL, data.table=FALSE, verbose=FALSE, mc.cores=1){   

            ## check BX exists
            if (check_index(self$bamdb_path, 'BX') != TRUE){
                stop("BX is not contained as an LMDB key in 'bamdb_path'. Please see documentation for details.")
            }

            if ((!is.null(barcodes)) & (!is.null(query))){
                stop("Both 'barcodes' and 'query' parameters cannot be used. Use method 'grab_bx()' by a character vector of BX barcodes, or a GRanges/data.table of a genomic region. Please see documentation for details.")
            }

            if ((is.null(barcodes)) & (is.null(query))){
                if (data.table==TRUE){
                    return(data.table())
                } else{
                    return(GRanges())
                }
            } else if ((!is.null(barcodes)) & (is.null(query))){

                if (!inherits(barcodes, "character")){
                    stop("Invalid barcodes. Input 'barcodes' must be a character vector. Must provide barcode with input BAM.")
                }

                barcodes = as.character(barcodes)

                if (length(barcodes)==1){

                    out = query_bam_index(self$bam_file, self$bamdb_path, "BX", barcodes)

                    out = as.data.table(out)
 
                    if (as.integer(dim(out)[1]) == 0){   ## query has no results, no rows returned
                        return(NA)
                    } else{
                        out[, BX := barcodes]  ## for one
        
                        if (any(nnix <<- out$cigar=='*')){
                            out$cigar[nnix] = NA
                        }
 
                        if (data.table == TRUE){
                            return(as.data.table(out)) 
                        } else{
                            return(parse_outputs(out))
                        }
                    }
                } else{
                    ## here use mclapply to loop through BX vector, and concatenate data.tables
                    query_list = NULL
                    multiple_barcodes = mclapply(1:length(barcodes), 
                        function(x){
                        output = query_bam_index(self$bam_file, self$bamdb_path, "BX", barcodes[x])
                        query_list = rbind(query_list, output)     ## create huge list of data.table'
                    }, mc.cores=mc.cores)

                    out = rbindlist(multiple_barcodes, fill=TRUE)
                    suppressWarnings(out[, BX := barcodes])  ## for one
        
                    if (any(nnix <<- out$cigar=='*')){
                        out$cigar[nnix] = NA
                    }

                    if (data.table == TRUE){
                        return(as.data.table(out)) 
                    } else{
                        return(parse_outputs(out))
                    }
                }

            } else if ((is.null(barcodes)) & (!is.null(query))){

                if (inherits(query, 'GRanges')){

                    if (is.null(query$BX)){

                        if (verbose){
                            message('Using read.bam to pull reads under GRanges query from BAM file and find their BX tags.')
                            now = Sys.time()
                        }

                        query_try = try(unlist(read.bam(self$bam_file,  gr = query, tag = c('BX','MD'), pairs.grl = FALSE)), silent=TRUE)   ### if no MD tag, bamUtils::read.bam() outputs an NA for this column
                        if (inherits(query_try, "try-error")){
                            if ('UCSC' %in% seqlevelsStyle(query)){
                                ## change from UCSC to Ensemble, e.g chr5 -> 5
                                message('Converting seqlevels of "query" GRanges from UCSC style to Ensembl style, e.g. from format chr# to #')
                                seqlevelsStyle(query) = 'Ensembl'
                                query_try = try(unlist(read.bam(self$bam_file,  gr = query, tag = c('BX','MD'), pairs.grl = FALSE))) 
                            } else if ('Ensembl' %in% seqlevelsStyle(query)){
                                ## change Ensembl to UCSC, e.g 5 -> chr5
                                message('Converting seqlevels of "query" GRanges from Ensembl style to UCSC style, e.g. from format # to chr#')
                                seqlevelsStyle(query) = 'UCSC'
                                query_try = try(unlist(read.bam(self$bam_file,  gr = query, tag = c('BX','MD'), pairs.grl = FALSE))) 
                            } else{
                                stop("Check GRanges input for 'query'. Appears to be formatted incorrectly.")
                            }
                            query_try = unlist(read.bam(self$bam_file,  gr = query, tag = c('BX','MD'), pairs.grl = FALSE))
                        }

                        if (verbose){
                            message(sprintf('Retrieved %s reads with %s unique barcodes:', length(query), length(unique(query$BX))))
                            print(Sys.time() - now)
                        }

                    query = query_try

                    } 
    
                    barcodes = query$BX  ## vector of BX's 
    
                } else if(inherits(query, 'data.frame') | inherits(query, 'data.table')){
             
                    if (is.null(query$BX)){

                        if (verbose){
                            message('Using read.bam to pull reads under data.table/data.frame query from BAM file and find their BX tags.')
                            now = Sys.time()
                        }

                        query = unlist(read.bam(self$bam_file,  gr = query, tag = c('BX','MD'), pairs.grl = FALSE))   ### if no MD tag, bamUtils::read.bam() outputs an NA for this column

                        if (verbose){
                            message(sprintf('Retrieved %s reads with %s unique barcodes:', length(query), length(unique(query$CB))))
                            print(Sys.time() - now)
                        }
                    } 

                    barcodes = query$BX  ## vector of BX's 

                } else{
                    stop("Invalid query. Input 'query' must be a data.table, data.frame or GRanges.")
                }

                if (length(query)==0){
                    return(NA)
                } else{

                    barcodes = as.character(barcodes)

                    ## here use mclapply to loop through BX vector, and concatenate data.tables
                    query_list = NULL
                    multiple_barcodes = mclapply(1:length(barcodes), 
                        function(x){
                        output = query_bam_index(self$bam_file, self$bamdb_path, "BX", barcodes[x])
                        query_list = rbind(query_list, output)     ## create huge list of data.table'
                    }, mc.cores=mc.cores)
    
                    out = rbindlist(multiple_barcodes, fill=TRUE)
                    suppressWarnings(out[, BX := barcodes])  ## for one
        
                    if (any(nnix <<- out$cigar=='*')){
                        out$cigar[nnix] = NA
                    }

                    if (data.table == TRUE){
                        return(as.data.table(out)) 
                    } else{
                        return(parse_outputs(out))
                    }
                }
            }
        }, 


        grab_cb = function(barcodes=NULL, query=NULL, data.table=FALSE, verbose=FALSE, mc.cores=1){   

            ## check CB exists
            if (check_index(self$bamdb_path, 'CB') != TRUE){
                stop("CB is not contained as an LMDB key in 'bamdb_path'. Please see documentation for details.")
            }

            if ((!is.null(barcodes)) & (!is.null(query))){
                stop("Both 'barcodes' and 'query' parameters cannot be used. Use method 'grab_cb()' by a character vector of CB barcodes, or a GRanges/data.table of a genomic region. Please see documentation for details.")
            }

            if ((is.null(barcodes)) & (is.null(query))){
                if (data.table==TRUE){
                    return(data.table())
                } else{
                    return(GRanges())
                }
            } else if ((!is.null(barcodes)) & (is.null(query))){

                if (!inherits(barcodes, "character")){
                    stop("Invalid barcodes. Input 'barcodes' must be a character vector. Must provide barcode with input BAM.")
                }

                barcodes = as.character(barcodes)

                if (length(barcodes)==1){

                    out = query_bam_index(self$bam_file, self$bamdb_path, "CB", barcodes)

                    out = as.data.table(out)

                    if (as.integer(dim(out)[1]) == 0){   ## query has no results, no rows returned
                        return(NA)
                    } else{

                        out[, CB := barcodes]  ## for one
        
                        if (any(nnix <<- out$cigar=='*')){
                            out$cigar[nnix] = NA
                        }

                        if (data.table == TRUE){
                            return(as.data.table(out)) 
                        } else{
                            return(parse_outputs(out))
                        }
                    }
                } else{
                    ## here use mclapply to loop through CB vector, and concatenate data.tables
                    query_list = NULL
                    multiple_barcodes = mclapply(1:length(barcodes), 
                        function(x){
                        output = query_bam_index(self$bam_file, self$bamdb_path, "CB", barcodes[x])
                        query_list = rbind(query_list, output)     ## create huge list of data.table'
                    }, mc.cores=mc.cores)

                    out = rbindlist(multiple_barcodes, fill=TRUE)
                    suppressWarnings(out[, CB := barcodes])  ## for one
        
                    if (any(nnix <<- out$cigar=='*')){
                        out$cigar[nnix] = NA
                    }

                    if (data.table == TRUE){
                        return(as.data.table(out)) 
                    } else{
                        return(parse_outputs(out))
                    }
                }

            } else if ((is.null(barcodes)) & (!is.null(query))){

                if (inherits(query, 'GRanges')){

                    if (is.null(query$CB)){

                        if (verbose){
                            message('Using read.bam to pull reads under GRanges query from BAM file and find their CB tags.')
                            now = Sys.time()
                        }

                        query_try = try(unlist(read.bam(self$bam_file,  gr = query, tag = c('CB','MD'), pairs.grl = FALSE)), silent=TRUE)   ### if no MD tag, bamUtils::read.bam() outputs an NA for this column
                        if (inherits(query_try, "try-error")){
                            if ('UCSC' %in% seqlevelsStyle(query)){
                                ## change from UCSC to Ensemble, e.g chr5 -> 5
                                message('Converting seqlevels of "query" GRanges from UCSC style to Ensembl style, e.g. from format chr# to #')
                                seqlevelsStyle(query) = 'Ensembl'
                                query_try = try(unlist(read.bam(self$bam_file,  gr = query, tag = c('CB','MD'), pairs.grl = FALSE))) 
                            } else if ('Ensembl' %in% seqlevelsStyle(query)){
                                ## change Ensembl to UCSC, e.g 5 -> chr5
                                message('Converting seqlevels of "query" GRanges from Ensembl style to UCSC style, e.g. from format # to chr#')
                                seqlevelsStyle(query) = 'UCSC'
                                query_try = try(unlist(read.bam(self$bam_file,  gr = query, tag = c('CB','MD'), pairs.grl = FALSE))) 
                            } else{
                                stop("Check GRanges input for 'query'. Appears to be formatted incorrectly.")
                            }
                            query_try = unlist(read.bam(self$bam_file,  gr = query, tag = c('CB','MD'), pairs.grl = FALSE))
                        }

                        if (verbose){
                            message(sprintf('Retrieved %s reads with %s unique barcodes:', length(query), length(unique(query$CB))))
                            print(Sys.time() - now)
                        }

                    query = query_try    

                    } 

                    barcodes = query$CB  ## vector of CB's 
    
                } else if(inherits(query, 'data.frame') | inherits(query, 'data.table')){
             
                    if (is.null(query$CB)){

                        if (verbose){
                            message('Using read.bam to pull reads under data.table/data.frame query from BAM file and find their CB tags.')
                            now = Sys.time()
                        }

                        query = unlist(read.bam(self$bam_file,  gr = query, tag = c('CB','MD'), pairs.grl = FALSE))  ### if no MD tag, bamUtils::read.bam() outputs an NA for this column

                        if (verbose){
                            message(sprintf('Retrieved %s reads with %s unique barcodes:', length(query), length(unique(query$CB))))
                            print(Sys.time() - now)
                        }
                    } 

                    barcodes = query$CB  ## vector of CB's 

                } else{
                    stop("Invalid query. Input 'query' must be a data.table, data.frame or GRanges.")
                }

                if (length(query)==0){
                    return(NA)
                } else{
                    barcodes = as.character(barcodes)

                    ## here use mclapply to loop through BX vector, and concatenate data.tables
                    query_list = NULL
                    multiple_barcodes = mclapply(1:length(barcodes), 
                        function(x){
                        output = query_bam_index(self$bam_file, self$bamdb_path, "CB", barcodes[x])
                        query_list = rbind(query_list, output)     ## create huge list of data.table'
                    }, mc.cores=mc.cores)

                    out = rbindlist(multiple_barcodes, fill=TRUE)
                    suppressWarnings(out[, BX := barcodes])  ## for one
        
                    if (any(nnix <<- out$cigar=='*')){
                        out$cigar[nnix] = NA
                    }

                    if (data.table == TRUE){
                        return(as.data.table(out)) 
                    } else{
                        return(parse_outputs(out))
                    }
                }
            }
        }, 

        grab_ub = function(barcodes=NULL, query=NULL, data.table=FALSE, verbose=FALSE, mc.cores=1){   

            ## check UB exists
            if (check_index(self$bamdb_path, 'UB') != TRUE){
                stop("UB is not contained as an LMDB key in 'bamdb_path'. Please see documentation for details.")
            }

            if ((!is.null(barcodes)) & (!is.null(query))){
                stop("Both 'barcodes' and 'query' parameters cannot be used. Use method 'grab_ub()' by a character vector of UB barcodes, or a GRanges/data.table of a genomic region. Please see documentation for details.")
            }

            if ((is.null(barcodes)) & (is.null(query))){
                if (data.table==TRUE){
                    return(data.table())
                } else{
                    return(GRanges())
                }
            } else if ((!is.null(barcodes)) & (is.null(query))){

                if (!inherits(barcodes, "character")){
                    stop("Invalid barcodes. Input 'barcodes' must be a character vector. Must provide barcode with input BAM.")
                }

                barcodes = as.character(barcodes)

                if (length(barcodes)==1){

                    out = query_bam_index(self$bam_file, self$bamdb_path, "UB", barcodes)

                    out = as.data.table(out)

                    if (as.integer(dim(out)[1]) == 0){   ## query has no results, no rows returned
                        return(NA)
                    } else{

                        out[, UB := barcodes]  ## for one
        
                        if (any(nnix <<- out$cigar=='*')){
                            out$cigar[nnix] = NA
                        }

                        if (data.table == TRUE){
                            return(as.data.table(out)) 
                        } else{
                            return(parse_outputs(out))
                        }
                    }
                } else{
                    ## here use mclapply to loop through UB vector, and concatenate data.tables
                    query_list = NULL
                    multiple_barcodes = mclapply(1:length(barcodes), 
                        function(x){
                        output = query_bam_index(self$bam_file, self$bamdb_path, "UB", barcodes[x])
                        query_list = rbind(query_list, output)     ## create huge list of data.table'
                    }, mc.cores=mc.cores)

                    out = rbindlist(multiple_barcodes, fill=TRUE)
                    suppressWarnings(out[, BX := barcodes])  ## for one
        
                    if (any(nnix <<- out$cigar=='*')){
                        out$cigar[nnix] = NA
                    }

                    if (data.table == TRUE){
                        return(as.data.table(out)) 
                    } else{
                        return(parse_outputs(out))
                    }
                }

            } else if ((is.null(barcodes)) & (!is.null(query))){

                if (inherits(query, 'GRanges')){

                    if (is.null(query$UB)){

                        if (verbose){
                            message('Using read.bam to pull reads under GRanges query from BAM file and find their UB tags.')
                            now = Sys.time()
                        }

                        query_try = try(unlist(read.bam(self$bam_file,  gr = query, tag = c('UB','MD'), pairs.grl = FALSE)), silent=TRUE)   ### if no MD tag, bamUtils::read.bam() outputs an NA for this column
                        if (inherits(query_try, "try-error")){
                            if ('UCSC' %in% seqlevelsStyle(query)){
                                ## change from UCSC to Ensemble, e.g chr5 -> 5
                                message('Converting seqlevels of "query" GRanges from UCSC style to Ensembl style, e.g. from format chr# to #')
                                seqlevelsStyle(query) = 'Ensembl'
                                query_try = try(unlist(read.bam(self$bam_file,  gr = query, tag = c('UB','MD'), pairs.grl = FALSE))) 
                            } else if ('Ensembl' %in% seqlevelsStyle(query)){
                                ## change Ensembl to UCSC, e.g 5 -> chr5
                                message('Converting seqlevels of "query" GRanges from Ensembl style to UCSC style, e.g. from format # to chr#')
                                seqlevelsStyle(query) = 'UCSC'
                                query_try = try(unlist(read.bam(self$bam_file,  gr = query, tag = c('UB','MD'), pairs.grl = FALSE))) 
                            } else{
                                stop("Check GRanges input for 'query'. Appears to be formatted incorrectly.")
                            }
                            query_try = unlist(read.bam(self$bam_file,  gr = query, tag = c('UB','MD'), pairs.grl = FALSE))
                        }

                        if (verbose){
                            message(sprintf('Retrieved %s reads with %s unique barcodes:', length(query), length(unique(query$UB))))
                            print(Sys.time() - now)
                        }

                    query = query_try

                    } 
    
                    barcodes = query$UB  ## vector of CB's 
    
    
                } else if(inherits(query, 'data.frame') | inherits(query, 'data.table')){
             
                    if (is.null(query$UB)){

                        if (verbose){
                            message('Using read.bam to pull reads under data.table/data.frame query from BAM file and find their BX tags.')
                            now = Sys.time()
                        }

                        query = unlist(read.bam(self$bam_file,  gr = query, tag = c('UB','MD'), pairs.grl = FALSE))   ### if no MD tag, bamUtils::read.bam() outputs an NA for this column

                        if (verbose){
                            message(sprintf('Retrieved %s reads with %s unique barcodes:', length(query), length(unique(query$UB))))
                            print(Sys.time() - now)
                        }
                    } 

                    barcodes = query$UB  ## vector of UB's 

                } else{
                    stop("Invalid query. Input 'query' must be a data.table, data.frame or GRanges.")
                }

                if (length(query)==0){
                    return(NA)
                } else{

                    barcodes = as.character(barcodes)

                    ## here use mclapply to loop through UB vector, and concatenate data.tables
                    query_list = NULL
                    multiple_barcodes = mclapply(1:length(barcodes), 
                        function(x){
                        output = query_bam_index(self$bam_file, self$bamdb_path, "UB", barcodes[x])
                        query_list = rbind(query_list, output)     ## create huge list of data.table'
                    }, mc.cores=mc.cores)

                    out = rbindlist(multiple_barcodes, fill=TRUE)
                    suppressWarnings(out[, BX := barcodes])  ## for one
        
                    if (any(nnix <<- out$cigar=='*')){
                        out$cigar[nnix] = NA
                    }

                    if (data.table == TRUE){
                        return(as.data.table(out)) 
                    } else{
                        return(parse_outputs(out))
                    }
                }
            }
        }, 

        fetch_by_tag = function(tag, tag_queries=NULL, query=NULL, data.table=FALSE, verbose=FALSE, mc.cores=1){     

            ## currently only one tag at a time supported
            if (!inherits(tag, "character")){
                stop("Invalid tag. Input 'tag' must be a single character string. Must provide valid BAM field.")
            }
            ## currently there's no infrastructure for multiple queries of a different type
            ## if (length(unique(tag)) != 1){
            ##    stop("Invalid tag. Currently, only one unique type of tag at a time supported. ")
            ##}

            if (length(tag) != 1){
                stop("Invalid tag input. Multiple tags at once is not currently supported with 'fetch_by_tag()'. Please see documentation for details.")
            }

            ## check tag exists
            if (check_index(self$bamdb_path, tag) != TRUE){
                stop("The input 'tag' is not contained as an LMDB key in 'bamdb_path'. Please see documentation for details.")
            }

            if ((!is.null(tag_queries)) & (!is.null(query))){
                stop("Both 'tag_queries' and 'query' parameters cannot be used. Use method 'fetch_by_tag()' by either a character vector of UB barcodes, or a GRanges/data.table of a genomic region. Please see documentation for details.")
            } 


            if ((is.null(tag_queries)) & (is.null(query))){
                if (data.table==TRUE){
                    return(data.table())
                } else{
                    return(GRanges())
                }
            } else if ((!is.null(tag_queries)) & (is.null(query))){

                ## I think this is untrue...
                ##
                ## if (!inherits(tag_queries, "character")){
                ##     stop("Invalid tag_queries. Input 'tag_queries' must be a character vector of a BAM field. Must provide valid BAM field from input BAM.")
                ## }

                tag_queries = as.character(tag_queries)

                if (length(tag_queries)==1){

                    out = query_bam_index(self$bam_file, self$bamdb_path, tag, tag_queries)

                    out = as.data.table(out)

                    if (as.integer(dim(out)[1]) == 0){   ## query has no results, no rows returned
                        return(NA)
                    } else{                    
                        out[, tag := tag_queries]  ## for one
        
                        if (any(nnix <<- out$cigar=='*')){
                            out$cigar[nnix] = NA
                        }

                        if (data.table == TRUE){
                            return(as.data.table(out)) 
                        } else{
                            return(parse_outputs(out))
                        }
                    }
                } else{
                    ## here use mclapply to loop through UB vector, and concatenate data.tables
                    query_list = NULL
                    multiple_barcodes = mclapply(1:length(tag_queries), 
                        function(x){
                        output = query_bam_index(self$bam_file, self$bamdb_path, tag, tag_queries[x])
                        query_list = rbind(query_list, output)     ## create huge list of data.table'
                    }, mc.cores=mc.cores)

                    out = rbindlist(multiple_barcodes, fill=TRUE)
                    suppressWarnings(out[, tag := tag_queries])  ## for one
        
                    if (any(nnix <<- out$cigar=='*')){
                        out$cigar[nnix] = NA
                    }

                    if (data.table == TRUE){
                        return(as.data.table(out)) 
                    } else{
                        return(parse_outputs(out))
                    }
                }

            } else if((is.null(tag_queries)) & (!is.null(query))){

                if (inherits(query, 'GRanges')){

                    if (is.null(query$tag)){

                        if (verbose){
                            message('Using read.bam to pull reads under GRanges query from BAM file and find their asocciated tags given the parameter "tag".')
                            now = Sys.time()
                        }

                        query_try = try(unlist(read.bam(self$bam_file,  gr = query, tag = c('MD'), pairs.grl = FALSE)), silent=TRUE)   ### if no MD tag, bamUtils::read.bam() outputs an NA for this column
                        if (inherits(query_try, "try-error")){
                            if ('UCSC' %in% seqlevelsStyle(query)){
                                ## change from UCSC to Ensemble, e.g chr5 -> 5
                                message('Converting seqlevels of "query" GRanges from UCSC style to Ensembl style, e.g. from format chr# to #')
                                seqlevelsStyle(query) = 'Ensembl'
                                query_try = try(unlist(read.bam(self$bam_file,  gr = query, tag = c('MD'), pairs.grl = FALSE))) 
                            } else if ('Ensembl' %in% seqlevelsStyle(query)){
                                ## change Ensembl to UCSC, e.g 5 -> chr5
                                message('Converting seqlevels of "query" GRanges from Ensembl style to UCSC style, e.g. from format # to chr#')
                                seqlevelsStyle(query) = 'UCSC'
                                query_try = try(unlist(read.bam(self$bam_file,  gr = query, tag = c('MD'), pairs.grl = FALSE))) 
                            } else{
                                stop("Check GRanges input for 'query'. Appears to be formatted incorrectly.")
                            }
                            query_try = unlist(read.bam(self$bam_file,  gr = query, tag = c('MD'), pairs.grl = FALSE))
                        }

                        if (verbose){
                            message(sprintf('Retrieved %s reads with %s unique barcodes:', length(query), length(unique(query$tag)))) ## vector of tag's ### ERROR 
                            print(Sys.time() - now)
                        }

                    query = query_try

                    } 
    
                    barcodes = query$tag  ## vector of tag's ### ERROR  
    
                } else if(inherits(query, 'data.frame') | inherits(query, 'data.table')){
             
                    if (is.null(query$tag)){

                        if (verbose){
                            message('Using read.bam to pull reads under data.table/data.frame query from BAM file and find their asocciated tags given the parameter "tag".')
                            now = Sys.time()
                        }

                        query = read.bam(self$bam_file, gr = dt2gr(query), tag = c('MD'), pairs.grl = FALSE)   ### if no MD tag, bamUtils::read.bam() outputs an NA for this column

                        if (verbose){
                            message(sprintf('Retrieved %s reads with %s unique tags:', length(query), length(unique(query$UB))))
                            print(Sys.time() - now)
                        }
                    } 

                    barcodes = query$tag  ## vector of tag's 

                } else{
                    stop("Invalid query. Input 'query' must be a data.table, data.frame or GRanges.")
                }

                if (length(query)==0){
                    return(NA)
                } else{

                    barcodes = as.character(barcodes)

                    ## here use mclapply to loop through tags, and concatenate data.tables
                    query_list = NULL
                    multiple_barcodes = mclapply(1:length(barcodes), 
                        function(x){
                        output = query_bam_index(self$bam_file, self$bamdb_path, tag, barcodes[x])
                        query_list = rbind(query_list, output)     ## create huge list of data.table'
                    }, mc.cores=mc.cores)

                    out = rbindlist(multiple_barcodes, fill=TRUE)
                    suppressWarnings(out[, tag := barcodes])  ## for one
        
                    if (any(nnix <<- out$cigar=='*')){
                        out$cigar[nnix] = NA
                    }

                    if (data.table == TRUE){
                        return(as.data.table(out)) 
                    } else{
                        return(parse_outputs(out))
                    }
                }
            }
        }
    )
)






