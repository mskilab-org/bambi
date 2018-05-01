
#' bxBam 
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
#' \code{$grab_cb()} grab CB tags (if exist)
#'
#' \code{$grab_ub()} grab CB tags (if exist)
#'
#'
#' @importFrom R6 R6Class
#' @name bxBam
#'
NULL

#' @export

bxBam = R6Class('bxBam',

    public = list(

        bam_file = NULL,
        bamdb_path = NULL,

        initialize = function(bam_file, bamdb_path){  

            if (!file.exists(bam_file) | is.null(bam_file)){
                stop("BAM file not found. A valid BAM for 'bam_file' must be provided.")
            }
            
            if (is.null(bamdb_path)){
                if (file.exists(bamdb_path <- gsub('.bam$', '_lmdb', bam_file))){
                    bamdb_path = bamdb_path 
                } else{
                    stop('Error: Pathname to bamdb LMDB subdirectory must be provided')
                }
            }

            self$bam_file = bam_file
            self$bamdb_path = bamdb_path  ## check that the BAM exists with...something
        },

        grab_cb = function(query, data.table = FALSE, verbose = FALSE, mc.cores = 1){     

            if (inherits(query, 'GRanges') | inherits(query, 'data.frame') | inherits(query, 'data.table')){

                if (is.null(query$CB)){

                    if (verbose){
                        message('Using read.bam to pull reads under GRanges query from bam file and find their CB tags')
                        now = Sys.time()
                    }

                    query = read.bam(bam,  gr = query, tag = c('MD'), pairs.grl = FALSE)   ### what if there are no MD tags? Use RSamtools to check for MD tags?

                    if (verbose){
                        message(sprintf('Retrieved %s reads with %s unique barcodes:', length(query), length(unique(query$CB))))
                        print(Sys.time() - now)
                    }
                } 

            query = query$CB ## vector of CB's 

            }

            query = setdiff(query, NA)

            
            get_cb = .C('print_cb_rows', self$bam, self$bamdb, query)
         
            ##if (get_cb == 1){  ## throw error
            ##    stop(paste0('Error: there are no CB indices in ', self$bamdb))   ## check C code---there should be several errors here; what if "one" CB tag is wrong?
            ##}
        
            if (verbose){
                now = Sys.time()
            }

            ## now, parse into data.table
            out = rbindlist(mclapply(get_cb, function(){
                out = dcast.data.table(get_cb, V1 ~ V2, value.var = 'V3')[, -1, with = FALSE]
                setnames(out, names(out), ifelse(names(out)=='CB', names(out), tolower(names(out))))
                return(out)
            }, mc.cores = mc.cores))

            if (nrow(out) > 0){
                out$pos = as.numeric(out$pos)
                out$mapq = as.integer(out$mapq)
            }    

            if (verbose){ 
                message('Built out table with ', nrow(out), ' rows')
                print(Sys.time() - now)
            }

            if (any(nnix <<- out$cigar=='*')){
                out$cigar[nnix] = NA
            }

            if (data.table == TRUE){
                return(as.data.table(out)) ### check format
            } 
            else{
                parse_outputs(out)
            }
        }, 
        
        grab_ub = function(query, data.table = FALSE, verbose = FALSE, mc.cores = 1){  

            if (inherits(query, 'GRanges') | inherits(query, 'data.frame') | inherits(query, 'data.table')){

                if (is.null(query$UB)){

                    if (verbose){
                        message('Using read.bam to pull reads under GRanges query from bam file and find their UB tags')
                        now = Sys.time()
                    }

                    query = read.bam(bam,  gr = query, tag = c('MD'), pairs.grl = FALSE)   ### what is there are no MD tags? Use RSamtools to check for MD tags?
                
                    if (verbose){
                        message(sprintf('Retrieved %s reads with %s unique barcodes:', length(query), length(unique(query$UB))))
                        print(Sys.time() - now)
                    }
                } 

            query = query$UB ## vector of UB's 

            }

            query = setdiff(query, NA)

            get_ub = .C('print_ub_rows', self$bam, self$bamdb, query)
         
            ##if (get_ub == 1){  ## throw error
            ##    stop(paste0("Error: there are no UB indices in ", self$bamdb))   ## check C code---there should be several errors here; what if "one" UB tag is wrong?
            ##}
        
            if (verbose){
                now = Sys.time()
            }

            ## now, parse into data.table
            ##  assuming get_ub is now a variable of tab-delimited C data
            out = as.data.table(get_ub, sep='\n')
            
            ##
            ## if get_ub can stream input from C code:
            ## out = rbindlist(mclapply(get_ub, function(){
            ##     out = dcast.data.table(get_ub, V1 ~ V2, value.var = "V3")[, -1, with = FALSE]
            ##     setnames(out, names(out), ifelse(names(out)=="UB", names(out), tolower(names(out))))
            ##     return(out)
            ## }, mc.cores = mc.cores))

            if (nrow(out) > 0){
                out$pos = as.numeric(out$pos)
                out$mapq = as.integer(out$mapq)
            }    

            if (verbose){ 
                message('Built out table with ', nrow(out), ' rows')
                print(Sys.time() - now)
            }

            if (any(nnix <<- out$cigar=='*')){
                out$cigar[nnix] = NA
            }

            if (data.table == TRUE){
                return(as.data.table(out)) ### check format
            } else{
                parse_outputs(out)
            }
        }, 

        grab_qname = function(query, data.table = FALSE, verbose = FALSE, mc.cores = 1){             
        
            if (inherits(query, 'GRanges') | inherits(query, 'data.frame') | inherits(query, 'data.table')){
            
                if (is.null(query$QNAME)){
             
                    if (verbose){
                        message('Using read.bam to pull reads under GRanges query from bam file and find their QNAME tags')
                        now = Sys.time()
                    }

                query = read.bam(bam,  gr = query, tag = c('MD'), pairs.grl = FALSE)  ### WHAT TO DO ABOUT BX??
                
                    if (verbose){                                                         ### what is there are no MD tags? Use RSamtools to check for MD tags?
                        message(sprintf('Retrieved %s reads with %s unique barcodes:', length(query), length(unique(query$UB))))
                        print(Sys.time() - now)
                    }

                } 

                query = query$QNAME ## vector of UB's 

            }

            query = setdiff(query, NA)

            get_qname = .C('print_qname_rows', self$bam, self$bamdb, query)
         
            ##if (get_qname == 1){  ## throw error
            ##    stop(paste0('ERROR: there are no UB indices in ', self$bamdb))   ## check C code---there should be several errors here; what if "one" UB tag is wrong?
            ##}
        
            if (verbose){
                now = Sys.time()
            }
            
            ## now, parse into data.table
            ##  assuming get_ub is now a variable of tab-delimited C data
            out = as.data.table(get_qname, sep='\n')
            

            ## now, parse into data.table
            out = rbindlist(mclapply(get_qname, function(){
                out = dcast.data.table(dat, V1 ~ V2, value.var = "V3")[, -1, with = FALSE]
                setnames(out, names(out), ifelse(names(out)=="QNAME", names(out), tolower(names(out))))
                return(out)
            }, mc.cores = mc.cores))

            if (nrow(out) > 0){
                out$pos = as.numeric(out$pos)
                out$mapq = as.integer(out$mapq)
            }    

            if (verbose){ 
                message('Built out table with ', nrow(out), ' rows')
                print(Sys.time() - now)
            }
  
            if (any(nnix <<- out$cigar=='*')){
                out$cigar[nnix] = NA
            }

            if (data.table == TRUE){
                return(as.data.table(out)) ### check format
            } 
            else{
                parse_outputs(out) 
            }
        },

        grab_tag = function(query, tag, data.table = FALSE, verbose = FALSE, mc.cores = 1){      

            if (inherits(query, 'GRanges') | inherits(query, 'data.frame') | inherits(query, 'data.table')){

                if (is.null(query$index)){
                    stop(paste0('Error: Input query must include tag ', index))
                }

                query = query$tag ## vector of CB's 
            }

            query = setdiff(query, NA)

            ### does this work with a vector of CBs? 
            ## for (cb in 1:length(query)){
                ## get_cb = .C("print_cb_rows", self$bam, self$bamdb, cb)      ### would need to rbind() all outputs
            ## }

            get_tag = .C('print_tag_rows', self$bam, self$bamdb, query, tag)
         
            ##if (get_tag == 1){  ## throw error
            ##    stop(paste0('Error: there are no ', tag, ' indices in ', self$bamdb))   ## check C code---there should be several errors here; what if "one" CB tag is wrong?
            ##}
        
            if (verbose){
                now = Sys.time()
            }

            ## now, parse into data.table
            out = rbindlist(mclapply(get_tag, function(){
                out = dcast.data.table(dat, V1 ~ V2, value.var = "V3")[, -1, with = FALSE]
                setnames(out, names(out), ifelse(names(out)==as.character(tag), names(out), tolower(names(out))))
                return(out)
            }, mc.cores = mc.cores))

            if (nrow(out)>0){
                out$pos = as.numeric(out$pos)
                out$mapq = as.integer(out$mapq)
            }    

            if (verbose){ 
                message('Built out table with ', nrow(out), ' rows')
                print(Sys.time() - now)
            }

            if (any(nnix <<- out$cigar=='*')){
                out$cigar[nnix] = NA
            }

            if (data.table == TRUE){
                return(as.data.table(out)) ### check format
            } 
            else{
                parse_outputs(out)
            }
        }
    )
)


