#############################################################################
## Evan Biederstedt
## New York Genome Center
## ebiederstedt@nygenome.org
##
## Marcin Imielinski
## The Broad Institute of MIT and Harvard / Cancer program.
## marcin@broadinstitute.org
## Weill-Cornell Medical College
## mai9037@med.cornell.edu
## New York Genome Center
## mimielinski@nygenome.org
##
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
## You should have received a copy of the GNU Lesser General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##############################################################################

#' @import rPython
#' @import data.table
#' @import Rsamtools
#' @import GenomicRanges
#' @import GenomicAlignments

#' @name countCigar
#' @title countCigar
#' @description
#'
#' Count bases in cigar string
#'
#' Counts the total number of bases, per cigar, that fall into D, I, M, S categories.
#' countCigar makes no distinction between, for instance 1S2M2S, 2S2M1S, or 3S2M
#' @param cigar character vector of cigar strings
#' @return a 4-column, length(cigar)-row matrix with the total counts for each type
#' @export
countCigar <- function(cigar) {
    
    cigar.vals <- unlist(strsplit(cigar, "\\d+"))
    cigar.lens <- strsplit(cigar, "[A-Z]")
    lens <- nchar(gsub('\\d+', '', cigar))
    lens[is.na(cigar)] <- 1
    
    cigar.lens <- as.numeric(unlist(cigar.lens))
    cigar.vals <- cigar.vals[cigar.vals != ""]
    repr       <- rep(seq_along(cigar), lens)
    dt         <- data.table(val=cigar.vals, lens=cigar.lens, group=repr, key="val")
    
    smr.d      <- dt["D",][, sum(lens), by=group]
    smr.i      <- dt["I",][, sum(lens), by=group]
    smr.m      <- dt["M",][, sum(lens), by=group]
    smr.s      <- dt["S",][, sum(lens), by=group]
    
    out <- matrix(nrow=length(cigar), ncol=4, 0)
    out[smr.d$group,1] <- smr.d$V1
    out[smr.i$group,2] <- smr.i$V1
    out[smr.m$group,3] <- smr.m$V1
    out[smr.s$group,4] <- smr.s$V1
    colnames(out) <- c('D','I','M','S')
    
    return(out)
}


#' @name bxbam-class
#' @title bxbam-class
#' @rdname bxbam-class
#' @description
#'
#' class::bxbam
#'
#' Class \code{bxbam} object to load bxbam and query GRanges
#'
#' @import rPython
#' @import data.table
#' @import Rsamtools
#' @import GenomicRanges
#' @import GenomicAlignments
#' @exportClass bxBam
#' @author Evan Biederstedt
setClass("bxBam", representation(.bxbamfile = 'character', .bamfile = 'BamFile', .sessionId = 'character', .sqllite = 'logical'))

setMethod('initialize', 'bxBam', function(.Object, bxbamfile = '', bamfile = '', tags = NULL, index_on = NULL, chunksize = 1e6, verbose = TRUE, overwrite = FALSE, nlimit = NULL, mc.cores = 1)
{
    sqllite = FALSE
    .Object@.sessionId <- paste0('session', runif(1))
    message(.Object@.bxbamfile)

    if (!file.exists(bamfile)) # will try and guess the bamfile name from the bxbam
    {
        if (!file.exists(bxbamfile))
            stop('Either bam or bxbam (.sqllite, .h5) file must be provided')
        
        bamfile = gsub('bxbamfile', 'bxbam', bxbamfile)
        if (!file.exists(bamfile))
            bamfile = paste0(bamfile, '.bam')

    }
        
    if (!is(bamfile, 'BamFile'))
        bamfile =  BamFile(bamfile)

    .Object@.bamfile = bamfile
    
    if (nchar(bxbamfile) == 0)
        bxbamfile = NULL
    
    if (is.null(bxbamfile)) ## if bxbam is not provided will create automatic name from bam path
    {       
        bxbamfile = gsub('.bam$', '.sqllite', bamfile)
        sqllite = TRUE           
    }
    else
        sqllite = grepl('.sqllite$', bxbamfile)
    
    .Object@.bxbamfile = bxbamfile
    
    if (!sqllite)
        python.exec(sprintf("sessions['%s'] = tables.open_file('%s').get_node('/bam_table/bam_fields')", .Object@.sessionId, .Object@.bxbamfile))
    else
    {
        if (!file.exists(bxbamfile) | overwrite)
        {
            if (file.exists(bxbamfile))
                {
                    message(paste('About to overwrite', bxbamfile, 'giving you a chance to think about it'))
                    Sys.sleep(1)
                }
            message('Creating .sqllite from .bam file')
            tags = union(c('MD', 'BX'), tags)
            index_on = union(c('qname', 'BX', 'rname', 'rnext'), index_on)                    
            bam2sqllite(.Object@.bxbamfile, Rsamtools::path(.Object@.bamfile), tags = tags, index_on = index_on, chunksize = chunksize, verbose = verbose, nlimit = nlimit, mc.cores = mc.cores)
        }
    }
        
    .Object@.sqllite = sqllite
    return(.Object)
})


#' @name bxBam
#' @title bxBam
#' @description
#' Initialize bxBam object specifying .bxbam file and optional .bam path
#' @export
#' @author Marcin Imielinski
bxBam = function(bxbamfile = '', bamfile = '', tags = NULL, index_on = NULL, chunksize = 1e6, verbose = TRUE, overwrite = FALSE, nlimit = NULL, mc.cores = 1)
    new('bxBam', bxbamfile = bxbamfile, bamfile = bamfile,
        tags = tags, index_on = index_on, chunksize = chunksize, verbose = verbose,
        overwrite = overwrite, nlimit = nlimit, mc.cores = mc.cores)


#' @name bam2sqllite
#' @title bam2sqllite
#' @description
#' Pipe output from samtools and make indexed sqllite table +/- drawing specific optional tags from bam file
#' and +/- creating optional indices
#' @export
#' @author Marcin Imielinski
bam2sqllite = function(sqllite_path, bam_path, tags = NULL, index_on = NULL, chunksize = 1e6, verbose = FALSE, nlimit = NULL, mc.cores = 1)
{
    if (!is.null(nlimit))
        verbose = TRUE
    
    if (verbose)
        message('Creating brand new SQLLite db in ', sqllite_path)

    create_str = sprintf("
CREATE TABLE reads (
 qname VARCHAR,
 flag INTEGER,
 rname VARCHAR,
 pos INTEGER,
 mapq INTEGER,
 cigar VARCHAR,
 rnext VARCHAR,
 pnext INTEGER,
 tlen INTEGER,
 seq VARCHAR,
 qual VARCHAR,
 %s
);", paste(tags, 'varchar', collapse = ','))
    
    system(paste('rm -rf', sqllite_path))
    mydb <- dbConnect(RSQLite::SQLite(), sqllite_path)
    index_cols = c('qname', 'BX', 'rname', 'rnext')
    dbExecute(mydb, create_str)
           
    if (verbose)
        message('Table created with indices on ', paste(index_cols, collapse = ', '), ' on additional tags ', paste(tags, collapse = ','),
                ' using SQL commands:', paste(c(create_str), collapse = ';\n'))
    
    p = pipe(paste('samtools view', bam_path), open = 'r')
    fields = c('qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual')
    tags = c('MD', 'BX')

  #   fandt = union(fields, tags)
  #    insert_frag = sprintf("INSERT INTO reads (%s) VALUES", 
  #                             paste(fandt, collapse = ','))


    begin = Sys.time()

    ## this DBI dbWithTransaction statement should cause the dbWriteTable
    ## statements to be inside a transaction
    DBI::dbWithTransaction(mydb,
    {
        while (length(lines <- readLines(p, n = chunksize))>0)
        {
            now = Sys.time()
            lchunks = split(lines, rep(1:mc.cores, ceiling(length(lines)/mc.cores)))
            chunk = rbindlist(mclapply(lchunks,
                                     .proclines, fields = fields, tags = tags, verbose = TRUE,
                                     mc.cores = mc.cores))
                                          
            if (verbose)
                message('\nprocessed ', nrow(chunk), ' from bam file, now writing to SQLlite')
            

            ## MUCH slower than dbWriteTable call below 
            ##     insert_statement = paste0(insert_frag, " ('",
            ##                             do.call(paste, c(as.list(chunk), list(sep = "','"))),
            ##                            "')")
            ##        bla = lapply(insert_statement, dbExecute, conn = mydb)

            now.write = Sys.time()
            dbWriteTable(mydb, "reads", chunk, append = TRUE, overwrite = FALSE)
            if (verbose)
            {
                nlines = dbGetQuery(mydb, 'SELECT COUNT(*) FROM reads')
                nsecs = as.numeric(difftime( Sys.time(), now, units = 'secs'))
                nsecs.write = as.numeric(difftime( Sys.time(), now.write, units = 'secs'))
                nsecs.begin = as.numeric(difftime( Sys.time(), begin, units = 'secs'))                
                message(prettyNum(nlines, big.mark = ','), ' records in table')
                print(Sys.time()-now)
                message(sprintf(' .. since beginning this chunk \n\t\t(parse + SQL write speed: %s records / second, \n\t\t SQL write speed: %s records / second)',
                                prettyNum(round(nrow(chunk) / nsecs, 2)),
                                prettyNum(round(nrow(chunk) / nsecs.write, 2))))
                print(Sys.time()-begin)
                message(sprintf(' .. since file creation\n\t\t(overall speed: %s records / second)',
                                prettyNum(round(nlines / nsecs.begin,2))))

                if (!is.null(nlimit))
                    if (nlines>nlimit)
                        break()
            }
        }
    })

    close(p)
    create_index_str = sapply(index_cols, function(x)
    {
        str = sprintf("CREATE INDEX %s ON reads(%s)", x,x)
        dbExecute(mydb, str)
        return(str)        
    })

    if (verbose)
        message('Table created with indices on ', paste(index_cols, collapse = ', '), ' on additional tags ', paste(tags, collapse = ','),
                ' using SQL commands:', paste(c(create_index_str), collapse = ';\n'))

                   
}

## process chunk of lines, used by bam2sqllite
.proclines = function(lines, fields, tags, verbose = FALSE)
{
    if (verbose)
        cat(".")
    linesp = strsplit(lines, '\t')
    chunk = as.data.table(do.call(rbind, lapply(linesp, function(x) x[1:11])))[, line := 1:length(V1)]
    m = munlist(lapply(linesp, function(x) x[-c(1:11)]))
    tagchunk = fread(paste(m[,3], collapse = '\n'), sep = ':')
    tagchunk[, line := as.numeric(m[,1])]
    tagchunk = dcast.data.table(tagchunk[V1 %in% tags, ], line ~ V1, value.var = 'V3')
    chunk = merge(chunk, tagchunk, by = 'line')[, -1, with = FALSE]
    setnames(chunk, 1:11, fields)
    chunk = chunk[, c(fields, tags), with = FALSE]
    if (verbose)
        cat("|")
    return(chunk)    
}

#' @name head
#' @title gets head of file
#' 
#' @exportMethod head
#' @export
#' @import RSQLite
#' @author Marcin Imielinski
setMethod('head', 'bxBam', function(x, n = 5)
{
    if (!.hasSlot(x, '.sqllite')) ## check for older version of sqllite
        sqllite = FALSE
    else
        sqllite = x@.sqllite

    if (sqllite)
    {
        mydb <- RSQLite::dbConnect(RSQLite::SQLite(), x@.bxbamfile)
        return( RSQLite::dbGetQuery(mydb, sprintf('SELECT * FROM reads LIMIT %s', n)))
    }
})


#' @name show
#' @title show
#' @description Display a \code{gTrack} object
#' @docType methods
#' @param object \code{gTrack} to display
#' @author Marcin Imielinski
setGeneric('reindex', function(.Object) standardGeneric('reindex'))
setMethod("reindex", "bxBam", function(.Object)
{    
    mydb <- RSQLite::dbConnect(RSQLite::SQLite(), .Object@.bxbamfile)
    DBI::dbExecute(mydb, 'REINDEX main.reads')
})

setValidity("bxBam", function(object){
    if (!file.exists(bxbamfile)) stop ("'bxbamfile' is missing; please set correct path to bxbam.h5 object")
})

#' @name show
#' @title show
#' @description Display a \code{gTrack} object
#' @docType methods
#' @param object \code{gTrack} to display
#' @author Marcin Imielinski
setMethod('show', 'bxBam', function(object)
{
    cat(sprintf('bxBam object stored in data file\n\t%s\nwith associated BAM file\n\t%s:\n', object@.bxbamfile, Rsamtools::path(object@.bamfile)))
    print(head(object))
})

#' @name get_bmates
#' @title get_bmates
#' @description
#' Grab all barcoded reads by bam field BX, return GRanges object
#'
#' @export
#' @author Evan Biederstedt
#' @author Marcin Imielinski
setGeneric('get_bmates', function(.Object, query, ...) standardGeneric('get_bmates'))
setMethod("get_bmates", "bxBam", function(.Object, query, verbose = FALSE){
    
    if (!.hasSlot(.Object, '.sqllite')) ## check for older version of bxbam
        sqllite = FALSE
    else
        sqllite = .Object@.sqllite
    
    if (inherits(query, 'GRanges') | inherits(query, 'data.frame'))
    {
        if (is.null(query$BX))
        {            
            if (verbose)
                {
                    message("BX field not found, will use read.bam to pull reads under query GRanges from bam file and find their bmates")
                    if (verbose)
                        now = Sys.time()
                }
            query = read.bam(.Object@.bamfile, gr = query, tag = c('BX', 'MD'), pairs.grl = FALSE)
            if (verbose)
                {
                    message('Retrieved reads:')
                    print(Sys.time()-now)
                }
        }    
        query = query$BX            
    }

    if (verbose)
        now = Sys.time()
    
    if (sqllite)
        {
            ## ahh how easy!
            mydb <- RSQLite::dbConnect(RSQLite::SQLite(), .Object@.bxbamfile)
            out = as.data.table(dbGetQuery(mydb, sprintf('SELECT * FROM reads WHERE BX in (%s)', paste0('"', query, '"', collapse = ','))))
        }
    else
        {            
            ## check if python session id exists and if not create
            python.exec( sprintf("
        if '%s' not in sessions.keys(): sessions['%s'] = tables.open_file('%s').get_node('/bam_table/bam_fields')",        
        .Object@.sessionId,
        .Object@.sessionId,
        .Object@.bxbamfile))
            

            if (any(nix <<- is.na(query)))
                query = query[!nix]

            if (length(query)==0)
                stop('Length 0 query, check input')

            qstring = paste(paste0('(BX==b"', query, '")'), collapse = "|")
            queryId = paste0('query', runif(1))

            query = sprintf("queries['%s'] = run_query(sessions['%s'], '%s')", queryId, .Object@.sessionId, qstring)
            if (verbose)
                message('Running query: ', query)
            python.exec(query)
            
            out = data.table(
                bx = python.get(sprintf("queries['%s'].BX.tolist()", queryId)),
                cigar = python.get(sprintf("queries['%s'].CIGAR.tolist()", queryId)),
                flag = python.get(sprintf("queries['%s'].FLAG.tolist()", queryId)),
                mapq = python.get(sprintf("queries['%s'].MAPQ.tolist()", queryId)),
                pnext = python.get(sprintf("queries['%s'].PNEXT.tolist()", queryId)),
                pos = python.get(sprintf("queries['%s'].POS.tolist()", queryId)),
                qname = python.get(sprintf("queries['%s'].QNAME.tolist()", queryId)),
                qual = python.get(sprintf("queries['%s'].QUAL.tolist()", queryId)),
                rname = python.get(sprintf("queries['%s'].RNAME.tolist()", queryId)),
                rnext = python.get(sprintf("queries['%s'].RNEXT.tolist()", queryId)),
                seq = python.get(sprintf("queries['%s'].SEQ.tolist()", queryId)),
                tlen = python.get(sprintf("queries['%s'].TLEN.tolist()", queryId)))
        }

    if (nrow(out)==0)
        return(out)
            
    if (verbose)
        print(Sys.time()-now)
    
    if (verbose)
        message('Built out table with ', nrow(out), ' rows')
    
    if (any(nnix <<- out$cigar=='*'))
        out$cigar[nnix] = NA

    if (nrow(out)==0)
        return(GRanges())
    
    cigs <- countCigar(out$cigar)
    out$pos2 <- out$pos + rowSums(cigs[, c("D", "M")], na.rm=T) - 1
    
    
    out$qwidth = nchar(out$seq)
    out$strand = bamflag(out$flag)[, "isMinusStrand"] == 1
    out$strand = ifelse(out$strand, "-", "+")

    
    unmapped = bamflag(out$flag)[,"isUnmappedQuery"] == 1
    if (any(unmapped))
    {
        out$pos[unmapped] = 1
        out$pos2[unmapped] = 0
        out$strand[unmapped] = "*"
    }
    
    bf = out$flag
    
    newdt <- data.table(pos = out$pos, pos2 = out$pos2, strand = out$strand, rname = out$rname)   # create data.table of start, end, strand, seqnames
    
    rr <- IRanges(newdt$pos, newdt$pos2)
    sf <- factor(newdt$strand, levels=c('+', '-', '*'))
    ff <- factor(newdt$rname, levels=unique(newdt$rname))
    gr.fields <- c("rname", "strand", "pos", "pos2")
    grobj <- GRanges(seqnames=ff, ranges=rr, strand=sf)
    
    vals = out[, setdiff(names(out), gr.fields), with=FALSE]
    values(grobj) <- vals
    return(grobj)
})



#' @name get_qmates
#' @title get_qmates
#' @description
#' Grab all barcoded reads by bam field QNAME, return GRanges object
#'
#' @export
#' @author Evan Biederstedt
#' @author Marcin Imielinski
setGeneric('get_qmates', function(.Object, query, ...) standardGeneric('get_qmates'))
setMethod("get_qmates", "bxBam", function(.Object, query, verbose = FALSE){

    if (!.hasSlot(.Object, '.sqllite')) ## check for older version of bxbam
        sqllite = FALSE
    else
        sqllite = .Object@.sqllite
    
    if (inherits(query, 'GRanges') | inherits(query, 'data.frame'))
    {
        if (is.null(query$BX))
        {            
            if (verbose)
            {
                message("BX field not found, will use read.bam to pull reads under query GRanges from bam file and find their bmates")
                if (verbose)
                    now = Sys.time()
            }
            query = read.bam(.Object@.bamfile, gr = query, tag = c('BX', 'MD'), pairs.grl = FALSE)
            if (verbose)
            {
                message('Retrieved reads:')
                print(Sys.time()-now)
            }
        }    
        query = query$qname
    }
    
    if (verbose)
        now = Sys.time()
    
    if (sqllite)
    {
        ## ahh how easy!
        mydb <- RSQLite::dbConnect(RSQLite::SQLite(), .Object@.bxbamfile)
        out = as.data.table(dbGetQuery(mydb, sprintf('SELECT * FROM reads WHERE qname in (%s)', paste0('"', query, '"', collapse = ','))))
    }
    else
        {
            ## check if python session id exists and if not create
            python.exec( sprintf("
        if '%s' not in sessions.keys(): sessions['%s'] = tables.open_file('%s').get_node('/bam_table/bam_fields')",        
        .Object@.sessionId,
        .Object@.sessionId,
        .Object@.bxbamfile))            

            if (inherits(query, 'GRanges') | inherits(query, 'data.frame'))
            {
                if (is.null(query$QNAME))
                {
                    warning("QNAME field not found, will use read.bam to pull reads under query GRanges from bam file and find their bmates")
                    query = read.bam(.Object@.bamfile, gr = query, tag = c('BX', 'MD'), pairs.grl = FALSE)    
                }      
                query = query$QNAME
            }
            
            if (any(nix <<- is.na(query)))
                query = query[!nix]

            if (length(query)==0)
                stop('Length 0 query, check input')

            
            qstring = paste(paste0('(QNAME==b"', query, '")'), collapse = "|")
            queryId = paste0('query', runif(1))
            python.exec(sprintf("queries['%s'] = run_query(sessions['%s'], '%s')", queryId, .Object@.sessionId, qstring))
            
            out = data.table(
                bx = python.get(sprintf("queries['%s'].BX.tolist()", queryId)),
                cigar = python.get(sprintf("queries['%s'].CIGAR.tolist()", queryId)),
                flag = python.get(sprintf("queries['%s'].FLAG.tolist()", queryId)),
                mapq = python.get(sprintf("queries['%s'].MAPQ.tolist()", queryId)),
                pnext = python.get(sprintf("queries['%s'].PNEXT.tolist()", queryId)),
                pos = python.get(sprintf("queries['%s'].POS.tolist()", queryId)),
                qname = python.get(sprintf("queries['%s'].QNAME.tolist()", queryId)),
                qual = python.get(sprintf("queries['%s'].QUAL.tolist()", queryId)),
                rname = python.get(sprintf("queries['%s'].RNAME.tolist()", queryId)),
                rnext = python.get(sprintf("queries['%s'].RNEXT.tolist()", queryId)),
                seq = python.get(sprintf("queries['%s'].SEQ.tolist()", queryId)),
                tlen = python.get(sprintf("queries['%s'].TLEN.tolist()", queryId)))
        }

    if (nrow(out)==0)
        return(out)    
    
    if (verbose)
        print(Sys.time()-now)
    
    if (any(nnix <<- out$cigar=='*'))
    out$cigar[nnix] = NA
    cigs <- countCigar(out$cigar)
    out$pos2 <- out$pos + rowSums(cigs[, c("D", "M")], na.rm=T) - 1
    
    out$qwidth = nchar(out$seq)
    out$strand = bamflag(out$flag)[, "isMinusStrand"] == 1
    out$strand = ifelse(out$strand, "-", "+")
    
    
    unmapped = bamflag(out$flag)[,"isUnmappedQuery"] == 1
    if (any(unmapped))
    {
        out$pos[unmapped] = 1
        out$pos2[unmapped] = 0
        out$strand[unmapped] = "*"
    }
    
    bf = out$flag
    
    newdt <- data.table(pos = out$pos, pos2 = out$pos2, strand = out$strand, rname = out$rname)   # create data.table of start, end, strand, seqnames
    
    rr <- IRanges(newdt$pos, newdt$pos2)
    sf <- factor(newdt$strand, levels=c('+', '-', '*'))
    ff <- factor(newdt$rname, levels=unique(newdt$rname))
    gr.fields <- c("rname", "strand", "pos", "pos2")
    grobj <- GRanges(seqnames=ff, ranges=rr, strand=sf)
    
    vals = out[, setdiff(names(out), gr.fields), with=FALSE]
    values(grobj) <- vals
    return(grobj)
})




#############################################################
#' @name munlist
#' @title munlist
#'
#' @description
#' unlists a list of vectors, matrices, data frames into a n x k matrix
#' whose first column specifies the list item index of the entry
#' and second column specifies the sublist item index of the entry
#' and the remaining columns specifies the value(s) of the vector
#' or matrices.
#'
#' force.cbind = T will force concatenation via 'cbind'
#' force.rbind = T will force concatenation via 'rbind'
#'
#' @param x list of vectors, matrices, or data frames
#' @param force.rbind logical flag to force concatenation via rbind (=FALSE), otherwise will guess
#' @param force.cbind logical flag to force concatenation via cbind (=FALSE), otherwise will guess
#' @param force.list logical flag to force concatenation via unlist (=FALSE), otherwise will guess
#' @return data.frame of concatenated input data with additional fields $ix and $iix specifying the list item and within-list index from which the given row originated from
#' @author Marcin Imielinski9
#' @export
#############################################################
munlist = function(x, force.rbind = F, force.cbind = F, force.list = F)
  {
    if (!any(c(force.list, force.cbind, force.rbind)))
      {
        if (any(sapply(x, function(y) is.null(dim(y)))))
          force.list = T
        if (length(unique(sapply(x, function(y) dim(y)[2]))) == 1)
          force.rbind = T
        if ((length(unique(sapply(x, function(y) dim(y)[1]))) == 1))
          force.cbind = T
      }
    else
      force.list = T

    if (force.list)
      return(cbind(ix = unlist(lapply(1:length(x), function(y) rep(y, length(x[[y]])))),
                   iix = unlist(lapply(1:length(x), function(y) if (length(x[[y]])>0) 1:length(x[[y]]) else NULL)),
                   unlist(x)))
    else if (force.rbind)
      return(cbind(ix = unlist(lapply(1:length(x), function(y) rep(y, nrow(x[[y]])))),
                   iix = unlist(lapply(1:length(x), function(y) if (nrow(x[[y]])>0) 1:nrow(x[[y]]) else NULL)),
                   do.call('rbind', x)))
    else if (force.cbind)
      return(t(rbind(ix = unlist(lapply(1:length(x), function(y) rep(y, ncol(x[[y]])))),
                     iix = unlist(lapply(1:length(x), function(y) if (ncol(x[[y]])>0) 1:ncol(x[[y]]) else NULL)),
                   do.call('cbind', x))))
  }
