#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h> // for optparse
#include <Rcpp.h>
#include <iostream>

#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/hts.h"

#include "include/bamdb.h"
#include "include/bam_sqlite.h"
#include "include/bam_lmdb.h"
#include "include/bam_api.h"
#include "include/bamdb_c.h"

using namespace Rcpp;
using namespace std;

/* Command to do successfully run.
barcodedReads('/gpfs/commons/groups/imielinski_lab/data/10X/bam3/HCC1143_BL_phased_possorted.bam','/gpfs/commons/groups/imielinski_lab/data/10X/bam3/HCC1143_BL_phased_possorted.db','GTGGTCGCAACGCTTA-1')
*/
