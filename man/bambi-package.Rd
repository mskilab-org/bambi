\name{bambi-package}
\alias{bambi-package}
\alias{bambi}
\docType{package}
\title{
  An R interface for querying indexed binary sequence data
}
\description{
  This package provides functions to query binary alignment
  \url{https://samtools.github.io/hts-specs/SAMv1.pdf} files that have been
  indexed with the bamdb software package.
}
\author{
Maintainer: Aaron Dornbrand-Lo
}
\keyword{ package }
\keyword{ interface }
\seealso{
  The bamdb library that this package extends currently lives at \url{https://github.com/D-Lo/bamdb}.
}
\examples{
  \dontrun{
  results = query_bam_index("sample_bam_file.bam", "/home/adornbrandlo/sample_bam_index", "BX", "CCGATAAGTGATGTGG-1")
  }
}
