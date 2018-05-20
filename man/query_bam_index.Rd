\name{query_bam_index}
\alias{query_bam_index}
\title{Fetch sequence rows matching a given index value}
\description{
  Wrapper around bamdb functions to query an indexed bam file. This requires
  a prebuilt index for the specific bam file that contains the desired column.
}
\usage{
  query_bam_index(bam_file_name, bam_index_path, index_name, key)
}
\arguments{
  \item{bam_file_name}{The path of the indexed sequence file}
  \item{bam_index_path}{The top level directory of the bamdb generated index}
  \item{index_name}{The name of indexed field to search in}
  \item{key}{The desired value of the indexed field}
}
\value{
  Returns a data frame with all of the matching sequnce rows on success. The 
  data frame will be empty of no results were found. If there is an error 
  with querying the index or file the result will be null.
}
