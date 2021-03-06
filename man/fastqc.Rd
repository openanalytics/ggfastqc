% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastqc-tables.R
\name{fastqc}
\alias{fastqc}
\title{Extract FastQC summary tables}
\usage{
fastqc(sample_info)
}
\arguments{
\item{sample_info}{Full path to file containing details about 
\code{samples} and their \code{paths}.}
}
\value{
An object of class \code{fastqc} which inherits from 
\code{"data.table"}, with two columns: \code{param} and \code{value}, where 
\code{value} is a list of \code{"data.table"}s.
}
\description{
\code{fastqc} function parses all the summary statistics of 
reports produced by \code{FastQC} tool and returns a \code{data.table} with 
two columns: \code{param} and \code{value}. 

Each row of the \code{value} column contains the data corresponding to that 
\code{param}, and is itself 
a \code{data.table}.
}
\details{
The file provided to \code{sample_info} argument should contain 
at least these three columns: 
\itemize{
  \item \code{sample} - contains the \code{sample} name.
  \item \code{pair} - in case of paired end reads, \code{1} or \code{2} 
  corresponding to first and second pair, and in case of single end reads, 
  \code{NA}.
  \item \code{path} - full path to the fastqc summary report (\code{.txt} 
  file) for each sample.
   
  If just the file name (\code{.txt}) is provided, it is assumed that the 
  file is in the same folder as the input file provided to 
  \code{sample_info} argument.
}

It can also optionally contain a \code{group} column. If present, the plots 
generated will take it into account and \code{color} / \code{facet} 
accordingly.
}
\examples{
path = system.file("tests/fastqc-sample", package="ggfastqc")
obj = fastqc(sample_info = file.path(path, "annotation.txt"))
}
\seealso{
\code{\link{plot_dup_stats}} \code{\link{plot_gc_stats}} 
\code{\link{plot_sequence_quality}} \code{\link{plot_total_sequence_stats}}
}

