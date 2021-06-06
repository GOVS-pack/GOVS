\name{scaCompute}
\alias{scaCompute}

\title{
SCA computation
}
\description{
Calculate hybrid special combining ability (SCA) based F1 phenotypic values.
}
\usage{
scaCompute(phe_df,which_male,which_female,trait,seqname)
}

\arguments{
  \item{phe_df}{
Phenotypic data frame, row represents F1 combination and includes the paternal information in columns.
}
  \item{which_male}{
The column index of paternal IDs.
}
  \item{which_female}{
The column index of maternal IDs.
}
  \item{trait}{
A character string to define which trait GCA will be computed, this function support two or more phenotypic GCA be computed  at the same time.
}
  \item{seqname}{
A character array regarding hybrid IDs
}
}
\value{
A data frame involves of target trait SCA.
}

\author{
Qian Cheng, Shuqin Jiang, Xiangfeng Wang
}

\seealso{
\code{\link{gcaCompute}}
}
\examples{
## simulation data
df <- data.frame(seqname = paste0("sample",1:500),
                 female_ID = sample(paste0("female",1:20), size = 500, replace = TRUE),
                 male_ID = sample(paste0("male",1:10), size = 500, replace = TRUE),
                 trait = rnorm(500, mean = 225, sd = 20))
plot(density(df$trait))

## compute SCA of hybrid lines
scaRes <- scaCompute(df,which_female = "female_ID",which_male = "male_ID",
                     trait = "trait",seqname = df$seqname)
}