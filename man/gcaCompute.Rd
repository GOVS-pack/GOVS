\name{gcaCompute}
\alias{gcaCompute}

\title{
GCA computation
}
\description{
Calculate parental general combining ability (GCA) based F1 phenotypic values.
}
\usage{
gcaCompute(phe_df,which,trait)
}

\arguments{
  \item{phe_df}{
Phenotypic data frame, row represents F1 combination and includes the paternal information in columns.
}
  \item{which}{
The column index of male or female to compute parternal or maternal GCA.
}
  \item{trait}{
A character string to define which trait GCA will be computed, this function support two or more phenotypic GCA be computed  at the same time.
}
}
\value{
A data frame involves of target trait GCA.
}

\author{
Qian Cheng, Shuqin Jiang, Xiangfeng Wang
}

\examples{
## simulation data
df <- data.frame(seqname = paste0("sample",1:500),
                 female_ID = sample(paste0("female",1:20), size = 500, replace = TRUE),
                 male_ID = sample(paste0("male",1:10), size = 500, replace = TRUE),
                 trait = rnorm(500, mean = 225, sd = 20))
plot(density(df$trait))

## compute GCA of female lines
gcaRes <- gcaCompute(df,which = "female_ID",trait = "trait")
}