\name{reviseFunc}
\alias{reviseFunc}

\title{
Data correct function
}
\description{
Scale two sets of data to a uniform distribution.
}
\usage{
reviseFunc(ori,aim,cut = 10,sample_names)
}

\arguments{
  \item{ori}{
  A numeric array, as reference for correction.
}
  \item{aim}{
  A numeric array, which is the object to implement correction.
}
  \item{cut}{
  Number of intervals to cut, default 10.
}
  \item{sample_names}{
  A character array consits of the names of 'aim'.
}
}

\author{
Qian Cheng, Shuqin Jiang, Xiangfeng Wang
}

\examples{
## Not run!
## simulation data
df <- data.frame(seqname = paste0("sample",1:500),
                 female_ID = sample(paste0("female",1:20), size = 500,
                 replace = TRUE),
                 male_ID = sample(paste0("male",1:10), size = 500,
                 replace = TRUE),
                 trait = rnorm(500, mean = 225, sd = 20))
plot(density(df$trait))

## correct data distribution
df1 <- data.frame(seqname = paste0("sample",1:500),
                  trait = rnorm(500, mean = 225, sd = 20))
df2 <- data.frame(seqname = paste0("sample",1:300),
                  trait = rnorm(300, mean = 170, sd = 30))

## comparison of original distribution
plot(density(df1$trait),xlim = c(30,300),xlab = "Value range",main = "")
lines(density(df2$trait),col = "red")
legend(30,0.02,legend = c("df1","df2"),col = c("black","red"),lty = c(1,1))

## scale distribution
correct_df2 <- reviseFunc(sample_names = df2$seqname,ori = df1$trait,
                          aim = df2$trait,cut = 500)

## comparison of corrected distribution and original distribution
plot(density(df1$trait),xlim = c(30,300),xlab = "Value range",main = "")
lines(density(df2$trait),col = "red")
lines(density(correct_df2),col = "blue")
legend(30,0.02,legend = c("df1","df2","correct_df2"),
       col = c("black","red","blue"),lty = c(1,1,1))
}
