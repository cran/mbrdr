\name{mps}
\alias{mps}
\title{Minneapolis School dataset}
\description{
The Minneapolis school dataset was collected to evaluate the performance of student The percentages of students in 63 Minneapolis schools in 1972. And, The dataset was reported in Star-Tribune in 1973.}

\usage{
data(mps)
}
\format{
A data frame of dimension is 63 x 15. Each row represents one elementary school. The first four columns correspond to  percentages of students in a school scoring above (A) and below (B) average on standardized fourth and sixth grade reading comprehension tests. Subtracting either pair of grade specific percentages from 100 gives the percentage of students scoring about average on the test. All the other variables are demographic informations for each school.
}

\details{
A4 = percentage of 4th graders scoring ABOVE average on a standard 4th grade vocabulary test in 1972.

B4 = percentage of 4th graders scoring BELOW average on a standard 4th grade vocabulary test in 1972.

A6 = percentage of 6th graders scoring BELOW average on a standard 6th grade comprehension test in 1972.

B6 = percentage of 6th graders scoring BELOW average on a standard 6th grade comprehension test in 1972.

AFDC = percentage of children receiving Aid to Families with Dependent Children

Attend = average percentage of childern in attendance during the year

B = percentage of children in the school not living with Both Parents

BthPts = percentage of children in the school living with Both Parents

Enrol = number of childeren enrolled in the school

HS = percent of adults in the school area who have completed high school

Minority = percent minority children in the area.

Mobility = percentage of children who started in a school, but did not finish there

Poverty = percentage of persons in the school area who are above the federal poverty levels

PTR = pupil-teacher ratio

School = names of school
}

\references{
Cook, R. D. and Setodji, C. M. (2003) A model-free test for reduced rank in multivariate regression. Journal of the American Statistical Association, 98, pp. 340-351.

JK. Yoo (2019) Unstructured principal fitted response reduction in multivariate regression. Journal of the Korean Statistical Society, 48, pp. 561-567.
}

\examples{
data(mps)
pairs(mps[,1:4])
}


