\name{hsq}
\alias{hsq}
\docType{data}
\title{Humor Styles}
\description{
The dataset was collected with an interactive online version of the Humor Styles Questionnaire (HSQ) which assesses four independent ways in which people express and appreciate humor (Martin et al. 2003): affiliative, defined as the benign uses of humor to enhance one's relationships with others; self-enhancing, indicating uses of humor to enhance the self; aggressive, the use of humor to enhance the self at the expense of others; self-defeating the use of humor to enhance relationships at the expense of oneself. The main part of the questionnaire consisted of 32 statements rated from 1 to 5 according to the respondents' level of agreement. Three more questions were included (age, gender and self-reported accuracy of answer). The number of respondents is 993, after removing the cases with missing values in the 32 statements.}

\usage{data("hsq")}
\format{
   A data frame with 993 observations on 35 variables. The first 32 variables are Likert-type statements with 5 response categories, ranging from 1 (strong agreement) to 5 (strong disagreement).
  \describe{
    \item{\code{AF1}}{I usually don't laugh or joke around much with other people}
    \item{\code{AF2}}{If I am feeling depressed, I can usually cheer myself up with humor}
    \item{\code{AF3}}{If someone makes a mistake, I will often tease them about it}
    \item{\code{AF4}}{I let people laugh at me or make fun at my expense more than I should}
    \item{\code{AF5}}{I don't have to work very hard at making other people laugh - I seem to be a naturally humorous person}
    \item{\code{AF6}}{Even when I'm by myself, I'm often amused by the absurdities of life}
    \item{\code{AF7}}{People are never offended or hurt by my sense of humor}
    \item{\code{AF8}}{I will often get carried away in putting myself down if it makes my family or friends laugh}
\item{\code{SE1}}{I rarely make other people laugh by telling funny stories about myself}
\item{\code{SE2}}{If I am feeling upset or unhappy I usually try to think of something funny about the situation to make myself feel better}
\item{\code{SE3}}{When telling jokes or saying funny things, I am usually not very concerned about how other people are taking it}
\item{\code{SE4}}{I often try to make people like or accept me more by saying something funny about my own weaknesses, blunders, or faults}
\item{\code{SE5}}{I laugh and joke a lot with my closest friends}
\item{\code{SE6}}{My humorous outlook on life keeps me from getting overly upset or depressed about things}
\item{\code{SE7}}{I do not like it when people use humor as a way of criticizing or putting someone down}
\item{\code{SE8}}{I don't often say funny things to put myself down}
\item{\code{AG1}}{I usually don't like to tell jokes or amuse people}
\item{\code{AG2}}{If I'm by myself and I'm feeling unhappy, I make an effort to think of something funny to cheer myself up}
\item{\code{AG3}}{Sometimes I think of something that is so funny that I can't stop myself from saying it, even if it is not appropriate for the situation}
\item{\code{AG4}}{I often go overboard in putting myself down when I am making jokes or trying to be funny}
\item{\code{AG5}}{I enjoy making people laugh}
\item{\code{AG6}}{If I am feeling sad or upset, I usually lose my sense of humor}
\item{\code{AG7}}{I never participate in laughing at others even if all my friends are doing it}
\item{\code{AG8}}{When I am with friends or family, I often seem to be the one that other people make fun of or joke about}
\item{\code{SD1}}{I don't often joke around with my friends}
\item{\code{SD2}}{It is my experience that thinking about some amusing aspect of a situation is often a very effective way of coping with problems}
\item{\code{SD3}}{If I don't like someone, I often use humor or teasing to put them down}
\item{\code{SD4}}{If I am having problems or feeling unhappy, I often cover it up by joking around, so that even my closest friends don't know how I really feel}
\item{\code{SD5}}{I usually can't think of witty things to say when I'm with other people}
\item{\code{SD6}}{I don't need to be with other people to feel amused - I can usually find things to laugh about even when I'm by myself}
\item{\code{SD7}}{Even if something is really funny to me, I will not laugh or joke about it if someone will be offended}
\item{\code{SD8}}{Letting others laugh at me is my way of keeping my friends and family in good spirits}
  }
}
\source{
Martin, R. A., Puhlik-Doris, P., Larsen, G., Gray, J., & Weir, K. (2003). Individual differences in uses of humor and their relation to psychological well-being: Development of the Humor Styles Questionnaire. \emph{Journal of Research in Personality}, \emph{37}(1), 48-75.}

\examples{
data(hsq)
}
\keyword{datasets}