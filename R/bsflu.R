##' bsflu
##'
##' An outbreak of influenza in an all-boys boarding school.
##'
##' Data are recorded from a 1978 flu outbreak in a closed population.
##' The variable \sQuote{B} refers to boys confined to bed on the corresponding day and \sQuote{C} to boys in convalescence,
##' i.e., not yet allowed back to class.
##' In total, 763 boys were at risk of infection and, over the course of the outbreak, 512 boys spent between 3 and 7 days away from class (either in bed or convalescent).
##' The index case was a boy who arrived at school from holiday six days before the next case.
##'
##' @name bsflu
##' @aliases bsflu
##' @rdname bsflu
##' @family datasets
##' @seealso \link{sir_models}
##'
##' @references
##' Anonymous (1978).
##' Influenza in a boarding school.
##' British Medical Journal 1:587
##'
##' @examples
##'
##' library(magrittr)
##' library(tidyr)
##' library(ggplot2)
##'
##' bsflu %>%
##'   gather(variable,value,-date,-day) %>%
##'   ggplot(aes(x=date,y=value,color=variable))+
##'   geom_line()+
##'   labs(y="number of boys",title="boarding school flu outbreak")+
##'   theme_bw()
##'
NULL
