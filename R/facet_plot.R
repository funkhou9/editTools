#' Turns two or more plots of edit_tables into a facet plot
#'
#' 
#' @param ... any number data.frames resulting from plot.edit_table(plot = FALSE)
#' @param names character vector for each sample provided in 
#' @param text_size numeric changing the appearence of text in the plot
#' @param group character providing the 'fill' aesthetic argument fo each facet
#' @return a gg object
#' @export
facet_plot <- function(...,
                       names = NULL,
                       text_size = 20,
                       percent_size = 10,
                       group = "Tissue") {
  
  comb_dat <- list(...) %>% 
                do.call(rbind, .)
  
  if (!is.null(names))
    levels(comb_dat$field) <- names
  
  if (group != "none") {
    g <- ggplot(comb_dat,
                aes(x = reorder(Event, -Freq),
                    y = Freq,
                    fill = reorder(group, Freq)))
    g <- g + geom_bar(stat = 'identity')
    
  } else {
    g <- ggplot(comb_dat,
                aes(x = reorder(Event, Freq),
                    y = Freq))
    g <- g + geom_bar(stat = "identity",
                      fill = I("darkgreen"))
  }
  
  
  
  g <- g + geom_text(aes(label = Total_prop),
                         na.rm = TRUE,
                         position = "stack",
                         hjust = 0.5,
                         vjust= -0.2,
                         size = percent_size)
  g <- g + ylab("Number of events")
  g <- g + xlab("Type of mismatch")
  g <- g + theme(axis.title.x = element_text(size = text_size),
                 axis.title.y = element_text(size = text_size),
                 axis.text.x = element_text(angle = 45,
                                            hjust = 1,
                                            size = text_size),
                 axis.text.y = element_text(size = text_size),
                 legend.text = element_text(size = text_size),
                 legend.title = element_text(size = text_size),
                 legend.key.height = unit(3, "line"),
                 strip.text.x = element_text(size = text_size),
                 text = element_text(family = "Arial"))
  g <- g + guides(fill = guide_legend(title = "Tissues"))
  g <- g + facet_wrap(~ field,
                      nrow = 1,
                      scale = "free_x")
  
  
  return (g)
}