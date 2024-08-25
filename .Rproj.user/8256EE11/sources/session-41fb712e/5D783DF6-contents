#' @title Epi-reaction_enrichment
#' @description Epi-reaction_enrichment to identify and plotting the significantly changed structural modification type of lipids in a disease.
#' @param result The result of the lipid_epimetabolite_matching function.
#' or differential features dataset: Please upload significant changed lipid epimetablites information you detected in a comparison of two conditions,
#' in the template format of "marker lipid epi-metabolites"
#' @param x_value Number of lipid metabolites identified by metabolic reaction (x value)
#' @param data2 EpiReactions: The information for the 94 epi-metabolic reactions that may occur in lipids included the reaction type,
#' reaction formula, the change of formula and mass after the reaction,
#' and the matching relationship between the lipid structure and reactions.
#'
#' @return list of plot and table
#' @export
#'
#' @examples
#' library(readxl)
#' library(Lipidepifind)
#' data_folder <- system.file("extdata", package = "Lipidepifind")
#' result <- read_xlsx(file.path(data_folder, "result.xlsx")) %>% as.data.table()
#' x_value <- read_xlsx(file.path(data_folder, "x value.xlsx")) %>% as.data.table()
#' data2 <- read_xlsx(file.path(data_folder, "Epireactions.xlsx")) %>% as.data.table()
#' epi_reaction_enrichment_list <- epi_reaction_enrichment(result,x_value,data2)
#' @references 1. Xinguang Liu, Ang Zhang, Yang Xie, Yanmin Shi, Jianya Yang, Yan Du, Jinyan Wu, Yitong Gao, Jiansheng Li.
#' Lipidepifind: A Workflow for the Characterization of Lipid Epi-metabolites as Diagnostic Biomarkers for Idiopathic Pulmonary Fibrosis, submitted.
epi_reaction_enrichment <- function(result,
                                    x_value,
                                    data2){
  p_value_df <- result %>%
    count(Reaction) %>%
    arrange(n) %>%
    inner_join(x_value, by = "Reaction") %>%
    mutate(
      p_value = 1 - phyper(n - 1, nrow(result), 3473, 174),
      log10_p_value = -log10(p_value),
      log10_p_value = if_else(is.infinite(log10_p_value),-1,log10_p_value)

    )

  plot <- ggplot(p_value_df, aes(n, Reaction, color = log10_p_value)) +
    geom_point() +
    labs(x = "Count", y = "Reaction", color = "-LOG10(Pvalue)") +
    scale_color_gradient(low = "blue", high = "red",
                         limits = c(min(p_value_df$log10_p_value),
                                    max(p_value_df$log10_p_value)))

  table <- p_value_df %>%
    inner_join(data2 %>% select(Reaction,`Enzyme/target`), by = "Reaction") %>%
    transmute(
      Reaction,
      Total = x,
      Hit = n,
      `P-value` = p_value,
      `Enzyme/target`
    )

  return(list(plot = plot, table = table))
}

