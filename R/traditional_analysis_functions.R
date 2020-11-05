#' Title
#'
#' @param ps_anspo_hc_rare
#'
#' @return
#' @export
#'
#' @examples
alpha_plot <- function(ps_anspo_hc_rare) {
  plot_richness(ps_anspo_hc_rare,
                              x= 'disease_group_for_study',
                              measures = c("Observed", "Shannon",
                                           "InvSimpson")) +
    geom_boxplot(aes(fill=disease_group_for_study)) + theme(legend.position="none")
}




#' Title
#'
#' @param ps_rare
#' @param ps_ord
#'
#' @return
#' @export
#'
#' @examples
beta_plot <- function(ps_rare, ps_ord) {
  p <- plot_ordination(ps_rare, ps_ord,
                        color =
                         "disease_group_for_study",
                       justDF = TRUE)

  beta_plot <- ggplot(p, aes_string(x="Axis.1",
                                    y="Axis.2",
                                    color="disease_group_for_study"#,
                                    #subject_id="subject_id"
                                    )) +
    geom_point()

  return(beta_plot)
}


#' Title
#'
#' @param ps_new
#' @param secondalpha
#'
#' @return
#' @export
#'
#' @examples
diff_analysis_lefse <- function(ps_new, secondalpha=0.05) {
  sample_data(ps_new)$disease_group_for_study <- factor(sample_data(ps_new)$disease_group_for_study)

  tx_tb <- data.frame(tax_table(ps_new))
  tx_tb <- data.frame(lapply(tx_tb, tidyr::replace_na, replace="NV"))

  rownames(tx_tb) <- rownames(data.frame(tax_table(ps_new)))
  tax_table(ps_new) <- tax_table(as.matrix(tx_tb))

  diffres <- MicrobiotaProcess::diff_analysis(obj=ps_new, classgroup="disease_group_for_study",
                           mlfun="lda",
                           filtermod="pvalue",
                           firstcomfun = "kruskal.test",
                           firstalpha=0.1,
                           strictmod=TRUE,
                           secondcomfun = "wilcox.test",
                           subclmin=3,
                           subclwilc=TRUE,
                           secondalpha=secondalpha,
                           lda=2)
  return(diffres)

}





#' Removes confusing multi level annotations from lefse analysis
#'
#' @param diffres - differential results analysis from diff_analysis_lefse
#'
#' @return
#' @export
#'
#' @examples
rewrite_annotation <- function(diffres){


  gsub_diffres <- function(diffres_taxda_column){

  remove_annotations <- c("s__", "g__", "f__", "c__", "NV_")


  lapply(remove_annotations, function(x){})

  }

  #lapply(diffres$)


  if(!is.null(diffres@taxda$Species)){
  diffres@taxda$Species <- gsub("s__", "", diffres@taxda$Species, fixed=TRUE)
  diffres@taxda$Species <- gsub("NV__", "", diffres@taxda$Species, fixed=TRUE)
  diffres@taxda$Species <- gsub("g__", "", diffres@taxda$Species, fixed=TRUE)
  }


  diffres@taxda$Genus <- gsub("g__", "", diffres@taxda$Genus, fixed=TRUE)
  diffres@taxda$Family <- gsub("f__", "", diffres@taxda$Family, fixed=TRUE)
  diffres@taxda$Class <- gsub("c__", "", diffres@taxda$Class, fixed=TRUE)

  diffres@taxda[diffres@secondvars$AnSpo$f,]

  diffres@mlres$f <- gsub("s__NV_", "", diffres@mlres$f)
  #diffres@mlres$f <- gsub("g__", "", diffres@mlres$f)
  #diffres@mlres$f <- gsub("s__", "", diffres@mlres$f)
  #diffres@mlres$f <- gsub("c__", "", diffres@mlres$f)

  diffres@mlres <- diffres@mlres %>% distinct(f, .keep_all=TRUE)


  return(diffres)
}

#' Filters lefse analysis to only contain one level of analysis
#'
#' @param diffres
#' @param level
#'
#' @return
#' @export
#'
#' @examples
keep_levels <- function(diffres, level = "Species") {

  vals <- c(Species = "^s__", Class = "^c__", Genus= "^g__")

  diffres@mlres <- diffres@mlres %>%
    dplyr::filter(grepl(vals[level], f))

  diffres@mlres$f <- gsub("s__NV_", "", diffres@mlres$f)
  #diffres@mlres$f <- gsub("g__", "", diffres@mlres$f)
  #diffres@mlres$f <- gsub("s__", "", diffres@mlres$f)
  #diffres@mlres$f <- gsub("c__", "", diffres@mlres$f)

  return(diffres)
}

make_manual_scale <- function(f){

  name_vec <- f

  name_vec <- gsub("g__", "", f)
  name_vec <- gsub("s__", "", name_vec)
  name_vec <- gsub("c__", "", name_vec)

  name_vec <- gsub("NV_[gscfo]__", "", name_vec, perl = TRUE)

  name_vec <- gsub("_", " ", name_vec)

  names(name_vec) <- f
  return(name_vec)

}


#' Returns a plot of lefse candidates - abundances separated by group
#'
#' @param diffres
#'
#' @return
#' @export
#'
#' @examples
lefse_plot <- function(diffres) {

  labs2 <- make_manual_scale(diffres@mlres$f)

  plotes_ab <- ggdiffbox(obj=diffres,
                         box_notch=FALSE,
                         addLDA=FALSE,
                         colorlist=c("#F8766D",
                                     "#00BFC4"),
                         l_xlabtext="relative abundance") +
    theme(legend.position="none") +
    scale_x_discrete(labels = labs2)
  return(plotes_ab)
}

#' Convenience function for extracting p-values from lefse analysis
#'
#' @param diffres
#'
#' @return
#' @export
#'
#' @examples
return_p_values <- function(diffres){

  diffres@kwres[diffres@kwres$f %in% diffres@mlres$f,]
}

#' Runs unifrac distance beta diversity analysis
#'
#' @param ps_rare
#' @param formula - a formula for the modeling. should be of the form
#' unifrac_dist ~ covariates. For example: #' unifrac_dist ~ sample_data(ps_rare)$disease_group_for_study
#' @param weighted - TRUE or FALSE. Whether to use weighted unifrac distance or not.
#' Default is unweighted.
#' @param seed - seed value for set.seed
#'
#'
#' @return
#' @export
#'
#' @examples
beta_diversity <- function(ps_rare, formula=NULL, weighted=FALSE, seed=100){
  set.seed(seed)

  unifrac_dist = phyloseq::distance(ps_rare, method="unifrac", weighted=weighted)

  if(is.null(formula)){
    formula <- unifrac_dist ~ sample_data(ps_rare)$disease_group_for_study
  }

  adonis(formula, permutations = 10000)

}

#' Returns alpha diversity measures and pvalues
#'
#' @param ps_rare
#'
#' @return
#' @export
#'
#' @examples
alpha_diversity <- function(ps_rare){
  alpha_measures_hc <- estimate_richness(ps_rare,
                                         measures = c("Observed",
                                                      "Shannon",
                                                      "InvSimpson"))

  test_frame <- data.frame(observed=alpha_measures_hc$Observed,
                           shannon = alpha_measures_hc$Shannon,
                           inv_simpson = alpha_measures_hc$InvSimpson,
                           disease_group_for_study = factor(sample_data(ps_rare)$disease_group_for_study))

  out_list <- list()

  out_list[[1]] <- t.test(observed ~ disease_group_for_study,
         data = test_frame) %>% broom::tidy() %>%
    mutate(method="observed")
  out_list[[2]] <- t.test(shannon ~ disease_group_for_study,
         data = test_frame) %>% broom::tidy() %>%
    mutate(method="shannon")
  out_list[[3]] <- t.test(inv_simpson ~ disease_group_for_study, data = test_frame) %>%
    broom::tidy() %>% mutate(method="inv_simpson")

  return(bind_rows(out_list))

}

#' Title
#'
#' @param ps_rare
#' @param subject_id
#'
#' @return
#' @export
#'
#' @examples
calc_alpha_diversity <- function(ps_rare, subject_id=subject_id){
  alpha_measures_hc <- estimate_richness(ps_rare,
                                         measures = c("Observed",
                                                      "Shannon",
                                                      "InvSimpson"))

  test_frame <- data.frame(subject_id = sample_data(ps_rare) %>% dplyr::pull({{subject_id}}),
                           observed=alpha_measures_hc$Observed,
                           shannon = alpha_measures_hc$Shannon,
                           inv_simpson = alpha_measures_hc$InvSimpson,
                           disease_group_for_study = factor(sample_data(ps_rare)$disease_group_for_study))

  return(test_frame)
}

