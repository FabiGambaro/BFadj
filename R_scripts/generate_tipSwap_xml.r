generate_xml_file <- function(file_name, metadata, input_dir, output_dir, template_file) {
  test <- read_tsv(paste0(input_dir, file_name), show_col_types = F)
  name.test <- gsub(".txt","",gsub("traits_", "", file_name))
  tab <- as.data.frame(merge(metadata, test, by = "traits", all = FALSE)) #merge with metadata
  
  xml <- scan(template_file, what = "", sep = "\n", quiet = TRUE, blank.lines.skip = FALSE)
  xml <- gsub("template", paste0(name.test, "_tipSwap"), xml)
  
  if (scheme_500) {
    xml <- gsub("template", paste0(name.test, "_tipSwap"), xml)
    xml <- gsub("chainLength=\"20000000\"", "chainLength=\"40000000\"",xml)
    xml <- gsub("logEvery=\"20000\"","logEvery=\"40000\"",xml)
  }
  
  output_file <- paste0(output_dir, name.test, "_tipSwap.xml")
  sink(file = output_file)
  
  for (n in 1:length(xml)) {
    cat(xml[n], "\n")
    
    if (grepl("<taxa id=\"taxa\">", xml[n])) {
      for (j in 1:nrow(tab)) {
        cat("\t\t<taxon id=\"", tab[j, "traits"], "\">\n", sep = "")
        cat("\t\t\t<date value=\"", tab[j,"samp.dates"], "\" direction=\"forwards\" units=\"years\"/>\n", sep = "")
        cat("\t\t\t<attr name=\"regions\">\n", sep = "")
        cat("\t\t\t\t", tab[j, "regions"], "\n", sep = "")
        cat("\t\t\t</attr>\n", sep = "")
        cat("\t\t</taxon>\n", sep = "")
      }
    }
    
    if (grepl("<alignment id=\"alignment\" dataType=\"nucleotide\">", xml[n])) {
      for (j in 1:nrow(tab)) {
        cat("\t\t<sequence>\n", sep = "")
        cat("\t\t\t<taxon idref=\"", tab[j, "traits"], "\"/>\n", sep = "")
        cat("\t\t\t\tNNNN\n", sep = "")
        cat("\t\t</sequence>\n", sep = "")
      }
    }
    
    if (grepl("<!-- insert number of states-->", xml[n])) {
      regions <- unique(tab$regions)
      for (k in 1:length(regions)) {
        cat("\t\t<state code=\"", regions[k], "\"/>\n", sep = "")
      }
    }
    
    if (grepl("<!-- Defining empirical tree distribution-->", xml[n])) {
      cat("\t<empiricalTreeDistributionModel id=\"treeModel\" fileName=\"", paste0(name.test, "_select100.trees"), "\">\n", sep = "")
      cat("\t\t<taxa idref=\"taxa\"/>\n", sep = "")
      cat("\t</empiricalTreeDistributionModel>\n", sep = "")
    }
    
    if (grepl("<!-- insert regions frequencies -->", xml[n])) {
      cat("\t\t\t<parameter id=\"regions.frequencies\" dimension=\"", length(regions), "\"/>\n", sep = "")
    }
    
    if (grepl("<!-- insert regions rates -->", xml[n])) {
      nb.transitions <- length(regions) * length(regions) - length(regions)
      cat("\t\t<parameter id=\"regions.rates\" dimension=\"", nb.transitions, "\" value=\"1.0\" lower=\"0.0\"/>\n", sep = "")
    }
    
    if (grepl("<!-- insert regions indicators -->", xml[n])) {
      cat("\t\t<parameter id=\"regions.indicators\" dimension=\"", nb.transitions, "\" value=\"1.0\"/>\n", sep = "")
    }
    
    if (grepl("<!-- insert regions root frequencies -->", xml[n])) {
      cat("\t\t\t\t<parameter id=\"regions.root.frequencies\" dimension=\"", length(regions), "\"/>\n", sep = "")
    }
    
    if (grepl("<!-- Add tipStateSwapOperator -->", xml[n])) {
      cat("\t<tipStateSwapOperator weight=\"", tipSwapOperator, "\" uniformRandomization=\"true\">\n", sep = "")
      cat("\t\t<ancestralTreeLikelihood idref=\"regions.treeLikelihood\"/>\n")
      cat("\t</tipStateSwapOperator>\n")
    }
    
    if (grepl("<!-- Add non-zero rates prior -->", xml[n])) {
      cat("\t\t<poissonPrior mean=\"", length(regions)-1, "\" offset=\"0.0\">\n", sep = "")
      cat("\t\t\t<statistic idref=\"regions.nonZeroRates\"/>\n")
      cat("\t\t</poissonPrior>")
    }
  }
  
  sink(NULL)
  message("Generated XML file:", output_file)
}