#' Age-Adjustment by Direct Method
#'
#' This function performs age-adjustment using the direct method. It allows for
#' adjustment with internal or external standard populations and supports custom age breaks.
#'
#' @param dataframe A data frame containing the data to be adjusted.
#' @param ageColumn The name of the column in `dataframe` that contains age groups.
#' @param deathColumn The name of the column in `dataframe` that contains the count of deaths.
#' @param variableColumn The name of the column in `dataframe` that contains the variable by which to stratify.
#' @param standardPopulation An optional vector of standard population sizes corresponding to age breaks; defaults to NULL.
#' @param externalSP Boolean, if TRUE, uses an external standard population provided via `standardPopulation`.
#' @param ageBreaks A vector of age breaks for categorizing the data; defaults to standard age categories.
#' @return A data list with age-adjusted rates.
#' @export
#' @examples
#' Age_Adjust_direct(data = myData, ageColumn = "age", deathColumn = "deaths",
#'                   variableColumn = "disease", standardPopulation = mySP, externalSP = TRUE,
#'                   ageBreaks = c(0, 5, 15, 25, 35, 45, 55, 65, 75, 85, Inf))


# Age Function
Age_Adjust_direct <- function(dataframe, ageColumn, deathColumn, variableColumn,
                              standardPopulation = NULL, externalSP = FALSE,
                              ageBreaks = c(0, 5, 15, 25, 35, 45, 55, 65, 75, 85, Inf)) {

  # Separate the data frame according to variableColumn
  datalist <- split(dataframe, dataframe[[variableColumn]])
  #return(datalist)

  # Assign age group
  Population_process_f <- function(data, ageBreaks) {
    data$AgeGroup <- cut(data[[ageColumn]], breaks = ageBreaks, right = FALSE)
    return(data)
  }

  datalist_pro <- lapply(datalist, function(data) Population_process_f(data, ageBreaks))

  #return(datalist_pro)

  # Calculate Freq for each age group in different population
  AgeGroups_count_f <- function(data) {
    return(as.data.frame(table(data$AgeGroup)))
  }

  AgeGroupCount <- lapply(datalist_pro, AgeGroups_count_f)

  AgeGroups_name_f <- function(data) {
    colnames(data)[colnames(data) == "Var1"] <- "AgeGroup"
    return(data)
  }

  AgeGroupCount <- lapply(AgeGroupCount, AgeGroups_name_f)

  #return(AgeGroupCount)

  # How many deaths for each age group in different population
  DeathCounts_f <- function(data, deathColumn) {
    aggregate(Deathstatus ~ AgeGroup, data = data, function(x) sum(x == 1))
  }

  DeathCount <- lapply(datalist_pro, DeathCounts_f)

  #return(DeathCount)

  # Death Rate for each age group in different population
  Age_merge_f <- function(data1, data2) {
    merge(data1, data2, by = "AgeGroup")
  }

  AgeFre_Deaths_merge <- Map(Age_merge_f, AgeGroupCount, DeathCount)

  #return(AgeFre_Deaths_merge)

  # Calculate Specific Death Rates for each populations
  DeathR_f <- function(data) {
    data$SpeRate <- data$Deathstatus / data$Freq
    return(data)
  }

  AgeFre_Deaths_merge_rate <- lapply(AgeFre_Deaths_merge, DeathR_f)

  #return(AgeFre_Deaths_merge_rate)

  if (externalSP && is.null(standardPopulation)) {
    stop("External standard population must be provided if externalSP is TRUE")
  }

  if (!externalSP && !is.null(standardPopulation)) {
    stop("externalSP is FALSE but a standardPopulation is provided.
         The provided standardPopulation will be ignored.
         Stand Population will be sum population.")
  }

  if (!externalSP) {
    # Sum Freq according to AgeGroup
    Freq_sum_f <- function(df_list) {
      combined_df <- do.call(rbind, df_list)
      SP <- aggregate(Freq ~ AgeGroup, data = combined_df, FUN=sum)
      return(SP)
    }

    SP_df <- Freq_sum_f(AgeFre_Deaths_merge_rate)
    #return(SP_df)

    # Calculate expect death for each populations
    Expect_rate_f <- function(df_list, df_sp) {
      result_list <- lapply(df_list, function(df_spe) {
        merged_df <- merge(df_spe, df_sp, by = "AgeGroup", suffixes = c("_spe", "_sp"))
        merged_df$Exp_deaths <- merged_df$SpeRate * merged_df$Freq_sp
        return(merged_df)
      })

      return(result_list)
    }

    Exp_rate_l <- Expect_rate_f(AgeFre_Deaths_merge_rate, SP_df)
    #return(Exp_rate_l)

    # Sum stand population and exp deaths for each population

    SP_total <- sum(SP_df$Freq)

    Exp_deaths_tot_f <- function(df) {
      Exp_deaths_tot <- sum(df$Exp_deaths)
      df$Exp_deaths_tot <- Exp_deaths_tot
      return(df)
    }

    Exp_rate_l <- lapply(Exp_rate_l, Exp_deaths_tot_f)

    Adj_rate_f <- function(df) {
      df$Adj_rate <- (df$Exp_deaths_tot / SP_total) * 1000
      return(df)
    }

    Exp_rate_l <- lapply(Exp_rate_l, Adj_rate_f)
    #return(Exp_rate_l)

    # Crude rate calculate
    Crude_f <- function(df) {
      df$Crude_rate <- sum(df$Deathstatus) / sum(df$Freq_spe)
      return(df)
    }
    Exp_rate_l <- lapply(Exp_rate_l, Crude_f)

    # Calculate Adj rate and Crude rate for age group
    Adj_age_f <- function(df) {
      df$Adj_rate_age <- df$Exp_deaths / df$Freq_sp
      return(df)
    }
    Exp_rate_l <- lapply(Exp_rate_l, Adj_age_f)

    Exp_rate_l$SP_total <- SP_total
    return(Exp_rate_l)
  }

  if (externalSP) {
    # Cut SP by age group
    result_df <- Population_process_f(standardPopulation, ageBreaks)
    #return(result_df)

    # Calculate Freq for each age group for SP_df
    SP_df <- AgeGroups_count_f(result_df)

    # Change name
    colnames(SP_df)[colnames(SP_df) == "Var1"] <- "AgeGroup"
    #return(SP_df)

    # Sum SP_df Freq
    SP_total <- sum(SP_df$Freq)
    #return(AgeFre_Deaths_merge_rate)

    # Calculate Exp Deaths for each population
    Expect_rate_f <- function(df_list, df_sp) {
      result_list <- lapply(df_list, function(df_spe) {
        merged_df <- merge(df_spe, df_sp, by = "AgeGroup", suffixes = c("_spe", "_sp"))
        merged_df$Exp_deaths <- merged_df$SpeRate * merged_df$Freq_sp
        return(merged_df)
      })

      return(result_list)
    }

    Exp_rate_l <- Expect_rate_f(AgeFre_Deaths_merge_rate, SP_df)
    #return(Exp_rate_l)

    Exp_deaths_tot_f <- function(df) {
      Exp_deaths_tot <- sum(df$Exp_deaths)
      df$Exp_deaths_tot <- Exp_deaths_tot
      return(df)
    }

    Exp_rate_l <- lapply(Exp_rate_l, Exp_deaths_tot_f)

    Adj_rate_f <- function(df) {
      df$Adj_rate <- (df$Exp_deaths_tot / SP_total) * 1000
      return(df)
    }

    Exp_rate_l <- lapply(Exp_rate_l, Adj_rate_f)
    #return(Exp_rate_l)

    # Crude rate calculate
    Crude_f <- function(df) {
      df$Crude_rate <- sum(df$Deathstatus) / sum(df$Freq_spe)
      return(df)
    }
    Exp_rate_l <- lapply(Exp_rate_l, Crude_f)

    # Calculate Adj rate and Crude rate for age group
    Adj_age_f <- function(df) {
      df$Adj_rate_age <- df$Exp_deaths / df$Freq_sp
      return(df)
    }
    Exp_rate_l <- lapply(Exp_rate_l, Adj_age_f)

    Exp_rate_l$SP_total <- SP_total

    return(Exp_rate_l)
  }
}

