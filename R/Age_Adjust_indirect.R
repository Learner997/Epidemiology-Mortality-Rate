#' Age-Adjustment by Indirect Method
#'
#' This function performs age-adjustment using the indirect method. It requires standard death rates
#' for comparison and supports custom age breaks. The indirect method is typically used when
#' the study population is small and the rates are unstable.
#'
#' @param dataframe A data frame containing the data to be adjusted.
#' @param ageColumn The name of the column in `dataframe` that contains age groups.
#' @param deathColumn The name of the column in `dataframe` that contains the count of deaths.
#' @param variableColumn The name of the column in `dataframe` that contains the variable by which to stratify.
#' @param SD_rate A numeric vector of standard death rates corresponding to specified age breaks.
#' These rates are typically taken from a larger, stable population.
#' @param ageBreaks A vector of age breaks for categorizing the data; defaults to standard age categories.
#' @return A data list with age-adjusted rates.
#' @export
#' @examples
#' Age_Adjust_indirect(data = myData, ageColumn = "age", deathColumn = "deaths",
#'                   variableColumn = "disease", SD_rate = c(50, 79, 105, 120, 250))


Age_Adjust_indirect <- function(dataframe, ageColumn, deathColumn, variableColumn,
                                SD_rate = c(50, 79, 105, 120, 250),
                                ageBreaks = c(0,15,45,55,65,Inf)) {

  # Separate the data frame according to variableColumn
  datalist <- split(dataframe, dataframe[[variableColumn]])
  #return(datalist)

  # Create a function to assign age group and use this function in data list
  Population_process_f <- function(data, ageBreaks) {
    data$AgeGroup <- cut(data[[ageColumn]], breaks = ageBreaks, right = FALSE)
    return(data)
  }

  datalist_pro <- lapply(datalist, function(data) Population_process_f(data, ageBreaks))
  #return(datalist_pro)

  # Find how many person in each age group and create a new list to store these information
  AgeGroups_count_f <- function(data) {
    return(as.data.frame(table(data$AgeGroup)))
  }

  AgeGroupCount <- lapply(datalist_pro, AgeGroups_count_f)
  #return(AgeGroupCount)

  # How many observed deaths in each age group
  Observed_deaths_f <- function(data, deathColumn) {
    aggregate(Deathstatus ~ AgeGroup, data = data, function(x) sum(x == 1))
  }

  Observed_deaths <- lapply(datalist_pro, Observed_deaths_f)
  #return(Observed_deaths)

  # Merge these two data frame, one with how many person in age group and one with how many observed deaths in each age group
  name_change <- function(data) {
    colnames(data)[colnames(data) == "Var1"] <- "AgeGroup"
    return(data)
  }
  AgeGroupCount <- lapply(AgeGroupCount, name_change)

  merge_f <- function(data1, data2) {
    merge(data1, data2, by = "AgeGroup")
  }

  Population_list <- Map(merge_f, Observed_deaths, AgeGroupCount)
  #return(Population_list)

  # Add Standard rates to data frame
  add_SD_rate <- function(df) {
    df$SD_rate <- SD_rate/100000
    return(df)
  }

  Population_list <- lapply(Population_list, add_SD_rate)
  #return(Population_list)

  # Except Deaths Calculate
  Exp_Deas_f <- function(df) {
    df$Expected_Deaths <- df$Freq * df$SD_rate
    return(df)
  }

  Population_list <- lapply(Population_list, Exp_Deas_f)
  #return(Population_list)

  # Calculate how many Observed cases in the population
  Obs_cas_f <- function(df) {
    df$Observed_cases <- sum(df$Deathstatus)
    return(df)
  }

  Population_list <- lapply(Population_list, Obs_cas_f)

  # Calculate how many Expected cases in the population
  Exp_cas_f <- function(df) {
    df$Expected_cases <- sum(df$Expected_Deaths)
    return(df)
  }

  Population_list <- lapply(Population_list, Exp_cas_f)

  # Calculate SMR for the population
  SMR_f <- function(df) {
    df$SMR_total <- df$Observed_cases / df$Expected_cases
    return(df)
  }

  Population_list <- lapply(Population_list, SMR_f)

  # Crude rate calculate
  Crude_f <- function(df) {
    df$Crude_rate <- sum(df$Deathstatus) / sum(df$Freq)
    return(df)
  }

  Population_list <- lapply(Population_list, Crude_f)

  # Calculate SMR and crude rate for each age group
  SMR_Crude_age_f <- function(df) {
    df$SMR_age <- df$Deathstatus / df$Expected_Deaths
    df$Crude_rate_age <- df$Deathstatus / df$Freq
    return(df)
  }

  Population_list <- lapply(Population_list, SMR_Crude_age_f)

  return(Population_list)
}

