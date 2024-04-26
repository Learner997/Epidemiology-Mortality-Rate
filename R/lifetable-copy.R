#'age-adjusted mortality life table
#'
#'@importFrom tidyr ggplot2
#'@param x a string of time
#'@param n a number of population
#'@param wx a number of withdraws for each age group
#'@param wt date of withdraw
#'@param dt date of death
#'@param dx a number of deaths for each age group
#'@param grp groups
#'@retrun a life table and a plot
#'@example example.R
#'@export

ensure_package <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

lifetable<-function(x=NULL,n=NULL,wx=NULL,wt=NULL,grp=NULL,dx=NULL,dt=NULL,data=NULL){

  ensure_package("tidyverse")

  env_data <- if (is.null(data)) parent.frame() else data
  if (!is.null(data)) {
    x <- eval(substitute(x), envir = env_data)
    n <- eval(substitute(n), envir = env_data)
    wx <- eval(substitute(wx), envir = env_data)
    dx <- eval(substitute(dx), envir = env_data)
    wt <- eval(substitute(wt), envir = env_data)
    dt <- eval(substitute(dt), envir = env_data)
    grp <- eval(substitute(grp), envir = env_data)
  }
  #group withdraw date and death date
  if(is.null(x)&is.null(n)){
    wt<-as.Date(wt,format="%m/%d/%y")
    dt<-as.Date(dt,format="%m/%d/%y")
    if (!is.null(grp)) {
      grp<-sort(c(as.Date(grp),Inf))
    }
    if (!is.null(wt)) {
      wg <- as.character(as.factor(cut(as.Date(wt),breaks=grp,include.lowest=TRUE)))
    }
    if (!is.null(dt)) {
      dg <- as.character(as.factor(cut(as.Date(dt),breaks=grp,include.lowest=TRUE)))
 }
    #new table
  lt<-data.frame(group=grp,
                 population=double(length=length(grp)),
                 withdraw=double(length=length(grp)),
                 death=double(length=length(grp)))
  lt$population[1]<-nrow(data)
  for(i in 1:nrow(lt)){
    lt$withdraw[i]<-sum(wg==lt$group[i] & !is.na(wg))
    lt$death[i]<-sum(!is.na(dg)&dg==lt$group[i])
  }
  for (i in 2:nrow(lt)){
    lt$population[i]<-lt$population[i-1]-lt$withdraw[i-1]-lt$death[i-1]
  }
  x<-lt$group
  n<-lt$population
  wx<-lt$withdraw
  dx<-lt$death
  #plot
  qtable<-data.frame(wt=data2$wt,dt=data2$dt)

  qtable$time<-qtable$wt
  qtable$time[!is.na(qtable$dt)]<-qtable$dt[!is.na(qtable$dt)]
  qtable$time<-sort(qtable$time)
  if(!is.na(qtable$wt[1])){
    qtable$risk[1]<- nrow(qtable)-0.5
  }else{
    qtable$risk[1]<- nrow(qtable)
  }
  if(!is.na(qtable$dt[1])){
    qtable$Q[1]<- 1-1/qtable&risk[1]
  }else{
    qtable$Q[1]<- 1
  }
  for(i in 2:nrow(qtable)){
    if(!is.na(qtable$wt[i])){
      qtable$risk[i]<- qtable$risk[i-1]-0.5
    }else{
      qtable$risk[i]<- qtable$risk[i-1]
    }
    if(!is.na(qtable$dt[i])){
      qtable$Q[i]<- qtable$Q[i-1]*(1-1/qtable$risk[i])
    }else{
      qtable$Q[i]<-qtable$Q[i-1]
    }
  }
  plot<-ggplot(qtable,aes(x=time,y=Q))+
    geom_step(aes(group = 1),color='blue')+
    #geom_vline(xintercept = qtable$time, linetype = "dashed", color = "orange") +
    labs(x='time',y='survival probability',
         title='survival curve')+
    theme_minimal()

  }
  #life table
  risk<-n-0.5*wx
  px<-dx/risk
  qx<-1-px
  table<-data.frame(`Time`=x,
                    `alive at beginning`=n,
                    `withdraw`=wx,
                    `death`=dx,
                    `at risk`=risk,
                    `proportion of survival`=qx,
                    `mortality`=px,
                    `cumulative survival`=cumprod(qx)
  )
  if(is.null(wt)&is.null(dt)){
    #plot
    ltable<-table%>%drop_na()
    ltable$Time <- factor(ltable$Time, levels = ltable$Time)
    plot<-ggplot(ltable,aes(x=Time,y=proportion.of.survival))+
      geom_line(aes(group = 1),color='blue')+
      geom_vline(xintercept = ltable$Time[ltable$Time != Inf], linetype = "dashed", color = "orange") +
    labs(x='time',y='survival probability',
         title='survival curve')+
    theme_minimal()
  }
  print(plot)
 return(table)
}














