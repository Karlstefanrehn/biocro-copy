## Very simple function for soil temperature

## The argument is a data frame with hour of the day and air temperature

stemp <- function(hr.temp, a=8, t0=0, lag=48){

    ## Sanity checks
    ## The object hr.temp should have hour in the first column and temperature in the second column 
    if(ncol(hr.temp) != 2) stop("number of columns should be equal to 2")
    
    ans <- numeric(nrow(hr.temp))

    for(i in 1:nrow(hr.temp)){

        if(i < lag){
            ans[i] <- hr.temp[i,1]
        }else{
            ## Collect the data for the last x hours, given by lag
            tmp <- hr.temp[(i-lag):i,2]
            mtmp <- mean(tmp, na.rm=TRUE)
            amp <- (max(tmp, na.rm=TRUE) - min(tmp, na.rm=TRUE))/a
            ans[i] <- mtmp + amp * -sin((hr.temp[i,1] - t0)/3.66)
        }

    }

    ans

}
