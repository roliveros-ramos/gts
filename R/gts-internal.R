
is_increasing = function(x) {
  return(all(diff(x)>0))
}

is_decreasing = function(x) {
  if(length(x)==1) return(FALSE)
  return(all(diff(x)<0))
}

is_monotonic = function(x) {
  return(is_increasing(x) | is_decreasing(x))
}

DateStamp = function(...) cat(..., "\t\t | ", date(), "\n\n")

logit = function(x) log(x/(1-x))
logistic = function(x) 1/(1 + exp(-x))

.strip = function(object, target) {
  tmp = 0*target + object
  return(as.numeric(tmp$x))
}

