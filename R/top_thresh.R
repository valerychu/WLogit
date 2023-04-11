top_thresh <-
function(vect, thresh){
    sorted_vect <- sort(vect,decreasing=TRUE)
    v=sorted_vect[thresh]
    ifelse(vect>=v,vect,v)
  }
