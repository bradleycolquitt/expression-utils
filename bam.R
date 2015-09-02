require(rbamtools)
pileup_5p = function(range) {
  coords = getCoords(range)    
  out_vector = vector("numeric", coords[3] - coords[2])
   while(!is.null(align)) {
     align = getNextAlign(range)
     pos = position(align)
     out_vector[pos] = out_vector[pos] + 1 
   }
   return(out_vector)
}