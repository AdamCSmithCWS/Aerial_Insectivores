transp.func <- function(cl = "", trans = 0.5){
transp <- 255*trans
 cl1 <- rgb(red = col2rgb(cl)[1,1],green = col2rgb(cl)[2,1], blue = col2rgb(cl)[3,1],alpha = transp,maxColorValue = 255)
 return(cl1)
}
