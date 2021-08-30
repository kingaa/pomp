\donttest{
  po <- sir()
  plot(po)
  coef(po)
  
  po <- sir2()
  plot(po)
  plot(simulate(window(po,end=3)))
  coef(po)
  
  po %>% as.data.frame() %>% head()
}
