# load grplassores

num_selected = c(0, 0, 0)
for (i in 1:length(grplassores)) {
  num_correctly_identified = grplassores[[i]]$number_correctly_selected
  print(num_correctly_identified)
  num_selected[num_correctly_identified+1] = num_selected[num_correctly_identified+1] + 1
}