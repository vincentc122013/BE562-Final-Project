scores <- c()
for (i in c(seq(1,10,by=1))){
  print(paste("dab",i))
  num = toString(i)
  extnd = '.nwk'
  
  str = 'output_newick_'
  fullname = paste(str,num,extnd,sep="")

  str2 = 'sub1_train_'
  fullname2 = paste(str2,num,extnd,sep="")
  
  a_tree = read.newick(fullname2)
  b_tree = read.newick(fullname)
  
  score = RF.dist(a_tree,b_tree,normalize=TRUE,rooted=TRUE)
  scores <-append(scores,score)
  output_vec_scores = cbind(scores)
}

