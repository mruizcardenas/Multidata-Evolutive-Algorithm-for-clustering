#Directorio de trabajo
setwd("C:/Users/Manuel Ruiz Cárdenas/Desktop/TFM/nuevos datos")
rm(list=ls())
######################################################################################################################################

#datasets usados
dataset1<-"caso1_5.csv"
dataset2<-"caso2_5.csv"
matriz_1<-read.csv(dataset1,sep=",")
matriz_2<-read.csv(dataset2,sep=",")
datalist<-list((matriz_1),(matriz_2))
######################################################################################################################################

#Paquetes y funciones externos

require("pheatmap")
heat<-function(matrix,name=" ",diag=FALSE){
  if(diag==FALSE){
    diag(matrix)<-0
  }
  pheatmap(matrix,color=colorRampPalette(c("black","red","yellow"))(1000),cluster_cols = FALSE, cluster_rows = FALSE,main =name)
  
}
medioids<-function(DATA,Z,I){
  return(apply(DATA[Z==I,],2,mean,na.rm=T))
}
cohesion<-function(Z,ncluster,DATA){
  SSE<-0
  for(i in 1:ncluster){
    medioid<-medioids(DATA,Z,i)
    clusterelements<-DATA[Z==i,]
    for(j in 1:nrow(clusterelements)){
      SSE<-SSE+(sum((medioid-clusterelements[j,])^2))
    }
    SSE<-SSE
    
  }
  return(SSE)
}
cohesion2<-function(Z,ncluster,DIST2){
  SSE<-0
  for(i in 1:ncluster){
    clusterelements<-1:nrow(DIST2)
    clusterelements<-clusterelements[Z==i]
    for(j in clusterelements){
      SSE<-SSE+(sum(DIST2[clusterelements,j])/length(clusterelements))
    }
    SSE<-SSE
    
  }
  return(SSE)
}
cohesion_multi_mean<-function(Z,ncluster,DATALIST){
  n<-length(DATALIST)
  cohesion_vector<-c()
  for(i in 1:n){
    DATA<-as.data.frame(DATALIST[i])
    cohesion_vector<-c(cohesion_vector,cohesion(Z,ncluster,DATA))
  }
  return(mean(cohesion_vector))
}
separacion<-function(Z,ncluster,DATA){
  SSB<-0
  meanDATA<-apply(DATA,2,mean,na.rm=T)
  for(i in 1:ncluster){
    SSB<-SSB+(sum(Z==i))*sum((medioids(DATA,Z,i)-meanDATA)^2)
  }
  return(SSB)
}
separacion_multi_mean<-function(Z,ncluster,DATALIST){
  n<-length(DATALIST)
  separacion_vector<-c()
  for(i in 1:n){
    DATA<-as.data.frame(DATALIST[i])
    separacion_vector<-c(separacion_vector,separacion(Z,ncluster,DATA))
  }
  return(mean(separacion_vector))
}
######################################################################################################################################

#Funciones de Evo_joint

parejasposibles<-function(size){
  if(size==2){
    return(matrix(c(1,2),1,2))
  }else{
    parejas<-matrix(c(0,0),1,2)
    for(i in 1:size){
      for(j in 1:size){
        if(j<i){
          parejas<-rbind(parejas,c(i,j))
        }
      }
    }
    return(parejas[-1,])
  }
}
hibridacion<-function(a,b,k,puntuacion_list=c(0.5,0.5)){
  l<-length(a)
  new<-rep(0,l)
  index<-1:l
  alpha<-puntuacion_list[1]/(puntuacion_list[1]+puntuacion_list[2])
  mantener_a_num<-floor(l*0.40*2*alpha)
  mantener_b_num<-floor(l*0.40*2*(1-alpha))
  mantener_a<-sample(index,mantener_a_num)
  index_<-setdiff(index,mantener_a)
  mantener_b<-sample(index_,mantener_b_num)
  index_<-setdiff(index_,mantener_b)
  new[mantener_a]<-a[mantener_a]
  new[mantener_b]<-b[mantener_b]
  new[index_]<-sample(k,length(index_),replace = TRUE)
  return(new)
}
sample_e<-function(x,s,prob){
  x1<-x[prob>0]
  p1<-prob[prob>0]
  return(sample(x1,s,prob = p1,replace=FALSE))
}
hibridacion_2<-function(a,b,k,puntuacion_list=c(0.5,0.5)){
  l<-length(a)
  new<-rep(0,l)
  MP<-matrix((1:(l+1)),1,l+1)
  for(i in 1:k){
    v<-rep(0.0001,l)
    index_i_a<-which(a==i)
    index_i_b<-which(b==i)
    size_i<-floor((length(index_i_b)+length(index_i_a))/(2))
    E<-table(c(rep(c(index_i_a),floor(1+(puntuacion_list[1]/puntuacion_list[2]))),rep(c(index_i_b),floor(1+(puntuacion_list[1]/puntuacion_list[2])) )))
    #E<-table(c(index_i_a,index_i_b))
    e<-as.numeric(names(E))
    P<-as.numeric(E)
    v[e]<-P
    v<-c(v,size_i)
    MP<-rbind(MP,v)
  }
  index<-c()
  for(j in 1:k){
    pc<-setdiff(1:l,index)
    fila<-as.numeric(MP[j+1,])
    h<-sample(pc,fila[l+1],prob = fila[pc])
    new[h]<-j
    index<-c(index,h)
  }
  cerosize<-sum(new==0)
  new[new==0]<-sample(k,cerosize,replace = TRUE)
  mutation_size<-floor(l*0.1)
  mutation_places<-sample(l,mutation_size,replace = FALSE)
  new[mutation_places]<-sample(k,mutation_size,replace = TRUE)
  return(new)
}
gen<-function(parejas,cluster_list,pl,k){
  l<-ncol(cluster_list)
  L<-matrix(0,1,l)
  for(i in 1:nrow(parejas)){
    a<-parejas[i,1]
    b<-parejas[i,2]
    p<-c(pl[a],pl[b])/(pl[a]+pl[b])
    new<-hibridacion_2(cluster_list[a,],cluster_list[b,],k,p)
    L<-rbind(L,new)  
  }
  L<-rbind(L,cluster_list)
  return(L[-1,])
}
criba<-function(cmatrix,puntuacion_list,k,num){
  n<-nrow(cmatrix)
  qc<-c()
  for(i in 1:n){
    qc<-c(qc,cohesion_multi_mean(cmatrix[i,],k,datalist))
  }
  oqc<-order(qc,decreasing = FALSE)
  punt<-(qc[oqc[1:num]])
  listareturn<-list(clustering=cmatrix[oqc[1:num],],puntuacion=punt)
  return(listareturn)
  #return(cmatrix[oqc[1:num],])
}
rbind_sinrep<-function(old,new){
  l<-ncol(old)
  if(is.null(nrow(old))){
    old<-matrix(old,1,l)
  }
  a<-nrow(old)
  if(is.null(nrow(new))){
    new<-matrix(new,1,l)
  }
  
  s<-nrow(new)
  resultado<-old
  for(i in 1:s){
    for(j in 1:a){
      if(sum(new[i,]==old[j,])==l){
        ins<-0
      }else{
        ins<-1
      }
    }
    if(ins==1){
      rbind(old,new[i,])
    }
  }
  return(resultado)
}
puntuaciones<-function(g_old=NULL,p_old=NULL,g_new,k,datalist){
  punt<-c()
  l_new<-nrow(g_new)
  #print(paste("g_old",as.numeric(dim(g_old)),"p_old",length(p_old),"g_new",as.numeric(dim(g_new))))
  l_old<-nrow(g_old)
  for(i in 1:l_new){
    f<-FALSE
    j<-1
    while(f==FALSE&&is.null(l_old)==FALSE&&j<=l_old){
      if(sum(g_old[j,]==g_new[i,])==ncol(g_old)){
        f<-TRUE
      }else{
        j<-j+1
      }
    }
    if(f==FALSE){
      punt<-c(punt,cohesion_multi_mean(g_new[i,],k,datalist))
    }else{
      punt<-c(punt,p_old[j])
    }
  }
  return(punt)
}
EVO_joint<-function(k,datalist,cluster_ini_list,iter){
  
  cluster_list<-cluster_ini_list
  evolutionofmethod<-c()
  i<-1
  num<-nrow(cluster_list)
  pl<-puntuaciones(g_new=cluster_list,k=k,datalist=datalist)
  time<-0
  while(i<=iter&&nrow(cluster_list)>1){
    t <- proc.time() 
    p<-parejasposibles(nrow(cluster_list))
    newgen<-gen(p,cluster_list,pl,k)
    if(i==iter){
      num<-1
    }
    pl<-puntuaciones(cluster_list,pl,newgen,k,datalist)
    #C<-criba(newandold,puntuacion_list,k,num)
    #cluster_list<-C$clustering
    #puntuacion_list<-C$puntuacion
    position<-order(pl)
    cluster_list<-newgen[position[1:num],]
    pl<-pl[position[1:num]]
    #evolutionofmethod<-c(evolutionofmethod,puntuacion_list[1])
    
    evolutionofmethod<-c(evolutionofmethod,pl[1])
    if(i>1&&evolutionofmethod[length(evolutionofmethod)]<evolutionofmethod[length(evolutionofmethod)-1]){
      print(paste("mejora del",(100*(1-(evolutionofmethod[length(evolutionofmethod)])/(evolutionofmethod[length(evolutionofmethod)-1]))),"%" ))
    }
    if(i>1&&evolutionofmethod[length(evolutionofmethod)]>evolutionofmethod[length(evolutionofmethod)-1]){
      print("empeora")
    }
    time<-time+proc.time()-t
    #print(paste("iteración",i,"quedan",(time*iter/i)%/%60,"minutos"))
    tiemporestante<-((time[3]*iter/i)-time[3])%/%60
    print(paste("iteración",i,",",tiemporestante,"minutos restantes"))
    
    i<-i+1
    
  }
  lista<-list(cluster=cluster_list,evo=evolutionofmethod)
  
  
  return(lista)
}
#####################################################################################################################################
c1<-c(rep(1,60),rep(2,20))
c2<-c(rep(1,20),rep(2,60))
c3<-c(rep(1,30),sample(2,20,replace = TRUE),rep(2,30))
#c1<-sample(2,80,replace = TRUE)
#c2<-sample(2,80,replace = TRUE)
cluster_ini<-rbind(c1,c2,c3)
heat(cluster_ini,diag = TRUE)
parejasposibles(2)
hibridacion(cluster_ini[1,],cluster_ini[2,],2)
heat(hibridacion(cluster_ini[1,],cluster_ini[2,],2),diag = TRUE)
a<-EVO_joint(2,datalist,cluster_ini,250)
plot(a$evo)
heat(t(a$cluster),diag = TRUE)
#write.table(cluster_ini,"cluster_ini_5.txt")
#write.table(a$cluster,"quinta  prueba.txt")
#write.table(a$evo,"quinta  prueba_2.txt")

