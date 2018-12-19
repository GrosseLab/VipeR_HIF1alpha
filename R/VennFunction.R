#-------------------------------------------
#' @title Venndiagramm Whrapper 
#' @author Claus Weinholdt
#' @description Venndiagramm for 2-Set
#' @param inList is list
#' @param name vector with names of sets
#' @param VennPlot is boolean to enable venn plot
#' @param VennOut is boolean to enable venn plot
#' @export
f.input.list <- function(inList,VennPlot=TRUE,VennOut=FALSE){
  
  N <- names(inList)
  
  switch( 
     as.character(length( N )),  
     "2"=f.input2(inList[[1]],inList[[2]],name=N,plotVENN=VennPlot),
     "3"=f.input3(inList[[1]],inList[[2]],inList[[3]],name=N,plotVENN=VennPlot,vennOut=VennOut),
     "4"=f.input4(inList[[1]],inList[[2]],inList[[3]],inList[[4]],name=N,plotVENN=VennPlot,vennOut=VennOut),
     "5"=f.input5(inList[[1]],inList[[2]],inList[[3]],inList[[4]],inList[[5]],name=N,plotVENN=VennPlot,vennOut=VennOut)
  )
  
  
}


#-------------------------------------------
#' @title Venndiagramm 2-Set
#' @author Claus Weinholdt
#' @description Venndiagramm for 2-Set
#' @param A is the 1. set
#' @param B is the 2. set
#' @param name vector with names of sets
#' @param plotVENN is boolean to enable venn plot
#' @export
f.input2 = function (p,q,name=c("A","B"),plotVENN=TRUE){
  require(gplots)
  input  <-list(A=p,B=q)
  names(input)<-name
  
  input <- lapply(input,function(x) 
    if(length(x)==1){
      if(nchar(x)==0){NULL }
      else{x}
    }else{x})
  
  if(sum(sapply(input,function(x) is.null(x) ) ) == length(input)){
    return(NULL)
  }
  
  #print(input)
  if(plotVENN){ venn(input) }
  
  i<-intersect(input[[1]],input[[2]])
  s1<-setdiff(input[[1]],input[[2]])
  s2<-setdiff(input[[2]],input[[1]])
  
  return(list(inter=i,diffAB=s1,diffBA=s2))
}

#-------------------------------------------
#' @title Venndiagramm 3-Set
#' @author Claus Weinholdt
#' @description Venndiagramm for 3-Set
#' @param A is the 1. set
#' @param B is the 2. set
#' @param C is the 3. set
#' @param name vector with names of sets
#' @param plotVENN is boolean to enable venn plot
#' @export
f.input3 = function (p,q,r,name=c("A","B","C"),plotVENN=TRUE,vennOut=FALSE){
  require(gplots)
  
  input  <-list(A=p,B=q,C=r)
  names(input)<-name
  #print(input)
  
  input <- lapply(input,function(x) 
    if(length(x)==1){
      if(nchar(x)==0){NULL }
      else{x}
    }else{x})
  
  if(sum(sapply(input,function(x) is.null(x) ) ) == length(input)){
    return(NULL)
  }
  
  if(plotVENN){ ve <- venn(input) }
  
  if(vennOut){ return(ve) 
  } else { return(intersect3(p,q,r)) }
}


#' @title Venndiagramm 4-Set
#' @description Venndiagramm for 4-Set
#' @author Claus Weinholdt
#' 
#' @export
f.input4 = function (p,q,r,s,t,name=c("A","B","C","D"),plotVENN=TRUE,vennOut=FALSE){
  require(gplots)
  
  input  <-list(A=p,B=q,C=r,D=s)
  names(input)<-name
  
  input <- lapply(input,function(x) 
    if(length(x)==1){
      if(nchar(x)==0){NULL }
      else{x}
    }else{x})
  
  if(sum(sapply(input,function(x) is.null(x) ) ) == length(input)){
    return(NULL)
  }
  
  #if(plotVENN){ venn(input,simplify=TRUE) }
  if(plotVENN){ venn(input,simplify=FALSE) }
  
  if(vennOut){ 
    ve <- venn(input,show.plot = F)
    return(ve)  
  } else { 
    return(intersect4( p,q,r,s)) 
  }
  
}

#' @title Venndiagramm 5-Set
#' @description Venndiagramm for 5-Set
#' @author Claus Weinholdt
#' @export
f.input5 = function (l1,l2,l3,l4,l5,name=c("A","B","C","D","E"),plotVENN=TRUE,vennOut=FALSE){
  require(gplots)
  
  input  <-list(A=l1,B=l2,C=l3,D=l4,E=l5)
  names(input)<-name
  
  input <- lapply(input,function(x) 
    if(length(x)==1){
      if(nchar(x)==0){NULL }
      else{x}
    }else{x})
  
  if(sum(sapply(input,function(x) is.null(x) ) ) == length(input)){
    return(NULL)
  }
  
  if(plotVENN){ venn(input,simplify=FALSE) }
  #print(input)
  
  
  if(vennOut){ 
    ve <- venn(input,show.plot = F)
    return(ve)  
  } else { 
    return(intersect5(l1,l2,l3,l4,l5)) 
  }
  
  
}
#' @title intersect 3-Set
#' @description intersect for 3-Set
#' @author Claus Weinholdt
#' @export
intersect3<-function(A,B,C){
  # O<-intersect(A,intersect(B,C))
  O <- Reduce(intersect,list(A,B,C))
  return(O)
}

#' @title intersect 4-Set
#' @description intersect for 4-Set
#' @author Claus Weinholdt
#' @export
intersect4<-function(A,B,C,D){
  # O<-intersect(A,intersect(B,intersect(C,D)))
  O <- Reduce(intersect,list(A,B,C,D))
  return(O)
}

#' @title intersect 5-Set
#' @description intersect for 5-Set
#' @author Claus Weinholdt
#' @export
intersect5<-function(A,B,C,D,E){
  # O<-intersect(A,intersect(B,intersect(C,intersect(D,E))))
  O <- Reduce(intersect,list(A,B,C,D,E))
  return(O)
}

#' @title union 3-Set
#' @description union for 3-Set
#' @author Claus Weinholdt
#' @export
union3<-function(A,B,C){
  # O<-union(A,union(B,C))
  O <- Reduce(union,list(A,B,C))
  return(O)
}

#' @title union 4-Set
#' @description union for 4-Set
#' @author Claus Weinholdt
#' @export
union4<-function(A,B,C,D){
  # O<-union(A,union(B,union(C,D)))
  O <- Reduce(union,list(A,B,C,D))
  return(O)
}

#' @title union 5-Set
#' @description union for 5-Set
#' @author Claus Weinholdt
#' @export
union5<-function(A,B,C,D,E){
  # O<-union(A,union(B,union(C,union(D,E))))
  O <- Reduce(intersect,list(A,B,C,D,E))
  return(O)
}

f.input5.pretty = function (p,q,r,s,e,name,VennName=""){
  input  <-list(A=p,B=q,C=r,D=s,E=e)
  names(input)<-name
  
  require("VennDiagram")
  venn.plot <- venn.diagram(
    x = input,
    filename = paste(sep="_",VennName,"Venn_5set_pretty.tiff"),
    col = "black",
    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    alpha = 0.50,
    cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
            1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
    cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    cat.cex = 1.5,
    cat.fontface = "bold",
    margin = 0.05
  );
  return(input)  
}

f.input4.pretty = function (p,q,r,s,name,VennName=""){
  input  <-list(A=p,B=q,C=r,D=s)
  names(input)<-name
  
  require("VennDiagram")
  venn.plot <- venn.diagram(
    x = input,
    filename = paste(sep="_",VennName,"Venn_4set_pretty.tiff"),
    col = "transparent",
    fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
    alpha = 0.50,
    label.col = c("orange", "white", "darkorchid4", "white", 
                  "white", "white", "white", "white", "darkblue", "white", 
                  "white", "white", "white", "darkgreen", "white"),
    cex = 1.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
    cat.cex = 1.5,
    cat.pos = 0,
    cat.dist = 0.07,
    cat.fontfamily = "serif",
    rotation.degree = 270,
    margin = 0.2
  );
}