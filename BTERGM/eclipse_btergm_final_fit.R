library(tergm)
library(network)
library(igraph)
library(intergraph)
library(btergm)

#------------------------------------------------------------------------------#
network_type = 'comments';

if (network_type == 'comments') {
  input_prefix = './eclipse/data/co-comment/Comments_modif50_';
  
} else if (network_type == 'changes') {
  input_prefix = './eclipse/data/co-change/Changes_modif50_';
}

#------------------------------------------------------------------------------#

get_network_list <- function(interval, start=1, end=50) {
 
   file_names_vector = seq(from=start, to=end, by=interval);
   networks_list <- list()
   
   for (x in file_names_vector) {
     file_name <- paste(input_prefix, x, '.graphml', sep='');
     net <- read_graph(file_name, format='graphml');
     net <- asNetwork(net);
     
     networks_list <- append(networks_list, list(net), after=length(networks_list));
   }
   
  
  #for (x in file_names_vector) {
    #file_name <- paste(input_prefix, x, '.net', sep='');
    #net <- read_graph(file_name, format='pajek');
    #net <- asNetwork(net);
    
    #networks_list <- append(networks_list, list(net), after=length(networks_list));
  #}
  
  return (networks_list);
}


preprocess <- function(networks_list) {
  new_list <- list()
  for (i in c(1:length(networks_list))) {
    print(networks_list[i])
    
    current_network <- networks_list[[i]];
    g <- asIgraph(current_network);

    #Remove isolated nodes
    
    Isolated <- which(degree(g) == 0);
    g <- igraph::delete.vertices(g, Isolated);
    
    current_network <- asNetwork(g)
    
    #finding the metrics of the changed network
    
    degree<-degree(g)
    closeness<-closeness(g)
    betweenness<-betweenness(g)
    pagerank<-page.rank(g)
    eigen_c<-eigen_centrality(g)
    clust_coeff<-transitivity(g, type="localundirected") 
    
    eigen_c<-eigen_c$vector
    pagerank<-pagerank$vector
    class(clust_coeff)
    clust_coeff[is.na(clust_coeff)] <- 0
    
    #adding the metrics to the nodes
    
    current_network %v% 'eigen_c' <- eigen_c
    current_network %v% 'degree' <- degree
    current_network %v% 'betweenness' <- betweenness
    current_network %v% 'clustcoeff' <- clust_coeff
    current_network %v% 'closeness' <- closeness
    current_network %v% 'pagerank' <- pagerank
    
    new_list <- append(new_list, list(current_network))
  }
  return (new_list);
}


x <- get_network_list(1, start = 5, end=50);
y <- preprocess(x)

#------------------------------------------------------------------------------#
#BTERGM fit
result <- btergm(y~nodecov("eigen_c")+nodecov("degree")+nodecov("betweenness")
                                     +nodecov("closeness")+nodecov("pagerank")+nodecov("clustcoeff"),
                                     parallel = "multicore", ncpus = 2, R=500, offset=TRUE,
)

summary(result)




