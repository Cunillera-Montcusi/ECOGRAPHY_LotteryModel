
### Cunillera-Montcusí et al. scripts for the functions used to run the lottery model to create the initial 
#   metacommunities based on the networks having 
#   different linkage distances. 
#   Also the script to simulate fake fires similar to the real one. 
#   Finally, the script to run the lottery model and recolonization after burning of specific nodes.

########################################################################################################################################################
########################################################################################################################################################
############################### METACOMMUNITY SPLITTING FUNCTION #######################################################################################
############################### (Manuscript Figure 2 FIRST STEP) #######################################################################################
########################################################################################################################################################

# Function to split the metacommunity between the SHARED and EXCLUSIVE for each species pool.
#   grafos= the list of the containing graphs 
#   metacom= the total metacommunity to be splitted
#   alphas= vector with all alpha values that has to correspond 
split.metacomm.alpha<-function(grafos, metacom, alphas=c(0.2)){
  N.pool<-length(grafos)  # identify number of graphs and therefore of dispersal perceptions
  S<-nrow(metacom)        # species richness
  gr.member<-rep(NA,S)    # Creates a vector of NAs as many times as species in the metacommunity 
  spp.i<-floor(nrow(metacom)/N.pool) # Split the metacommunity species into the number of dispersal perceptions
  for(gg in 1:N.pool){ # for each dispersal perception do the following actions
    gr.member[(1+(gg-1)*spp.i):(gg*spp.i)]<-gg # creates the dummy variable. It is a matrix with as many columns as dispersal perceptions. For each column a proportional 
    #number of species (rows) are filled with the corresponding number of the N.pool. 
    # Doing this, we create a dummy value that links metacommunity species with the considered dispersal perceptions 
  }
  if(length(which(is.na(gr.member)))>0)  gr.member[which(is.na(gr.member))]<-N.pool # CONDITION: If gr.member has no more NAs add the number of the N.pool 
  gr.member<-dummy(gr.member) # convert the gr.member into a dummy variable (1 or 0)
  colnames(gr.member)<-1:N.pool # Set the colnames of the gr.member with the corresponding graph number 
  
  metacom.e<-NULL  # Generate a null object for the EXCLUSIVE part  
  metacom.s<-NULL  # Generate a null object for the SHARED part 
  for(me in 1:N.pool){ # do the following actions for each graph
    ii<-which(gr.member[,me]==1) # create a vector containing all the species with the same dispersal perception 
    for(i in 1:length(ii)){ # Do the following actions for all the species having different dispersal perceptions 
      metacom.e<-rbind(metacom.e,round(metacom[ii[i],]*(1-alphas[me])))# splits the whole metacommunity abundance according to 1-aplha value creating the EXCLUSIVE part 
      metacom.s<-rbind(metacom.s,round(metacom[ii[i],]*alphas[me]))    # splits the whole metacommunity abundance according to 1-aplha value creating the SHARED part 
    }
  }
  out<-list( "shared"=metacom.s, "exclusive"=metacom.e) # Generates the output 
  out
}


########################################################################################################################################################
########################################################################################################################################################
######################## GENERATION OF METACOMMUNITIES FROM BEFORE THE WILDFIRE  (LOTTERY MODEL) #######################################################
##################################### (Manuscript Figure 2 SECOND STEP) ################################################################################
########################################################################################################################################################

### The following functions are used to generate the initial metacommunities based on each network (graph) before the wildfires with accompainig functions
## this first script is obtained and modified from Borthagaray et al. (2014). Ecography (doi: 10.1111/j.1600-0587.2013.00366.x)

### N0 is the metacommunity abundance per species (it could be a staritng condition as uniform distribution or the output of previous simulations)
### metacom.s and metacomm.e are the original matrices. starting with: metacom<-matrix(rep(No,N.m),length(No),N.m)# generates a matrix spp x community (N.m= N nodos)
### out.t<-matrix(c(1,0,0),1,3)
### n: is the number of individuals dead and replaced in each local communities in a single iteration
### grafos: is a list of graphs, each graph represents a taxon-dependent metacommunity network
### alphas: is the fraction of individuals (resources) shared among species pools (e.g. niche overlap); 
### (1-aplpha) is the fraction of individuals (resources) exclusive of the species pool j

# metacom<-matrix(rep(No,N.m),length(No),N.m)
# out.t<-matrix(c(1,0,0),1,3)
# split.metacomm.alpha(grafos = redes,metacom = entrada[[2]], alphas = c(0.1,.5,.9))->mm
# entrada<-list("contrlo.S"=out.t, metacom.s=mm[[1]], metacom.e=mm[[2]], metacom=metacom)
meta.comm_neutral_asimetrica<-function(grafos, it, metacom.e=entrada[[4]], metacom.s=entrada[[3]], m.in, 
                                       m.out, n=50, out.t=entrada[[1]], cada=cada, alphas=c(0.2,0.2)){
  require(igraph)
  require(dummies)
  #  if(length(redes)!=length(alphas)) return("number of networks and number of alphas have to be equal")  
  #  par(mfrow=c(1,1))
  N.pool<-length(grafos)                                     # identify number of species pools (i.e., species perceptions)
  n.e<-round(n*(1-alphas))                                   # numbers of individuals to kill and replace internal to the species group (no shared)
  n.s<-n-n.e                                                 # numbers of individuals shared with other species
  n.s<-sum(n.s)
  S<-nrow(metacom.s)
  gr.member<-rep(NA,S)                           # Dummy variables with species group memebership of species
  spp.i<-floor(nrow(metacom.s)/N.pool)           # Florr rounds numbers downwards ("a la baixa")
  for(gg in 1:N.pool){
    gr.member[(1+(gg-1)*spp.i):(gg*spp.i)]<-gg # Sets species dispersal ability to equal parts of the metacommunity pool
  }
  if(length(which(is.na(gr.member)))>0)  gr.member[which(is.na(gr.member))]<-N.pool # Sets dispersal ability (the maximum) to the last species of the pool 
  gr.member<-dummy(gr.member) # Makes the gr.members variable a dummy (a factor-like).
  colnames(gr.member)<-1:N.pool
  N.m<-length(degree(grafos[[1]]))                                      # identifies the number of local communities (metacommunity perceptions)
  puntos<-seq(1,it,cada) # vector to calculate some metrics (alfa, CV.alfa and beta) to check stavility once 
  # plot(c(-1,-1),ylim=c(0,sum(metacom.e[1:round(nrow(metacom.e)/N.pool,0),])+sum(metacom.s[1:round(nrow(metacom.e)/N.pool,0),])), xlim=c(1,it), xlab="iteration", ylab="groups abundances")
  #__________________ Iterations starting here: 
  # Here start iterations #   
  for (i in 1:it){
    cat("va en  ", i, "\n")
    orden<-sample(N.m)                                                  # randomize the sequence of commmunities to update 
    for(j in orden){                                                    # update local communities abundances
      submeta.in<-rep(0,S)
      submeta.out<-rep(0,S)
      # Lines for generate the neutral assembly for the EXCLUSIVE PART:
      for (g in 1:N.pool){
        id.neighbor.in<-neighborhood(grafos[[g]],1,j, mode="in")[[1]][-1] #  identifies which are the neighbors of the community to update (remove the first node that is itself)
        if(length(id.neighbor.in)>1)  submeta.in.temp<- apply((metacom.s[,id.neighbor.in]+metacom.e[,id.neighbor.in]),1,sum) # if have neighbors, sum the abudnance lper species among all community neighbors (this will be the species abundances for immigration)
        if(length(id.neighbor.in)==1) submeta.in.temp<-       (metacom.s[,id.neighbor.in]+metacom.e[,id.neighbor.in]) # if it only have one neighbors the same but not sum
        if(length(id.neighbor.in)==0) submeta.in.temp<-                     metacom.s[,j]+metacom.e[,j] # if no neighbours, your value is yourself
        submeta.in<-submeta.in+submeta.in.temp*gr.member[,g]
        
        # the same from above but with the out part
        id.neighbor.out<-neighborhood(grafos[[g]],1,j, mode="out")[[1]][-1] # identifies which are the neighbors of the community to update (remove the first node that is itself)
        if(length(id.neighbor.out)>1)  submeta.out.temp<- apply((metacom.s[,id.neighbor.out] + metacom.e[,id.neighbor.in]),1,sum) # if have neighbors, sum the abudnance lper species among all community neighbors (this will be the species abundances for immigration)
        if(length(id.neighbor.out)==1) submeta.out.temp<-        metacom.s[,id.neighbor.out] + metacom.e[,id.neighbor.in]
        if(length(id.neighbor.out)==0) submeta.out.temp<-                      metacom.s[,j] + metacom.e[,j]
        submeta.out<-submeta.out+submeta.out.temp*gr.member[,g]
        
        ii<-which(gr.member[,g]==1) # create an ID vector with the species of an specific dispersal ability 
        
        if(sum(metacom.e[ii,j])<n.e[g]) n.e[g]<-sum(metacom.e[ii, j])-1 # CONDITION: If the abundance in the site and for the selected species is less than n.e. change n.e. for a smaller number 
        if(n.e[g]>0 & sum(metacom.e[ii,j])>0) { # CONDITION: If all these values are greater than 0, carry the following lines.
          metacom.e[ii,j]<-neutral.2(mcomm.update=as.vector(metacom.e[ii,j]), mcomm.immigrant.in = submeta.in[ii], 
                                     mcomm.immigrant.out = submeta.out[ii],  m.in, m.out, n.e[g]) # Do the function neutral.2 to update abundances following neutral drift 
        }
      }
      # Lines for generate the neutral assembly for the SHARED PART:
      if(sum(metacom.s[,j])<n.s)   n.s<-sum(metacom.s[,j])-1 # CONDITION: If the abundance in the site and for the selected species is less than n.e. change n.e. for a smaller number
      metacom.s[,j]<-neutral.2(mcomm.update=as.vector(metacom.s[,j]), mcomm.immigrant.in=as.vector(submeta.in),
                               mcomm.immigrant.out = as.vector(submeta.out),  m.in, m.out, n.s)# Do the function neutral.2 to update abundances following neutral drift
    } 
    
    if(any(puntos==i)==T){ # When the puntos vector coincides with the iteration number do the following functions. 
      # The following functions are made to keep a control the stability  
      
      #controla<-NULL
      #for(control in 1:N.pool){
      #  iii<-which(gr.member[,control]==1)
      #  controla<-c(controla,(sum((metacom.e[iii,]))))
      #  points(controla[length(controla)]~i, pch=19, cex=.6, col="red")  
      #  points(controla[length(controla)]~i, pch=19, cex=.8, col=control)  
      #  controla<-c(controla,(sum((metacom.s[iii,]))))
      #  points(controla[length(controla)]~i, pch=19, cex=.6, col="blue")  
      #  points(controla[length(controla)]~i, pch=17, cex=.8, col=control)
      #}
      #image(t(ifelse(metacom.e>0,1,0)))
      #image(t(ifelse(metacom.s>0,1,0)))
      
      metacom<-metacom.e + metacom.s # Creates the global metacommunity joining the shared and the exclusive parts
      sm<-mean(rich(metacom))/nrow(metacom) # Calculates mean alfa richness for the metacom 
      des<-sd(rich(metacom)) # desviation of the alfa richness  
      cvv<-des/sm # variation coeficient of the alfa richness
      bet<-mean(beta_M(metacom)) # beta diversity of the metacom 
      out.t<-rbind(out.t,c(sm,cvv,bet)) # bind all these three values
    }
  }
  colnames(out.t)<-c("alfa","CV.alfa","beta") # change the colnames of the output
  out<-list("transient"=as.matrix(out.t),  "Species.pool"=gr.member, 
            "Metacomunity Shared"=as.matrix(metacom.s),"Metacomunity exclusive"=as.matrix(metacom.e),
            "Metacomunidad"=as.matrix(metacom.s+metacom.e))
  out # Attach the out.t to the first list 
}

#######
neutral.2<-function(mcomm.update, mcomm.immigrant.in, mcomm.immigrant.out, m.in, m.out, n){ # Function to carry the neutral drift on the metacommunity 
  # (separatedly for the SHARED and the EXCLUDED)
  for(rr in 1:n){ # n is the number of individuals eliminated from the metacommunity
    rem<-rmultinom(1, 1,mcomm.update)     # select individuals using a multinomial probability
    mcomm.update<-mcomm.update-rem        # the selected individual is removed
  }
  born<-rmultinom(1,n,
                  mcomm.immigrant.in/sum(mcomm.immigrant.in)*m.in+
                    mcomm.immigrant.out/sum(mcomm.immigrant.out)*m.out+
                    mcomm.update/sum(mcomm.update)*(1-m.in-m.out)) # replace from neighbour (Mmetacom*m) o the same community (mcomm*(1-m)) 
  mcomm.update<-mcomm.update+born # add the new individuals (from neighbouring or born within the community)
  mcomm.update
}
################################ end ############################

#################################
# Function specific metrics
#################################
rich<-function(M){ # Calculates the alfa richness from the metacommunity
  apply(ifelse(M>0,1,0),2,sum)
}
###############################
beta_mean<-function(M){ 
  require(vegan)
  quantile(vegdist(M,method="jaccard",binary = T), c(0.025,.5,.975),na.rm=TRUE)
}
############################
##Estimates the matrix of dissimilarity 
beta_M<-function(M){# Calculates the beta richness from the metacommunity
  require(vegan)
  vegdist(t(M),method="jaccard",binary = T)->B ####DISIMILITUD
  as.matrix(B)->B
  B
}

#______________________________________________________________________________________________________________________________________________________#

########################################################################################################################################################
########################################################################################################################################################
######################## FUNCTION THAT GENERATES FAKE FIRES ############################################################################################
######################### (Manuscript Figure 2 THIRD STEP) #############################################################################################
########################################################################################################################################################

#   This function will simulate a "wildfire" which basically will be the creation of a data.frame with a list of all ponds and then several 
# columns with "1" or "0" that correspond to being burned or not. Each column will correspond to all wildfire intensities (0.05,0.01, etc.).
#   In addition, an extra column is created with the corresponding size of the wildfire (jj).

# XY= UTM coordinates from the studied network that will to be burned. In our case: The albera network 
# fire0 is the "id" or the vector of "ids" for the rows with the pond where the fire will starts
# niveles-> How many "cuts" there will be between the 0.05 and 0.95 on the distance that will affect ponds. In our case twenty.
# efectividad-> How many "cuts" there will be between 0 and 1 indicating the percentage of affected ponds within the simulated area. In our case twenty.

rangos<-function(XY, fire0, niveles, efectividad){         # Creation of the function
  DD<-as.matrix(dist(XY,method = "euclidian"))[, fire0]    # Generation of the matrix of distances in relation to FIRE0 pond (where fire starts)
  cortes.fuegos<-c(0,quantile(DD, seq(0.05,.95,,niveles))) # Vector creating the distances that will be used according to "NIVELES" 
  IDs<-NULL                                                # IDs is a NULL thing that will be filled 
  jj<-0                                                    #
  ef.t<-seq(0,1,,efectividad+1)[-1]                        # Vector creating the values of efficiency of burning of the pond within the area
  for(i in cortes.fuegos[2:length(cortes.fuegos)]){        # For each of the areas ("cortes.fuego") run this
    jj<-jj+1                                               # Adds 1 at the jj which will be the "counter" for each area
    ii<-which(DD<i & DD>cortes.fuegos[jj])                 # Selects all the IDs of the ponds within the selected distances 
    out.t<-cbind(jj,round(i,2),ii)                         # Generates an out.t-> binding the before created vectors
    for(prop in ef.t){                                     # For each of the proportions ("ef.t") runs the following actions
      n.patch<-round(length(ii)*prop,0)                      # Multiply each proportion for the number of ponds "ii" creating the number of "n.patch"
      affected<-sample(c(rep(0,length(ii)-n.patch),rep(1,n.patch)), replace=FALSE) # Randomly adds a 0 to all ponds-n.patch and adds a 1 to n.patch  
      out.t<-cbind(out.t,affected)                           # Adds to the out.t the a column with 0 and 1 for each pond id
    }                                                      # 
    IDs<-rbind(IDs,out.t)                                  # Binds everything at IDs
  }                                                        #
  colnames(IDs)[4:ncol(IDs)]<-round(ef.t,2)                # Change the names of the columns
  IDs                                                      #
}

#_______________________________________________________________________________________________________________________________________________________#

########################################################################################################################################################
########################################################################################################################################################
##################### FUNCTION OF BURNING AND RECOLONIZATION (LOTTERY MODEL) OF AFFECTED WATER BODIES ##################################################
############################### (Manuscript Figure 2 FOURTH STEP) ######################################################################################
########################################################################################################################################################

# This function must be used with the "rangos" function output. 
# It basically do two actions:
#   1: Burns the specified ponds in the "rangos" output for each level of Distance (jj) and Intensity (every column of the output)  
#   2: Recolonize these burned ponds following neutral drift with the function fill.n.2 
#
# M is the output of the "rangos" function. It is a data.frame with the names of all ponds and columns with "1" or "0" that indicate if the ponds are burned or not.
# grafos is a list of metacommunity graphs defined by the linkage distance associated to each dispersal group
# out[[5]] is the fifth element in the output of "meta.comm_neutral_asimetrica()"
# out[[3]] is the shared fraction of the metacommunity
# out[[4]] is the exclusive fraction of the metacommunity
# migra is the migration rate (same as the meta.comm_neutral_asimetrica())
# it is the number of iterations
# alphas are the alphas values (same than in meta.comm_neutral_asimetrica())
# dummy.pool is the second element in the output of "meta.comm_neutral_asimetrica()" that corresponds to the different dispersal groups.
# lista.n.grupos is a list with the same number of elements than dispersal distances. It is filled within the function.
burning_recolonization_function.2<-function(M, metacom=out[[5]], metacom.s=out[[3]], metacom.e=out[[4]], 
                                            migra, grafos, it, alphas, dummy.pool= out[[2]],lista.n.grupos){               
  out.spp<-matrix(NA,nrow=length(4:ncol(M)),ncol=length(unique(M[,1]))) # Creates a matrix with rows having the % efficiency and columns having the burned "areas" 
  colnames(out.spp)<-sort(unique(M[,1])) # Set the names of the output columns according to the different burned areas (M[,1])
  rownames(out.spp)<-colnames(M)[4:ncol(M)] # Set the names of the output rowa according to the different burning intensitites (the columns of M matrix)
  out.it<-out.spp # Copi out.spp format to the out.it
  out.sd.spp<-out.spp # Copi out.spp format to the out.it
  out.it.sd<-out.spp # Copi out.spp format to the out.it
  
  lista.iterations <- lista.n.grupos
  lista.sd.iterations <- lista.n.grupos
  lista.sd.n.grupos <- lista.n.grupos
  
  for(gg in 1:length(alphas)){ # for each one of the alphas (one by each dispersal ability) 
    lista.n.grupos[[gg]]<-out.spp #copy the out.spp structure into the elements of the list lista.n.grupos
    lista.iterations[[gg]] <-out.it
    lista.sd.n.grupos[[gg]] <- out.sd.spp
    lista.sd.iterations[[gg]] <- out.it.sd
  }
  
  fila<-0                                                     # Sets fila=0 to build the final matrix
  for(i in 4:ncol(M)){                                        # For each column of M data.frame do the following actions (each column corresponds to all fire intensities)
    fila<-fila+1                                              # Add +1 to each fila (increase the row number in every iteration)
    columna<-0                                                # Sets column=0 to build the final matrix
    for(j in sort(unique(M[,1]))){  # For the column of M data.frame do the following actions (ion this column  all fire areas are contained as numbers (from 1 to the specified "rangos" num.)
      columna<-columna+1                                      # Add +1 to each columna (increase the column number in every iteration)
      cat("Intensity ", colnames(M)[i], "Distance ", j, "\n") # output to follow iterations progress 
      ii<-M[which(M[,1]<=j & M[,i]==1),3]                     # for each pond that have the value "j" and is burned in the column i=1 do the following actions:
      #return(ii)
      metacom[,ii]<-0                                         # brings each burned pond to 0 for the whole metacommunity
      metacom.s[,ii]<-0                                       # brings each burned pond to 0 for the shared part of the metacommunity 
      metacom.e[,ii]<-0                                       # brings each burned pond to 0 for the exclusive part of the metacommunity
      #return(metacom)
      FILL<-fill.n.2(grafos=grafos, # grafos is the list of graphs                                           
                     J.metacomm=metacom, # the whole metacom
                     J.metacomm.s=metacom.s, # the shared part of the whole metacom
                     J.metacomm.e=metacom.e, # the exclusive part of the whole metacom
                     dummy.pools= dummy.pool,  m=migra, alphas=alphas, it=it) # Do the fill.n.2 for the metacommunity (most items are inherated from the initial inputs of the functions)
      #return(FILL)
      E<-FILL[[1]]  # Retains the whole metacommunity (both shared and exclusive parts) of the fill.n.2 function 
      lista.n.grupos[[1]][fila,columna]<-mean(apply(ifelse(as.matrix(E[which(dummy.pool[,1]==1),ii])>0,1,0),2,sum))# For the first dispersal ability, calculates mean richness after colonizaiton
      lista.n.grupos[[2]][fila,columna]<-mean(apply(ifelse(as.matrix(E[which(dummy.pool[,2]==1),ii])>0,1,0),2,sum))# For the second dispersal ability, calculates mean richness after colonizaiton
      lista.n.grupos[[3]][fila,columna]<-mean(apply(ifelse(as.matrix(E[which(dummy.pool[,3]==1),ii])>0,1,0),2,sum))# For the third dispersal ability, calculates mean richness after colonizaiton
      lista.n.grupos[[4]][fila,columna]<-mean(apply(ifelse(as.matrix(E[which(dummy.pool[,4]==1),ii])>0,1,0),2,sum))# For the fourth dispersal ability, calculates mean richness after colonizaiton
      lista.n.grupos[[5]][fila,columna]<-mean(apply(ifelse(as.matrix(E[which(dummy.pool[,5]==1),ii])>0,1,0),2,sum))# For the fifth dispersal ability, calculates mean richness after colonizaiton
      lista.n.grupos[[6]][fila,columna]<-mean(apply(ifelse(as.matrix(E[which(dummy.pool[,6]==1),ii])>0,1,0),2,sum))# For the sixth dispersal ability, calculates mean richness after colonizaiton
      lista.n.grupos[[7]][fila,columna]<-mean(apply(ifelse(as.matrix(E[which(dummy.pool[,7]==1),ii])>0,1,0),2,sum))# For the seventh dispersal ability, calculates mean richness after colonizaiton
      lista.n.grupos[[8]][fila,columna]<-mean(apply(ifelse(as.matrix(E[which(dummy.pool[,8]==1),ii])>0,1,0),2,sum))# For the eightth dispersal ability, calculates mean richness after colonizaiton
      lista.n.grupos[[9]][fila,columna]<-mean(apply(ifelse(as.matrix(E[which(dummy.pool[,9]==1),ii])>0,1,0),2,sum))# For the nineth dispersal ability, calculates mean richness after colonizaiton
      
      lista.iterations[[1]][fila,columna]<-mean(FILL[[2]][1,][-which(FILL[[2]][1,]==0)]) 
      lista.iterations[[2]][fila,columna]<-mean(FILL[[2]][2,][-which(FILL[[2]][2,]==0)]) 
      lista.iterations[[3]][fila,columna]<-mean(FILL[[2]][3,][-which(FILL[[2]][3,]==0)]) 
      lista.iterations[[4]][fila,columna]<-mean(FILL[[2]][4,][-which(FILL[[2]][4,]==0)]) 
      lista.iterations[[5]][fila,columna]<-mean(FILL[[2]][5,][-which(FILL[[2]][5,]==0)]) 
      lista.iterations[[6]][fila,columna]<-mean(FILL[[2]][6,][-which(FILL[[2]][6,]==0)]) 
      lista.iterations[[7]][fila,columna]<-mean(FILL[[2]][7,][-which(FILL[[2]][7,]==0)]) 
      lista.iterations[[8]][fila,columna]<-mean(FILL[[2]][8,][-which(FILL[[2]][8,]==0)]) 
      lista.iterations[[9]][fila,columna]<-mean(FILL[[2]][9,][-which(FILL[[2]][9,]==0)]) 
      
      lista.sd.n.grupos[[1]][fila,columna]<-sd(apply(ifelse(as.matrix(E[which(dummy.pool[,1]==1),ii])>0,1,0),2,sum))# For the first dispersal ability, calculates mean richness after colonizaiton
      lista.sd.n.grupos[[2]][fila,columna]<-sd(apply(ifelse(as.matrix(E[which(dummy.pool[,2]==1),ii])>0,1,0),2,sum))# For the second dispersal ability, calculates mean richness after colonizaiton
      lista.sd.n.grupos[[3]][fila,columna]<-sd(apply(ifelse(as.matrix(E[which(dummy.pool[,3]==1),ii])>0,1,0),2,sum))# For the third dispersal ability, calculates mean richness after colonizaiton
      lista.sd.n.grupos[[4]][fila,columna]<-sd(apply(ifelse(as.matrix(E[which(dummy.pool[,4]==1),ii])>0,1,0),2,sum))# For the fourth dispersal ability, calculates mean richness after colonizaiton
      lista.sd.n.grupos[[5]][fila,columna]<-sd(apply(ifelse(as.matrix(E[which(dummy.pool[,5]==1),ii])>0,1,0),2,sum))# For the fifth dispersal ability, calculates mean richness after colonizaiton
      lista.sd.n.grupos[[6]][fila,columna]<-sd(apply(ifelse(as.matrix(E[which(dummy.pool[,6]==1),ii])>0,1,0),2,sum))# For the sixth dispersal ability, calculates mean richness after colonizaiton
      lista.sd.n.grupos[[7]][fila,columna]<-sd(apply(ifelse(as.matrix(E[which(dummy.pool[,7]==1),ii])>0,1,0),2,sum))# For the seventh dispersal ability, calculates mean richness after colonizaiton
      lista.sd.n.grupos[[8]][fila,columna]<-sd(apply(ifelse(as.matrix(E[which(dummy.pool[,8]==1),ii])>0,1,0),2,sum))# For the eightth dispersal ability, calculates mean richness after colonizaiton
      lista.sd.n.grupos[[9]][fila,columna]<-sd(apply(ifelse(as.matrix(E[which(dummy.pool[,9]==1),ii])>0,1,0),2,sum))# For the nineth dispersal ability, calculates mean richness after colonizaiton
      
      lista.sd.iterations[[1]][fila,columna]<-sd(FILL[[2]][1,][-which(FILL[[2]][1,]==0)]) 
      lista.sd.iterations[[2]][fila,columna]<-sd(FILL[[2]][2,][-which(FILL[[2]][2,]==0)]) 
      lista.sd.iterations[[3]][fila,columna]<-sd(FILL[[2]][3,][-which(FILL[[2]][3,]==0)]) 
      lista.sd.iterations[[4]][fila,columna]<-sd(FILL[[2]][4,][-which(FILL[[2]][4,]==0)]) 
      lista.sd.iterations[[5]][fila,columna]<-sd(FILL[[2]][5,][-which(FILL[[2]][5,]==0)]) 
      lista.sd.iterations[[6]][fila,columna]<-sd(FILL[[2]][6,][-which(FILL[[2]][6,]==0)]) 
      lista.sd.iterations[[7]][fila,columna]<-sd(FILL[[2]][7,][-which(FILL[[2]][7,]==0)]) 
      lista.sd.iterations[[8]][fila,columna]<-sd(FILL[[2]][8,][-which(FILL[[2]][8,]==0)]) 
      lista.sd.iterations[[9]][fila,columna]<-sd(FILL[[2]][9,][-which(FILL[[2]][9,]==0)])
      
      #for (ttt in 1:length(alphas)) {
      #if(length(which(FILL[[2]][ttt,]>0))>0) # Creates the output (out.it) containing the number of iterations to fill the burned ponds 
      #  lista.iterations[[ttt]][fila,columna]<-mean(FILL[[2]][ttt,][-which(FILL[[2]][ttt,]==0)]) 
      #}
      
      #return(FILL)    
      #      for (g in 1:length(alphas)){
      #       iii<-which(dummy.pool[,g]==1)
      #        S.mean<-mean(apply(ifelse(as.matrix(E[iii,ii])>0,1,0),2,sum))   # Creates the S.mean by summ the final richness for affected ponds and calculate de mean
      #        lista.n.grupos[[g]][fila,columna]<-S.mean/length(iii)           # creates the out with the columns and the rows corresponding to columna and fila values
      #        if(length(which(FILL[[2]]>0))>0)
      #         out.it[[g]][fila,columna]<-mean(FILL[[2]][-which(FILL[[2]]==0)])
      #      }
    }
  }
  # Rownames
  out<-list("Prop.spp.disponibles"=lista.n.grupos, 
            "Iteraciones_para_fill"=lista.iterations,
            "sd.spp.disponibles"=lista.sd.n.grupos,
            "sd.iteraciones.fill"=lista.sd.iterations) # Creates the out of the whole function 
  out
}

######### function to refill the emptied ponds by the wildfire accordingly to their non-burned neighbours used at the end of burning_recolonization_function function
# The filling process is based on neutral drift and considers the burned ponds neighbours.
# grafos is a list with graphs associated to the dispersal distances to consider
# alphas is a vector with the fraction of individuals in each dispersal group associated to shared (alpha) or exclusive (1-alpha) resources (see above)
fill.n.2<-function(grafos, J.metacomm,J.metacomm.s,J.metacomm.e, dummy.pools,  m, alphas, it){
  require(igraph)        # Package needed to work with graphs
  J<-sum(J.metacomm[,1]) #  J is equal for all local comunity, therefore this value is the same for all columnes (i.e., column)
  N.pool<-length(grafos) # N.pool is the number of graphs, thus the number of considered dispersal perceptions
  out.n.it<-matrix(0,nrow=N.pool,ncol=ncol(J.metacomm)) # out.n.it is a matrix with the as many columns as the number of ponds 
  out<-list("Metacomunidad"=J.metacomm, 
            "n.it.fill"=out.n.it,
            "Metacomunidad.s"=J.metacomm.s,
            "Metacomunidad.e"=J.metacomm.e,
            "Dummy.Spp.Pool"=dummy.pools) # generates the out of the function, which will be used to run the fill function to fill burned ponds through neutral drift  
  
  N.local.per.group <- matrix(NA, nrow =length(alphas), ncol=ncol(J.metacomm))
  
  for(k in 1:it){       # Number of times that all this process will be repeated. 
    for (eee in 1:N.pool) {
      N.local.per.group[eee,] <- apply(as.matrix(out[[1]][which(out[[5]][,eee]==1),]),2,sum)
    }
    
    #return(N.local.per.group)
    N.local <- as.matrix(ifelse(N.local.per.group>0,1,0))
    #return(N.local)
    id.burn <- which(apply(N.local,2,sum)<N.pool)   #quines rises quan trobis aixo!!!!! VIVA LED ZEPPELIN
    if(length(id.burn)==0) return(out)  # CONDITION: If there is no pond with 0 abundance return the same J.metacomm
    
    for(i in id.burn){                       # Repeat the following steps for each burned pond
      pool.i.shared<-rep(0,nrow(out[[1]])) # shared pool of available species  
      total.neigh<-NULL # a null object that will be used 
      
      # ITERATION COUNTER    
      for(ooo in 1:N.pool){ # Do the following actions for all dispersal abilities
        # CONDITION: If the abundance of the dispersal group (ooo) in the burned pond (i) is equal to zero implies that a new iteration 
        #will be needed to reach the current pond (i) 
        if(N.local.per.group[ooo,i]==0) out[[2]][ooo,i]<-out[[2]][ooo,i]+1  # Adds a +1 each time that these action is done for the i pond
      }
      
      # REFILLING PROCESS following neutral drift for the EXCLUSIVE PART OF THE METACOMMUNITY
      for (g in 1:N.pool){   # Repeat the following actions for each N.pool value (each dispersal ability)
        
        ii<-which(out[[5]][,g]==1) # A vector with the "id" of species with an specific dispersal ability 
        #J.e<-round(J*(1-alphas[g]),0) 
        J.e<-round((J/N.pool)*(1-alphas[g]),0) # Metacommunity size of the exclusive part of the metacommunity 
        id.neigh<-neighborhood(grafos[[g]],order=1,mode = "all", node=i)[[1]][-1] # Creates the vector with the id of burned ponds neighbours for this dispersal distance
        disponibles.e<-length(which(is.na(match(id.neigh,id.burn))==T)) # calculate the number of non-burned neighbours available (connected to the burned pond)
        #return(id.neigh)
        if(disponibles.e>0){ #CONDITION:  If the value of non-burned neigh is greater than 1 keeps filling them
          
          if(length(id.neigh)>1)  pool.i<-apply(out[[1]][,id.neigh],1,sum)  #CONDITION: If more than 1 neighbour 
          # creates a vector with species abundance of nonburned ponds connected to the pond to fill
          if(length(id.neigh)==1) pool.i<-out[[1]][,id.neigh]               #CONDITION: If only 1 neighbour 
          # creates a vector with species abundance of the nonburned pond connected to the pond to fill
          if(length(id.neigh)>0)  { # CONDITION: If the number of neighbours is greater than 0 do the following actions
            
            if(sum(pool.i[ii])==0) out[[4]][ii,i] <- 0 #CONDITION: If pool.i is 0 (neighbours with 0 abundance (i.e., burned)) maintain to 0 the abundance of the exclusiva part 
            # of the pond to fill. 
            if(sum(pool.i[ii])>0) # CONDITION: If pool.i is greater than zero (some neighbours with individuals) do the following function
              
              out[[4]][ii,i]<-fill(metacomm = pool.i[ii],
                                   J = J.e,
                                   m.local = out[[4]][ii,i],
                                   m=m) # Do the function fill to each burned pond. Add the output to the EXCLUSIVE part out.
            # 1: pool.i as metacommunity
            # 2: J as the final abundance size
            # 3: community of the burned pond being filled
            # 4: m: migration rate
            
            pool.i.shared[ii]<-pool.i.shared[ii]+ pool.i[ii] # It accumulates all the pools in one vector to built a global pool for the following shared metacommunity
            total.neigh<-c(total.neigh,id.neigh) # It accumulates all the neighbours IDs in one vector in order to built a global pool of neighbours for the following shared metacommunity
          }
        }
      }
      
      # REFILLING PROCESS following neutral drift for the SHARED PART OF THE METACOMMUNITY      
      J.s<-round(sum((J/N.pool)*alphas),0) # Metacommunity size of the exclusive part of the metacommunity 
      if(length(total.neigh)==0) out # CONDITION: there are no neighbours at all (completely disconnected burned ponds) return out directly
      
      if(length(total.neigh)>0) # CONDITION: if there are some neighbours do the following actions
        for (gg in 1:N.pool){ # For each dispersal ability do the following functions 
          ee<-which(out[[5]][,gg]==1) # Create a vector with the species IDs that correspond to each dispersal ability
          
          if(sum(pool.i.shared[ee])==0) out[[3]][ee,i] <- 0 # CONDITION: If the global pool (all pools joined together) is 0 (no organisms in the whole pool) 
          # maintain the exlusive part of the pond to 0    
          if(sum(pool.i.shared[ee])>0) # CONDITION: If the global pool is greater than zero do the following actions 
            
            out[[3]][ee,i]<-fill(metacomm = pool.i.shared[ee],
                                 J = round(J.s),
                                 m.local = out[[3]][ee,i], # 
                                 m = m) # Do the function fill to each burned pond. Add the output to the SHARED part out.  
          
        }
      out[[1]][,i]<-out[[3]][,i]+out[[4]][,i] # join the two fractions of the community to the whole metacommunity
      #}  
    }
  }
  out
}

######### funcion fill used at the end of the function fill.

# Function to do the filling based on neutral drift of each pond. It uses de species pool (metacom) based on the non-burned or recolonized neighbours of the selected 
#pond and accounting also with the selected pond metacommunity (m.local) used within the fill.n.2 function.
# J is the size of the community that has to be filled
# m is the migration rate
# the emptied ponds by the wildfire accordingly to their non-burned neighbours 
fill<-function(metacomm, J, m.local,m){
  m.local<-(m.local + rmultinom(1,1,metacomm)) # creates a vector randomly chosing species from the metacom
  for(i in 2:J){                                # Repeats the function as many times as the total abundance of the community
    pool<-m.local*(1-m)+metacomm*m  #   Generates a vector joining spp from the local com (m.local) and from the metacommunity and thus
    # the "pool" is a vector with the proportinal (probability) of sample individuals from each species.
    m.local<-m.local+rmultinom(1,1,pool) #  Adds to the m.local with a randomly chosen spp from the pool (the vector created after)
    # because the structure of rmultinom-> Then, Pool operates as a vector of probabilities
  }
  m.local
}

#______________________________________________________________________________________________________________________________________________________#


######################################################################################################################## 
######################################## END ###########################################################################
######################################################################################################################## 



