 

    
    mutSOMA<-function(pedigree.data, p0aa, p0cc, p0tt, p0gg, Nstarts, out.dir, out.name)
    {
    
      allow.neg.intercept="no"
      
    
            ##### Defining the divergence function
            divergence <- function(pedigree, p0aa, p0cc, p0tt, p0gg, param)
            {
              
              ## Initializing parameters
              # c("AA", "CC", "TT", "GG", "AC", "AT", "AG", "CA", "CT", "CG", "TA", "TC", "TG", "GA", "GC", "GT")
              PrAA <- p0aa; PrCC <- p0cc; PrTT <- p0tt; PrGG <- p0gg; PrAC <- p0ac; PrAT <- p0at; PrAG <- p0ag; PrCA <- p0ca; PrCT <- p0ct
              PrCG <- p0cg; PrTA <- p0ta; PrTC <- p0tc; PrTG <- p0tg; PrGA <- p0ga; PrGC <- p0gc; PrGT <- p0gt
                              
              g <- param[1]
             
              
              ## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
              ## Keeping this in this format because I may want to introduce weights for the "unobserved" initial heterozygotes
              svGzero    <- c(PrAA,
                              PrCC,
                              PrTT,
                              PrGG,
                              PrAC,
                              PrAT,
                              PrAG,
                              PrCA,
                              PrCT,
                              PrCG,
                              PrTA,
                              PrTC,
                              PrTG,
                              PrGA,
                              PrGC,
                              PrGT)
              
              ## Defining the generation (or transition) matrix for the mitotic case
              Tmat<-matrix(NA, nrow=16, ncol=16)
              
              Ta<-matrix(c(
                (1-g)^2, 1/9*g^2, 1/9*g^2, 1/9*g^2, 
                1/9*g^2, (1-g)^2, 1/9*g^2, 1/9*g^2, 
                1/9*g^2, 1/9*g^2, (1-g)^2, 1/9*g^2, 
                1/9*g^2, 1/9*g^2, 1/9*g^2, (1-g)^2), nrow=4, byrow=TRUE)
              
              
              #              1                2            3          4            5           6              7           8            9            10           11           12
              Tb<-matrix(c(
                1/3*(1-g)*g,   1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,     1/9*g^2, 
                1/3*(1-g)*g,   1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,
                1/9*g^2,       1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,     1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g,
                1/9*g^2,       1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g), 
                nrow=4, byrow=TRUE)
              
              
              Tc<-matrix(c(
                1/3*(1-g)*g,  1/3*(1-g)*g, 1/9*g^2,     1/9*g^2, 
                1/3*(1-g)*g,  1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,
                1/3*(1-g)*g,  1/9*g^2,     1/9*g^2,     1/3*(1-g)*g,
                1/3*(1-g)*g,  1/3*(1-g)*g, 1/9*g^2,     1/9*g^2, 
                1/9*g^2,      1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*g^2, 
                1/9*g^2,      1/3*(1-g)*g, 1/9*g^2,     1/3*(1-g)*g,
                1/3*(1-g)*g,  1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,
                1/9*g^2,      1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*g^2,
                1/9*g^2,      1/9*g^2,     1/3*(1-g)*g, 1/3*(1-g)*g,
                1/3*(1-g)*g,  1/9*g^2,     1/9*g^2,     1/3*(1-g)*g,
                1/9*g^2,      1/3*(1-g)*g, 1/9*g^2,     1/3*(1-g)*g,
                1/9*g^2,      1/9*g^2,     1/3*(1-g)*g, 1/3*(1-g)*g), nrow=12, byrow=TRUE)
              
              #              1                2            3           4            5           6            7            8           9            10           11           12
              Td<-matrix(c(
                (1-g)^2,     1/3*(1-g)*g,  1/3*(1-g)*g,  1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,
                1/3*(1-g)*g, (1-g)^2,      1/3*(1-g)*g,  1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g,
                1/3*(1-g)*g, 1/3*(1-g)*g,  (1-g)^2,      1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,
                1/9*(g)^2,   1/9*(g)^2,    1/9*(g)^2,    (1-g)^2,     1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,
                1/9*(g)^2,   1/3*(1-g)*g,  1/9*(g)^2,    1/3*(1-g)*g, (1-g)^2,     1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g,
                1/9*(g)^2,   1/9*(g)^2,    1/3*(1-g)*g,  1/3*(1-g)*g, 1/3*(1-g)*g, (1-g)^2,     1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,
                1/9*(g)^2,   1/9*(g)^2,    1/9*(g)^2,    1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   (1-g)^2,     1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,
                1/3*(1-g)*g, 1/9*(g)^2,    1/9*(g)^2,    1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, (1-g)^2,     1/3*(1-g)*g, 1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,
                1/9*(g)^2,   1/9*(g)^2,    1/3*(1-g)*g,  1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, (1-g)^2,     1/9*(g)^2,   1/9*(g)^2, 1/9*(g)^2,
                1/9*(g)^2,   1/9*(g)^2,    1/9*(g)^2,    1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   (1-g)^2,     1/3*(1-g)*g, 1/3*(1-g)*g,
                1/3*(1-g)*g, 1/9*(g)^2,    1/9*(g)^2,    1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g,  1/9*(g)^2,  1/3*(1-g)*g, (1-g)^2,     1/3*(1-g)*g, 
                1/9*(g)^2,   1/3*(1-g)*g,  1/9*(g)^2,    1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/3*(1-g)*g, (1-g)^2), nrow=12, byrow=TRUE)
              
              
              ## Transition matrix
              Tmat[1:4, 1:4]<-Ta
              Tmat[1:4, 5:16]<-Tb
              Tmat[5:16, 1:4]<-Tc
              Tmat[5:16, 5:16]<-Td
              Genmatrix <- Tmat
          
              ## Coding the divergence between samples i and j
              Deffects<-matrix(as.numeric(c("0", "1", "1", "1", "0.5", "0.5", "0.5", "0.5", "1", "1", "0.5", "1", "1", "0.5", "1", "1", "1", "0", "1", "1", "0.5", "1", "1", 
                                           "0.5", "0.5", "0.5", "1", "0.5", "1", "1", "0.5", "1", "1", "1", "0", "1", "1", "0.5", "1", "1", "0.5", "1", "0.5", "0.5", "0.5", 
                                           "1", "1", "0.5", "1", "1", "1", "0", "1", "1", "0.5", "1", "1", "0.5", "1", "1", "0.5", "0.5", "0.5", "0.5", "0.5", "0.5", "1", 
                                           "1", "0", "0.5", "0.5", "1", "1", "1", "1", "0.5", "1", "1", "0.5", "1", "0.5", "1", "0.5", "1", "0.5", "0", "0.5", "1", "0.5", 
                                           "1", "1", "1", "1", "1", "1", "0.5", "0.5", "1", "1", "0.5", "0.5", "0.5", "0", "1", "1", "0.5", "1", "1", "0.5", "1", "1", 
                                           "1", "0.5", "0.5", "1", "1", "1", "1", "1", "0", "0.5", "0.5", "0.5", "1", "1", "0.5", "1", "1", "1", "0.5", "0.5", "1", "1", 
                                           "0.5", "1", "0.5", "0", "0.5", "1", "1", "1", "1", "1", "0.5", "1", "0.5", "1", "0.5", "1", "1", "0.5", "0.5", "0.5", "0", "1", 
                                           "1", "0.5", "1", "1", "1", "0.5", "1", "0.5", "1", "1", "1", "1", "0.5", "1", "1", "0", "0.5", "0.5", "0.5", "1", "1", "1", 
                                           "0.5", "0.5", "1", "0.5", "1", "1", "1", "1", "1", "0.5", "0", "0.5", "1", "0.5", "1", "1", "1", "0.5", "0.5", "1", "1", "0.5", 
                                           "1", "1", "0.5", "0.5", "0.5", "0", "1", "1", "1", "0.5", "1", "1", "0.5", "1", "1", "1", "0.5", "1", "1", "0.5", "1", "1", 
                                           "0", "0.5", "0.5", "1", "0.5", "1", "0.5", "0.5", "1", "1", "1", "1", "1", "1", "0.5", "1", "0.5", "0", "0.5", "1", "1", "0.5", 
                                           "0.5", "1", "0.5", "1", "1", "0.5", "1", "1", "1", "1", "0.5", "0.5", "0")), nrow=16, ncol=16, byrow=TRUE)
              
             
              
              ## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
              Dt1t2<-NULL
              
              ## Matrix of all possible initial states
              tvec.mat<-diag(16)
              
        
                  for (p in 1:nrow(pedigree))
                  {
                   
                    svt0<-NULL
                    svt1<-list()
                    svt2<-list()
                    
                    svt0<-t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
                    
                        for (i in 1:nrow(tvec.mat))
                        {
                          
                            svt1[i]<-list(t(tvec.mat[1,])%*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1])))
                            svt2[i]<-list(t(tvec.mat[1,]) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1])))
                        }
                    
                   
                        
                    DivProbt1t2<-NULL
                    
                      for (j in 1:nrow(tvec.mat))
                      {
        
                          t1in<-svt1[[j]]
                          t2in<-svt2[[j]]
                            
                          jointPROB     <- expand.grid(a = t1in, b = t2in)
                          jointPROB     <- jointPROB[,1] * jointPROB[,2] 
                          jointPROBt1t2 <- matrix(jointPROB, 16, 16, byrow=F)
                          DivProbt1t2[j]   <- svt0[j]*sum(Deffects * jointPROBt1t2)
                      }
                      
                      
                
                ## Total (weighted) divergence 
                Dt1t2[p]<- sum(DivProbt1t2)
                
                
              }
              
              divout<-list(Dt1t2)
              
              return(divout)
              
            }
            
            
            ###### Defining the Least Square function to be minimized
            ###### Note the equilibrium constraint, which can be made as small as desired.
            
            LSE_intercept<-function(param_int) 
            {
              sum((pedigree[,4] - param_int[2] - divergence(pedigree, p0aa, p0cc, p0tt, p0gg, param_int[1])[[1]])^2)	
            }
            
            
            
            ###### Calculating the initial proportions 
            p0aa <- p0aa
            p0cc <- p0cc
            p0tt <- p0tt
            p0gg <- p0gg
            p0ac <- 0
            p0at <- 0
            p0ag <- 0
            p0ca <- 0
            p0ct <- 0
            p0cg <- 0
            p0ta <- 0
            p0tc <- 0
            p0tg <- 0
            p0ga <- 0
            p0gc <- 0
            p0gt <- 0
            
            if(sum(c(p0aa, p0cc, p0tt, p0gg), na.rm =T) != 1) 
            {stop("The initial state probabilities don't sum to 1")}
            
          
            ##### Initializing
            optim.method<-"Nelder-Mead"
            final<-NULL
            counter<-0
            opt.out<-NULL
            pedigree<-pedigree.data
            
            
            for (s in 1:Nstarts)
            {
              
              ## Draw random starting values
              g.start  <-10^(runif(1, log10(10^-11), log10(10^-3)))
              intercept.start <-runif(1,0,max(pedigree[,4]))
              param_int0 = c(g.start, intercept.start)
              
              ## Initializing
              counter<-counter+1
              
              cat("Progress: ", counter/Nstarts, "\n")
              
              
              opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
              opt.out <-cbind(opt.out, g.start, intercept.start)
              final<-rbind(final, opt.out)
              
              
            } # End of Nstarts loop
            
            colnames(final)[1:2]<-c("gamma", "intercept")
            
            
            
            final<-final[order(final[,"value"]),]
            
            if (allow.neg.intercept == "yes") 
            { index.1<-which(final["gamma"] > 0 & final["convcode"] == 0)}
            
            if (allow.neg.intercept == "no")
            {index.1<-which(final["gamma"] > 0 & final["intercept"] > 0 & final["convcode"] == 0)}
            
      
            index.2<-setdiff(1:nrow(final), index.1)
            final.1<-final[index.1,]
            final.2<-final[index.2,]
            
            
          
    
            
            ##### Calculting the predicted values based on the 'best' model (i.e. that with the lowest least square)
            PrAA <- p0aa; PrCC <- p0cc; PrTT <- p0tt; PrGG <- p0gg; PrAC <- p0ac; PrAT <- p0at; PrAG <- p0ag; PrCA <- p0ca; PrCT <- p0ct
            PrCG <- p0cg; PrTA <- p0ta; PrTC <- p0tc; PrTG <- p0tg; PrGA <- p0ga; PrGC <- p0gc; PrGT <- p0gt
            
            g  <- final.1[1, "gamma"]
            intercept<-final.1[1,"intercept"]
            
            
           
            svGzero    <- c(PrAA,
                            PrCC,
                            PrTT,
                            PrGG,
                            PrAC,
                            PrAT,
                            PrAG,
                            PrCA,
                            PrCT,
                            PrCG,
                            PrTA,
                            PrTC,
                            PrTG,
                            PrGA,
                            PrGC,
                            PrGT)
            
            ## Defining the generation (or transition) matrix for the mitotic case
            Tmat<-matrix(NA, nrow=16, ncol=16)
            
            Ta<-matrix(c(
              (1-g)^2, 1/9*g^2, 1/9*g^2, 1/9*g^2, 
              1/9*g^2, (1-g)^2, 1/9*g^2, 1/9*g^2, 
              1/9*g^2, 1/9*g^2, (1-g)^2, 1/9*g^2, 
              1/9*g^2, 1/9*g^2, 1/9*g^2, (1-g)^2), nrow=4, byrow=TRUE)
            
            
            #              1                2            3          4            5           6              7           8            9            10           11           12
            Tb<-matrix(c(
              1/3*(1-g)*g,   1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,     1/9*g^2, 
              1/3*(1-g)*g,   1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,
              1/9*g^2,       1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,     1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g,
              1/9*g^2,       1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g), 
              nrow=4, byrow=TRUE)
            
            
            Tc<-matrix(c(
              1/3*(1-g)*g,  1/3*(1-g)*g, 1/9*g^2,     1/9*g^2, 
              1/3*(1-g)*g,  1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,
              1/3*(1-g)*g,  1/9*g^2,     1/9*g^2,     1/3*(1-g)*g,
              1/3*(1-g)*g,  1/3*(1-g)*g, 1/9*g^2,     1/9*g^2, 
              1/9*g^2,      1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*g^2, 
              1/9*g^2,      1/3*(1-g)*g, 1/9*g^2,     1/3*(1-g)*g,
              1/3*(1-g)*g,  1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,
              1/9*g^2,      1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*g^2,
              1/9*g^2,      1/9*g^2,     1/3*(1-g)*g, 1/3*(1-g)*g,
              1/3*(1-g)*g,  1/9*g^2,     1/9*g^2,     1/3*(1-g)*g,
              1/9*g^2,      1/3*(1-g)*g, 1/9*g^2,     1/3*(1-g)*g,
              1/9*g^2,      1/9*g^2,     1/3*(1-g)*g, 1/3*(1-g)*g), nrow=12, byrow=TRUE)
            
            #              1                2            3           4            5           6            7            8           9            10           11           12
            Td<-matrix(c(
              (1-g)^2,     1/3*(1-g)*g,  1/3*(1-g)*g,  1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,
              1/3*(1-g)*g, (1-g)^2,      1/3*(1-g)*g,  1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g,
              1/3*(1-g)*g, 1/3*(1-g)*g,  (1-g)^2,      1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,
              1/9*(g)^2,   1/9*(g)^2,    1/9*(g)^2,    (1-g)^2,     1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,
              1/9*(g)^2,   1/3*(1-g)*g,  1/9*(g)^2,    1/3*(1-g)*g, (1-g)^2,     1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g,
              1/9*(g)^2,   1/9*(g)^2,    1/3*(1-g)*g,  1/3*(1-g)*g, 1/3*(1-g)*g, (1-g)^2,     1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,
              1/9*(g)^2,   1/9*(g)^2,    1/9*(g)^2,    1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   (1-g)^2,     1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,
              1/3*(1-g)*g, 1/9*(g)^2,    1/9*(g)^2,    1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, (1-g)^2,     1/3*(1-g)*g, 1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,
              1/9*(g)^2,   1/9*(g)^2,    1/3*(1-g)*g,  1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, (1-g)^2,     1/9*(g)^2,   1/9*(g)^2, 1/9*(g)^2,
              1/9*(g)^2,   1/9*(g)^2,    1/9*(g)^2,    1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   (1-g)^2,     1/3*(1-g)*g, 1/3*(1-g)*g,
              1/3*(1-g)*g, 1/9*(g)^2,    1/9*(g)^2,    1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g,  1/9*(g)^2,  1/3*(1-g)*g, (1-g)^2,     1/3*(1-g)*g, 
              1/9*(g)^2,   1/3*(1-g)*g,  1/9*(g)^2,    1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/3*(1-g)*g, (1-g)^2), nrow=12, byrow=TRUE)
            
            
            ## Transition matrix
            Tmat[1:4, 1:4]<-Ta
            Tmat[1:4, 5:16]<-Tb
            Tmat[5:16, 1:4]<-Tc
            Tmat[5:16, 5:16]<-Td
            Genmatrix <- Tmat
            
            ## Coding the divergence between samples i and j
            Deffects<-matrix(as.numeric(c("0", "1", "1", "1", "0.5", "0.5", "0.5", "0.5", "1", "1", "0.5", "1", "1", "0.5", "1", "1", "1", "0", "1", "1", "0.5", "1", "1", 
                                          "0.5", "0.5", "0.5", "1", "0.5", "1", "1", "0.5", "1", "1", "1", "0", "1", "1", "0.5", "1", "1", "0.5", "1", "0.5", "0.5", "0.5", 
                                          "1", "1", "0.5", "1", "1", "1", "0", "1", "1", "0.5", "1", "1", "0.5", "1", "1", "0.5", "0.5", "0.5", "0.5", "0.5", "0.5", "1", 
                                          "1", "0", "0.5", "0.5", "1", "1", "1", "1", "0.5", "1", "1", "0.5", "1", "0.5", "1", "0.5", "1", "0.5", "0", "0.5", "1", "0.5", 
                                          "1", "1", "1", "1", "1", "1", "0.5", "0.5", "1", "1", "0.5", "0.5", "0.5", "0", "1", "1", "0.5", "1", "1", "0.5", "1", "1", 
                                          "1", "0.5", "0.5", "1", "1", "1", "1", "1", "0", "0.5", "0.5", "0.5", "1", "1", "0.5", "1", "1", "1", "0.5", "0.5", "1", "1", 
                                          "0.5", "1", "0.5", "0", "0.5", "1", "1", "1", "1", "1", "0.5", "1", "0.5", "1", "0.5", "1", "1", "0.5", "0.5", "0.5", "0", "1", 
                                          "1", "0.5", "1", "1", "1", "0.5", "1", "0.5", "1", "1", "1", "1", "0.5", "1", "1", "0", "0.5", "0.5", "0.5", "1", "1", "1", 
                                          "0.5", "0.5", "1", "0.5", "1", "1", "1", "1", "1", "0.5", "0", "0.5", "1", "0.5", "1", "1", "1", "0.5", "0.5", "1", "1", "0.5", 
                                          "1", "1", "0.5", "0.5", "0.5", "0", "1", "1", "1", "0.5", "1", "1", "0.5", "1", "1", "1", "0.5", "1", "1", "0.5", "1", "1", 
                                          "0", "0.5", "0.5", "1", "0.5", "1", "0.5", "0.5", "1", "1", "1", "1", "1", "1", "0.5", "1", "0.5", "0", "0.5", "1", "1", "0.5", 
                                          "0.5", "1", "0.5", "1", "1", "0.5", "1", "1", "1", "1", "0.5", "0.5", "0")), nrow=16, ncol=16, byrow=TRUE)
            
            
            
            ## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
            Dt1t2<-NULL
            
            ## Matrix of all possible initial states
            tvec.mat<-diag(16)
            
            
            for (p in 1:nrow(pedigree))
            {
              
              svt0<-NULL
              svt1<-list()
              svt2<-list()
              
              svt0<-t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
              
              for (i in 1:nrow(tvec.mat))
              {
                
                svt1[i]<-list(t(tvec.mat[1,])%*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1])))
                svt2[i]<-list(t(tvec.mat[1,]) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1])))
              }
              
              
              
              DivProbt1t2<-NULL
              
              for (j in 1:nrow(tvec.mat))
              {
                
                t1in<-svt1[[j]]
                t2in<-svt2[[j]]
                
                jointPROB     <- expand.grid(a = t1in, b = t2in)
                jointPROB     <- jointPROB[,1] * jointPROB[,2] 
                jointPROBt1t2 <- matrix(jointPROB, 16, 16, byrow=F)
                DivProbt1t2[j]   <- svt0[j]*sum(Deffects * jointPROBt1t2)
              }
              
              
              
              ## Total (weighted) divergence 
              Dt1t2[p]<- sum(DivProbt1t2)
              
              
            }
           
        
            ## Calculating the least square part
            Residual<-(pedigree[,4] - intercept - Dt1t2)
            
            ##### Augmenting pedigree
            delta.t<-pedigree[,2] + pedigree[,3] - 2*pedigree[,1]
            #pedigree<-cbind(pedigree,delta.t)
            pedigree<-cbind(pedigree, delta.t, Dt1t2 + intercept, Residual)
            colnames(pedigree)[c(4,5,6,7)]<-c("div.obs", "delta.t","div.pred", "residual")
            
            
            ##### Making info about settings
            info<-c("p0aa", "p0cc", "p0tt", "p0gg", "Nstarts", "optim.method")
            info2<-c(p0aa, p0cc, p0tt, p0gg, Nstarts, optim.method)
            info.out<-data.frame(info, info2)
            colnames(info.out)<-c("Para", "Setting")
            
            
            
            
            
            ###### Generating theoretical fit 
            
            ## Reading in pedigree
            obs<-pedigree[,"div.obs"]
            dtime<-pedigree[,"delta.t"]
            
            ## Reading in parameter estimates
            est <-final.1
            g <-as.numeric(est[1,1])
            intercept<-as.numeric(est[1,2])
            
            ## Reading initial state vector
            settings<-info.out
            PrAA<-p0aa<-as.numeric(as.character(settings[1,2]))
            PrCC<-p0cc<-as.numeric(as.character(settings[2,2]))
            PrTT<-p0tt<-as.numeric(as.character(settings[3,2]))
            PrGG<-p0gg<-as.numeric(as.character(settings[4,2]))
            PrAC <- PrAT <- PrAG <-PrCA <-PrCT <- PrCG <- PrTA <- PrTC <- PrTG <- PrGA <-PrGC <- PrGT <- 0
            time1<- seq(1,max(c(pedigree[,2], pedigree[,3])))
            time2<- seq(1,max(c(pedigree[,2], pedigree[,3])))
            time.out<-expand.grid(time1,time2)
            #time0<- rep(min(pedigree[,1]), nrow(time.out))
            time0<- rep(0, nrow(time.out))
            pedigree.new<-as.matrix(cbind(time0,time.out))
            pedigree.new<-cbind(pedigree.new, c(pedigree.new[,2] + pedigree.new[,3] - 2*pedigree.new[,1]))
            pedigree.new<-pedigree.new[!duplicated(pedigree.new[,4]), ]
            pedigree.new<-pedigree.new[,1:3]
           
            
            ## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
            Dt1t2<-NULL
            
            ## Matrix of all possible initial states
            tvec.mat<-diag(16)
            
            
            for (p in 1:nrow(pedigree.new))
            {
              
              svt0<-NULL
              svt1<-list()
              svt2<-list()
              
              svt0<-t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,1]))
              
              for (i in 1:nrow(tvec.mat))
              {
                
                svt1[i]<-list(t(tvec.mat[1,])%*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1])))
                svt2[i]<-list(t(tvec.mat[1,]) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1])))
              }
              
              
              
              DivProbt1t2<-NULL
              
              for (j in 1:nrow(tvec.mat))
              {
                
                t1in<-svt1[[j]]
                t2in<-svt2[[j]]
                
                jointPROB     <- expand.grid(a = t1in, b = t2in)
                jointPROB     <- jointPROB[,1] * jointPROB[,2] 
                jointPROBt1t2 <- matrix(jointPROB, 16, 16, byrow=F)
                DivProbt1t2[j]   <- svt0[j]*sum(Deffects * jointPROBt1t2)
              }
              
              
              
              ## Total (weighted) divergence 
              Dt1t2[p]<- sum(DivProbt1t2)
              
              
            }
            
           
            pedigree.new<-cbind(pedigree.new, Dt1t2+intercept, c(pedigree.new[,2] + pedigree.new[,3] - 2*pedigree.new[,1]))
            colnames(pedigree.new)<-c("time0", "time1", "time2", "div.sim", "delta.t")
            pedigree.new<-pedigree.new[order(pedigree.new[,5]),]
            
            
            model<-"mutSOMA.R"
            
            abfree.out<-list(final.1, final.2, pedigree, info.out, model, pedigree.new)
            names(abfree.out)<-c("estimates", "estimates.flagged", "pedigree", "settings", "model", "for.fit.plot")
            
            
            
            ## Ouputting result datasets
            dput(abfree.out, paste(out.dir, out.name, ".Rdata", sep=""))
            return(abfree.out)
            
            
      } #End of function
    
  
   
  
   
    