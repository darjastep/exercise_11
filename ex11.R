rm(list=ls())
setwd("V:/MPA-PRG/exercise_11")
load('HMM1.Rdata')
load('HMM2.Rdata')
load('HMM3.Rdata')

Viterbi <- function(obs, HMM) {    
        N <- HMM$N                 # hidden states.. "H" "L"
        M <- HMM$M                 # emission characters.. "A" "C" "G" "T"
        A <- HMM$A                 # transition matrix
        B <- HMM$B                 # emission matrix
        pi <- HMM$pi               # initial probabilities
        
        
        print("DEBUG:")
        print(obs)
        print(class(obs))
        T_len <- length(obs)  # number of observations
        K <- length(N)    #number of states
        
        #matice
        #delta = matice s nejvetsimi pravdepodobnostmi
        delta <- matrix(-Inf, nrow = K, ncol = T_len)
        
        #psi = pamatuje si trasu
        psi   <- matrix(0, nrow = K, ncol = T_len)
        
        
        for (i in 1:K) {
          obs_index <- match(as.character(obs[1]), M)
          delta[i, 1] <- pi[i] + B[i, obs_index]   # log(a*b) = log(a)+log(b)
          psi[i, 1] <- 0
        }
        
        for (time in 2:T_len) { # pro kazdy "čas"/observation
          obs_index <- match(as.character(obs[time]), M)
          
          for (j in 1:K) { #pro kazdy stav
            trans_probs <- delta[, time-1] + A[, j]   
            psi[j, time] <- which.max(trans_probs)
            delta[j, time] <- max(trans_probs) + B[j, obs_index]
          }
        }
          states <- numeric(T_len)
          states[T_len] <- which.max(delta[, T_len])
          
         
          for (time in (T_len-1):1) {
            states[time] <- psi[states[time+1], time+1]
          }
          
          # return results
          list(
            delta = delta,
            states = N[states]
          )
            
            
        
        
}

library(Biostrings)
s <- AAString("TGA")

Viterbi(s, HMM1)
