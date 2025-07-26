P_rw = matrix(0, 8, 8)
P_imh = matrix(0, 8, 8)
post = exp( c(0, 63.9847, -2.7647,  207.669,90.4603, 88.6873,148.945, 204.904 ))

post = post/sum(post)

for (i in 1:8){
  binary_vector <- 1* as.logical(intToBits(i-1))[1:3]
  print(binary_vector)
}


for (i in 1:8){
  binary_vector <- 1* as.logical(intToBits(i-1))[1:3]
  # integer_value <- binary_vector[1]*1 + binary_vector[2]*2 + binary_vector[3]*4 + 1
  for (j in 1:3){
    next_binary = binary_vector
    next_binary[j] = 1 - binary_vector[j]
    next_integer_value = next_binary[1]*1 + next_binary[2]*2 + next_binary[3]*4 + 1
    P_rw[i, next_integer_value] = min(1, post[next_integer_value]/post[i]) /3
  }
  P_rw[i,i] = 1 - sum(P_rw[i,])
  
}
min(svd(P_rw)$d)




for (i in 1:8){
  binary_vector <- 1* as.logical(intToBits(i-1))[1:3]
  # integer_value <- binary_vector[1]*1 + binary_vector[2]*2 + binary_vector[3]*4 + 1
  current_post_value = post[i]
  next_post_value = numeric(3)
  for(j in 1:3){
    next_binary = binary_vector
    next_binary[j] = 1 - binary_vector[j]
    next_integer_value = next_binary[1]*1 + next_binary[2]*2 + next_binary[3]*4 + 1
    next_post_value[j] = min( max(post[next_integer_value] / post[i], 3), 9)
  }
  # c(pi(100), pi(010), pi(001))
  
  
  
  for (j in 1:3){
    next_binary = binary_vector
    next_binary[j] = 1 - binary_vector[j]
    next_integer_value = next_binary[1]*1 + next_binary[2]*2 + next_binary[3]*4 + 1
    # pi(100)
    next_next_post_value = numeric(3)
    #pi(000, 110, 101)
    for (k in 1:3){
      next_next_binary = next_binary
      next_next_binary[k] = 1 - next_binary[k]
      next_next_integer_value = next_next_binary[1]*1 + next_next_binary[2]*2 + next_next_binary[3]*4 + 1
      next_next_post_value[k] = min( max(post[next_next_integer_value] / post[next_integer_value], 3), 9)
    }
    
    P_imh[i, next_integer_value] = 
      next_post_value[j]/sum(next_post_value) * # Kernel value
      
      min(1, 
          post[next_integer_value]/post[i] * 
            
            sum(next_post_value) / sum(next_next_post_value)  * min( max(current_post_value / post[next_integer_value], 3), 9) / next_post_value[j] 
          
          )   
              
  }
  P_imh[i,i] = 1 - sum(P_imh[i,])
}



mu = c(1,rep(0,7)) %*% P_rw
nu = c(1,rep(0,7)) %*% P_imh

for (i in 1:100){
  mu = mu %*% P_rw
  nu = nu %*% P_imh
}

round(mu, 6)
round(post, 6)
round(nu, 6)
eigen(P_rw)$values
eigen(P_imh)$values
1.00000000-0.41770491

#data <- read.table("results/snr_0_method_single_init_mid_clip_0.txt", header = FALSE, row.names = NULL)  # Adjust file name and path as needed
data <- read.table("results/informed_hit_time_init_mid.txt", header = TRUE)  # Adjust file name and path as needed
hit_time <- data[, 2]
data <- read.table("results/informed_hit_iter_init_mid.txt", header = TRUE)  # Adjust file name and path as needed
hit_iter <- data[, 2]
data <- read.table("results/informed_wall_time_init_mid.txt", header = TRUE)  # Adjust file name and path as needed
wall_time <- data[, 2]
# Print the vector
print(median(hit_time))
print(median(hit_iter))
print(median(wall_time))

data <- read.table("results/informed_hit_time_init_bad.txt", header = TRUE)  # Adjust file name and path as needed
hit_time <- data[, 2]
data <- read.table("results/informed_hit_iter_init_bad.txt", header = TRUE)  # Adjust file name and path as needed
hit_iter <- data[, 2]
data <- read.table("results/informed_wall_time_init_bad.txt", header = TRUE)  # Adjust file name and path as needed
wall_time <- data[, 2]
# Print the vector
print(median(hit_time))
print(median(hit_iter))
print(median(wall_time[1:80]))
