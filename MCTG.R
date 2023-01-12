# NAME CHE BAIHUI
# MATRIC NO. G2202211L

setwd("D:/mydoc/NTU Analytics MSc/MH6311 Stochastic Processes for Data Science")

library(qdapDictionaries)

# Read the text
text <- readChar("example_text.txt",nchars=20000)

# Break the text file into characters
chars <- unlist (strsplit(gsub ("[^a-z]", "_", tolower(text)), ""))

# 10% of the text
chars_short <- chars[1:as.integer(round(length(chars)*0.1,0))]

# 50% of the text
chars_half <- chars[1:as.integer(round(length(chars)*0.5,0))]

# 100% of the text
chars_full <- chars

# Display lengths of characters resp.
length(chars_short);length(chars_half); length(chars_full)


######################################################
# Function to generate a transition matrix
######################################################

kth_order_matrix <- function(transition_order, characters) {
  # Extract unique characters set
  uniq_char <- unique(characters)
  # transition order
  trans_order <- transition_order
  # prepare for the existing state vector
  exist_state <- rep(NA, length(characters)-trans_order)
  # put elements in the states bucket
  for (i in 1:(length(characters)-trans_order)) {
    exist_state[i] <- paste(characters[i:(i+trans_order-1)], collapse='')
  }
  # Create a zero matrix P ready for a transition matrix
  uniq_str <- unique(exist_state)
  P <- matrix(nrow=length(uniq_str), ncol=length(uniq_char), 0, 
              dimnames=list(uniq_str, uniq_char))
  # put (k+1)-gram elements in another states bucket
  exist_state2 <- rep(NA, length(characters)-trans_order-1)
  for (i in 1:(length(characters)-trans_order-1)) {
    exist_state2[i] <- paste(characters[i:(i+trans_order)], collapse='')
  }
  # Summarize the frequency of an k-gram element followed by one letter/symbol
  for (i in 1:(length(exist_state2))) {
    subchar <- unlist(strsplit(exist_state2[i],""))
    matrow <- paste(subchar[1:trans_order], collapse='')
    matcol <- subchar[trans_order+1]
    P[matrow, matcol] <- P[matrow, matcol] + 1
  }
  # Compute the relative frequency for each row of the matrix
  for (i in 1:length(uniq_str)) {
    P[i,] <- P[i,] / sum(P[i,])
  }
  
  # Generate the initial state
  pi <- table(exist_state) / length(exist_state)
  state0 <- sample(rownames(P), size=1, prob=pi)
  # Return transition matrix, probability distribution, initial state
  rtn_lst <- list("tran_mat" = P, "prob_dtb" = pi, "init_state" = state0)
  return(rtn_lst)
}

#####################################################################
# Ouput word accuracy (implemented in the text generation function)
#####################################################################

# Set up a function to exclude "" in the characters vector
my_is.null <- function(char_vec) {
  tf_box <- rep(NA, length(char_vec))
  for (i in 1:length(char_vec)) {
    if (char_vec[i] == "") {
      tf_box[i] <- TRUE
    } else {tf_box[i] <- FALSE}
  }
  return(tf_box)
}

# Set up a function to compute the word accuracy
# Package "qdapDictionaries" needed
word_acc_func <- function(char_vec) {
  is.word_box <- rep(NA, length(char_vec))
  for (i in 1:length(char_vec)) {
    if (char_vec[i] %in% GradyAugmented) {
      is.word_box[i] <- TRUE
    } else {is.word_box[i] <- FALSE}
  }
  # Total real words count
  word_count <- sum(is.word_box)
  # Total letters combinations count
  total_count <- length(is.word_box)
  # Real word ratio
  word_acc_ratio <- word_count  / total_count
  
  rtn_lst <- list("word_count"=word_count, "total_count"=total_count,
                  "word_accuracy"=word_acc_ratio)
  return(rtn_lst)
}

# Function to output word accuracy of a generated text
txt_word_acc <- function(text_generated) {
  # Separate text by "_"
  txt_prc <- unlist(strsplit(text_generated, "_"))
  # Remove "" and keep real words
  txt_prc2 <- txt_prc[!my_is.null(txt_prc)]
  # Compute the word accuracy, total letters combinations, total real words
  out_acc_group <- word_acc_func(txt_prc2)
  out_acc <- out_acc_group$word_accuracy
  total_words <- out_acc_group$total_count
  total_real_words <- out_acc_group$word_count
  
  acc_rtn_lst <- list("word_count"=total_real_words, "total_count"=total_words,
                      "word_accuracy"=out_acc)
  return(acc_rtn_lst)
}

######################################################
# Function to generate a text (simple version)
######################################################

kth_order_txt_simple <- function(transition_order, characters, 
                                 text_length=length(characters)) {
  
  # Start timing of the whole process
  total_start_time <- Sys.time()
  
  # Extract unique characters set
  uniq_char <- unique(characters)
  
  # Summarize transition matrix, probability distribution
  # Collect the first transition matrix making time
  mat_start_time <- Sys.time()
  kth_order_mat <- kth_order_matrix(transition_order, characters)
  mat_end_time <- Sys.time()
  mat_con_time <- mat_end_time - mat_start_time
  P <- kth_order_mat$tran_mat
  
  # If the state exists but it is the end of the text,
  # then the relative frequency would be NaN.
  
  # Generally it always happens when order > 1 if it happens when order = 2.
  # It could happen as well when order = 1,
  # e.g. the only "." is at the end of the text when not replaced with "_".
  # Remove this row.
  # In this step removing this row does not affect the existing state
  # probability distribution,
  # and makes the next steps for dealing with non-existent state in 
  # the matrix more conveniently.
  if (anyNA(P[nrow(P),])) {
    P <- P[-c(nrow(P)),]
  }
  
  # Generate an initial state
  state0 <- kth_order_mat$init_state
  
  # Prepare for a text container
  txt_gen <- rep(NA, text_length)
  
  # Create a container to collect the time generating per character
  time_bucket <- rep(0, text_length)
  
  # Generate the first k characters in the container
  for (i in 1:transition_order) {
    start_time <- Sys.time()
    txt_gen[i] <- unlist(strsplit(state0, ""))[i]
    end_time <- Sys.time()
    time_bucket[i] <- end_time - start_time    
  }
  
  # Generate characters in the rest space
  for (i in 1:(text_length-transition_order)) {
    
    # Start timing
    start_time <- Sys.time()    
    
    # Get the current state ready
    curr_state <- paste(txt_gen[i:(i+transition_order-1)], collapse='')

    # When the k-gram state does not exist in the transition matrix row indices,
    # we just generate the character by the initial distribution
    con1 <- curr_state %in% rownames(P) == FALSE
    
    if (con1) {
      pi <- kth_order_matrix(1, characters)$prob_dtb
      txt_gen[i+transition_order] <- sample(uniq_char, size=1, prob=pi)
    } else { # otherwise generate the character by the transition matrix
      txt_gen[i+transition_order] <- sample(uniq_char, size=1, prob=P[curr_state,])
    }
    # End timing    
    end_time <- Sys.time()
    # Put the time generating a character into the container
    time_bucket[i+transition_order] <- end_time - start_time
  }
  
  # Connect characters to a text
  txt_gen <- paste(txt_gen, collapse='')
  
  # End timing of the whole process
  total_end_time <- Sys.time()
  
  # Compute timing of the whole process
  total_time <- total_end_time - total_start_time
  
  # Compute the word accuracy of the text generated
  # Summarize the total real words generated
  # Summarize the total letters combinations generated
  text_gen_word_acc <- txt_word_acc(txt_gen)
  gen_word_accuracy <- text_gen_word_acc$word_accuracy
  gen_real_word_count <- text_gen_word_acc$word_count
  gen_word_count <- text_gen_word_acc$total_count
  
  # Return text generated, word accuracy of the generated text, 
  # first transition matrix making time,
  # generating time per character, total processing time
  gnt_rtn_lst <- list("text_gen"= txt_gen, "text_gen_word_acc"=gen_word_accuracy, 
                      "combined_letters_count"=gen_word_count, "real_word_count"=gen_real_word_count,
                      "mat_make_time"=mat_con_time, "time_per_char"= time_bucket, 
                      "total_time"=total_time)
  return(gnt_rtn_lst)
  
}



######################################################
# Function to generate a text (less simple version)
######################################################

kth_order_txt <- function(transition_order, characters, 
                          text_length=length(characters)) {
  
  # Start timing of the whole process
  total_start_time <- Sys.time()
  
  # Extract unique characters set
  uniq_char <- unique(characters)
  
  # Summarize transition matrix, probability distribution
  # Collect the first transition matrix making time
  mat_start_time <- Sys.time()
  kth_order_mat <- kth_order_matrix(transition_order, characters)
  mat_end_time <- Sys.time()
  mat_con_time <- mat_end_time - mat_start_time
  P <- kth_order_mat$tran_mat
  
  # If the state exists but it is the end of the text,
  # then the relative frequency would be NaN.
  
  # Generally it always happens when order > 1 if it happens when order = 2.
  # It could happen as well when order = 1,
  # e.g. the only "." is at the end of the text when not replaced with "_".
  # Remove this row.
  # In this step removing this row does not affect the existing state
  # probability distribution,
  # and makes the next steps for dealing with non-existent state in 
  # the matrix more conveniently.
  if (anyNA(P[nrow(P),])) {
    P <- P[-c(nrow(P)),]
  }
  
  # Generate an initial state
  state0 <- kth_order_mat$init_state
  
  # Prepare for a text container
  txt_gen <- rep(NA, text_length)
  
  # Create a container to collect the time generating per character
  time_bucket <- rep(0, text_length)
  
  # Create a container to collect loop times q
  loop_times <- rep(NA, text_length)
  for (i in 1:transition_order) {
    loop_times[i] <- 0
  }
  
  # Generate the first k characters in the container
  for (i in 1:transition_order) {
    start_time <- Sys.time()
    txt_gen[i] <- unlist(strsplit(state0, ""))[i]
    end_time <- Sys.time()
    time_bucket[i] <- end_time - start_time
  }
  
  # Generate characters in the rest space
  for (i in 1:(text_length-transition_order)) {
    # Start timing
    start_time <- Sys.time()
    
    # Get the current state ready
    curr_state <- paste(txt_gen[i:(i+transition_order-1)], collapse='')
    
    # When the k-gram state does not exist in the transition matrix P rows
    con1 <- curr_state %in% rownames(P) == FALSE
    
    # Set the variable to collect the while loop time below per generation
    q <- 0
    
    # If the first condition is satisfied ...
    if (con1) {
      
      # Assign transition order, transition matrix to new variables
      k <- transition_order
      new_P <- P
      
      # Set the 2nd condition
      # Judging whether the current state is in the 
      # transition matrix new_P rows
      con2 <- curr_state %in% rownames(new_P) == FALSE
      
      # Do the loop when the second condition and the transition order > 1
      # When the current state is not in the transition matrix new_P,
      # then remove the first letter/symbol of the current state
      # to make a new current state,
      # and use the transition matrix with one order lower,
      # and see if the new current state is in the new transition matrix,
      # if not, do the same as above recursively,
      # until transition order = 1
      while ((con2) & (k>1)) {
        # print("Meeting transition order reduction")
        k <- k - 1
        q <- q + 1
        new_mat <- kth_order_matrix(k, characters)
        new_P <- new_mat$P
        if (anyNA(new_P[nrow(new_P),])) {new_P <- new_P[-c(nrow(new_P)),]}
        # Create the new current state by removing the first character of the
        # previous current state
        curr_state <- paste(unlist(strsplit(curr_state,""))[-1], collapse = '')
      }
      # When the transition order is one,
      # jumps out of the loop above and do the if statement
      # If the transition order is one, but the current state is
      # still not in the transition matrix,
      # then use the initial distribution to randomly generate
      # a new character
      # Generally this does not happen based on the design of
      # our transition matrix ...
      if ((k == 1) & (con2)){
        pi <- kth_order_matrix(1, characters)$prob_dtb
        txt_gen[i+transition_order] <- sample(uniq_char, size=1, prob=pi)
      } else { # If not, then use the transition matrix with a specific order to generate a character
        txt_gen[i+transition_order] <- sample(uniq_char, size=1, prob=new_P[curr_state,])
      }
    } else { # This is the choice when the first condition is not satisfied
      txt_gen[i+transition_order] <- sample(uniq_char, size=1, prob=P[curr_state,])
    }
    # End timing
    end_time <- Sys.time()
    # Put the time generating a character into the container
    time_bucket[i+transition_order] <- end_time - start_time
    # Put the loop times generating a character into the container
    loop_times[i+transition_order] <- q
  }
  
  # Connect characters to a text
  txt_gen <- paste(txt_gen, collapse='')
  
  # End timing of the whole process
  total_end_time <- Sys.time()
  
  # Compute timing of the whole process
  total_time <- total_end_time - total_start_time
  
  # Compute the word accuracy of the text generated
  # Summarize the total real words generated
  # Summarize the total letters combinations generated
  text_gen_word_acc <- txt_word_acc(txt_gen)
  gen_word_accuracy <- text_gen_word_acc$word_accuracy
  gen_real_word_count <- text_gen_word_acc$word_count
  gen_word_count <- text_gen_word_acc$total_count
  
  
  # Return text generated, word accuracy of the generated text, 
  # first transition matrix making time,
  # generating time per character, looping time per generation,
  # total processing time
  gnt_rtn_lst <- list("text_gen"= txt_gen, "text_gen_word_acc"=gen_word_accuracy, 
                      "combined_letters_count"=gen_word_count, "real_word_count"=gen_real_word_count, 
                      "mat_make_time"=mat_con_time, "time_per_char"= time_bucket, 
                      "loop_times"=loop_times, "total_time"=total_time)
  
  return(gnt_rtn_lst)
}



################################################
# Use kth_order_txt_simple to generate a text
################################################

mytxt <- kth_order_txt_simple(4, chars_full, 1200)

# text generated
text_gen <- mytxt$text_gen

# Word accuracy of the generated text
wrd_acc <- mytxt$text_gen_word_acc

# Total letters combination generated
total_words <- mytxt$combined_letters_count

# Total real words generated
total_real_words <- mytxt$real_word_count

# First transition Matrix Making time
time_matrix <- mytxt$mat_make_time

# Generating time per character
time_per_char_gen <- mytxt$time_per_char

# Total generating time of the whole process
total_time <- mytxt$total_time

# Plot the characters generating time series
x_char <- c(1:length(time_per_char_gen))
plot(x_char, time_per_char_gen, type="l", lwd=1.6, col="blue", xlab="Generating Order", 
     ylab="Generating Time (sec)", main="Generating Time per Character",
     cex.main=1.5, cex.lab=1.36); grid()


##########################################
# Use kth_order_txt to generate a text
##########################################

mytxt <- kth_order_txt(4, chars_full, 1200)

# text generated
text_gen <- mytxt$text_gen

# Word accuracy of the generated text
wrd_acc <- mytxt$text_gen_word_acc

# Total letters combination generated
total_words <- mytxt$combined_letters_count

# Total real words generated
total_real_words <- mytxt$real_word_count

# Generating time per character
time_per_char_gen <- mytxt$time_per_char

# First transition Matrix Making time
time_matrix <- mytxt$mat_make_time

# Generating loops per character
loop_per_char_gen <- mytxt$loop_times

# Total generating time of the whole process
total_time <- mytxt$total_time

# Plot the characters generating time series
x_char <- c(1:length(time_per_char_gen))
plot(x_char, time_per_char_gen, type="l", lwd=1.6, col="blue", xlab="Generating Order", 
     ylab="Generating Time (sec)", main="Generating Time per Character",
     cex.main=1.5, cex.lab=1.36); grid()


# Plot the character generating loops series
x_char <- c(1:length(time_per_char_gen))
plot(x_char, loop_per_char_gen, type="l", lwd=1.6, col="blue", xlab="Generating Order", 
     ylab="Generating Loops", main="Generating Loops per Character",
     cex.main=1.5, cex.lab=1.36)



##############################
# Experiments Outcome Storage
##############################

expr_simple_gen <- function(transition_order, characters, text_length, expr_time) {
  
  total_words_group <- rep(NA, expr_time)
  total_real_words_group <- rep(NA, expr_time)
  times_per_char_gen <- list()
  times_matrix <- rep(NA, expr_time)
  total_times <- rep(NA, expr_time)
  
  for (i in (1:expr_time)) {
    mytxt <- kth_order_txt_simple(transition_order, characters, text_length)
    
    # Total letters combination generated
    total_words_group[i] <- mytxt$combined_letters_count
    
    # Total real words generated
    total_real_words_group[i] <- mytxt$real_word_count
    
    # Generating time per character
    times_per_char_gen[[i]] <- mytxt$time_per_char
    
    # First transition Matrix Making time
    times_matrix[i] <- mytxt$mat_make_time
    
    # Total generating time of the whole process
    total_times[i] <- mytxt$total_time
  }
  
  expr_lst <- list("Average word accuracy" = sum(total_real_words_group) / sum(total_words_group),
                   "Time per character generation (AVG)" = mean(unlist(times_per_char_gen)),
                   "Time per character generation (VAR)" = var(unlist(times_per_char_gen)),
                   "Initial transition matrix making time (AVG)" = mean(times_matrix),
                   "Total text generation time (AVG)"=mean(total_times),
                   "Total text generation time (VAR)"=var(total_times))
  return(expr_lst)
}

expr_complete_gen <- function(transition_order, characters, text_length, expr_time) {

  total_words_group <- rep(NA, expr_time)
  total_real_words_group <- rep(NA, expr_time)
  times_per_char_gen <- list()
  times_matrix <- rep(NA, expr_time)
  loops_per_char_gen <- list()
  total_times <- rep(NA, expr_time)
  
  for (i in (1:expr_time)) {
    mytxt <- kth_order_txt(transition_order, characters, text_length)
    
    # Total letters combination generated
    total_words_group[i] <- mytxt$combined_letters_count
    
    # Total real words generated
    total_real_words_group[i] <- mytxt$real_word_count
    
    # Generating time per character
    times_per_char_gen[[i]] <- mytxt$time_per_char
    
    # First transition Matrix Making time
    times_matrix[i] <- mytxt$mat_make_time
    
    # Generating loops per character
    loops_per_char_gen[[i]] <- mytxt$loop_times
    
    # Total generating time of the whole process
    total_times[i] <- mytxt$total_time
  }

  expr_lst <- list("Average word accuracy" = sum(total_real_words_group) / sum(total_words_group),
                   "Time per character generation (AVG)" = mean(unlist(times_per_char_gen)),
                   "Time per character generation (VAR)" = var(unlist(times_per_char_gen)),
                   "Loops per character generation (AVG)" = mean(unlist(loops_per_char_gen)),
                   "Loops per character generation (VAR)" = var(unlist(loops_per_char_gen)),
                   "Initial transition matrix making time (AVG)" = mean(times_matrix),
                   "Total text generation time (AVG)"=mean(total_times),
                   "Total text generation time (VAR)"=var(total_times))
  return(expr_lst)
}



# Group 1
expr_simple_gen(1, chars_short, 1200, 30)

expr_simple_gen(4, chars_short, 1200, 30)

expr_simple_gen(8, chars_short, 1200, 30)


# Group 2
expr_simple_gen(1, chars_half, 1200, 30)

expr_simple_gen(4, chars_half, 1200, 30)

expr_simple_gen(8, chars_half, 1200, 30)


# Group 3
expr_simple_gen(1, chars_full, 1200, 30)

expr_simple_gen(4, chars_full, 1200, 30)

expr_simple_gen(8, chars_full, 1200, 30)


# Group 4
expr_complete_gen(1, chars_short, 1200, 30)

expr_complete_gen(4, chars_short, 1200, 30)

expr_complete_gen(8, chars_short, 1200, 30)


# Group 5
expr_complete_gen(1, chars_half, 1200, 30)

expr_complete_gen(4, chars_half, 1200, 30)

expr_complete_gen(8, chars_half, 1200, 30)


# Group 6
expr_complete_gen(1, chars_full, 1200, 30)

expr_complete_gen(4, chars_full, 1200, 30)

expr_complete_gen(8, chars_full, 1200, 30)

