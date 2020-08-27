#Quick script to compare outputs that may be rounded.
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(Xmisc))
parser <- ArgumentParser$new()
parser$add_description("Script to quickly compare to files and see if they are sufficiently close, since there may be rounding error.")
parser$add_argument("--original", type = 'character', help = "Path to the original file")
parser$add_argument("--new", type = "character", help = "Path to file to compare it with")
parser$add_argument("--round", type = "numeric", help = "How many dps to round to.", default = 5)
parser$add_argument("--sum", type = "character", action = "store_true", help = "Specify a file to add to the new")
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$helpme()
args <- parser$get_args()
original <- read_tsv(args$original) %>% arrange(IID)
new <- read_tsv(args$new) %>% arrange(IID)

print("The dimensions of the data:")
paste("New", paste(dim(new)))
paste("Original", paste(dim(original)))
print("Do the IDs match?")
print(sum(new$IID !=  original$IID) == 0) #should be 0
row_names <- original$IID
#Okay, let's round everything to 3 dps and see what we get....\
original_r <- as.matrix(round(original %>% select(-IID), digits = args$round))
if(length(args$sum) > 0)
{
    add_in <- read_tsv(args$sum) %>% arrange(IID)
    add_in_m <- as.matrix(add_in %>% select(-IID))
    to_add <- as.matrix(new %>% select(-IID))
    new_r <- round(add_in_m + to_add, digits = args$round)
} else
 {
new_r <- as.matrix(round(new %>% select(-IID), digits = args$round))
}

diff <- abs(original_r - new_r) 
print(paste("Max error we see", max(diff)))

if(max(diff) == 0) {
print("R claims these rounded datasets have no difference. Is this true?")
print(paste("Test: 1 is true, anything else is false:", 1- sum(original_r != new_r)))
}
#Max expected error:
err <- 10^(-(args$round))
print(paste("Rounding error", err))
max_th_error <- err * dim(diff)[1] * dim(diff)[2]
max_true_error <- max(diff)
max_practical_error <- max_true_error * sum(diff != 0)
total_true_error <- sum(diff)
print(paste("Max theoretical error possible:", max_th_error))
print(paste("Max possible error in dataset:", max_true_error, "x", sum(diff != 0), " = ", max_practical_error ))
print(paste("True error", total_true_error))
if(max_th_error >= total_true_error && max_practical_error >= total_true_error)
{
    print("The difference between these results cannot be confidently assigned to anything other than rounding error.")
} else {
    print("This data likely came from a pipeline error. Take a closer look")
}
