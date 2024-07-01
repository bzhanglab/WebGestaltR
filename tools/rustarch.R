# Function to determine Rust target using rustup default
get_rust_target <- function() {
    # Get the output of `rustup default`
    rustup_output <- system("rustup default", intern = TRUE)

    # Extract the target from the output
    target <- sub(" .*", "", rustup_output)
    target_split <- strsplit(target, "-")
    target <- paste(target_split[[1]][-1], collapse = "-")

    # Return the determined target
    return(target)
}

# Get the Rust target
rust_target <- get_rust_target()

# Print the result
cat(rust_target)
