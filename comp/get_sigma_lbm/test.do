* Clear workspace
clear all
set more off
mata mata clear

mata:
    // Create the same test distance matrix as in R
    distmat = (
        0.0000, 0.6124, 0.2171, 0.4829, 0.8331 \ 
        0.6124, 0.0000, 0.5251, 0.3798, 0.7481 \
        0.2171, 0.5251, 0.0000, 0.3091, 0.6509 \
        0.4829, 0.3798, 0.3091, 0.0000, 0.3553 \
        0.8331, 0.7481, 0.6509, 0.3553, 0.0000
    )
    
    // Print the input matrix
    printf("Input distance matrix:\n")
    distmat
    
    // Run get_sigma_lbm on the test matrix
    sigma_lbm = get_sigma_lbm(distmat)
    
    // Print the result
    printf("\nComputed LBM covariance matrix (Stata implementation):\n")
    sigma_lbm
    
    // Print specific values for easier comparison
    printf("\nKey values for comparison:\n")
    printf("First element [1,1]: %g\n", sigma_lbm[1,1])
    printf("Element [2,3]: %g\n", sigma_lbm[2,3])
    printf("Element [4,5]: %g\n", sigma_lbm[4,5])
    printf("Matrix trace: %g\n", trace(sigma_lbm))
end