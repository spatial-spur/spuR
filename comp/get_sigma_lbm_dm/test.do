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
    
    // Call get_sigma_lbm_dm
    sigma_lbm_dm = get_sigma_lbm_dm(distmat)
    
    // Print the result
    printf("\nDemeaned sigma_lbm matrix (Stata implementation):\n")
    sigma_lbm_dm
    
    // Check properties of demeaned matrix
    printf("\nMatrix properties check:\n")
    printf("Matrix dimension: %d x %d\n", rows(sigma_lbm_dm), cols(sigma_lbm_dm))
    printf("Row sums (should be ~0):\n")
    rowsum(sigma_lbm_dm)'
    printf("Column sums (should be ~0):\n")
    colsum(sigma_lbm_dm)
    printf("Trace: %g\n", trace(sigma_lbm_dm))
    
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
    
    // Call get_sigma_lbm_dm
    sigma_lbm_dm = get_sigma_lbm_dm(distmat)
    
    // Print the result
    printf("\nDemeaned sigma_lbm matrix (Stata implementation):\n")
    sigma_lbm_dm
    
    // Check properties of demeaned matrix
    printf("\nMatrix properties check:\n")
    printf("Matrix dimension: %d x %d\n", rows(sigma_lbm_dm), cols(sigma_lbm_dm))
    printf("Row sums (should be ~0):\n")
    rowsum(sigma_lbm_dm)'
    printf("Column sums (should be ~0):\n")
    colsum(sigma_lbm_dm)
    printf("Trace: %g\n", trace(sigma_lbm_dm))
    
    // Print specific elements for comparison
    printf("\nSpecific elements for comparison:\n")
    printf("[1,1]: %g\n", sigma_lbm_dm[1,1])
    printf("[2,3]: %g\n", sigma_lbm_dm[2

end