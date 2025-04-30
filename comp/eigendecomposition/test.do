* Clear workspace
clear all
set more off
mata mata clear

mata:
    // Create the same test distance matrix
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
    
    // Test 1: Get demeaned sigma_lbm matrix
    printf("\nTest 1: get_sigma_lbm_dm function\n")
    sigma_lbm_dm = get_sigma_lbm_dm(distmat)
    printf("Result dimension: %d x %d\n", rows(sigma_lbm_dm), cols(sigma_lbm_dm))
    printf("Matrix trace: %g\n", trace(sigma_lbm_dm))
    printf("Row sums (should be ~0):")
    rowsum(sigma_lbm_dm)'
    printf("Full matrix:\n")
    sigma_lbm_dm
    
    // Test 3: Check eigendecomposition (matching R)
    printf("\nTest 3: Check eigendecomposition\n")
    V = .
    d = .
    symeigensystem(sigma_lbm_dm, V, d)
    d = d'
    
    // Sort in descending order
    ii = order(d, -1)
    eval = d[ii]
    evec = V[.,ii]
    
    // Apply threshold
    small = 1e-10
    ii = eval :> small
    eval = select(eval, ii)
    evec = select(evec, ii')
    
    // Print eigenvalues for comparison
    printf("Number of eigenvalues above threshold: %d\n", length(eval))
    printf("Eigenvalues (descending):")
    eval
    printf("First eigenvector:")
    evec[.,1]'
    printf("Last eigenvector:")
    evec[.,cols(evec)]'
    
    // Test 4: Create and verify GLS transformation matrix
    dsi = 1 :/ sqrt(eval)
    Dsi = diag(dsi)
    H = evec*Dsi*evec'
    
    printf("\nTest 4: Final transformation matrix properties\n")
    printf("Matrix shape: %d x %d\n", rows(H), cols(H))
    printf("Row sums:")
    rowsum(H)[1..5]'
    printf("First element [1,1]: %g\n", H[1,1])
    printf("Element [2,3]: %g\n", H[2,3])
    printf("Element [4,5]: %g\n", H[4,5])
end