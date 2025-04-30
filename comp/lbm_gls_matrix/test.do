* Clear workspace
clear all
set more off
mata mata clear

* Create test dataset with spatial coordinates
input str15 city lat lon
"New York"     40.7128 -74.0060
"Los Angeles"  34.0522 -118.2437
"Chicago"      41.8781 -87.6298
"Houston"      29.7604 -95.3698
"Miami"        25.7617 -80.1918
end

* Display input data
list city lat lon

* Create touse variable for all observations
gen touse = 1

* Test in Mata
mata:
    // Get coordinates
    s = st_data(., ("lat", "lon"), "touse")
    
    // Test 1: Using Euclidean distances
    latlongflag = 0
    
    H_euclidean = lbm_gls_matrix(s)
    
    printf("\nTest 1: lbm_gls_matrix with Euclidean distances\n")
    printf("Matrix dimension: %d x %d\n", rows(H_euclidean), cols(H_euclidean))
    printf("Matrix trace: %g\n", trace(H_euclidean))  // Changed format to %g
    printf("Full matrix:\n")
    H_euclidean  // Print the full matrix
    
    // Test 2: Using lat/long distances
    latlongflag = 1
    
    H_latlong = lbm_gls_matrix(s)
    
    printf("\nTest 2: lbm_gls_matrix with lat/long distances\n")
    printf("Matrix dimension: %d x %d\n", rows(H_latlong), cols(H_latlong))
    printf("Matrix trace: %g\n", trace(H_latlong))  // Changed format to %g
    printf("Full matrix:\n")
    H_latlong  // Print the full matrix
end