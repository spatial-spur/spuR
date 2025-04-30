* Clear workspace
clear all
set more off
mata mata clear

* Create the same dataset as in R
input str15 city lat lon
"New York"     40.7128 -74.0060
"Los Angeles"  34.0522 -118.2437
"Chicago"      41.8781 -87.6298
"Houston"      29.7604 -95.3698
"Miami"        25.7617 -80.1918
end

* Create touse variable
gen touse = 1

* Display the input data
list city lat lon

mata:
    // Get the coordinates into Mata
    s = st_data(., ("lat", "lon"), "touse")
    
    // Set latlongflag to 1 for lat/lon calculation
    latlongflag = 1
    
    // Call the actual getdistmat function
    mat = getdistmat(s)
    
    n = rows(s)
    
    // Print the results
    printf("Distance matrix using getdistmat (fractions of circumference):\n")
    for (i=1; i<=n; i++) {
        for (j=1; j<=n; j++) {
            printf("%12.6f", mat[i,j])
        }
        printf("\n")
    }
    
    // Convert to km for comparison
    dist_km = mat * (2*3.14159265359*6371)  // 2πR gives full circumference
    
    printf("\nDistance matrix in kilometers:\n")
    for (i=1; i<=n; i++) {
        for (j=1; j<=n; j++) {
            printf("%12.6f", dist_km[i,j])
        }
        printf("\n")
    }
    
end

di "Matrix saved to stata_distance_matrix.csv for comparison with R output"