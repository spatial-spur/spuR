
* Clear workspace
clear all
set more off
mata mata clear

* Create the same test dataset
input x y
1.5 0.5
2.3 1.2
3.7 1.8
4.2 2.5
5.8 3.1
end

* Create a touse variable for all observations
gen touse_all = 1

* Create a touse variable for subsetting
gen touse_sub = 0
replace touse_sub = 1 in 1
replace touse_sub = 1 in 3
replace touse_sub = 1 in 4

* Define the transform function in mata
mata:
real vector transform(string scalar varname, real matrix H, string scalar touse, string scalar transformation)
{
    real vector y, hy
    
    y = st_data(., varname, touse)
    
    hy = H * y
    
    if (transformation == "lbmgls") {
        hy = hy:-mean(hy)
    }
    
    return(hy)
}
end

* Create the test matrices in mata
mata:
// 1. Identity matrix
H_identity = I(5)

// 2. Simple averaging matrix
H_avg = J(5,5,0)
for(i=1; i<=5; i++) {
    if(i == 1) {
        H_avg[i, 1..2] = (0.7, 0.3)
    } 
    else if(i == 5) {
        H_avg[i, 4..5] = (0.3, 0.7)
    } 
    else {
        H_avg[i, (i-1)..(i+1)] = (0.2, 0.6, 0.2)
    }
}

// 3. Symmetric matrix
H_sym = (
    2.1, -0.3, -0.5, -0.2, -0.1 \ 
    -0.3, 1.9, -0.4, -0.1, -0.1 \ 
    -0.5, -0.4, 2.5, -0.6, -0.3 \ 
    -0.2, -0.1, -0.6, 1.8, -0.5 \ 
    -0.1, -0.1, -0.3, -0.5, 2.0  
)

// Run tests
printf("Test 1: Identity matrix, 'x' variable, no demeaning\n")
result1 = transform("x", H_identity, "touse_all", "iso")
result1

printf("\nTest 2: Averaging matrix, 'x' variable, no demeaning\n")
result2 = transform("x", H_avg, "touse_all", "iso")
result2

printf("\nTest 3: Symmetric matrix, 'y' variable, with demeaning (lbmgls)\n")
result3 = transform("y", H_sym, "touse_all", "lbmgls")
result3
printf("Mean of result3: %g\n", mean(result3))

printf("\nTest 4: Identity matrix, 'x' variable, with subsetting\n")
// Extract submatrix for the subset
idx = selectindex(st_data(., "touse_sub") :== 1)
H_sub = H_identity[idx,idx]
result4 = transform("x", H_sub, "touse_sub", "iso")
result4

// Save to text file
fh = fopen("transform_stata_results.txt", "w")
fwrite(fh, "Stata transform() test results:" + char(10))
fwrite(fh, "Test 1: " + strofreal(result1', "%8.6f") + char(10))
fwrite(fh, "Test 2: " + strofreal(result2', "%8.6f") + char(10))
fwrite(fh, "Test 3: " + strofreal(result3', "%8.6f") + char(10))
fwrite(fh, "Test 4: " + strofreal(result4', "%8.6f") + char(10))
fclose(fh)

end