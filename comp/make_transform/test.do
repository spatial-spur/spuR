* Clear workspace
clear all
set more off
mata mata clear

* Create sample data with the same values as your R script
input str15 city s_1 s_2
"New York"     40.7128 -74.0060
"Los Angeles"  34.0522 -118.2437
"Chicago"      41.8781 -87.6298
"Houston"      29.7604 -95.3698
"Phoenix"      33.4484 -112.0740
"Philadelphia" 39.9526 -75.1652
end

* Display the input data
list

* Create a touse variable
gen touse = 1

* Test 1: LBM-GLS transformation
di as text "==== Testing LBM-GLS transformation ===="
mata: s = get_s_matrix("touse", "latlong")
mata: H = make_transform(s, "lbmgls", -1, 0)
mata: rows(H)
mata: cols(H)
mata: H[1..3,1..3]
mata: trace(H)
* same as R