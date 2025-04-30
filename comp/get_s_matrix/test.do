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

* Create a touse variable (this is what get_s_matrix expects)
gen touse = 1

* Display the input data
list

* Test with latlong = FALSE (Euclidean distance)
di as text "==== Testing with latlong = FALSE ===="
mata: s = get_s_matrix("touse", "")
mata: s
mata: latlongflag

/*

result:
                  1              2
    +-------------------------------+
  1 |   40.71279907    -74.0059967  |
  2 |   34.05220032   -118.2436981  |
  3 |   41.87810135   -87.62979889  |
  4 |   29.76040077   -95.36979675  |
  5 |   33.44839859   -112.0739975  |
  6 |   39.95259857   -75.16519928  |
    +-------------------------------

#/