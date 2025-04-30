* Clear workspace
clear all
set more off

* Create sample data
input str15 city cluster str5 region
"New York"      1 "East"
"Los Angeles"   2 "West"
"Chicago"       1 "East"
"Houston"       2 "West"
"Phoenix"       2 "West"
"Philadelphia"  1 "East"
end

di as text "==== Test 1: All rows ===="
gen touse_all = 1
mata: result1 = get_cluster_matrix("cluster", "touse_all")
mata: result1
mata: rows(result1)
mata: cols(result1)
