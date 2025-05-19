// clpam.ado
// ---------------------------------------------------------
// This program implements PAM (Partitioning Around Medoids) for unsupervised clustering.
// It is used to partition a dataset into a specified number of clusters based on a pairwise distance matrix.

// Jun 27 2017 13:12:10
// Copyright 2017 Brendan Halpin

// pam: Partitioning Around Medoids

// Newvar: variable holding the cluster solution
// Mandatory Options:
// - DISTMAT: pairwise distance matrix
// - IDVAR: ID variable to verify sort order of data and DISTMAT match
// - MEDOIDS: Either the number of clusters to fit, or a binary variable indicating which cases are medoids
// Optional options:
// - MANY or GA: mutually exclusive, if both are absent start depends on medoids() (if number, random)
//   - MANY: Start from multiple random starting points, select best
//   - GA: use a genetic algorithm to find an approximate best solution
// - NONVORONOI: use slightly slower non-Voronoi algorithm 

program clpam, rclass
syntax newvarlist(max=1), DISTmat(string) IDvar(varname) MEDoids(string) [ MANY GA NONVORonoi]

// Determine whether medoids option is a variable or a number
local medoidstype "undetermined"
capture su `medoids'
if (_rc == 0) {
  // If summary statistics succeed, then `medoids' is an existing variable.
  qui su `medoids'
  di "Medoids (n=`r(sum)') taken from `medoids'"
  local medoidstype "var"
}
else {
  // If summary fails, check if `medoids' can be interpreted as a number.
  if (_rc==111) {
    di "medoids var doesn't exist"
  }
  else {
    if (real("`medoids'")+0==`medoids') {
      di "Random starting medoids (Nk=`=floor(`medoids')')"
      local medoidstype "num"
    }
    else {
      di "Something is wrong with medoid variable"
    }
  }
}

// Process the NONVORonoi option: set flag (default is 0, i.e., use Voronoi algorithm).
local nonv 0
if ("`nonvoronoi'"!="") {
  local nonv 1
}

// Check that the distance matrix specified in DISTmat exists.
qui matlist `distmat'[1,1]

// Ensure that the observations are uniquely identified and sorted 
// according to the variable specified in IDvar.
isid `idvar', sort

// Initialize some Mata global variables for memoization used in the genetic algorithm routines.
mata: ga_memo = .
mata: ga_hit = .
mata: ga_insert = .

// Load the distance matrix into Mata.
mata: pam_temp_distmat = st_matrix("`distmat'")

// Determine which search strategy to use.
// Options: default (vanilla), MANY (multiple random starts), or GA (genetic algorithm).
if (("`many'"!="") & ("`ga'"!="")) {
  di as error "Only one of MANY and GA may be selected"
  exit
}
local searchtype 1
if ("`many'"!="") {
  local searchtype 2
}
if ("`ga'"!="") {
  local searchtype 3
}

// Call the appropriate Mata function based on searchtype and medoidstype.
if (("`searchtype'"=="1") | ("`searchtype'"=="")) {
  if ("`medoidstype'"=="var") {
    // If initial medoids are provided as a variable.
    mata: pam_temp_medoids = pam_vanilla(pam_temp_distmat, medoids_from_var("`medoids'"),`nonv')
  }
  else {
    // If the medoids option is given as a number (number of clusters).
    mata: pam_temp_medoids = pam_vanilla(pam_temp_distmat, random_medoids(`medoids', pam_temp_distmat), `nonv')
  }
}
else if ("`searchtype'" == "2") {
  di "Trying multiple starting points"
  if ("`medoidstype'"=="var") {
    di as error "PAM with multiple starting points (option MANY) is not compatible with specifying starting medoids"
    exit
  }
  else {
    // Use many random starts and select the best solution.
    mata: pam_temp_medoids = pam_vanilla(pam_temp_distmat, random_medoids_many(`medoids', pam_temp_distmat), `nonv')
  }
}
else if ("`searchtype'"=="3") {
  if ("`medoidstype'"=="var") {
    di as error "PAM with option GA is not compatible with specifying starting medoids"
    exit
  }
  else {
    // Use genetic algorithm to search for an approximate best solution.
    mata: pam_temp_medoids = pam_ga(pam_temp_distmat, `medoids', 100, `nonv')[1,]
  }
}

// Return the medoids
tempname ret_med
mata: st_matrix("`ret_med'",pam_temp_medoids)
return matrix medoids = `ret_med'

// Compute the group assignment:
// Each observation is assigned to the cluster (group) of the closest medoid.
mata: pam_temp_group = getgroup(pam_temp_distmat, pam_temp_medoids)

// Write the group variable back to Stata.
getmata `varlist'=pam_temp_group

// Clean up temporary Mata variables.
mata: mata drop pam_temp_distmat pam_temp_medoids pam_temp_group ga_memo ga_hit ga_insert

end // program pam

// ------------------------------
// Mata Functions used by clpam.ado
// ------------------------------

// Function: getgroup
// Purpose: Given a distance matrix and a vector of medoid indices,
//          assign each observation to the nearest medoid (i.e., determine the cluster group).
capture mata mata drop getgroup()
mata
real matrix getgroup (real matrix distmat, real matrix medoids) {
  // For each observation, find the minimum distance to any medoid.
  nearestvalue = rowmin(distmat[.,medoids])
  // Identify which medoid attains the minimum distance.
  nearestcase  = distmat[.,medoids] :== nearestvalue
  // Assign a group number based on the medoid index.
  group = rowmax(nearestcase:*transposeonly(range(1,length(medoids),1)))
  return(group)
}
end

// Function: medoids_from_var
// Purpose: Extract medoid indices from a variable passed in from Stata.
capture mata mata drop medoids_from_var()
mata
real matrix medoids_from_var(string varname) {
  // Create a view of the variable data.
  V = .
  st_view(V,.,varname)
  // Return the vector of medoid indices.
  return(select(range(1,length(V),1), V))
}
end


// Function: reassign_medoids
// Purpose: Given a set of medoids, reassign them by computing the best medoid for each cluster
//          using the Voronoi cell approach.
capture mata mata drop reassign_medoids()
mata
real matrix function reassign_medoids (matrix medoids, matrix distmat) {
  // Given a set of medoids, identify sets of elements closest to each
  // medoid, then return the set of true medoids of these groups
  // N: total number of observations; Nk: number of medoids.
  N = rows(distmat)
  Nk = length(medoids)
  // Obtain the group assignment based on current medoids.
  group = getgroup(distmat,medoids)
  // For each cluster, find the observation that minimizes the total within-cluster distance.
  newmedoids = J(Nk,1,.)
  for (i=1;i<=Nk;i++) {
    // Select the indices belonging to cluster i.
    selector = select(range(1,N,1), (group :== i))
    clusterdistmat = distmat[selector, selector]

    clustersize = rows(clusterdistmat)
    // Compute the average distance within the cluster.
    centroiddist = mean(clusterdistmat)
    mindist = min(centroiddist)
    // Choose the observation that minimizes the centroid distance.
    for (row = 1; row <= clustersize; row++) {
      if (mindist == centroiddist[row]) {
        newmedoids[i] = selector[row] // reverse select
      }
    }
  }
  newmedoids = uniqrows(newmedoids)
  return(newmedoids[order(newmedoids,1)])
}
end

// Function: reassign_medoids_non_voronoi
// Purpose: Similar to reassign_medoids but uses a non-Voronoi approach.
//          It selects a random candidate that improves the medoid cost.
capture mata mata drop reassign_medoids_non_voronoi()
mata
real matrix function reassign_medoids_non_voronoi (matrix medoids, matrix distmat) {
  // Given a set of medoids, identify sets of elements closest to each
  // medoid, then return the set of revised medoids of these groups
  // Non-voronoi: take first random medoid that improves
  N = rows(distmat)
  Nk = length(medoids)
  group = getgroup(distmat,medoids)
  // For each group find an improved medoid
  newmedoids = medoids //J(Nk,1,.)
  for (i=1;i<=Nk;i++) {
    // Obtain all indices in cluster i.
    selector = select(range(1,N,1), (group :== i))
    clusterdistmat = distmat[selector,]

    centroiddist = mean(clusterdistmat)
    maxdist = max(centroiddist)
    newdist = maxdist+1;
    // Try up to a fixed number of random candidates.
    for (j=1; j<=min((N,100)); j++) { // note break below
      newmedoid = runiformint(1,1,1,N) // 1  + floor(runiform(1,1)*N)
      newdist = centroiddist[newmedoid]
      if (newdist<centroiddist[medoids[i]]) {
        newmedoids[i] = newmedoid
        break
      }
    }
  }
  // If there is a change and the new medoids are valid, update medoids.
  if ((medoids != newmedoids) & (min(distmat[newmedoids,newmedoids] + I(length(newmedoids))) != 0)) {
    medoids = newmedoids
  }
  return(medoids[order(medoids,1)])
}
end


// Function: random_medoids
// Purpose: Generate an initial random set of medoids given a specified number Nk.
capture mata mata drop random_medoids()
mata
real matrix function random_medoids(real Nk, real matrix distmat) {
  N = rows(distmat)
  // Initialise with random medoids (Nk)
  medoids = uniqrows(floor(runiform(Nk,1):*N) :+ 1)
  // Ensure that there are no duplicates and that the selected medoids are valid.
  while ( (length(medoids)<Nk) | (min(distmat[medoids,medoids] + I(length(medoids))) == 0)) {
    medoids = uniqrows(floor(runiform(Nk,1):*N) :+ 1)
  }
  return(medoids)
}
end

// Function: random_medoids_many
// Purpose: Generate multiple sets of random medoids for use with the MANY option.
capture mata mata drop random_medoids_many()
mata
real matrix random_medoids_many(real Nk, real matrix distmat) {
  many = 100
  meds = J(many,Nk,.)
  for (i=1;i<=many; i++) {
    meds[i,] = transposeonly(random_medoids(Nk,distmat))
  }
  return(meds)
}
end

// Function: pam_sum_within_clusters
// Purpose: Compute the sum of within-cluster distances (a measure of clustering quality).
//          This function uses memoization to speed up repeated calculations.
capture mata mata drop pam_sum_within_clusters()
mata
real scalar function pam_sum_within_clusters (real matrix distmat, real matrix medoids) {
  external ga_memo, ga_hit, ga_insert
  N = rows(distmat)
  Nk = length(medoids)
  group = getgroup(distmat,medoids)
  
  if (asarray_contains(ga_memo, transposeonly(medoids))) {
    // If this set of medoids has been evaluated before, retrieve the stored value.
    SS = asarray(ga_memo, transposeonly(medoids))
    ga_hit++
  } else {
    // Otherwise, compute the sum of within-cluster distances for each group.
    SS = 0
    for (i=1;i<=Nk;i++) {
      // Select the cluster
      selector = select(range(1,N,1), (group :== i))
      clusterdistmat = distmat[selector, selector]
      SStemp =sum(clusterdistmat)/rows(clusterdistmat)
      SS = SS + SStemp
    }
    // Store the computed value for memoization.
    asarray(ga_memo, transposeonly(medoids), SS)
    ga_insert++
  }
  return(SS)
}
end

// Function: pam_iterate
// Purpose: Iteratively reassign medoids until convergence (i.e., the medoid set does not change).
//          This function implements a hill-climbing procedure.
capture mata mata drop pam_iterate()
mata
real matrix pam_iterate (real matrix distmat, real matrix medoids, real scalar nonvor) {
  // With a proposed set of medoids, keep reassigning them until stability
  // Effectively a hill-climb function
  Nk=rows(medoids)
  oldmed = J(Nk,1,0)
  i=0
  // Continue iterating until the medoids do not change.
  while ((oldmed != medoids)) {
    i++
    oldmed = medoids
    // If defective, randomise
    if ( (length(medoids)<Nk) | (min(distmat[medoids,medoids] + I(length(medoids))) == 0)) {
      medoids = random_medoids(Nk, distmat)
    }
    if (nonvor) {
      medoids = reassign_medoids_non_voronoi(medoids, distmat)
    }
    else {
      medoids = reassign_medoids(medoids, distmat)
    }
  }
  return(medoids)
}
end

// Function: pam_ga
// Purpose: Use a genetic algorithm to search for an approximate optimal set of medoids.
//          This function creates a population of candidate medoid sets, evaluates their fitness,
//          and uses crossover and mutation to evolve better solutions.
capture mata mata drop pam_ga()
mata
real matrix function pam_ga (real matrix distmat, real scalar Nk, real scalar npop, real scalar nonv) {
  external ga_memo, ga_hit, ga_insert

  // Memoise the fitness function
  // Is worth about 10-12% in time.
  // Used by pam_sum_within_clusters()
  ga_memo = asarray_create("real",Nk)
  ga_hit = 0
  ga_insert = 0

  npop = 2*floor(npop/2) // even
  nparents = floor(npop/2) // change from 50% to intensify selection
  nsurv = nparents
  nnew = floor(npop*0.20)
  ntoconv = ceil(npop*0.1)
  // Use top 50% as core; these will reproduce to create new candidates.

  // Create a random population, each row a medoid set
  population = J(npop,Nk,.)
  sigmadist = J(npop,1,.)
  for (i=1; i<=npop; i++) {
    population[i,.] = transposeonly(random_medoids(Nk, distmat))
//    sigmadist[i] = pam_sum_within_clusters(distmat, transposeonly(population[i,.]))
  }

  // Main iteration
  iter = 0
  while (iter==0 |(max(sigmadist[1..ntoconv]) != min(sigmadist[1..ntoconv]))) {
    newpop = J(npop,Nk,.)
    iter++
    // Improve each candidate by applying the PAM hill-climbing procedure.
    for (i=1; i<=npop; i++) {
      medoids = transposeonly(population[i,.])
      newpop[i,.] = transposeonly(pam_iterate(distmat,medoids,nonv))
    }

    // Calculate fitness (lower is better)
    for (i=1; i<=npop; i++) {
      medoids = transposeonly(newpop[i,.])
      sigmadist[i] = pam_sum_within_clusters(distmat, medoids)
    }
    
    // Order the new medoid sets per fitness (best (lowest) first)
    newpop = newpop[order(sigmadist,1),][1..nparents,]
    sigmadist = sigmadist[order(sigmadist,1)]
    "Iter "+strofreal(iter,"%-3.0f")+": � dist: min "+strofreal(min(sigmadist[1..ntoconv]))+"; mean: "+strofreal(mean(sigmadist[1..ntoconv]))+"; max: "+strofreal(max(sigmadist[1..ntoconv]))
    displayflush()
  
    // Keep top surviving candidates.
    population[1..nsurv,] = newpop[1..nsurv,]

    // Crossover: combine pairs of parent candidates to generate new candidates.
    for (i=1+nsurv; i<=npop-nnew; i++) {
      j = 1+floor((1-runiform(1,1)^0.1)*nparents)
      k = 1+floor((1-runiform(1,1)^0.1)*nparents)
      l = runiformint(1,1,1,Nk-1)
      population[i,1..l] = newpop[j,1..l]
      population[i,l+1..Nk] = newpop[k,l+1..Nk]
    }
    // Mutation: replace a portion of the population with new random candidates.
    for (i=1+npop-nnew; i<=npop; i++) {
      if (runiform(1,1)<0.1) {
        population[i,.] = transposeonly(random_medoids(Nk, distmat))
      }
    }
  }
  sprintf("Memoisation: hits: %1.0f; inserts: %1.0f.", ga_hit,ga_insert)
  return(population)
}
end


// Function: pam_vanilla
// Purpose: The basic (vanilla) implementation of PAM. This function iteratively improves
//          the initial medoid set until convergence using a hill-climbing approach.
capture mata mata drop pam_vanilla()
mata
real matrix function pam_vanilla(real matrix distmat, real medoids, real scalar nonv) {
  external ga_memo, ga_hit, ga_insert

  // If medoids is a matrix, each row is a medoid set
  // If it's a vector, it's a col
  many_meds = min((rows(medoids),cols(medoids)))
  if (many_meds>1) {
    // When multiple candidate medoid sets are provided.
    Nk = length(transposeonly(medoids[1,]))
  }
  else {
    // When medoids is a single vector.
    Nk = length(medoids)
  }

  // Memoise the fitness function
  // Is worth about 10-12% in time.
  // Used by pam_sum_within_clusters()
  ga_memo = asarray_create("real",Nk)
  ga_hit = 0
  ga_insert = 0

  if (many_meds>1) {
    results = J(many_meds,length(medoids[1,]),.)
    scores  = J(many_meds,1,.)
    for (i=1; i<=many_meds; i++) {
      results[i,1..Nk] = transposeonly(pam_iterate(distmat,transposeonly(medoids[i,]),nonv))
      scores[i] = pam_sum_within_clusters(distmat, transposeonly(medoids[i,]))
    }
    result=results[order(scores,1),][1,]
  }
  else {
    result = pam_iterate(distmat,medoids,nonv)
  }
  return(result)
}
end
