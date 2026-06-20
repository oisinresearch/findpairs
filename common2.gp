\\ =========================================================================
\\ 1. Initialization and Helper Functions
\\ =========================================================================

addhelp(process_relations, "process_relations(f1, f2, R1, R2, B, n): \nTakes pre-computed shared relations and field polynomials, constructs ideal factor bases/valuations, calculates left kernels, and returns the stacked block matrix Mout.\n\nReturns: [R1, R2, ri1, ri2, M1, M2, Mout]");

\\ Factor Base Constructor
construct_FB(f, B_bound) = {
  my(K = nfinit(f));
  my(primes = []);
  forprime(p=2, B_bound,
    my(dec = idealprimedec(K, p));
    for(i=1, #dec,
      primes = concat(primes, [idealhnf(K, dec[i])]);
    )
  );
  return([K, primes]);
}

\\ =========================================================================
\\ 2. Main Execution Function
\\ =========================================================================
process_relations(f1, f2, R1, R2, B, n) = {
  my(num_common = #R1);
  
  if(num_common != #R2,
    error("Error: R1 and R2 must contain the same number of relations.");
  );
  if(num_common == 0,
    print("Error: Empty relation vectors provided.");
    return([]);
  );

  print("Step 1: Constructing Ideal Factor Bases up to B = ", B, "...");
  my(FB_data1 = construct_FB(f1, B));
  my(K1 = FB_data1[1]); my(FB1 = FB_data1[2]);

  my(FB_data2 = construct_FB(f2, B));
  my(K2 = FB_data2[1]); my(FB2 = FB_data2[2]);

  my(num_p1 = #FB1);
  my(num_p2 = #FB2);
  print("   -> Field 1 Factor Base Size: ", num_p1);
  print("   -> Field 2 Factor Base Size: ", num_p2);

  \\ Maps for fast O(1) index lookup
  my(map1 = Map()); for(j=1, num_p1, mapput(map1, Str(FB1[j]), j));
  my(map2 = Map()); for(j=1, num_p2, mapput(map2, Str(FB2[j]), j));

  print("Step 2: Calculating Ideal Valuations & Residues for relations...");
  my(M1 = matrix(num_common, num_p1));
  my(M2 = matrix(num_common, num_p2));
  my(ri1 = vector(num_common));
  my(ri2 = vector(num_common));

  for(i=1, num_common,
    my(elA = R1[i]);
    ri1[i] = Mod(abs(nfeltnorm(K1, elA)), n);
    my(facA = idealfactor(K1, elA));
    for(r=1, matsize(facA)[1],
      my(P_hnf = idealhnf(K1, facA[r,1]));
      my(key = Str(P_hnf));
      if(mapisdefined(map1, key),
        my(idx = mapget(map1, key));
        M1[i, idx] = facA[r,2];
      )
    );

    my(elB = R2[i]);
    ri2[i] = Mod(abs(nfeltnorm(K2, elB)), n);
    my(facB = idealfactor(K2, elB));
    for(r=1, matsize(facB)[1],
      my(P_hnf = idealhnf(K2, facB[r,1]));
      my(key = Str(P_hnf));
      if(mapisdefined(map2, key),
        my(idx = mapget(map2, key));
        M2[i, idx] = facB[r,2];
      )
    );
  );

  print("Step 3: Calculating Left Kernels...");
  \\ matker(M~)~ finds the left kernel (rows are the kernel vectors)
  my(k1 = matker(M1~)~);
  my(k2 = matker(M2~)~);

  my(d1 = matsize(k1)[1]);
  my(d2 = matsize(k2)[1]);
  print("   -> Field 1 Kernel Dimensions: ", d1, " x ", num_common);
  print("   -> Field 2 Kernel Dimensions: ", d2, " x ", num_common);

  print("Step 4: Stacking Kernels Vertically into Mout...");
  my(Mout);
  \\ Use block-matrix construction natively via matconcat to prevent loops
  if(d1 == 0 && d2 == 0,
    Mout = matrix(0, num_common);
  , d1 == 0,
    Mout = k2
  , d2 == 0,
    Mout = k1
  ,
    Mout = matconcat([k1; k2])
  );

  print("\n=================== DONE ===================");
  print("Returning block data: [R1, R2, ri1, ri2, M1, M2, Mout]");

  return([R1, R2, ri1, ri2, M1, M2, Mout]);
}

