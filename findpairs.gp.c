/*-*- compile-command: "cc -c -o findpairs.gp.o -g -O3 -Wall -fno-strict-aliasing -ffp-contract=off -fPIC -I\"/usr/include/x86_64-linux-gnu\" findpairs.gp.c && cc -o findpairs.gp.so -shared -g -O3 -Wall -fno-strict-aliasing -ffp-contract=off -fPIC -Wl,-shared -Wl,-z,relro findpairs.gp.o -lc -lm -L/usr/lib/x86_64-linux-gnu -lpari"; -*-*/
#include <pari/pari.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
/*
install("init_findpairs","v","init_findpairs","./findpairs.gp.so");
install("issmooth","lD0,G,D0,G,","issmooth","./findpairs.gp.so");
install("fact_str","D0,G,","fact_str","./findpairs.gp.so");
install("process_one_field","D0,G,D0,G,D0,G,D0,G,D0,G,D0,G,","process_one_field","./findpairs.gp.so");
install("findpairs","GD0,G,D0,G,D0,G,D0,G,D0,G,D0,G,D0,G,D0,G,p","findpairs","./findpairs.gp.so", "findpairs(num_fields, target, B, N0, start_D, mlim, threads, A):\nSearches for matching smooth norms across fundamental discriminants.\n  num_fields : Number of fundamental discriminants to process\n  target     : Target number of smooth norms to find per field\n  B          : Smoothness bound (max prime factor)\n  N0         : Minimum acceptable norm magnitude\n  start_D    : Starting point for fundamental discriminants\n  mlim       : Timeout in minutes per field\n  threads    : Number of parallel worker threads\n  A          : Sieve region parameter (defines [-A/2, A/2) x [1, 2A])");
install("anon_0","D0,G,GGGGGG","anon_0","./findpairs.gp.so");
install("anon_1","D0,G,GG","anon_1","./findpairs.gp.so");
*/
void init_findpairs(void);
long issmooth(GEN n, GEN bound);
GEN fact_str(GEN n);
GEN process_one_field(GEN D, GEN target, GEN B, GEN N0, GEN mlim, GEN A_sieve_in);
GEN findpairs(GEN num_fields, GEN target, GEN B, GEN N0, GEN start_D, GEN mlim, GEN t, GEN A_sieve_in, long prec);
GEN anon_0(GEN i, GEN D_list, GEN target, GEN B, GEN N0, GEN mlim, GEN A_sieve_in);
GEN anon_1(GEN i, GEN num_fields, GEN field_norms);
/*End of prototype*/

void
init_findpairs(void)	  /* void */
{
  /* Export functions required by the parallel workers */
  //export(strtofunction("issmooth"));
  //export(strtofunction("process_one_field"));
  return;
}

/* ========================================================================= */
/* 1. Define Helper Functions */
/* ========================================================================= */

long
issmooth(GEN n, GEN bound)
{
  GEN fa;
  if (gcmpgs(n, 1) <= 0)
    return 1;
  fa = factor(n);
  return gcmp(gcoeff(fa, glength(gtrans(fa)), 1), bound) <= 0;
}

GEN
fact_str(GEN n)
{
  GEN v;
  GEN p1;	  /* vec */
  GEN s;
  long l2;
  /* Handle negative norms by prepending a minus sign recursively */
  if (gcmpgs(n, 0) < 0)
  {
    GEN abs_str = fact_str(gneg(n));
    GEN p_neg = cgetg(3, t_VEC);
    gel(p_neg, 1) = strtoGENstr("-");
    gel(p_neg, 2) = abs_str;
    return gconcat1(p_neg);
  }
  v = factor(n);
  {
    long l3, l4;
    p1 = cgetg(3, t_MAT);
    for (l4 = 1; l4 <= 2; ++l4)
    {
      gel(p1, l4) = cgetg(1, t_COL);
      for (l3 = 1; l3 <= 0; ++l3)
        gcoeff(p1, l3, l4) = gen_0;
    }
  }
  if (gequal(v, p1))
    return strtoGENstr("1");
  s = strtoGENstr("");
  l2 = glength(gel(v, 1));
  {
    long i;
    for (i = 1; i <= l2; ++i)
    {
      GEN p, e, term;
      p = gcopy(gcoeff(v, i, 1));
      e = gcopy(gcoeff(v, i, 2));
      term = gcopy(GENtoGENstr(p));
      if (gcmpgs(e, 1) > 0)
      {
        GEN p5;	  /* vec */
        p5 = cgetg(4, t_VEC);
        gel(p5, 1) = gcopy(term);
        gel(p5, 2) = strtoGENstr("^");
        gel(p5, 3) = gcopy(GENtoGENstr(e));
        term = gconcat1(p5);
      }
      if (gequal(s, strtoGENstr("")))
        s = term;
      else
      {
        GEN p6;	  /* vec */
        p6 = cgetg(4, t_VEC);
        gel(p6, 1) = gcopy(s);
        gel(p6, 2) = strtoGENstr(" * ");
        gel(p6, 3) = gcopy(term);
        s = gconcat1(p6);
      }
    }
  }
  return s;
}

/* Heavy lifting worker function */
GEN
process_one_field(GEN D, GEN target, GEN B, GEN N0, GEN mlim, GEN A_sieve_in)	  /* vec */
{
  pari_sp av = avma; // Save current stack pointer
  GEN b, a = gen_1, c, f, m, count = gen_0, start_time, timed_out = gen_0, FB_list, FB;
  long l1;
  GEN norms_vec, mat_res;
  GEN p2, p3;	  /* vec */
  /* 1. Construct optimized polynomial f(x) = ax^2 + bx + c */
  /* Ensures a=1, b and c are O(sqrt(D)), and b^2 - 4ac = D */
  b = sqrtint(D);
  if (!gequal(gmodgs(b, 2), gmodgs(D, 2)))
    b = gsubgs(b, 1);
  c = gdivgs(gsub(gsqr(b), D), 4);
  f = gadd(gadd(gmul(a, gsqr(pol_x(0))), gmul(b, pol_x(0))), c);
  m = mkmap();
  start_time = getwalltime();
  /* 2. Compute the Factor Base */
  /* Collects prime ideals (p, r) where p <= B and r is a root of f mod p */
  FB_list = mklist();
  {
    forprime_t iter;	  /* forprime */
    u_forprime_init(&iter, 2, gtos(B));;
    {
      long p;
      while ((p = u_forprime_next(&iter)))
      {
        GEN roots;
        long l4;
        roots = polrootsmod(f, stoi(p));
        l4 = glength(roots);
        {
          long k;
          for (k = 1; k <= l4; ++k)
          {
            GEN p5;	  /* vec */
            p5 = cgetg(3, t_VEC);
            gel(p5, 1) = stoi(p);
            gel(p5, 2) = lift(gel(roots, k));
            listput0(FB_list, p5, 0);
          }
        }
      }
    }
  }
  FB = gtovec(FB_list);
  /* ======================================================================= */
  /* 3. Franke-Kleinjung Efficient 2D Lattice Sieve (Outer Loop) */
  /* ======================================================================= */
  long A_sieve = itos(A_sieve_in);
  long J_max = 2 * A_sieve;
  long total_elements = A_sieve * J_max;

  /* Allocate the raw sieve array strictly in native C */
  uint8_t* sieve_array = (uint8_t*) calloc(total_elements, sizeof(uint8_t));

  l1 = glength(FB);
  {
    long i;
    for (i = 1; i <= l1; ++i)
    {
      GEN p_gen, r_gen;
      long p_val, r_val;

      /* Early timeout check to prevent long-tail hangs on difficult fields */
      if (gcmp(gsub(getwalltime(), start_time), gmulgs(mlim, 60000)) > 0)
      {
        timed_out = gen_1;
        break;
      }

      /* Extract elements to native C types for the inner loop */
      p_gen = gel(gel(FB, i), 1);
      r_gen = gel(gel(FB, i), 2);
      p_val = itos(p_gen);
      r_val = itos(r_gen);

      /* ------------------------------------------------------------------- */
      /* Inner Lattice Sieve Routine (Franke-Kleinjung) */
      /* ------------------------------------------------------------------- */
      long i0 = -p_val, j0 = 0;
      long i1 = r_val,  j1 = 1;

      while (1) {
          if (labs(i1) == 1 || labs(i1) < A_sieve) break;
          long a_k = labs(i0) / labs(i1);
          long next_i = i0 + a_k * i1;
          long next_j = j0 + a_k * j1;
          i0 = i1; j0 = j1;
          i1 = next_i; j1 = next_j;
      }

      long a_adj = 0;
      if (labs(i1) > 0) {
          a_adj = (labs(i0) - A_sieve) / labs(i1) + 1;
      }

      long vec_a_i, vec_a_j, vec_c_i, vec_c_j;
      if (i1 > 0) {
          vec_a_i = i0 + a_adj * i1; vec_a_j = j0 + a_adj * j1;
          vec_c_i = i1;              vec_c_j = j1;
      } else {
          vec_a_i = i1;              vec_a_j = j1;
          vec_c_i = i0 + a_adj * i1; vec_c_j = j0 + a_adj * j1;
      }

      long b0 = -vec_a_i;
      long b1 = A_sieve - vec_c_i;

      long step_a = vec_a_j * A_sieve + vec_a_i;
      long step_c = vec_c_j * A_sieve + vec_c_i;

      uint8_t log_p = (uint8_t)round(log2((double)p_val));
      if (log_p == 0) log_p = 1;

      long curr_j = 1;
      long curr_i = r_val % p_val;
      if (curr_i < 0) curr_i += p_val;

      while (curr_i >= A_sieve / 2) { curr_i += vec_a_i; curr_j += vec_a_j; }
      while (curr_i < -A_sieve / 2) { curr_i += vec_c_i; curr_j += vec_c_j; }

      long x = curr_j * A_sieve + (curr_i + A_sieve / 2);
      long max_x = A_sieve * (J_max + 1);

      /* Tight inner loop: pure C operations updating the native array */
      while (x < max_x) {
          long i_idx = x % A_sieve;
          if (x >= A_sieve && x < max_x) {
              sieve_array[x - A_sieve] += log_p;
          }

          if (i_idx >= b1) {
              x += step_a;
          } else if (i_idx < b0) {
              x += step_c;
          } else {
              x += step_a + step_c; /* The missing diagonal step! */
          }
      }
    }
  }

  /* ------------------------------------------------------------------- */
  /* Post-Sieve Scan & Norm Evaluation */
  /* ------------------------------------------------------------------- */
  if (gequal0(timed_out)) {
    double log_N0 = log2(gtodouble(N0));
    uint8_t threshold = (uint8_t)(log_N0 * 0.75);

    long idx;
    for (idx = 0; idx < total_elements; ++idx) {
      if (sieve_array[idx] >= threshold) {
        long cx = (idx % A_sieve) - A_sieve / 2;
        long cy = (idx / A_sieve) + 1;

        if (cx == 0 && cy == 0) continue;
        if (ugcd((unsigned long)labs(cx), (unsigned long)labs(cy)) != 1) continue;

        /* Evaluate exact norm using polresultant:
           The norm of (b*theta - a) where (a, b) corresponds to (cx, cy)
           is the resultant of f(x) and (cy*x - cx).
        */
        GEN cx_g = stoi(cx);
        GEN cy_g = stoi(cy);
        
        /* Linear form: L(x) = cy*x - cx */
        GEN L = gsub(gmul(cy_g, pol_x(0)), cx_g);
        
        /* Resultant of the defining polynomial f(x) and the linear form L(x).
           f is defined globally in this context as f = x^2 + b*x + c.
        */
        GEN val_g = polresultant0(f, L, -1, 0);
        GEN abs_val = gabs(val_g, 0);

        /* Compare magnitude against N0, but test smoothness on the absolute value */
        if (gcmp(abs_val, N0) >= 0) {
          if (issmooth(abs_val, B)) {
            /* The sieve identifies roots for (cx - cy*theta). 
               Reconstruct as [cx, -cy]~ for GP compatibility.
            */
            GEN cy_neg = stoi(-cy);
            GEN vec = cgetg(3, t_COL);
            gel(vec, 1) = cx_g;
            gel(vec, 2) = cy_neg;

            /* Store using the signed val_g as the map key */
            mapput(m, val_g, vec);

            count = gaddgs(count, 1);
            if (gcmp(count, target) >= 0) break;
          }
        }
      }
    }
  }

  free(sieve_array);
  /* ======================================================================= */
  if (!gequal0(timed_out))
    pari_printf("   [Timeout] Truncated D = %Ps after %Ps min (found %Ps norms)\n", D, mlim, count);
  else
    pari_printf("   [Progress] Finished processing D = %Ps\n", D);
  norms_vec = cgetg(1, t_VEC);
  {
    long l6, l7;
    p2 = cgetg(3, t_MAT);
    for (l7 = 1; l7 <= 2; ++l7)
    {
      gel(p2, l7) = cgetg(1, t_COL);
      for (l6 = 1; l6 <= 0; ++l6)
        gcoeff(p2, l6, l7) = gen_0;
    }
  }
  mat_res = p2;
  if (gcmpgs(count, 0) > 0)
  {
    mat_res = gtomat(m);
    norms_vec = vecsort0(gtrans(gel(mat_res, 1)), NULL, 0);
  }
  p3 = cgetg(5, t_VEC);
  gel(p3, 1) = gcopy(D);
  gel(p3, 2) = gcopy(norms_vec);
  gel(p3, 3) = gcopy(mat_res);
  gel(p3, 4) = gcopy(f);
  return gerepileupto(av, p3);
}

GEN
anon_0(GEN i, GEN D_list, GEN target, GEN B, GEN N0, GEN mlim, GEN A_sieve_in)
{
  return process_one_field(gel(D_list, gtos(i)), target, B, N0, mlim, A_sieve_in);
}

GEN
anon_1(GEN i, GEN num_fields, GEN field_norms)	  /* vec */
{
  GEN loc_max = gen_m1, loc_best_j = gen_0, p1;
  GEN p2;	  /* vec */
  p1 = gaddgs(i, 1);
  {
    GEN j;
    for (j = p1; gcmp(j, num_fields) <= 0; j = gaddgs(j, 1))
    {
      GEN common;
      common = stoi(lg(setintersect(gel(field_norms, gtos(i)), gel(field_norms, gtos(j))))-1);
      if (gcmp(common, loc_max) > 0)
      {
        loc_max = common;
        loc_best_j = j;
      }
    }
  }
  p2 = cgetg(4, t_VEC);
  gel(p2, 1) = /* The last evaluated expression is returned by the thread */
  gcopy(/* The last evaluated expression is returned by the thread */
  loc_max);
  gel(p2, 2) = gcopy(i);
  gel(p2, 3) = gcopy(loc_best_j);
  return p2;
}

/* ========================================================================= */
/* 2. Main Encapsulated Routine */
/* ========================================================================= */

GEN
findpairs(GEN num_fields, GEN target, GEN B, GEN N0, GEN start_D, GEN mlim, GEN t, GEN A_sieve_in, long prec)
{
  pari_sp av = avma; /* <-- FIX: Save stack state */
  GEN D_list, d, field_data;
  GEN p1;	  /* vec */
  GEN field_norms;
  GEN p2;	  /* vec */
  GEN field_maps;
  GEN p3;	  /* vec */
  GEN local_bests, p4;
  GEN p5;	  /* vec */
  GEN max_common = gen_m1, best_pair;
  GEN p6;	  /* vec */
  GEN p7, idxA, idxB, all_common, display_count;
  D_list = cgetg(1, t_VEC);
  d = gsubgs(start_D, 1);
  default0("nbthreads", GENtostr_unquoted(t));
  /* Generate fields dynamically based on configuration */
  while (gcmpsg(glength(D_list), num_fields) < 0)
  {
    d = gaddgs(d, 1);
    if (isfundamental(d))
      D_list = gconcat(D_list, d);
  }
  pari_printf("Step 1: Generated %Ps fundamental discriminants.\n", num_fields);
  pari_printf("Step 2: Enumerating elements in parallel...\n");
  {
    long l8;
    p1 = cgetg(gtos(num_fields)+1, t_VEC);
    for (l8 = 1; gcmpsg(l8, num_fields) <= 0; ++l8)
      gel(p1, l8) = gen_0;
  }
  field_data = p1;
  /* Parallel execution block  */
  {
    parfor_t iter;	  /* parfor */
    parfor_init(&iter, gen_1, num_fields, strtoclosure("anon_0", 6, D_list,
	  target, B, N0, mlim, A_sieve_in));
    {
      GEN i, res, p9;
      while ((p9 = parfor_next(&iter)))
      {
        res = gcopy(gel(p9, 2));
        i = gcopy(gel(p9, 1));
        gel(field_data, gtos(i)) = gcopy(res);
      }
    }
  }
  /* Unpack parallel results dynamically */
  {
    long l10;
    p2 = cgetg(gtos(num_fields)+1, t_VEC);
    for (l10 = 1; gcmpsg(l10, num_fields) <= 0; ++l10)
      gel(p2, l10) = gen_0;
  }
  field_norms = p2;
  {
    long l11;
    p3 = cgetg(gtos(num_fields)+1, t_VEC);
    for (l11 = 1; gcmpsg(l11, num_fields) <= 0; ++l11)
      gel(p3, l11) = gen_0;
  }
  field_maps = p3;
  {
    GEN i;
    for (i = gen_1; gcmp(i, num_fields) <= 0; i = gaddgs(i, 1))
    {
      GEN p12 = gen_0;
      gel(field_norms, gtos(i)) = gcopy(gel(gel(field_data, gtos(i)), 2));
      if (glength(gel(gel(field_data, gtos(i)), 2)) > 0)
        p12 = listinit(gtomap(gel(gel(field_data, gtos(i)), 3))); /* <-- Restored to original */
      else
        p12 = mkmap();
      gel(field_maps, gtos(i)) = p12;
    }
  }
  pari_printf("Successfully collected smooth norms for all %Ps fields.\n", num_fields);
  /* ======================================================================= */
  /* 3. Combinatorial Intersection Pair Matcher (Parallelized) */
  /* ======================================================================= */
  pari_printf("Step 3: Calculating intersections between all pairs in parallel...\n");
  /* Pre-allocate a vector to store the best matches from each thread */
  p4 = gsubgs(num_fields, 1);
  {
    long l13;
    p5 = cgetg(gtos(p4)+1, t_VEC);
    for (l13 = 1; gcmpsg(l13, p4) <= 0; ++l13)
      gel(p5, l13) = gen_0;
  }
  local_bests = p5;
  /* Parallel execution: map the inner loop to available threads */
  {
    parfor_t iter;	  /* parfor */
    parfor_init(&iter, gen_1, gsubgs(num_fields, 1), strtoclosure("anon_1", 2, num_fields, field_norms));
    {
      GEN res, p14;
      while ((p14 = parfor_next(&iter)))
      {
        res = gcopy(gel(p14, 2));
        gel(local_bests, gtos(gel(res, 2))) = gcopy(res);
      }
    }
  }
  p6 = cgetg(3, t_VEC);
  gel(p6, 1) = gen_0;
  gel(p6, 2) = gen_0;
  best_pair = p6;
  p7 = gsubgs(num_fields, 1);
  {
    GEN k;
    for (k = gen_1; gcmp(k, p7) <= 0; k = gaddgs(k, 1))
    {
      if (gcmp(gel(gel(local_bests, gtos(k)), 1), max_common) > 0)
      {
        GEN p15;	  /* vec */
        max_common = gcopy(gel(gel(local_bests, gtos(k)), 1));
        p15 = cgetg(3, t_VEC);
        gel(p15, 1) = gcopy(gel(gel(local_bests, gtos(k)), 2));
        gel(p15, 2) = gcopy(gel(gel(local_bests, gtos(k)), 3));
        best_pair = p15;
      }
    }
  }
  idxA = gcopy(gel(best_pair, 1));
  idxB = gcopy(gel(best_pair, 2));
  all_common = setintersect(gel(field_norms, gtos(idxA)), gel(field_norms, gtos(idxB)));
  display_count = stoi(minss(10, glength(all_common)));
  GEN dataA = gel(field_data, gtos(idxA));
  GEN dataB = gel(field_data, gtos(idxB));
  long A_val = itos(A_sieve_in);
  /* ======================================================================= */
  /* 4. Output Presentation Block */
  /* ======================================================================= */
  pari_printf("\n=================== RESULTS ===================\n");
  pari_printf("Sieve Region (a, b): [-%ld, %ld) x [1, %ld]\n",
              A_val / 2, A_val / 2, 2 * A_val);
  pari_printf("Field A: Index %Ps | Discriminant = %Ps\n", idxA, gel(D_list, gtos(idxA)));
  pari_printf("   -> Polynomial: %Ps\n", gel(dataA, 4));
  pari_printf("   -> Total Smooth Norms Found: %ld\n", glength(gel(field_norms, gtos(idxA))));
  pari_printf("Field B: Index %Ps | Discriminant = %Ps\n", idxB, gel(D_list, gtos(idxB)));
  pari_printf("   -> Polynomial: %Ps\n", gel(dataB, 4));
  pari_printf("   -> Total Smooth Norms Found: %ld\n", glength(gel(field_norms, gtos(idxB))));
  pari_printf("Shared Smooth Norms in Common: %Ps\n", max_common);
  pari_printf("===============================================\n");
  if (glength(all_common) == 0)
  {
    pari_printf("\nNo shared smooth norms found among the processed pairs.\n");
    return gerepilecopy(av, cgetg(1, t_VEC)); /* <-- FIX: Safe stack return */
  }
  pari_printf("\nTop %Ps Shared Norm Elements:\n", display_count);
  {
    GEN k;
    for (k = gen_1; gcmp(k, display_count) <= 0; k = gaddgs(k, 1))
    {
      GEN n_val, elemA_vec, elemB_vec, elemA_alg, elemB_alg;
      //n_val = gcopy(gel(all_common, gtos(k)));
	  n_val = gcopy(gel(all_common, glength(all_common) - gtos(k) + 1));
      elemA_vec = mapget(gel(field_maps, gtos(idxA)), n_val);
      elemB_vec = mapget(gel(field_maps, gtos(idxB)), n_val);
      elemA_alg = gadd(gel(elemA_vec, 1), gmul(gel(elemA_vec, 2), pol_x(0)));
      elemB_alg = gadd(gel(elemB_vec, 1), gmul(gel(elemB_vec, 2), pol_x(0)));
      /* Formatting string concatenation appropriately */
      pari_printf("Field A element %Ps, Norm %Ps\n", GENtoGENstr(elemA_alg), fact_str(n_val));
      pari_printf("Field B element %Ps, Norm %Ps\n", GENtoGENstr(elemB_alg), fact_str(n_val));
      pari_printf("\n");
    }
  }

  /* ======================================================================= */
  /* 5. Pack and Return R1, R2, f1, f2 */
  /* ======================================================================= */
  GEN R1 = cgetg(glength(all_common) + 1, t_VEC);
  GEN R2 = cgetg(glength(all_common) + 1, t_VEC);
  long k_idx;
  for (k_idx = 1; k_idx <= glength(all_common); ++k_idx)
  {
    GEN n_val = gel(all_common, k_idx);
    GEN elemA_vec = mapget(gel(field_maps, gtos(idxA)), n_val);
    GEN elemB_vec = mapget(gel(field_maps, gtos(idxB)), n_val);
    
    /* Reconstruct algebraic form: cx - cy*x */
    GEN elemA_alg = gadd(gel(elemA_vec, 1), gmul(gel(elemA_vec, 2), pol_x(0)));
    GEN elemB_alg = gadd(gel(elemB_vec, 1), gmul(gel(elemB_vec, 2), pol_x(0)));
    
    gel(R1, k_idx) = elemA_alg;
    gel(R2, k_idx) = elemB_alg;
  }

  GEN f1 = gel(dataA, 4);
  GEN f2 = gel(dataB, 4);

  GEN res = mkvec4(R1, R2, f1, f2);
  
  /* <-- FIX: Deep copy the return vector and reset the stack cleanly */
  return gerepilecopy(av, res); 
}

