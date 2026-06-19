install("init_findpairs","v","init_findpairs","./findpairs.gp.so");
install("issmooth","lD0,G,D0,G,","issmooth","./findpairs.gp.so");
install("fact_str","D0,G,","fact_str","./findpairs.gp.so");
install("process_one_field","D0,G,D0,G,D0,G,D0,G,D0,G,D0,G,","process_one_field","./findpairs.gp.so");
install("findpairs","vD0,G,D0,G,D0,G,D0,G,D0,G,D0,G,D0,G,D0,G,p","findpairs","./findpairs.gp.so");

\\ Use addhelp to attach your description to the symbol
addhelp(findpairs, "findpairs(num_fields, target, B, N0, start_D, mlim, threads, A):\nSearches for matching smooth norms across fundamental discriminants.\n  num_fields : Number of fundamental discriminants to process\n  target     : Target number of smooth norms to find per field\n  B          : Smoothness bound (max prime factor)\n  N0         : Minimum acceptable norm magnitude\n  start_D    : Starting point for fundamental discriminants\n  mlim       : Timeout in minutes per field\n  threads    : Number of parallel worker threads\n  A          : Sieve region parameter (defines [-A/2, A/2) x [1, 2A])");

install("anon_0","D0,G,GGGGGG","anon_0","./findpairs.gp.so");
install("anon_1","D0,G,GG","anon_1","./findpairs.gp.so");

default(parisizemax, 4*1024^3);
default(threadsize, "20M");  \\ Increased thread size to prevent worker thread crashes!
