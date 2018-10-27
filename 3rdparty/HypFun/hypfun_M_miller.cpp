// 
// Computes the hypergeometric function \mathbf{M}(a+k; c+k; z), for k a large
// positive integer, using \mathbf{M}(a; c; z) and Miller's algorithm, with
// n a number much larger than k.
// 
// Copyright John W. Pearson 2014

Tp
hypfun_M_recplusplus_miller(a, c, z, n, k)
{
  // Initialise f and v in Miller's algorithm
  f = zeros(k+1,1);
  v = zeros(n+1,1);
  a1 = zeros(n-1,1);
  c1 = zeros(n-1,1);

  // Define coefficients in recurrence relation
  for (auto i = 1; i <= n - 1; ++i)
    {
      a1[i] = -1 / ((a + i) * z);
      c1[i] = (c + i - z - 1) / ((a + i) * z);
    }

  // Input minimal solution with k=0
  f1 = hypergeom(a, c, z) / gamfun(c);
  v[end] = 0;
  v[end-1] = 1;

  // Compute recurrence backwards
  for (auto i = 2; i <= n; ++i)
    v[n + 1 - i] = -(v[n + 3 - i] + c1[n + 1 - i] * v[n + 2 - i]) / a1[n + 1 - i];

  // Apply last line of Miller's algorithm
  for (auto i = 1; i <= k+1; ++i)
    f[i] = f1 / v[1] * v[i] * gamfun(c + i - 1);

  // Return solution
  return f.back() / gamfun(c + k);
}
