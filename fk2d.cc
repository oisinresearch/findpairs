#include <iostream>
#include <cmath>
#include <cstdint>
#include <cstdlib>

using std::abs;
using std::round;
using std::log2;
using std::cout;

struct Point {
    int32_t i;
    int32_t j;
};

void franke_kleinjung_sieve(int32_t p, int32_t r, int32_t A, uint8_t* sieve_array) {
    // --- Step 1: Basis Generation (Proposition 1) ---
    int32_t i0 = -p, j0 = 0;
    int32_t i1 = r,  j1 = 1;

    while (true) {
        if (abs(i1) == 1 || abs(i1) < A) {
            break;
        }
        int32_t a_k = abs(i0) / abs(i1);
        int32_t next_i = i0 + a_k * i1;
        int32_t next_j = j0 + a_k * j1;

        i0 = i1; j0 = j1;
        i1 = next_i; j1 = next_j;
    }

    int32_t a_adj = 0;
    if (abs(i1) > 0) {
        a_adj = (abs(i0) - A) / abs(i1) + 1;
    }

    Point vec_a, vec_c;
    if (i1 > 0) {
        vec_a = { i0 + a_adj * i1, j0 + a_adj * j1 };
        vec_c = { i1, j1 };
    } else {
        vec_a = { i1, j1 };
        vec_c = { i0 + a_adj * i1, j0 + a_adj * j1 };
    }

    // --- Step 2: Setup Loop Parameters (Remark 1) ---
    int32_t b0 = -vec_a.i; 
    int32_t b1 = A - vec_c.i; 

    int32_t step_a = vec_a.j * A + vec_a.i;
    int32_t step_c = vec_c.j * A + vec_c.i;

    uint8_t log_p = static_cast<uint8_t>(round(log2(p)));
    if (log_p == 0) log_p = 1;

    // --- Step 3: Fixed Initial Alignment ---
    int32_t J_max = 2 * A;
    int32_t max_x = A * (J_max + 1);

    // Calculate a valid baseline lattice point on row j = 1
    int32_t curr_j = 1;
    int32_t curr_i = r % p;
    if (curr_i < 0) curr_i += p;

    // Walk the point into the horizontal bounds [-A/2, A/2 - 1]
    while (curr_i >= A / 2) {
        curr_i += vec_a.i; 
        curr_j += vec_a.j; 
    }
    while (curr_i < -A / 2) {
        curr_i += vec_c.i; 
        curr_j += vec_c.j;
    }

    // Map safely into the flat 1D index tracking variable 'x'
    int32_t x = curr_j * A + (curr_i + A / 2);

    // --- Step 4: Inner Sieving Loop ---
    int32_t hits = 0;
    while (x < max_x) {
        // FIXED: Added missing int32_t type declaration here
        int32_t i = x % A;

        // Sieve bounds check: 1 <= j <= 2A
        if (x >= A && x <= A * J_max) {
            sieve_array[x - A] += log_p; 
            hits++;
        }

        if (i >= b1) x += step_a;
        if (i < b0)  x += step_c;
    }
    cout << "Lattice points found: " << hits << " (Expected ~518-519)\n";
}

int main() {
    int32_t A = 512;
    int32_t p = 1009;
    int32_t r = 921;

    int32_t total_elements = A * (2 * A);
    uint8_t* sieve_array = new uint8_t[total_elements]();

    franke_kleinjung_sieve(p, r, A, sieve_array);

    delete[] sieve_array;
    return 0;
}

