//! This module implements functionalities to run Number Theoretic Transform
//! computations.

const std = @import("std");
const testing = std.testing;

const utils = @import("utils.zig");

/// Number Theoretic Transform based on the 2n-th root of unity.
/// This implementation is using the polynomial quotient ring Z_q[x]/(x^n + 1)
/// where q is the coefficient modulus and n is the degree of the cyclotomic
/// polynomial.
/// The NTT algorithm is based on the paper "Low-Cost and Area-Efficient FPGA
/// Implementations of Lattice-Based Cryptography" by Aysu et al.
// See: https://schaumont.dyn.wpi.edu/schaum/pdf/papers/2013hostb.pdf
pub fn NTT(comptime q: i64, comptime n: i64) type {
    return struct {
        const Self = @This();

        /// Coefficient modulus.
        comptime q: i64 = q,
        /// Degree of cyclotomic polynomial.
        comptime n: i64 = n,
        /// Inverse of degree of cyclotomic polynomial.
        n_inverse: i64,
        /// Powers of psi from psi^0 to psi^(n - 1).
        psi_powers: [n]i64,
        /// Powers of psi^-1 = ipsi from ipsi^0 to ipsi^(n - 1).
        psi_inverse_powers: [n]i64,

        /// Initializes a new NTT instance. Ensures that the cyclotomic
        /// polynomial is a power of 2.
        pub fn init() !Self {
            // Check if degree of cyclotomic polynomial is a power of 2.
            if (!utils.isPowerOfTwo(n)) {
                return error.InvalidDegree;
            }

            // Compute psi, the 2n-th root of unity and its inverse.
            const psi = try utils.findRootOfUnity(2 * n, q);
            const psi_inverse = try utils.modInv(psi, q);

            // Compute the inverse of the degree of the cyclotomic polynomial.
            const n_inverse = try utils.modInv(n, q);

            // Compute powers of psi as well as powers of psi^-1.
            var psi_powers = [_]i64{1} ** n;
            var psi_inverse_powers = [_]i64{1} ** n;

            for (1..n) |i| {
                psi_powers[i] = @mod(psi_powers[i - 1] * psi, q);
                psi_inverse_powers[i] = @mod(psi_inverse_powers[i - 1] * psi_inverse, q);
            }

            return Self{
                .n_inverse = n_inverse,
                .psi_powers = psi_powers,
                .psi_inverse_powers = psi_inverse_powers,
            };
        }

        /// Runs a forward pass of NTT with the given coefficients.
        /// Note that the slice of coefficients is mutated in-place.
        pub fn fwd(self: Self, coefficients: []i64) ![]i64 {
            // Length of coefficients must equal the degree of the cyclotomic polynomial.
            if (coefficients.len != self.n) {
                return error.InvalidLength;
            }

            for (0..coefficients.len) |i| {
                coefficients[i] = @mod(coefficients[i] * self.psi_powers[i], self.q);
            }

            return self.ntt(coefficients, self.psi_powers);
        }

        /// Runs iNTT (Inverse NTT) with the given coefficients.
        /// Note that the slice of coefficients is mutated in-place.
        pub fn inv(self: Self, coefficients: []i64) ![]i64 {
            // Length of coefficients must equal the degree of the cyclotomic polynomial.
            if (coefficients.len != self.n) {
                return error.InvalidLength;
            }

            var result = try self.ntt(coefficients, self.psi_inverse_powers);

            for (0..coefficients.len) |i| {
                result[i] = @mod(result[i] * self.psi_inverse_powers[i] * self.n_inverse, self.q);
            }

            return result;
        }

        /// Runs an iterative version of NTT with the given coefficients and twiddles
        /// which are the powers of the roots of unity (i.e. powers of psi / powers
        /// of psi^-1).
        /// Note that the slice of coefficients is mutated in-place.
        fn ntt(self: Self, coefficients: []i64, twiddles: [n]i64) ![]i64 {
            // Length of coefficients and twiddles must be the same.
            if (coefficients.len != twiddles.len) {
                return error.InvalidLength;
            }

            const log2_n = std.math.log2_int(usize, @intCast(self.n));
            var result = try utils.bitReverseSlice(i64, self.n, coefficients);

            for (0..log2_n) |i| {
                const in_1 = @as(usize, 1) << @intCast(i); // 2^i
                const in_2 = @as(usize, 1) << @intCast(i + 1); // 2^(i + 1)
                const in_3 = @as(usize, @intCast(self.n)) >> @intCast(i + 1); // n >> (i + 1) = n / 2^(i + 1)

                for (0..in_1) |j| {
                    for (0..in_3) |t| {
                        const index_even = (t * in_2) + j; // (t * 2^(i + 1)) + j
                        const index_odd = index_even + in_1; // (t * 2^(i + 1)) + j + 2^i

                        const c = result[index_even];
                        const d = result[index_odd];

                        const twiddle_index = @as(usize, j) << @intCast(1 + log2_n - (i + 1)); // j << (1 + log2(n) - (i + 1))
                        const twiddle_factor = @mod(twiddles[twiddle_index] * d, self.q);

                        const butterfly_plus = @mod(c + twiddle_factor, self.q);
                        const butterfly_minus = @mod(c - twiddle_factor, self.q);

                        result[index_even] = butterfly_plus;
                        result[index_odd] = butterfly_minus;
                    }
                }
            }

            return result;
        }
    };
}

test "ntt - init" {
    {
        const q = 7681;
        const n = 5;

        const expected = error.InvalidDegree;
        const result = NTT(q, n).init();

        try testing.expectError(expected, result);
    }
}

test "ntt - fwd" {
    {
        const q = 7681;
        const n = 4;

        const ntt = try NTT(q, n).init();

        var coefficients = [_]i64{ 1, 2, 3, 4 };
        const expected = [_]i64{ 1467, 2807, 3471, 7621 };

        const result = try ntt.fwd(&coefficients);

        try testing.expectEqualSlices(i64, &expected, result);
    }

    {
        const q = 7681;
        const n = 4;

        const ntt = try NTT(q, n).init();

        var coefficients = [_]i64{ 5, 6, 7, 8 };
        const expected = [_]i64{ 2489, 7489, 6478, 6607 };

        const result = try ntt.fwd(&coefficients);

        try testing.expectEqualSlices(i64, &expected, result);
    }

    {
        const q = 7681;
        const n = 4;

        const ntt = try NTT(q, n).init();

        var coefficients = [_]i64{ 1, 2, 3 };

        const expected = error.InvalidLength;
        const result = ntt.fwd(&coefficients);

        try testing.expectError(expected, result);
    }
}

test "ntt - inv" {
    {
        const q = 7681;
        const n = 4;

        const ntt = try NTT(q, n).init();

        var coefficients = [_]i64{ 1467, 2807, 3471, 7621 };
        const expected = [_]i64{ 1, 2, 3, 4 };

        const result = try ntt.inv(&coefficients);

        try testing.expectEqualSlices(i64, &expected, result);
    }

    {
        const q = 7681;
        const n = 4;

        const ntt = try NTT(q, n).init();

        var coefficients = [_]i64{ 2489, 7489, 6478, 6607 };
        const expected = [_]i64{ 5, 6, 7, 8 };

        const result = try ntt.inv(&coefficients);

        try testing.expectEqualSlices(i64, &expected, result);
    }

    {
        const q = 7681;
        const n = 4;

        const ntt = try NTT(q, n).init();

        var coefficients = [_]i64{ 1, 2, 3 };

        const expected = error.InvalidLength;
        const result = ntt.inv(&coefficients);

        try testing.expectError(expected, result);
    }
}

test "ntt - ntt" {
    {
        const q = 7681;
        const n = 4;

        const ntt = try NTT(q, n).init();

        const twiddles = ntt.psi_powers;
        var coefficients = [_]i64{ 1, 2, 3, 4 };
        const expected = [_]i64{ 10, 913, 7679, 6764 };

        const result = try ntt.ntt(&coefficients, twiddles);

        try testing.expectEqualSlices(i64, &expected, result);
    }

    {
        const q = 7681;
        const n = 4;

        const ntt = try NTT(q, n).init();

        const twiddles = ntt.psi_powers;
        var coefficients = [_]i64{ 1, 2 };

        const expected = error.InvalidLength;
        const result = ntt.ntt(&coefficients, twiddles);

        try testing.expectError(expected, result);
    }
}

test "ntt - convolution" {
    {
        const q = 7681;
        const n = 4;

        const ntt = try NTT(q, n).init();

        var coefficients_1 = [_]i64{ 1, 2, 3, 4 };
        var coefficients_2 = [_]i64{ 5, 6, 7, 8 };
        // Note that x^n + 1 = x^4 + 1 which means that x^4 = -1
        //   (1 + 2x + 3x^2 + 4x^3) * (5 + 6x + 7x^2 + 8x^3)
        // = 5 + 6x + 7x^2 + 8x^3 +
        //   10x + 12x^2 + 14x^3 + 16x^4 +
        //   15x^2 + 18x^3 + 21x^4 + 24x^5 +
        //   20x^3 + 24x^4 + 28x^5 + 32x^6
        // = 5 + 16x + 34x^2 + 60x^3 + 61x^4 + 52x^5 + 32x^6
        // = 5 + 16x + 34x^2 + 60x^3 + 61(x^4) + x(52x^4) + x^2(32x^4)
        // = 5 + 16x + 34x^2 + 60x^3 + (61 * -1) + x(52 * -1) + x^2(32 * -1)
        // = 5 + 16x + 34x^2 + 60x^3 - 61 - 52x - 32x^2
        // = -56 - 36x + 2x^2 + 60x^3
        const expected = [_]i64{ 7625, 7645, 2, 60 }; // = { -56, -36, 2, 60 }

        const fwd_1 = try ntt.fwd(&coefficients_1);
        const fwd_2 = try ntt.fwd(&coefficients_2);

        var interim = [_]i64{0} ** n;
        for (0..n) |i| {
            interim[i] = fwd_1[i] * fwd_2[i];
        }

        const result = try ntt.inv(&interim);

        try testing.expectEqualSlices(i64, &expected, result);
    }
}
