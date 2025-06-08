const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;

const utils = @import("utils.zig");

/// Instance of a Number Theoretic Transform based on 2n-th root of unity.
/// This implementation is using the polynomial quotient ring Z_q[x]/(x^n + 1)
/// where q is the coefficient modulus and n is the degree of the cyclotomic
/// polynomial.
/// The NTT algorithm is based on the paper "Low-Cost and Area-Efficient FPGA
/// Implementations of Lattice-Based Cryptography" by Aysu et al.
// See: https://schaumont.dyn.wpi.edu/schaum/pdf/papers/2013hostb.pdf
pub const NTT = struct {
    const Self = @This();

    /// Coefficient modulus.
    q: i64,
    /// Degree of cyclotomic polynomial.
    n: i64,
    /// Inverse of degree of cyclotomic polynomial.
    n_inverse: i64,
    /// Powers of psi from psi^0 to psi^(n - 1).
    psi_powers: []i64,
    /// Powers of psi^-1 = ipsi from ipsi^0 to ipsi^(n - 1).
    psi_inverse_powers: []i64,
    /// Allocator used for internal memory allocations.
    allocator: Allocator,

    /// Initializes a new NTT instance with the given coefficient modulus and
    /// cyclotomic polynomial degree (which MUST be a power of 2).
    pub fn init(allocator: Allocator, q: i64, n: i64) !Self {
        // Check if degree of cyclotomic polynomial is a power of 2.
        if (!utils.isPowerOfTwo(n)) {
            return error.InvalidDegree;
        }

        // Compute psi, the 2n-th root of unity and its inverse.
        const psi = try utils.findRootOfUnity(allocator, 2 * n, q);
        const psi_inverse = try utils.modInv(psi, q);

        // Compute the inverse of the degree of the cyclotomic polynomial.
        const n_inverse = try utils.modInv(n, q);

        // Compute powers of psi.
        var powers = try std.ArrayList(i64).initCapacity(allocator, @intCast(n));
        defer powers.deinit();
        try powers.append(1);
        for (1..@intCast(n)) |_| {
            const power = @mod(powers.getLast() * psi, q);
            try powers.append(power);
        }
        const psi_powers = try powers.toOwnedSlice();

        // Compute powers of psi^-1.
        var inverse_powers = try std.ArrayList(i64).initCapacity(allocator, @intCast(n));
        defer inverse_powers.deinit();
        try inverse_powers.append(1);
        for (1..@intCast(n)) |_| {
            const power = @mod(inverse_powers.getLast() * psi_inverse, q);
            try inverse_powers.append(power);
        }
        const psi_inverse_powers = try inverse_powers.toOwnedSlice();

        return NTT{
            .q = q,
            .n = n,
            .n_inverse = n_inverse,
            .psi_powers = psi_powers,
            .psi_inverse_powers = psi_inverse_powers,
            .allocator = allocator,
        };
    }

    /// Release all allocated memory.
    pub fn deinit(self: Self) void {
        self.allocator.free(self.psi_powers);
        self.allocator.free(self.psi_inverse_powers);
    }

    /// Runs a forward pass of NTT with the given coefficients.
    /// The caller owns the returned memory.
    pub fn fwd(self: Self, coefficients: []const i64) ![]const i64 {
        // Length of coefficients must equal the degree of the cyclotomic polynomial.
        if (coefficients.len != self.n) {
            return error.InvalidLength;
        }

        var input = try self.allocator.alloc(i64, coefficients.len);
        defer self.allocator.free(input);

        for (0..coefficients.len) |i| {
            input[i] = @mod(coefficients[i] * self.psi_powers[i], self.q);
        }

        return self.ntt(input, self.psi_powers);
    }

    /// Runs INTT (Inverse NTT) with the given coefficients.
    /// The caller owns the returned memory.
    pub fn inv(self: Self, coefficients: []const i64) ![]const i64 {
        // Length of coefficients must equal the degree of the cyclotomic polynomial.
        if (coefficients.len != self.n) {
            return error.InvalidLength;
        }

        const unscaled_result = try self.ntt(coefficients, self.psi_inverse_powers);
        defer self.allocator.free(unscaled_result);

        var result = try self.allocator.alloc(i64, coefficients.len);

        for (0..coefficients.len) |i| {
            result[i] = @mod(unscaled_result[i] * self.psi_inverse_powers[i] * self.n_inverse, self.q);
        }

        return result;
    }

    /// Runs an iterative version of NTT with the given coefficients and twiddles.
    /// The caller owns the returned memory.
    fn ntt(self: Self, coefficients: []const i64, twiddles: []const i64) ![]const i64 {
        // Length of coefficients and twiddles must be the same.
        if (coefficients.len != twiddles.len) {
            return error.InvalidLength;
        }

        const log2_n = std.math.log2_int(usize, @intCast(self.n));
        const reversed = try utils.bitReverseSlice(self.allocator, coefficients);
        defer self.allocator.free(reversed);

        var result = try self.allocator.dupe(i64, reversed);

        for (0..log2_n) |i| {
            var temp_twiddle: i64 = 1;
            var final_twiddle: i64 = 1;
            const twiddle = twiddles[i + 1];

            const in_1 = @as(usize, 1) << @intCast(i); // 2^i
            const in_2 = @as(usize, 1) << @intCast(i + 1); // 2^(i + 1)
            const in_3 = @as(usize, @intCast(self.n)) / in_2; // n / (2^(i + 1))

            for (0..in_1) |j| {
                for (0..in_3) |t| {
                    const index_1 = (t * in_2) + j; // (t * 2^(i + 1)) + j
                    const index_2 = index_1 + in_1; // (t * 2^(i + 1)) + j + 2^i

                    const c = result[index_1];
                    const d = result[index_2];

                    const butterfly_plus = @mod(c + (final_twiddle * d), self.q);
                    const butterfly_minus = @mod(c - (final_twiddle * d), self.q);

                    result[index_1] = butterfly_plus;
                    result[index_2] = butterfly_minus;

                    temp_twiddle = @mod(temp_twiddle * twiddle, self.q);
                }
                final_twiddle = temp_twiddle;
            }
        }

        return result;
    }
};

test "ntt - init" {
    const allocator = testing.allocator;

    {
        const q = 7681;
        const n = 5;

        const expected = error.InvalidDegree;
        const result = NTT.init(allocator, q, n);

        try testing.expectError(expected, result);
    }
}

test "ntt - fwd" {
    const allocator = testing.allocator;

    {
        const q = 7681;
        const n = 4;

        const ntt = try NTT.init(allocator, q, n);
        defer ntt.deinit();

        const coefficients = [_]i64{ 1, 2, 3, 4 };
        const expected = [_]i64{ 1467, 2807, 3471, 7621 };

        const result = try ntt.fwd(&coefficients);
        defer allocator.free(result);

        try testing.expectEqualSlices(i64, &expected, result);
    }

    {
        const q = 7681;
        const n = 4;

        const ntt = try NTT.init(allocator, q, n);
        defer ntt.deinit();

        const coefficients = [_]i64{ 5, 6, 7, 8 };
        const expected = [_]i64{ 2489, 7489, 6478, 6607 };

        const result = try ntt.fwd(&coefficients);
        defer allocator.free(result);

        try testing.expectEqualSlices(i64, &expected, result);
    }

    {
        const q = 7681;
        const n = 4;

        const ntt = try NTT.init(allocator, q, n);
        defer ntt.deinit();

        const coefficients = [_]i64{ 1, 2, 3 };

        const expected = error.InvalidLength;
        const result = ntt.fwd(&coefficients);

        try testing.expectError(expected, result);
    }
}

test "ntt - inv" {
    const allocator = testing.allocator;

    {
        const q = 7681;
        const n = 4;

        const ntt = try NTT.init(allocator, q, n);
        defer ntt.deinit();

        const coefficients = [_]i64{ 1467, 2807, 3471, 7621 };
        const expected = [_]i64{ 1, 2, 3, 4 };

        const result = try ntt.inv(&coefficients);
        defer allocator.free(result);

        try testing.expectEqualSlices(i64, &expected, result);
    }

    {
        const q = 7681;
        const n = 4;

        const ntt = try NTT.init(allocator, q, n);
        defer ntt.deinit();

        const coefficients = [_]i64{ 2489, 7489, 6478, 6607 };
        const expected = [_]i64{ 5, 6, 7, 8 };

        const result = try ntt.inv(&coefficients);
        defer allocator.free(result);

        try testing.expectEqualSlices(i64, &expected, result);
    }

    {
        const q = 7681;
        const n = 4;

        const ntt = try NTT.init(allocator, q, n);
        defer ntt.deinit();

        const coefficients = [_]i64{ 1, 2, 3 };

        const expected = error.InvalidLength;
        const result = ntt.inv(&coefficients);

        try testing.expectError(expected, result);
    }
}

test "ntt - ntt" {
    const allocator = testing.allocator;

    {
        const q = 7681;
        const n = 4;

        const ntt = try NTT.init(allocator, q, n);
        defer ntt.deinit();

        const twiddles = ntt.psi_powers;
        const coefficients = [_]i64{ 1, 2, 3, 4 };
        const expected = [_]i64{ 10, 913, 7679, 6764 };

        const result = try ntt.ntt(&coefficients, twiddles);
        defer allocator.free(result);

        try testing.expectEqualSlices(i64, &expected, result);
    }

    {
        const q = 7681;
        const n = 4;

        const ntt = try NTT.init(allocator, q, n);
        defer ntt.deinit();

        const twiddles = ntt.psi_powers;
        const coefficients = [_]i64{ 1, 2 };

        const expected = error.InvalidLength;
        const result = ntt.ntt(&coefficients, twiddles);

        try testing.expectError(expected, result);
    }
}
