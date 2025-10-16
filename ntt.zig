const std = @import("std");
const Allocator = std.mem.Allocator;
const testing = std.testing;
const math = std.math;

const utils = @import("src/utils.zig");

const Mode = enum {
    fwd,
    inv,
};

fn ntt(allocator: Allocator, coefficients: []const i64, q: i64, mode: Mode) ![]const i64 {
    const n = coefficients.len;

    // Ensure that length of coefficients is a power of 2.
    if (!utils.isPowerOfTwo(@intCast(n))) {
        return error.InvalidLength;
    }

    // Compute psi, the 2n-th root of unity as well as its inverse.
    const psi = try utils.findRootOfUnity(allocator, 2 * @as(i64, @intCast(n)), q);
    const psi_inv = try utils.modInv(psi, q);

    // Compute the inverse of the degree of the cyclotomic polynomial.
    const n_inv = try utils.modInv(@as(i64, @intCast(n)), q);

    // Compute powers of psi as well as powers of psy^-1.
    var powers = try std.ArrayList(i64).initCapacity(allocator, n);
    var powers_inv = try std.ArrayList(i64).initCapacity(allocator, n);
    defer powers.deinit();
    defer powers_inv.deinit();

    try powers.append(1);
    try powers_inv.append(1);

    for (0..n) |_| {
        const power = @mod(powers.getLast() * psi, q);
        const power_inv = @mod(powers_inv.getLast() * psi_inv, q);
        try powers.append(power);
        try powers_inv.append(power_inv);
    }

    const psi_powers = try powers.toOwnedSlice();
    const psi_powers_inv = try powers_inv.toOwnedSlice();
    defer allocator.free(psi_powers);
    defer allocator.free(psi_powers_inv);

    // Prepare input for NTT.
    const input = switch (mode) {
        .fwd => blk: {
            var result = try allocator.alloc(i64, n);
            for (0..n) |i| {
                result[i] = @mod(coefficients[i] * psi_powers[i], q);
            }
            break :blk result;
        },
        .inv => try allocator.dupe(i64, coefficients),
    };
    defer allocator.free(input);

    // Run iterative NTT.
    const log2_n = math.log2_int(usize, @intCast(n));

    const reversed = try utils.bitReverseSlice(i64, allocator, input);
    defer allocator.free(reversed);

    var result = try allocator.dupe(i64, reversed);

    for (0..log2_n) |i| {
        const in_1 = @as(usize, 1) << @intCast(i); // 2^i
        const in_2 = @as(usize, 1) << @intCast(i + 1); // 2^(i + 1)
        const in_3 = @as(usize, n) >> @intCast(i + 1); // n >> (i + 1) = n / 2^(i + 1)

        for (0..in_1) |j| {
            for (0..in_3) |t| {
                const index_even = (t * in_2) + j; // (t * 2^(i + 1)) + j
                const index_odd = index_even + in_1; // (t * 2^(i + 1)) + j + 2^i

                const c = result[index_even];
                const d = result[index_odd];

                const rou_index = @as(usize, j) << @intCast(1 + log2_n - (i + 1)); // j << (1 + log2(n) - (i + 1))

                const psi_factor = switch (mode) {
                    .fwd => @mod(psi_powers[rou_index] * d, q),
                    .inv => @mod(psi_powers_inv[rou_index] * d, q),
                };

                const butterfly_plus = @mod(c + psi_factor, q);
                const butterfly_minus = @mod(c - psi_factor, q);

                result[index_even] = butterfly_plus;
                result[index_odd] = butterfly_minus;
            }
        }
    }

    // Rescale if iNTT is run.
    if (mode == Mode.inv) {
        for (0..n) |i| {
            result[i] = @mod(result[i] * psi_powers_inv[i] * n_inv, q);
        }
    }

    return result;
}

test "ntt" {
    const allocator = testing.allocator;

    const q = 7681;

    {
        const coefficients = [_]i64{ 1, 2, 3, 4 };
        const expected = [_]i64{ 1467, 2807, 3471, 7621 };

        const result = try ntt(allocator, &coefficients, q, Mode.fwd);
        defer allocator.free(result);

        try testing.expectEqualSlices(i64, &expected, result);
    }

    {
        const coefficients = [_]i64{ 5, 6, 7, 8 };
        const expected = [_]i64{ 2489, 7489, 6478, 6607 };

        const result = try ntt(allocator, &coefficients, q, Mode.fwd);
        defer allocator.free(result);

        try testing.expectEqualSlices(i64, &expected, result);
    }

    {
        const coefficients = [_]i64{ 1467, 2807, 3471, 7621 };
        const expected = [_]i64{ 1, 2, 3, 4 };

        const result = try ntt(allocator, &coefficients, q, Mode.inv);
        defer allocator.free(result);

        try testing.expectEqualSlices(i64, &expected, result);
    }

    {
        const coefficients = [_]i64{ 2489, 7489, 6478, 6607 };
        const expected = [_]i64{ 5, 6, 7, 8 };

        const result = try ntt(allocator, &coefficients, q, Mode.inv);
        defer allocator.free(result);

        try testing.expectEqualSlices(i64, &expected, result);
    }

    {
        const coefficients_1 = [_]i64{ 1, 2, 3, 4 };
        const coefficients_2 = [_]i64{ 5, 6, 7, 8 };
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

        const n = coefficients_1.len;
        const fwd_1 = try ntt(allocator, &coefficients_1, q, Mode.fwd);
        const fwd_2 = try ntt(allocator, &coefficients_2, q, Mode.fwd);
        defer allocator.free(fwd_1);
        defer allocator.free(fwd_2);

        var interim = [_]i64{0} ** n;
        for (0..n) |i| {
            interim[i] = fwd_1[i] * fwd_2[i];
        }

        const result = try ntt(allocator, &interim, q, Mode.inv);
        defer allocator.free(result);

        try testing.expectEqualSlices(i64, &expected, result);
    }
}
