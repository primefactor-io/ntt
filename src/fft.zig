const std = @import("std");
const math = std.math;
const testing = std.testing;
const Allocator = std.mem.Allocator;
const pi = math.pi;
const Complex = math.Complex;

const utils = @import("utils.zig");
const expectEqualComplexSlices = @import("testing.zig").expectEqualComplexSlices;

/// Instance of a Fast Fourier Transform based on n-th root of unity.
/// The FFT algorithm is based on the paper "Low-Cost and Area-Efficient FPGA
/// Implementations of Lattice-Based Cryptography" by Aysu et al.
// See: https://schaumont.dyn.wpi.edu/schaum/pdf/papers/2013hostb.pdf
pub const FFT = struct {
    const Self = @This();

    /// Length of FFT vector.
    n: i64,
    /// Powers of omega from omega^0 to omega^(n - 1).
    omega_powers: []Complex(f64),
    /// Powers of omega^-1 = iomega from imoega^0 to iomega^(n - 1).
    omega_inverse_powers: []Complex(f64),
    /// Allocator used for internal memory allocations.
    allocator: Allocator,

    /// Initializes a new FFT instance with the given FFT vector length n.
    pub fn init(allocator: Allocator, n: i64) !Self {
        // Check if vector length n is a power of 2.
        if (!utils.isPowerOfTwo(n)) {
            return error.InvalidLength;
        }

        // Compute powers of omega as well as powers of omega^-1.
        var powers = try std.ArrayList(Complex(f64)).initCapacity(allocator, @intCast(n));
        var powers_inverse = try std.ArrayList(Complex(f64)).initCapacity(allocator, @intCast(n));
        defer powers.deinit();
        defer powers_inverse.deinit();

        for (0..@intCast(n)) |i| {
            const angle = 2 * pi * @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(n));
            const power = Complex(f64).init(@cos(angle), @sin(angle));
            const power_inverse = Complex(f64).init(@cos(-angle), @sin(-angle));
            try powers.append(power);
            try powers_inverse.append(power_inverse);
        }
        const omega_powers = try powers.toOwnedSlice();
        const omega_inverse_powers = try powers_inverse.toOwnedSlice();

        return FFT{
            .n = n,
            .omega_powers = omega_powers,
            .omega_inverse_powers = omega_inverse_powers,
            .allocator = allocator,
        };
    }

    /// Release all allocated memory.
    pub fn deinit(self: Self) void {
        self.allocator.free(self.omega_powers);
        self.allocator.free(self.omega_inverse_powers);
    }

    /// Runs a forward pass of FFT with the given coefficients.
    /// The caller owns the returned memory.
    pub fn fwd(self: Self, coefficients: []const Complex(f64)) ![]const Complex(f64) {
        // Length of coefficients must equal n.
        if (coefficients.len != self.n) {
            return error.InvalidLength;
        }

        return self.fft(coefficients, self.omega_powers);
    }

    /// Runs IFFT (Inverse FFT) with the given coefficients.
    /// The caller owns the returned memory.
    pub fn inv(self: Self, coefficients: []const Complex(f64)) ![]const Complex(f64) {
        // Length of coefficients must equal n.
        if (coefficients.len != self.n) {
            return error.InvalidLength;
        }

        const unscaled_result = try self.fft(coefficients, self.omega_inverse_powers);
        defer self.allocator.free(unscaled_result);

        var result = try self.allocator.alloc(Complex(f64), coefficients.len);

        const n = Complex(f64).init(@floatFromInt(self.n), 0);
        for (0..coefficients.len) |i| {
            result[i] = unscaled_result[i].div(n);
        }

        return result;
    }

    /// Runs an iterative version of FFT with the given coefficients and twiddles
    /// which are the powers of the roots of unity (i.e. powers of omega / powers
    /// of omega^-1).
    /// The caller owns the returned memory.
    fn fft(self: Self, coefficients: []const Complex(f64), twiddles: []const Complex(f64)) ![]const Complex(f64) {
        // Length of coefficients and twiddles must be the same.
        if (coefficients.len != twiddles.len) {
            return error.InvalidLength;
        }

        const log2_n = std.math.log2_int(usize, @intCast(self.n));
        const reversed = try utils.bitReverseSlice(Complex(f64), self.allocator, coefficients);
        defer self.allocator.free(reversed);

        var result = try self.allocator.dupe(Complex(f64), reversed);

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

                    const twiddle_index = @as(usize, (j * @as(usize, @intCast(self.n)))) >> @intCast(i + 1); // (j * n) >> i + 1 = (j * n) / 2^(i + 1)
                    const twiddle_factor = twiddles[twiddle_index].mul(d);

                    const butterfly_plus = c.add(twiddle_factor);
                    const butterfly_minus = c.sub(twiddle_factor);

                    result[index_even] = butterfly_plus;
                    result[index_odd] = butterfly_minus;
                }
            }
        }

        return result;
    }
};

test "fft - init" {
    const allocator = testing.allocator;

    {
        const n = 5;

        const expected = error.InvalidLength;
        const result = FFT.init(allocator, n);

        try testing.expectError(expected, result);
    }
}

test "fft - fwd" {
    const allocator = testing.allocator;

    {
        const n = 4;

        const fft = try FFT.init(allocator, n);
        defer fft.deinit();

        const coefficients = [_]Complex(f64){
            Complex(f64).init(1, 0),
            Complex(f64).init(2, 0),
            Complex(f64).init(3, 0),
            Complex(f64).init(4, 0),
        };
        const expected = [_]Complex(f64){
            Complex(f64).init(10, 0), //
            Complex(f64).init(-2, -2),
            Complex(f64).init(-2, 0),
            Complex(f64).init(-2, 2),
        };

        const result = try fft.fwd(&coefficients);
        defer allocator.free(result);

        try expectEqualComplexSlices(f64, &expected, result);
    }

    {
        const n = 4;

        const fft = try FFT.init(allocator, n);
        defer fft.deinit();

        const coefficients = [_]Complex(f64){
            Complex(f64).init(5, 0),
            Complex(f64).init(6, 0),
            Complex(f64).init(7, 0),
            Complex(f64).init(8, 0),
        };
        const expected = [_]Complex(f64){
            Complex(f64).init(26, 0), //
            Complex(f64).init(-2, -2),
            Complex(f64).init(-2, 0),
            Complex(f64).init(-2, 2),
        };

        const result = try fft.fwd(&coefficients);
        defer allocator.free(result);

        try expectEqualComplexSlices(f64, &expected, result);
    }

    {
        const n = 8;

        const fft = try FFT.init(allocator, n);
        defer fft.deinit();

        const coefficients = [_]Complex(f64){
            Complex(f64).init(1, 0), //
            Complex(f64).init(1, 0),
            Complex(f64).init(1, 0),
            Complex(f64).init(1, 0),
            Complex(f64).init(0, 0),
            Complex(f64).init(0, 0),
            Complex(f64).init(0, 0),
            Complex(f64).init(0, 0),
        };
        const expected = [_]Complex(f64){
            Complex(f64).init(4, 0), //
            Complex(f64).init(1, 2.41421),
            Complex(f64).init(0, 0),
            Complex(f64).init(1, 0.41421),
            Complex(f64).init(0, 0),
            Complex(f64).init(1, -0.41421),
            Complex(f64).init(0, 0),
            Complex(f64).init(1, -2.41421),
        };

        const result = try fft.fwd(&coefficients);
        defer allocator.free(result);

        try expectEqualComplexSlices(f64, &expected, result);
    }

    {
        const n = 4;

        const fft = try FFT.init(allocator, n);
        defer fft.deinit();

        const coefficients = [_]Complex(f64){
            Complex(f64).init(1, 0), //
            Complex(f64).init(2, 0),
            Complex(f64).init(3, 0),
        };

        const expected = error.InvalidLength;
        const result = fft.fwd(&coefficients);

        try testing.expectError(expected, result);
    }
}

test "fft - inv" {
    const allocator = testing.allocator;

    {
        const n = 4;

        const fft = try FFT.init(allocator, n);
        defer fft.deinit();

        const coefficients = [_]Complex(f64){
            Complex(f64).init(10, 0), //
            Complex(f64).init(-2, -2),
            Complex(f64).init(-2, 0),
            Complex(f64).init(-2, 2),
        };
        const expected = [_]Complex(f64){
            Complex(f64).init(1, 0),
            Complex(f64).init(2, 0),
            Complex(f64).init(3, 0),
            Complex(f64).init(4, 0),
        };

        const result = try fft.inv(&coefficients);
        defer allocator.free(result);

        try expectEqualComplexSlices(f64, &expected, result);
    }

    {
        const n = 4;

        const fft = try FFT.init(allocator, n);
        defer fft.deinit();

        const coefficients = [_]Complex(f64){
            Complex(f64).init(26, 0), //
            Complex(f64).init(-2, -2),
            Complex(f64).init(-2, 0),
            Complex(f64).init(-2, 2),
        };
        const expected = [_]Complex(f64){
            Complex(f64).init(5, 0),
            Complex(f64).init(6, 0),
            Complex(f64).init(7, 0),
            Complex(f64).init(8, 0),
        };

        const result = try fft.inv(&coefficients);
        defer allocator.free(result);

        try expectEqualComplexSlices(f64, &expected, result);
    }

    {
        const n = 8;

        const fft = try FFT.init(allocator, n);
        defer fft.deinit();

        const coefficients = [_]Complex(f64){
            Complex(f64).init(4, 0), //
            Complex(f64).init(1, 2.41421),
            Complex(f64).init(0, 0),
            Complex(f64).init(1, 0.41421),
            Complex(f64).init(0, 0),
            Complex(f64).init(1, -0.41421),
            Complex(f64).init(0, 0),
            Complex(f64).init(1, -2.41421),
        };
        const expected = [_]Complex(f64){
            Complex(f64).init(1, 0), //
            Complex(f64).init(1, 0),
            Complex(f64).init(1, 0),
            Complex(f64).init(1, 0),
            Complex(f64).init(0, 0),
            Complex(f64).init(0, 0),
            Complex(f64).init(0, 0),
            Complex(f64).init(0, 0),
        };

        const result = try fft.inv(&coefficients);
        defer allocator.free(result);

        try expectEqualComplexSlices(f64, &expected, result);
    }

    {
        const n = 4;

        const fft = try FFT.init(allocator, n);
        defer fft.deinit();

        const coefficients = [_]Complex(f64){
            Complex(f64).init(1, 0), //
            Complex(f64).init(2, 0),
            Complex(f64).init(3, 0),
        };

        const expected = error.InvalidLength;
        const result = fft.inv(&coefficients);

        try testing.expectError(expected, result);
    }
}

test "fft - fft" {
    const allocator = testing.allocator;

    {
        const n = 4;

        const fft = try FFT.init(allocator, n);
        defer fft.deinit();

        const twiddles = fft.omega_powers;
        const coefficients = [_]Complex(f64){
            Complex(f64).init(1, 0),
            Complex(f64).init(2, 0),
            Complex(f64).init(3, 0),
            Complex(f64).init(4, 0),
        };
        const expected = [_]Complex(f64){
            Complex(f64).init(10, 0), //
            Complex(f64).init(-2, -2),
            Complex(f64).init(-2, 0),
            Complex(f64).init(-2, 2),
        };

        const result = try fft.fft(&coefficients, twiddles);
        defer allocator.free(result);

        try expectEqualComplexSlices(f64, &expected, result);
    }

    {
        const n = 4;

        const fft = try FFT.init(allocator, n);
        defer fft.deinit();

        const twiddles = fft.omega_powers;
        const coefficients = [_]Complex(f64){
            Complex(f64).init(1, 0), //
            Complex(f64).init(2, 0),
        };

        const expected = error.InvalidLength;
        const result = fft.fft(&coefficients, twiddles);

        try testing.expectError(expected, result);
    }
}

test "fft - convolution" {
    const allocator = testing.allocator;

    {
        const n = 4;

        const fft = try FFT.init(allocator, n);
        defer fft.deinit();

        const coefficients_1 = [_]Complex(f64){
            Complex(f64).init(3, 0), //
            Complex(f64).init(2, 0),
            Complex(f64).init(0, 0),
            Complex(f64).init(0, 0),
        };
        const coefficients_2 = [_]Complex(f64){
            Complex(f64).init(1, 0), //
            Complex(f64).init(5, 0),
            Complex(f64).init(0, 0),
            Complex(f64).init(0, 0),
        };
        //   (3 + 2x + 0x^2 + 0x^3) * (1 + 5x + 0x^2 + 0x^3)
        // = 3 + 15x + 2x + 10x^2
        // = 3 + 17x + 10x^2
        const expected = [_]Complex(f64){
            Complex(f64).init(3, 0), //
            Complex(f64).init(17, 0),
            Complex(f64).init(10, 0),
            Complex(f64).init(0, 0),
        };

        const fwd_1 = try fft.fwd(&coefficients_1);
        const fwd_2 = try fft.fwd(&coefficients_2);
        defer allocator.free(fwd_1);
        defer allocator.free(fwd_2);

        var interim = [_]Complex(f64){Complex(f64).init(0, 0)} ** n;
        for (0..n) |i| {
            interim[i] = fwd_1[i].mul(fwd_2[i]);
        }

        const result = try fft.inv(&interim);
        defer allocator.free(result);

        try expectEqualComplexSlices(f64, &expected, result);

        // Modify result to turn complex number into integer.
        const expected_casted = [_]i64{ 3, 17, 10, 0 };
        var result_casted = try allocator.alloc(i64, n);
        defer allocator.free(result_casted);
        for (0..n) |i| {
            result_casted[i] = @intFromFloat(result[i].re + 0.5);
        }

        try testing.expectEqualSlices(i64, &expected_casted, result_casted);
    }
}
