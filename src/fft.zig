//! This module implements functionalities to run Fast Fourier Transform
//! computations.

const std = @import("std");
const math = std.math;
const testing = std.testing;
const pi = math.pi;
const Complex = math.Complex;

const utils = @import("utils.zig");
const expectEqualComplexSlices = @import("testing.zig").expectEqualComplexSlices;

/// Fast Fourier Transform based on n-th root of unity.
/// The FFT algorithm is based on the paper "Low-Cost and Area-Efficient FPGA
/// Implementations of Lattice-Based Cryptography" by Aysu et al.
// See: https://schaumont.dyn.wpi.edu/schaum/pdf/papers/2013hostb.pdf
pub fn FFT(comptime n: i64) type {
    return struct {
        const Self = @This();

        /// Length of FFT vector.
        comptime n: i64 = n,
        /// Powers of omega from omega^0 to omega^(n - 1).
        omega_powers: [n]Complex(f64),
        /// Powers of omega^-1 = iomega from imoega^0 to iomega^(n - 1).
        omega_inverse_powers: [n]Complex(f64),

        /// Initializes a new FFT instance. Ensures that the FFT vector length
        /// is a power of 2.
        pub fn init() !Self {
            // Check if FFT vector length n is a power of 2.
            if (!utils.isPowerOfTwo(n)) {
                return error.InvalidLength;
            }

            // Compute powers of omega as well as powers of omega^-1.
            var omega_powers = [_]Complex(f64){Complex(f64).init(0, 0)} ** n;
            var omega_inverse_powers = [_]Complex(f64){Complex(f64).init(0, 0)} ** n;

            for (0..n) |i| {
                const angle = 2 * pi * @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(n));
                omega_powers[i] = Complex(f64).init(@cos(angle), @sin(angle));
                omega_inverse_powers[i] = Complex(f64).init(@cos(-angle), @sin(-angle));
            }

            return Self{
                .n = n,
                .omega_powers = omega_powers,
                .omega_inverse_powers = omega_inverse_powers,
            };
        }

        /// Runs a forward pass of FFT with the given coefficients.
        /// Note that the slice of coefficients is mutated in-place.
        pub fn fwd(self: Self, coefficients: []Complex(f64)) ![]Complex(f64) {
            // Length of coefficients must equal n.
            if (coefficients.len != self.n) {
                return error.InvalidLength;
            }

            return self.fft(coefficients, self.omega_powers);
        }

        /// Runs iFFT (Inverse FFT) with the given coefficients.
        /// Note that the slice of coefficients is mutated in-place.
        pub fn inv(self: Self, coefficients: []Complex(f64)) ![]Complex(f64) {
            // Length of coefficients must equal n.
            if (coefficients.len != self.n) {
                return error.InvalidLength;
            }

            var result = try self.fft(coefficients, self.omega_inverse_powers);

            const n_complex = Complex(f64).init(@floatFromInt(self.n), 0);
            for (0..coefficients.len) |i| {
                result[i] = result[i].div(n_complex);
            }

            return result;
        }

        /// Runs an iterative version of FFT with the given coefficients and twiddles
        /// which are the powers of the roots of unity (i.e. powers of omega / powers
        /// of omega^-1).
        /// Note that the slice of coefficients is mutated in-place.
        fn fft(self: Self, coefficients: []Complex(f64), twiddles: [n]Complex(f64)) ![]Complex(f64) {
            // Length of coefficients and twiddles must be the same.
            if (coefficients.len != twiddles.len) {
                return error.InvalidLength;
            }

            const log2_n = std.math.log2_int(usize, @intCast(self.n));
            var result = try utils.bitReverseSlice(Complex(f64), self.n, coefficients);

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
}

test "fft - init" {
    {
        const n = 5;

        const expected = error.InvalidLength;
        const result = FFT(n).init();

        try testing.expectError(expected, result);
    }
}

test "fft - fwd" {
    {
        const n = 4;

        const fft = try FFT(n).init();

        var coefficients = [_]Complex(f64){
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

        try expectEqualComplexSlices(f64, &expected, result);
    }

    {
        const n = 4;

        const fft = try FFT(n).init();

        var coefficients = [_]Complex(f64){
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

        try expectEqualComplexSlices(f64, &expected, result);
    }

    {
        const n = 8;

        const fft = try FFT(n).init();

        var coefficients = [_]Complex(f64){
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

        try expectEqualComplexSlices(f64, &expected, result);
    }

    {
        const n = 4;

        const fft = try FFT(n).init();

        var coefficients = [_]Complex(f64){
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
    {
        const n = 4;

        const fft = try FFT(n).init();

        var coefficients = [_]Complex(f64){
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

        try expectEqualComplexSlices(f64, &expected, result);
    }

    {
        const n = 4;

        const fft = try FFT(n).init();

        var coefficients = [_]Complex(f64){
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

        try expectEqualComplexSlices(f64, &expected, result);
    }

    {
        const n = 8;

        const fft = try FFT(n).init();

        var coefficients = [_]Complex(f64){
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

        try expectEqualComplexSlices(f64, &expected, result);
    }

    {
        const n = 4;

        const fft = try FFT(n).init();

        var coefficients = [_]Complex(f64){
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
    {
        const n = 4;

        const fft = try FFT(n).init();

        const twiddles = fft.omega_powers;
        var coefficients = [_]Complex(f64){
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

        try expectEqualComplexSlices(f64, &expected, result);
    }

    {
        const n = 4;

        const fft = try FFT(n).init();

        const twiddles = fft.omega_powers;
        var coefficients = [_]Complex(f64){
            Complex(f64).init(1, 0), //
            Complex(f64).init(2, 0),
        };

        const expected = error.InvalidLength;
        const result = fft.fft(&coefficients, twiddles);

        try testing.expectError(expected, result);
    }
}

test "fft - convolution" {
    {
        const n = 4;

        const fft = try FFT(n).init();

        var coefficients_1 = [_]Complex(f64){
            Complex(f64).init(3, 0), //
            Complex(f64).init(2, 0),
            Complex(f64).init(0, 0),
            Complex(f64).init(0, 0),
        };
        var coefficients_2 = [_]Complex(f64){
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

        var interim = [_]Complex(f64){Complex(f64).init(0, 0)} ** n;
        for (0..n) |i| {
            interim[i] = fwd_1[i].mul(fwd_2[i]);
        }

        const result = try fft.inv(&interim);

        try expectEqualComplexSlices(f64, &expected, result);

        // Modify result to turn complex number into integer.
        const expected_casted = [_]i64{ 3, 17, 10, 0 };
        var result_casted = [_]i64{0} ** n;
        for (0..n) |i| {
            result_casted[i] = @intFromFloat(result[i].re + 0.5);
        }

        try testing.expectEqualSlices(i64, &expected_casted, &result_casted);
    }
}
