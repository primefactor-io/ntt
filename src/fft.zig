const std = @import("std");
const math = std.math;
const testing = std.testing;
const Allocator = std.mem.Allocator;
const pi = math.pi;
const Complex = math.Complex;

const utils = @import("utils.zig");

// TODO: Add comment
pub const FFT = struct {
    const Self = @This();

    /// Length of FFT vector.
    length: i64,
    /// Powers of omega from omega^0 to omega^(length - 1).
    omega_powers: []Complex(f64),
    /// Powers of omega^-1 = iomega from imoega^0 to iomega^(length - 1).
    omega_inverse_powers: []Complex(f64),
    /// Allocator used for internal memory allocations.
    allocator: Allocator,

    /// Initializes a new FFT instance with the given FFT vector length.
    pub fn init(allocator: Allocator, length: i64) !Self {
        // Check if length is a power of 2.
        if (!utils.isPowerOfTwo(length)) {
            return error.InvalidLength;
        }

        // Compute powers of omega as well as powers of omega^-1.
        var powers = try std.ArrayList(Complex(f64)).initCapacity(allocator, @intCast(length));
        var inverse_powers = try std.ArrayList(Complex(f64)).initCapacity(allocator, @intCast(length));
        defer powers.deinit();
        defer inverse_powers.deinit();

        for (0..@intCast(length)) |i| {
            const angle = 2 * pi * @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(length));
            const power = Complex(f64).init(@cos(angle), @sin(angle));
            const inverse_power = Complex(f64).init(@cos(-angle), @sin(-angle));
            try powers.append(power);
            try inverse_powers.append(inverse_power);
        }
        const omega_powers = try powers.toOwnedSlice();
        const omega_inverse_powers = try inverse_powers.toOwnedSlice();

        return FFT{
            .length = length,
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

    /// Runs an iterative version with the given coefficients and twiddles.
    /// The caller owns the returned memory.
    fn fft(self: Self, coefficients: []const Complex(f64), twiddles: []const Complex(f64)) ![]const Complex(f64) {
        // Length of coefficients and twiddles must be the same.
        if (coefficients.len != twiddles.len) {
            return error.InvalidLength;
        }

        const log2_n = std.math.log2_int(usize, @intCast(self.length));
        const reversed = try utils.bitReverseSlice(Complex(f64), self.allocator, coefficients);
        defer self.allocator.free(reversed);

        var result = try self.allocator.dupe(Complex(f64), reversed);

        // TODO: See: https://rosettacode.org/wiki/Fast_Fourier_transform#C#
        // var n: usize = 2;
        // while (n <= coefficients.len) : (n <<= 1) {
        //     var i: usize = 0;
        //     while (i < coefficients.len) : (i += n) {
        //         var k: usize = 0;
        //         while (k < n / 2) : (k += 1) {
        //             const evenIdx = i + k;
        //             const oddIdx = evenIdx + (n / 2);

        //             const even = result[evenIdx];
        //             const odd = result[oddIdx];

        //             const term = -2 * pi * @as(f64, @floatFromInt(k)) / @as(f64, @floatFromInt(n));
        //             const exp = Complex(f64).init(@cos(term), @sin(term)).mul(odd);

        //             result[evenIdx] = even.add(exp);
        //             result[oddIdx] = even.sub(exp);
        //         }
        //     }
        // }

        for (0..log2_n) |i| {
            var temp_twiddle = Complex(f64).init(1, 0);
            var final_twiddle = Complex(f64).init(1, 0);
            const twiddle = twiddles[i];

            const in_1 = @as(usize, 1) << @intCast(i); // 2^i
            const in_2 = @as(usize, 1) << @intCast(i + 1); // 2^(i + 1)
            const in_3 = @as(usize, @intCast(self.length)) / in_2; // length / (2^(i + 1))

            for (0..in_1) |j| {
                for (0..in_3) |t| {
                    const index_1 = (t * in_2) + j; // (t * 2^(i + 1)) + j
                    const index_2 = index_1 + in_1; // (t * 2^(i + 1)) + j + 2^i

                    const c = result[index_1];
                    const d = result[index_2];

                    const butterfly_plus = c.add(final_twiddle.mul(d));
                    const butterfly_minus = c.sub(final_twiddle.mul(d));

                    result[index_1] = butterfly_plus;
                    result[index_2] = butterfly_minus;

                    temp_twiddle = temp_twiddle.mul(twiddle);
                }
                final_twiddle = temp_twiddle;
            }
        }

        return result;
    }
};

test "fft - init" {
    const allocator = testing.allocator;

    {
        const length = 5;

        const expected = error.InvalidLength;
        const result = FFT.init(allocator, length);

        try testing.expectError(expected, result);
    }
}

test "fft - fft" {
    const allocator = testing.allocator;

    {
        const length = 8;

        const fft = try FFT.init(allocator, length);
        defer fft.deinit();

        const twiddles = fft.omega_powers;
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
            Complex(f64).init(1, -2.41421),
            Complex(f64).init(0, 0),
            Complex(f64).init(1, -0.41421),
            Complex(f64).init(0, 0),
            Complex(f64).init(1, 0.41421),
            Complex(f64).init(0, 0),
            Complex(f64).init(1, 2.41421),
        };

        const result = try fft.fft(&coefficients, twiddles);
        defer allocator.free(result);

        try testing.expectEqualSlices(Complex(f64), &expected, result);
        // try expectSimilarComplexSlices(f64, &expected, result);
    }

    // {
    //     const length = 4;

    //     const fft = try FFT.init(allocator, length);
    //     defer fft.deinit();

    //     const twiddles = fft.omega_powers;
    //     const coefficients = [_]Complex(f64){
    //         Complex(f64).init(0, 0), //
    //         Complex(f64).init(1, 0),
    //         Complex(f64).init(4, 0),
    //         Complex(f64).init(5, 0),
    //     };
    //     const expected = [_]Complex(f64){
    //         Complex(f64).init(10, 0), //
    //         Complex(f64).init(-4, -4),
    //         Complex(f64).init(-2, 0),
    //         Complex(f64).init(-4, 4),
    //     };

    //     const result = try fft.fft(&coefficients, twiddles);
    //     defer allocator.free(result);

    //     // try testing.expectEqualSlices(Complex(f64), &expected, result);
    //     try expectSimilarComplexSlices(f64, &expected, result);
    // }

    {
        const length = 4;

        const fft = try FFT.init(allocator, length);
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
    // TODO: See: Antti Laaksonen - Guide to Competitive Programming - 3rd Edition
    // const allocator = testing.allocator;

    // {
    //     const length = 4;

    //     const fft = try FFT.init(allocator, length);
    //     defer fft.deinit();

    //     const coefficients_1 = [_]f64{3, 2, 0, 0};
    //     const coefficients_2 = [_]f64{1, 5, 0, 0};

    //     const expected = [_]f64{3, 17, 10, 0}; // (int)(p[i].real()+0.5)

    // }
}

fn expectSimilarComplexSlices(comptime T: type, a: []const Complex(T), b: []const Complex(T)) !void {
    if (a.len != b.len) {
        return error.InvalidLength;
    }

    const threshold = 0.00001;

    for (0..a.len) |i| {
        if (@abs(a[i].re - b[i].re) > threshold or (@abs(a[i].im - b[i].im) > threshold)) {
            return error.DifferentSlices;
        }
    }
}
