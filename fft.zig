const std = @import("std");
const Allocator = std.mem.Allocator;
const testing = std.testing;
const math = std.math;
const pi = math.pi;
const Complex = math.Complex;

const utils = @import("src/utils.zig");

const Mode = enum {
    fwd,
    inv,
};

fn fft(allocator: Allocator, numbers: []const Complex(f64), mode: Mode) ![]const Complex(f64) {
    const n = numbers.len;

    // Ensure that length of numbers is a power of 2.
    if (!utils.isPowerOfTwo(@intCast(n))) {
        return error.InvalidLength;
    }

    // Compute powers of omega as well as powers of omega^-1.
    var powers = try std.ArrayList(Complex(f64)).initCapacity(allocator, n);
    var powers_inv = try std.ArrayList(Complex(f64)).initCapacity(allocator, n);
    defer powers.deinit();
    defer powers_inv.deinit();

    for (0..n) |i| {
        const angle = 2 * pi * @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(n));
        const power = Complex(f64).init(@cos(angle), @sin(angle));
        const power_inv = Complex(f64).init(@cos(-angle), @sin(-angle));
        try powers.append(power);
        try powers_inv.append(power_inv);
    }

    const omega_powers = try powers.toOwnedSlice();
    const omega_powers_inv = try powers_inv.toOwnedSlice();
    defer allocator.free(omega_powers);
    defer allocator.free(omega_powers_inv);

    // Run iterative FFT.
    const log2_n = math.log2_int(usize, @intCast(n));

    const reversed = try utils.bitReverseSlice(Complex(f64), allocator, numbers);
    defer allocator.free(reversed);

    var result = try allocator.dupe(Complex(f64), reversed);

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

                const rou_index = @as(usize, (j * n)) >> @intCast(i + 1); // (j * n) >> i + 1 = (j * n) / 2^(i + 1)

                const omega_factor = switch (mode) {
                    .fwd => omega_powers[rou_index].mul(d),
                    .inv => omega_powers_inv[rou_index].mul(d),
                };

                const butterfly_plus = c.add(omega_factor);
                const butterfly_minus = c.sub(omega_factor);

                result[index_even] = butterfly_plus;
                result[index_odd] = butterfly_minus;
            }
        }
    }

    // Rescale if iFFT is run.
    if (mode == Mode.inv) {
        for (0..n) |i| {
            result[i] = result[i].div(Complex(f64).init(@floatFromInt(n), 0));
        }
    }

    return result;
}

test "fft" {
    const allocator = testing.allocator;

    // See: https://github.com/sarojaerabelli/py-fhe/blob/master/tests/test_ntt.py#L29-L31
    {
        const numbers = [_]Complex(f64){
            Complex(f64).init(0, 0), //
            Complex(f64).init(1, 0),
            Complex(f64).init(4, 0),
            Complex(f64).init(5, 0),
        };
        const expected = [_]Complex(f64){
            Complex(f64).init(10, 0), //
            Complex(f64).init(-4, -4),
            Complex(f64).init(-2, 0),
            Complex(f64).init(-4, 4),
        };

        const result = try fft(allocator, &numbers, Mode.fwd);
        defer allocator.free(result);

        // try testing.expectEqualSlices(Complex(f64), &expected, result);
        try expectSimilarComplexSlices(f64, &expected, result);
    }

    // See: https://rosettacode.org/wiki/Fast_Fourier_transform#C#
    {
        const numbers = [_]Complex(f64){
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

        const result = try fft(allocator, &numbers, Mode.fwd);
        defer allocator.free(result);

        // try testing.expectEqualSlices(Complex(f64), &expected, result);
        try expectSimilarComplexSlices(f64, &expected, result);
    }

    // See: https://rosettacode.org/wiki/Fast_Fourier_transform#C#
    {
        const numbers = [_]Complex(f64){
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

        const result = try fft(allocator, &numbers, Mode.inv);
        defer allocator.free(result);

        // try testing.expectEqualSlices(Complex(f64), &expected, result);
        try expectSimilarComplexSlices(f64, &expected, result);
    }

    // See: Antti Laaksonen - Guide to Competitive Programming - 3rd Edition
    {
        const numbers_1 = [_]Complex(f64){
            Complex(f64).init(3, 0), //
            Complex(f64).init(2, 0),
            Complex(f64).init(0, 0),
            Complex(f64).init(0, 0),
        };
        const numbers_2 = [_]Complex(f64){
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

        const n = numbers_1.len;
        const fwd_1 = try fft(allocator, &numbers_1, Mode.fwd);
        const fwd_2 = try fft(allocator, &numbers_2, Mode.fwd);
        defer allocator.free(fwd_1);
        defer allocator.free(fwd_2);

        var interim = [_]Complex(f64){Complex(f64).init(0, 0)} ** n;
        for (0..n) |i| {
            interim[i] = fwd_1[i].mul(fwd_2[i]);
        }

        const result = try fft(allocator, &interim, Mode.inv);
        defer allocator.free(result);

        // try testing.expectEqualSlices(Complex(f64), &expected, result);
        try expectSimilarComplexSlices(f64, &expected, result);

        const expected_rescaled = [_]i64{
            3,
            17,
            10,
            0,
        };
        var result_rescaled = try allocator.alloc(i64, n);
        defer allocator.free(result_rescaled);
        for (0..n) |i| {
            result_rescaled[i] = @intFromFloat(result[i].re + 0.5);
        }

        try testing.expectEqualSlices(i64, &expected_rescaled, result_rescaled);
    }
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
