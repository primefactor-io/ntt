const std = @import("std");
const Allocator = std.mem.Allocator;
const testing = std.testing;
const math = std.math;
const pi = math.pi;
const Complex = math.Complex;

const utils = @import("src/utils.zig");

// TODO: pass-in twiddles?!
//  if so, length of number and twiddles must be the same.
fn fft(allocator: Allocator, numbers: []const Complex(f64)) ![]const Complex(f64) {
    const n = numbers.len;

    if (!utils.isPowerOfTwo(@intCast(n))) {
        return error.InvalidLength;
    }

    const lgn = math.log2_int(usize, @intCast(n));
    const reversed = try utils.bitReverseSlice(Complex(f64), allocator, numbers);
    defer allocator.free(reversed);

    var result = try allocator.dupe(Complex(f64), reversed);

    var s: usize = 1;
    while (s <= lgn) : (s += 1) {
        const m = math.pow(usize, 2, s); // 2^s
        // const m = @as(usize, 1) << @intCast(s); // 2^s

        const angle = 2 * pi / @as(f64, @floatFromInt(m));
        const omega_m = Complex(f64).init(@cos(angle), @sin(angle));

        var k: usize = 0;
        while (k < n) : (k += m) {
            var omega = Complex(f64).init(1, 0);

            var j: usize = 0;
            while (j < m / 2) : (j += 1) {
                const index_1 = k + j;
                const index_2 = k + j + (m / 2);

                const t = omega.mul(result[index_2]);
                const u = result[index_1];

                result[index_1] = u.add(t);
                result[index_2] = u.sub(t);

                omega = omega.mul(omega_m);
            }
        }
    }

    return result;
}

test "fft" {
    const allocator = testing.allocator;

    {
        // See: https://github.com/sarojaerabelli/py-fhe/blob/master/tests/test_ntt.py#L29-L31
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

        const result = try fft(allocator, &numbers);
        defer allocator.free(result);

        // try testing.expectEqualSlices(Complex(f64), &expected, result);
        try expectSimilarComplexSlices(f64, &expected, result);
    }

    {
        // See: https://rosettacode.org/wiki/Fast_Fourier_transform#C#
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
            Complex(f64).init(1, -2.41421),
            Complex(f64).init(0, 0),
            Complex(f64).init(1, -0.41421),
            Complex(f64).init(0, 0),
            Complex(f64).init(1, 0.41421),
            Complex(f64).init(0, 0),
            Complex(f64).init(1, 2.41421),
        };
        // const expected = [_]Complex(f64){
        //     Complex(f64).init(4, 0), //
        //     Complex(f64).init(1, 2.41421),
        //     Complex(f64).init(0, 0),
        //     Complex(f64).init(1, 0.41421),
        //     Complex(f64).init(0, 0),
        //     Complex(f64).init(1, -0.41421),
        //     Complex(f64).init(0, 0),
        //     Complex(f64).init(1, -2.41421),
        // };

        const result = try fft(allocator, &numbers);
        defer allocator.free(result);

        try testing.expectEqualSlices(Complex(f64), &expected, result);
        // try expectSimilarComplexSlices(f64, &expected, result);
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
