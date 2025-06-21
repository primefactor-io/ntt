//! This module implements testing utilities.

const std = @import("std");
const math = std.math;
const Complex = math.Complex;

/// Checks if two slices with complex numbers are approximately equal.
pub fn expectEqualComplexSlices(comptime T: type, a: []const Complex(T), b: []const Complex(T)) !void {
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
