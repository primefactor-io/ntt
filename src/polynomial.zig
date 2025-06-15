const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;

const NTT = @import("ntt.zig").NTT;
const utils = @import("utils.zig");

// TODO: Add doc comment.
pub const Polynomial = struct {
    const Self = @This();

    /// Coefficient modulus.
    q: i64,
    /// Degree of cyclotomic polynomial.
    n: i64,
    /// Coefficients of the polynomial.
    coefficients: []const i64,
    /// Allocator used for internal memory allocations.
    allocator: Allocator,

    /// Initializes a new Polynomial instance with the given coefficient modulus and
    /// cyclotomic polynomial degree (which MUST be a power of 2).
    /// The coefficients are copied and therefore the memory used is managed by
    /// the instance.
    pub fn init(allocator: Allocator, coefficients: []const i64, q: i64) !Self {
        // Derive degree from coefficient slice.
        const n: i64 = @intCast(coefficients.len);

        // Check if degree of cyclotomic polynomial is a power of 2.
        if (!utils.isPowerOfTwo(n)) {
            return error.InvalidDegree;
        }

        const coeffs = try allocator.dupe(i64, coefficients);

        return Polynomial{ .q = q, .n = n, .coefficients = coeffs, .allocator = allocator };
    }

    /// Release all allocated memory.
    pub fn deinit(self: Self) void {
        self.allocator.free(self.coefficients);
    }

    // TODO:
    pub fn add(self: Self, other: Self) !Self {
        try self.checkEquality(other);

        const result = try self.allocator.alloc(i64, @intCast(self.n));
        defer self.allocator.free(result);

        for (0..@intCast(self.n)) |idx| {
            result[idx] = @mod(self.coefficients[idx] + other.coefficients[idx], self.q);
        }

        return Polynomial.init(self.allocator, result, self.q);
    }

    // TODO:
    pub fn subtract(self: Self, other: Self) !Self {
        _ = other;

        return self;
    }

    // TODO:
    pub fn multiply(self: Self, other: Self) !Self {
        _ = other;
        // TODO: Create NTT.

        return self;
    }

    /// Performs modular reduction on all of the polynomial's coefficients.
    /// The caller owns the returned memory.
    pub fn mod(self: Self) ![]const i64 {
        const result = try self.allocator.alloc(i64, @intCast(self.n));

        for (self.coefficients, 0..) |coefficient, idx| {
            result[idx] = @mod(coefficient, self.q);
        }

        return result;
    }

    // TODO:
    fn checkEquality(self: Self, other: Self) !void {
        if (self.q != other.q) {
            return error.UnequalModulus;
        } else if (self.n != other.n) {
            return error.UnequalDegree;
        }
    }
};

test "polynomial - init" {
    const allocator = testing.allocator;

    {
        const q = 7681;
        const coefficients = [_]i64{ 1, 2, 3, 4, 5 };

        const expected = error.InvalidDegree;
        const result = Polynomial.init(allocator, &coefficients, q);

        try testing.expectError(expected, result);
    }
}

test "polynomial - add" {
    const allocator = testing.allocator;

    {
        const q = 7681;

        const coefficients_1 = [_]i64{ 1, 2, 3, 4 };
        const poly_1 = try Polynomial.init(allocator, &coefficients_1, q);
        defer poly_1.deinit();

        const coefficients_2 = [_]i64{ 5, 6, 7, 8 };
        const poly_2 = try Polynomial.init(allocator, &coefficients_2, q);
        defer poly_2.deinit();

        const expected = [_]i64{ 6, 8, 10, 12 };
        const result = try poly_1.add(poly_2);
        defer result.deinit();

        try testing.expectEqualSlices(i64, &expected, result.coefficients);
    }

    {
        const q = 7681;

        const coefficients_1 = [_]i64{ 1, 2, 3, 4 };
        const poly_1 = try Polynomial.init(allocator, &coefficients_1, q);
        defer poly_1.deinit();

        const coefficients_2 = [_]i64{ 7681, 7682, 7683, 7684 };
        const poly_2 = try Polynomial.init(allocator, &coefficients_2, q);
        defer poly_2.deinit();

        const expected = [_]i64{ 1, 3, 5, 7 };
        const result = try poly_1.add(poly_2);
        defer result.deinit();

        try testing.expectEqualSlices(i64, &expected, result.coefficients);
    }

    {
        const q_1 = 7681;
        const coefficients_1 = [_]i64{ 1, 2, 3, 4 };
        const poly_1 = try Polynomial.init(allocator, &coefficients_1, q_1);
        defer poly_1.deinit();

        const q_2 = 7682;
        const coefficients_2 = [_]i64{ 1, 2, 3, 4 };
        const poly_2 = try Polynomial.init(allocator, &coefficients_2, q_2);
        defer poly_2.deinit();

        const expected = error.UnequalModulus;
        const result = poly_1.add(poly_2);

        try testing.expectError(expected, result);
    }

    {
        const q = 7681;

        const coefficients_1 = [_]i64{ 1, 2, 3, 4 };
        const poly_1 = try Polynomial.init(allocator, &coefficients_1, q);
        defer poly_1.deinit();

        const coefficients_2 = [_]i64{ 1, 2, 3, 4, 5, 6, 7, 8 };
        const poly_2 = try Polynomial.init(allocator, &coefficients_2, q);
        defer poly_2.deinit();

        const expected = error.UnequalDegree;
        const result = poly_1.add(poly_2);

        try testing.expectError(expected, result);
    }
}

test "polynomial - mod" {
    const allocator = testing.allocator;

    {
        const q = 7681;
        const coefficients = [_]i64{ 1, 2, 3, 4 };

        const poly = try Polynomial.init(allocator, &coefficients, q);
        defer poly.deinit();

        const expected = [_]i64{ 1, 2, 3, 4 };
        const result = try poly.mod();
        defer allocator.free(result);

        try testing.expectEqualSlices(i64, &expected, result);
    }

    {
        const q = 7681;
        const coefficients = [_]i64{ 7681, 7682, 7683, 7684 };

        const poly = try Polynomial.init(allocator, &coefficients, q);
        defer poly.deinit();

        const expected = [_]i64{ 0, 1, 2, 3 };
        const result = try poly.mod();
        defer allocator.free(result);

        try testing.expectEqualSlices(i64, &expected, result);
    }
}

test "polynomial - checkEquality" {
    const allocator = testing.allocator;

    {
        const q_1 = 7681;
        const coefficients_1 = [_]i64{ 1, 2, 3, 4 };
        const poly_1 = try Polynomial.init(allocator, &coefficients_1, q_1);
        defer poly_1.deinit();

        const q_2 = 7682;
        const coefficients_2 = [_]i64{ 1, 2, 3, 4 };
        const poly_2 = try Polynomial.init(allocator, &coefficients_2, q_2);
        defer poly_2.deinit();

        const expected = error.UnequalModulus;
        const result = poly_1.checkEquality(poly_2);

        try testing.expectError(expected, result);
    }

    {
        const q = 7681;

        const coefficients_1 = [_]i64{ 1, 2, 3, 4 };
        const poly_1 = try Polynomial.init(allocator, &coefficients_1, q);
        defer poly_1.deinit();

        const coefficients_2 = [_]i64{ 1, 2, 3, 4, 5, 6, 7, 8 };
        const poly_2 = try Polynomial.init(allocator, &coefficients_2, q);
        defer poly_2.deinit();

        const expected = error.UnequalDegree;
        const result = poly_1.checkEquality(poly_2);

        try testing.expectError(expected, result);
    }
}
