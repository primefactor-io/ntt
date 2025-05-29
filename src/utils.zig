const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;

/// Finds a root of unity with the given order n (i.e. an nth root of unity) in
/// the given prime modulus m.
/// Note that the modulus MUST be prime.
pub fn findRootOfUnity(allocator: Allocator, n: i64, m: i64) !i64 {
    // Ensure that order n | m - 1.
    if (@mod(m - 1, n) != 0) {
        return error.InvalidOrder;
    }

    const generator = try findPrimitiveRoot(allocator, m);

    const base = generator;
    const exponent = @divFloor(m - 1, n);

    const result = try pow(base, exponent, m);
    if (result == 1) {
        return findRootOfUnity(allocator, n, m);
    }

    return result;
}

/// Finds a primitive root (a generator) in the given prime modulus.
/// Note that the modulus MUST be prime.
// See: https://cp-algorithms.com/algebra/primitive-root.html#implementation
pub fn findPrimitiveRoot(allocator: Allocator, m: i64) !i64 {
    var fact = std.ArrayList(i64).init(allocator);
    defer fact.deinit();

    const phi = m - 1;
    var n = phi;

    var i: i64 = 2;
    while (i * i <= n) : (i += 1) {
        if (@mod(n, i) == 0) {
            try fact.append(i);
            while (@mod(n, i) == 0) {
                n = @divFloor(n, i);
            }
        }
    }

    if (n > 1) {
        try fact.append(n);
    }

    var result: i64 = 2;
    while (result <= m) : (result += 1) {
        var ok = true;

        for (fact.items) |f| {
            const base = result;
            const exponent = @divFloor(phi, f);

            if (try pow(base, exponent, m) == 1) {
                ok = false;
                break;
            }
        }

        if (ok) {
            return result;
        }
    }

    return error.PrimitiveRootNotFound;
}

/// Returns (gcd, x, y) such that a * x + b * y = gcd = gcd(a, b).
pub fn egcd(a: i64, b: i64) struct { gcd: i64, x: i64, y: i64 } {
    if (a == 0) {
        const gcd = b;
        const x = 0;
        const y = 1;
        return .{ .gcd = gcd, .x = x, .y = y };
    }
    const result = egcd(@mod(b, a), a);
    const gcd = result.gcd;
    const x = result.x;
    const y = result.y;

    return .{ .gcd = gcd, .x = y - (@divFloor(b, a) * x), .y = x };
}

/// Calculates x such that (x * a) % m == 1.
pub fn modInv(a: i64, m: i64) !i64 {
    const result = egcd(a, m);
    if (result.gcd != 1) {
        return error.InvalidGCD;
    }
    return @mod(result.x, m);
}

/// Calculates a^b mod m.
/// Note that b can be either >= 0 or -1. If it's -1, then the modular
/// multiplicative inverse is calculated.
// See: https://cp-algorithms.com/algebra/primitive-root.html#implementation
pub fn pow(a: i64, b: i64, m: i64) !i64 {
    if (b == -1) {
        return modInv(a, m);
    } else if (b < -1) {
        return error.InvalidExponent;
    }

    var base = a;
    var exponent = b;
    const modulus = m;
    var result: i64 = 1;

    while (exponent != 0) {
        if (exponent & 1 == 1) {
            result = @mod(result * base, modulus);
            exponent -= 1;
        } else {
            base = @mod(base * base, modulus);
            exponent >>= 1;
        }
    }

    return result;
}

test "findRootOfUnity - core" {
    var n: i64 = undefined;
    var m: i64 = undefined;
    var expected: i64 = undefined;
    var result: i64 = undefined;

    n = 2;
    m = 5;
    expected = 4;
    result = try findRootOfUnity(testing.allocator, n, m);
    try testing.expectEqual(expected, result);

    n = 3;
    m = 7;
    expected = 2;
    result = try findRootOfUnity(testing.allocator, n, m);

    n = 5;
    m = 11;
    expected = 4;
    result = try findRootOfUnity(testing.allocator, n, m);
}

test "findRootOfUnity - invalid order" {
    const n: i64 = 3;
    const m: i64 = 11;

    const expected = error.InvalidOrder;
    const result = findRootOfUnity(testing.allocator, n, m);

    try testing.expectError(expected, result);
}

test "findPrimitiveRoot" {
    var m: i64 = undefined;
    var expected: i64 = undefined;
    var result: i64 = undefined;

    m = 5;
    expected = 2;
    result = try findPrimitiveRoot(testing.allocator, m);
    try testing.expectEqual(expected, result);

    m = 7;
    expected = 3;
    result = try findPrimitiveRoot(testing.allocator, m);
    try testing.expectEqual(expected, result);

    m = 11;
    expected = 2;
    result = try findPrimitiveRoot(testing.allocator, m);
    try testing.expectEqual(expected, result);
}

test "egcd" {
    const a = 38;
    const b = 97;

    const result = egcd(a, b);
    const gcd = result.gcd;
    const x = result.x;
    const y = result.y;

    try testing.expectEqual(1, gcd);
    try testing.expectEqual(23, x);
    try testing.expectEqual(-9, y);
}

test "modInv" {
    const a = 38;
    const m = 97;

    const result = try modInv(a, m);

    try testing.expectEqual(23, result);
}

test "pow - b > 0" {
    const a = 4;
    const b = 3;
    const m = 5;

    const result = try pow(a, b, m);

    try testing.expectEqual(4, result);
}

test "pow - b = 0" {
    const a = 4;
    const b = 0;
    const m = 5;

    const result = try pow(a, b, m);

    try testing.expectEqual(1, result);
}

test "pow - b = -1" {
    const a = 38;
    const b = -1;
    const m = 97;

    const result = try pow(a, b, m);

    try testing.expectEqual(23, result);
}

test "pow - b < -1" {
    const a = 4;
    const b = -2;
    const m = 5;

    const result = pow(a, b, m);

    try testing.expectError(error.InvalidExponent, result);
}
