const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;

/// Finds a root of unity with the given order n (i.e. an nth root of unity) in
/// the given prime modulus m.
pub fn findRootOfUnity(allocator: Allocator, n: i64, m: i64) !i64 {
    if (!try isPrime(m)) {
        return error.NotPrime;
    }

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
// See: https://cp-algorithms.com/algebra/primitive-root.html#implementation
pub fn findPrimitiveRoot(allocator: Allocator, m: i64) !i64 {
    if (!try isPrime(m)) {
        return error.NotPrime;
    }

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

/// Checks if a given number is (probably) prime using the Miller-Rabin
/// primality test algorithm.
pub fn isPrime(number: i64) !bool {
    if (number < 4) {
        return number == 2 or number == 3;
    }

    // Check if number is even.
    if (@mod(number, 2) == 0) {
        return false;
    }

    var exponent = number - 1;
    while (@mod(exponent, 2) == 0) {
        exponent = @divFloor(exponent, 2);
    }

    // See: https://zig.guide/standard-library/random-numbers/
    var prng = std.Random.DefaultPrng.init(blk: {
        var seed: u64 = undefined;
        try std.posix.getrandom(std.mem.asBytes(&seed));
        break :blk seed;
    });
    const rand = prng.random();

    // Perform 1_000 trials.
    for (0..1_000) |_| {
        const rand_val: i64 = rand.intRangeAtMost(i64, 1, number - 1);
        var new_exponent = exponent;

        var power = try pow(rand_val, new_exponent, number);

        while (new_exponent != number - 1 and power != 1 and power != number - 1) {
            power = @mod((power * power), number);
            new_exponent *= 2;
        }

        if (power != number - 1 and @mod(new_exponent, 2) == 0) {
            return false;
        }
    }

    return true;
}

/// Checks if a given number is a power of two.
// See: https://stackoverflow.com/a/600306
pub fn isPowerOfTwo(number: i64) bool {
    if (number <= 1) {
        return false;
    }

    return number & (number - 1) == 0;
}

/// Reverses the bits of the number using the specified bit width.
/// The value 9 = 0b1001 with a desired width of 6 would be 0b001001 and its
/// reverse would therefore be 0b100100 = 36.
pub fn bitReverseNumber(number: i64, width: i64) i64 {
    var x = number;
    var result: i64 = 0;

    for (0..@intCast(width)) |_| {
        result = (result << 1) | (x & 1);
        x >>= 1;
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

test "findRootOfUnity - m not prime" {
    const n: i64 = 3;
    const m: i64 = 10;

    const expected = error.NotPrime;
    const result = findRootOfUnity(testing.allocator, n, m);

    try testing.expectError(expected, result);
}

test "findPrimitiveRoot - core" {
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

test "findPrimitiveRoot - m not prime" {
    const m = 10;

    const expected = error.NotPrime;
    const result = findPrimitiveRoot(testing.allocator, m);

    try testing.expectError(expected, result);
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

test "isPrime" {
    var number: i64 = undefined;
    var result: bool = undefined;

    number = 1;
    result = try isPrime(number);
    try testing.expectEqual(false, result);

    number = 2;
    result = try isPrime(number);
    try testing.expectEqual(true, result);

    number = 3;
    result = try isPrime(number);
    try testing.expectEqual(true, result);

    number = 4;
    result = try isPrime(number);
    try testing.expectEqual(false, result);

    number = 5;
    result = try isPrime(number);
    try testing.expectEqual(true, result);

    number = 6;
    result = try isPrime(number);
    try testing.expectEqual(false, result);

    number = 7;
    result = try isPrime(number);
    try testing.expectEqual(true, result);

    number = 8;
    result = try isPrime(number);
    try testing.expectEqual(false, result);

    number = 9;
    result = try isPrime(number);
    try testing.expectEqual(false, result);

    number = 10;
    result = try isPrime(number);
    try testing.expectEqual(false, result);

    number = 39877;
    result = try isPrime(number);
    try testing.expectEqual(true, result);

    number = 74903;
    result = try isPrime(number);
    try testing.expectEqual(true, result);
}

test "isPowerOfTwo" {
    var number: i64 = undefined;
    var result: bool = undefined;

    number = 0;
    result = isPowerOfTwo(number);
    try testing.expectEqual(false, result);

    number = 1;
    result = isPowerOfTwo(number);
    try testing.expectEqual(false, result);

    number = 2;
    result = isPowerOfTwo(number);
    try testing.expectEqual(true, result);

    number = 256;
    result = isPowerOfTwo(number);
    try testing.expectEqual(true, result);

    number = 2048;
    result = isPowerOfTwo(number);
    try testing.expectEqual(true, result);

    number = 1234;
    result = isPowerOfTwo(number);
    try testing.expectEqual(false, result);
}

test "bitReverseNumber" {
    {
        const number = 100; // 0b1100100
        const width = 7;

        const expected = 19; // 0b0010011
        const result = bitReverseNumber(number, width);

        try testing.expectEqual(expected, result);
    }

    {
        const number = 100; // 0b0001100100
        const width = 10;

        const expected = 152; // 0b0010011000
        const result = bitReverseNumber(number, width);

        try testing.expectEqual(expected, result);
    }

    {
        const number = 100;
        const width = 12;

        const number_1 = number; // 0b000001100100
        const expected_1 = 608; // 0b001001100000
        const result_1 = bitReverseNumber(number_1, width);
        try testing.expectEqual(expected_1, result_1);

        const number_2 = result_1; // 0b001001100000
        const expected_2 = number; // 0b000001100100
        const result_2 = bitReverseNumber(number_2, width);
        try testing.expectEqual(expected_2, result_2);
    }
}
