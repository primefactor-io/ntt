const std = @import("std");
const testing = std.testing;

fn mul(a: i32, b: i32) i32 {
    return a * b;
}

test "mul" {
    try testing.expectEqual(6, mul(2, 3));
}
