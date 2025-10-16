const std = @import("std");

pub fn main() !void {
    // --- FixedBufferAllocator ---
    // var buffer: [1024]u8 = undefined;
    var buffer: [1024]u8 = @splat(0);

    var fba = std.heap.FixedBufferAllocator.init(&buffer);
    const allocator = fba.allocator();

    const memory = try allocator.alloc(u8, 100);
    defer allocator.free(memory);

    for (memory, 0..) |*item, idx| {
        item.* = @intCast(idx);
    }

    for (memory[0..10], 0..) |value, idx| {
        std.debug.print("memory[{d}] = {d}\n", .{ idx, value });
    }

    std.debug.print("{any}\n", .{buffer});

    // --- SIMD ---
    const v1 = @Vector(4, u32){ 4, 12, 37, 9 };
    const v2 = @Vector(4, u32){ 10, 22, 5, 12 };
    const v3 = v1 + v2;

    std.debug.print("{any}\n", .{v3});

    // --- Comptime String ---
    // const string =
    //     \\ Hello
    //     \\
    //     \\ World
    // ;
    // const string: [*:0]u8 = "uint n = " ++ "42";
    const shader_code = @embedFile("src/ntt.wgsl");
    const string: [*:0]const u8 = //
        "const foo : i32 = 24;\n" ++
        "const bar : i32 = 42;\n" ++
        "\n" ++ shader_code;

    std.debug.print("{s}\n", .{string});
}
