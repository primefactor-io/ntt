//! This module is the root module which exports functionalities of all other
//! modules.

const std = @import("std");
const testing = std.testing;

pub const ntt = @import("ntt.zig");
pub const fft = @import("fft.zig");
pub const utils = @import("utils.zig");

pub const gpu = @import("gpu.zig");
pub const ggpu = @import("gpu/gpu.zig");

test {
    @import("std").testing.refAllDecls(@This());
}
